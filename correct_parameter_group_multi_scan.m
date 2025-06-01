function [corrected_params_group, success_correction] = ...
    correct_parameter_group_multi_scan(mpc_with_error, line_to_correct_idx, ...
                                       multi_scan_measurements_z, ...
                                       initial_states_multi_scan, R_variances_vec, baseMVA)
%CORRECT_PARAMETER_GROUP_MULTI_SCAN
%  Performs multi-scan augmented-state estimation to correct R and X
%  of a specific line (line_to_correct_idx).
%
%  Inputs:
%    - mpc_with_error: MATPOWER struct with the *erroneous* R,X in the branch matrix
%    - line_to_correct_idx: scalar index of the line whose R,X we want to correct
%    - multi_scan_measurements_z: matrix of size (nz_actual x num_scans)
%    - initial_states_multi_scan: each column is an initial guess of [V;theta(deg)]
%        for that scan (assumes the first nb rows are V, the next nb rows are angle in deg)
%    - R_variances_vec: vector of measurement variances (length nz_actual)
%    - baseMVA: system base MVA
%
%  Outputs:
%    - corrected_params_group = [R_est; X_est]
%    - success_correction = 1 if converged, 0 otherwise

    fprintf('    Entering correct_parameter_group_multi_scan for line %d...\n', line_to_correct_idx);

    %% 1) Basic Setup
    nb = size(mpc_with_error.bus, 1);
    nl = size(mpc_with_error.branch, 1);
    s  = size(multi_scan_measurements_z, 2);  % Number of scans
    success_correction = 0;

    % Current param guesses (start from the erroneous ones in mpc_with_error)
    p_g_current = [
        mpc_with_error.branch(line_to_correct_idx, 3);  % R
        mpc_with_error.branch(line_to_correct_idx, 4);  % X
    ];
    corrected_params_group = p_g_current;  % default return value

    % Build the per-scan measurement covariance (R_single_scan), then make block diagonal
    nz_single = size(multi_scan_measurements_z, 1);  % = 3*nb + 4*nl in your case
    R_single_scan = spdiags(R_variances_vec, 0, nz_single, nz_single);
    R_inv_single  = inv(R_single_scan);
    R_v_inv       = kron(speye(s), R_inv_single);  %#ok<NASGU> (if you need the big block)

    % We remove the reference bus angle => each scan has (2*nb - 1) states.
    % The augmented dimension = s*(2*nb -1) for states + 2 for [R,X].
    n_states_per_scan = 2*nb - 1;
    N = s*n_states_per_scan + 2;  % dimension of the augmented vector

    % Prepare initial guess for the augmented vector x_v:
    % We'll store (theta minus ref, V for all buses) for each scan + param(2).
    x_states_all_scans_current = zeros(s*n_states_per_scan, 1);

    % Identify the reference bus
    ref_bus = find(mpc_with_error.bus(:,2) == 3);
    if isempty(ref_bus)
        warning('No type-3 (slack) bus found; defaulting ref_bus=1');
        ref_bus = 1;
    end

    %% 2) Build the initial augmented state x_v_current
    for k_scan = 1:s
        % initial_states_multi_scan(:,k_scan) => first nb are V, next nb are angle (deg)
        v_init     = initial_states_multi_scan(1:nb, k_scan);
        a_init_deg = initial_states_multi_scan(nb+1 : 2*nb, k_scan);
        a_init_rad_full = a_init_deg * (pi/180);

        % Remove reference angle
        a_sub = a_init_rad_full;
        a_sub(ref_bus) = [];  % drop that angle => (nb-1) remain

        % Our reduced state for this scan is [a_sub; v_init], length=2*nb -1
        x_k_scan = [a_sub; v_init];

        % Insert into the big augmented vector
        idx_start = (k_scan-1)*n_states_per_scan + 1;
        idx_end   = k_scan*n_states_per_scan;
        x_states_all_scans_current(idx_start : idx_end) = x_k_scan;
    end

    % Append the 2 line parameters [R, X] at the end
    x_v_current = [x_states_all_scans_current; p_g_current];

    %% 3) Iterative ASE to correct R and X
    max_iter_corr = 20;
    tol_corr      = 5e-3;

    fprintf('    Starting iterative correction for line %d with %d scans...\n', ...
            line_to_correct_idx, s);

    for iter = 1:max_iter_corr

        % Build the big gain matrix G_v (N x N) and the big residual vector (N x 1)
        G_v = sparse(N, N);
        RHS = zeros(N,1);

        for k_scan = 1:s
            %% 3a) Extract sub-vector of states for this scan
            idx_start = (k_scan-1)*n_states_per_scan + 1;
            idx_end   = k_scan*n_states_per_scan;
            x_scan_k  = x_v_current(idx_start : idx_end);  % length (2*nb -1)

            %% 3b) Expand it back to full dimension => (nb angles, nb volts)
            [theta_full, V_full] = expand_state_vector(x_scan_k, nb, ref_bus);

            %% 3c) Overwrite the line_to_correct_idx in a copy of mpc_with_error
            mpc_iter = mpc_with_error;
            mpc_iter.branch(line_to_correct_idx, 3) = x_v_current(s*n_states_per_scan + 1); % R
            mpc_iter.branch(line_to_correct_idx, 4) = x_v_current(s*n_states_per_scan + 2); % X

            %% 3d) Compute predicted measurements h(x,k) and mismatch
            hx_k = calculate_hx(mpc_iter, theta_full, V_full);   % (nz_single x 1)
            z_k  = multi_scan_measurements_z(:, k_scan);         % measured vector
            delta_z_k = z_k - hx_k;                              % mismatch

            %% 3e) Build Jacobian w.r.t. states => H_x(q) using makeJaco
            [Ybus_iter, Yf_iter, Yt_iter] = makeYbus(mpc_iter);
            fbus_iter = mpc_iter.branch(:,1);
            tbus_iter = mpc_iter.branch(:,2);

            % The "full" state is [theta_full; V_full]
            x_full_scan = [theta_full; V_full];
            Vc_full     = V_full .* exp(1j*theta_full);

            % makeJaco returns a (3*nb + 4*nl) x (2*nb) matrix => we then remove ref angle column
            [J_full, ~, ~] = makeJaco(x_full_scan, Ybus_iter, Yf_iter, Yt_iter, ...
                                      nb, nl, fbus_iter, tbus_iter, Vc_full);

            % Remove the reference angle's column => yields size (3*nb + 4*nl) x (2*nb - 1)
            J_full(:, ref_bus) = [];

            H_x_k = J_full;  % rename for clarity

            %% 3f) Build the Jacobian w.r.t. the line parameters => H_p(q)
            %     (only partials w.r.t. [R,X] of that line)
            H_p_k = calculate_param_jacobian_for_line(mpc_iter, line_to_correct_idx, ...
                                                      [theta_full; V_full], ...
                                                      Ybus_iter, Yf_iter, Yt_iter);

            %% 3g) Accumulate in the big G_v and RHS, weighting by R_inv_single
            Wblock  = R_inv_single;   % size (nz_single x nz_single)

            % Precompute for normal equations
            Ht_W    = H_x_k' * Wblock;   % ((2*nb -1) x nz_single)
            Hp_t_W  = H_p_k' * Wblock;   % (2 x nz_single)

            Ht_W_dz    = Ht_W  * delta_z_k;  
            Hp_t_W_dz  = Hp_t_W * delta_z_k;

            % Indices for the state part & param part in the augmented vector
            row_s  = idx_start : idx_end;  
            col_s  = row_s;
            row_p  = s*n_states_per_scan + (1:2);
            col_p  = row_p;

            % (A) G_v(state,state)
            G_v(row_s, col_s) = G_v(row_s, col_s) + (Ht_W * H_x_k);

            % (B) G_v(state,param)
            G_v(row_s, col_p) = G_v(row_s, col_p) + (Ht_W * H_p_k);

            % (C) G_v(param,state)
            G_v(row_p, col_s) = G_v(row_p, col_s) + (Hp_t_W * H_x_k);

            % (D) G_v(param,param)
            G_v(row_p, col_p) = G_v(row_p, col_p) + (Hp_t_W * H_p_k);

            % (E) RHS(state)
            RHS(row_s) = RHS(row_s) + Ht_W_dz;

            % (F) RHS(param)
            RHS(row_p) = RHS(row_p) + Hp_t_W_dz;
        end % end loop over scans

        %% 3h) Solve for delta_x_v
        if rcond(full(G_v)) < 1e-14
            fprintf('    G_v is ill-conditioned at iteration %d. Aborting correction.\n', iter);
            success_correction = 0;
            return;
        end

        delta_x_v = G_v \ RHS;
        x_v_current = x_v_current + 0.2 * delta_x_v;  % update

        %% 3i) Check convergence
        max_update = max(abs(delta_x_v));
        p_g_current = x_v_current(s*n_states_per_scan + (1:2));  % last 2 entries => [R,X]
        fprintf('    Iter %d: max update=%.2e, R_est=%.4f, X_est=%.4f\n', ...
                iter, max_update, p_g_current(1), p_g_current(2));

        if max_update < tol_corr
            success_correction    = 1;
            corrected_params_group = p_g_current;
            fprintf('    Correction converged in %d iterations.\n', iter);
            break;
        end

        if iter == max_iter_corr
            success_correction    = 0;
            corrected_params_group = p_g_current;
            fprintf('    Correction did NOT converge by iteration %d. Returning last estimate.\n', iter);
        end
    end

    fprintf('    Exiting correct_parameter_group_multi_scan. success=%d\n', success_correction);
end