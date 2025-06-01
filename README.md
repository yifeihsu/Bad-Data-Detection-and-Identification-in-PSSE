# Bad Data Detection and Identification in PSSE

This repository provides MATLAB examples for state estimation and bad data analysis in power systems.  The code relies on the [MATPOWER](https://matpower.org/) toolbox and includes routines for detecting parameter and topology errors as well as correcting erroneous line data.

## Features

- Weighted least squares (WLS) state estimator (`SE.m`)
- Lagrangian multiplier methods for error detection (`LagrangianM.m`, `LagrangianMtopo.m`)
- Weighted least absolute value estimator (`WLAV.m`) using YALMIP/Gurobi
- Parameter correction and multi-scan identification (`main_pe_correction.m`, `correct_parameter_group_multi_scan.m`)
- Example workflows for IEEE test systems (`Lag14Test.m`, `multi_ge_118_angles.m`, etc.)
- Utility for exporting MATPOWER cases to OpenDSS (`exportToOpenDSS.m`)

## Requirements

- MATLAB with the MATPOWER package available on the path
- [YALMIP](https://yalmip.github.io/) and a compatible solver (e.g., Gurobi) for the WLAV example
- Example data such as `loadwotime.xlsx` and the MATPOWER case files (e.g., `case9`, `case118`)

## Usage

Add this repository to your MATLAB path and run any of the demonstration scripts. Typical entry points include:

```matlab
traditional_process      % Basic NLM pipeline
Lag14Test                % Parameter and topology detection on IEEE 14 bus
multi_ge_118_angles      % Multi-error example on IEEE 118 bus
```

Several helper functions compute Jacobians, generate measurements, and apply corrections.  Refer to the comments in each script for details on the workflow.

## File Overview

- `LagrangianM.m` – Parameter error detection via the Lagrangian multiplier method
- `LagrangianMtopo.m` – Topology error detection
- `WLAV.m` – Weighted least absolute value estimator
- `main_pe_correction.m` – Example of parameter error correction
- `test_grouped_residuals*.m` – Demonstrations of grouped residual analysis
- `exportToOpenDSS.m` – Convert MATPOWER cases to OpenDSS format

## License

This project is licensed under the Apache 2.0 license. See `LICENSE.md` for details.
