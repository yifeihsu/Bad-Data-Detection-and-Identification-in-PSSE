# GCN-BDI
1. traditional_process: 传统NLM pipeline.
   - SE: 考虑diverge情况下的。
2. resi_para_error: 参数错误的NLM pipeline。
3. LagrangianM: 参数NLM法。
4. LagrangianMtopo: 拓扑NLM法。
5. maxM_topo与maxindex: 配合拓扑NLM法确定最大的乘子（移除NaN项）。
6. modelselection: 选择不同运行模式。
7. Lagtest4: 测试量测样本的参数错误辨识准确率（有错误/无错误/是否错误识别）。