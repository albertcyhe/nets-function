NET‑F vs NET‑M（ΔFM）脑转移特异性失衡：从假设到推论（含可复现指引）

假设（Hypothesis）
- 我们提出 ΔFM = z(F) − z(M) 作为量化 NETs“功能（F）/标记（M）”失衡的指标，其中：
  - NET‑F：功能端，选用 ELANE/PRTN3/CTSG（避免上游 NOX 等产生 ROS 的组件混入）。
  - NET‑M：标记端，以 PADI4 为锚，并在中性粒富集样本中做共表达扩展（v2，互斥于 F）。
- 预期：ΔFM 低（F<M）不是普遍存在，而是“脑转移微环境”的特异性状态，并与 THBS1/ECM/血管生成等外生长相关轴耦合。

设计（Design）
- Bulk（主验证）：
  - 配对转移（优先）：GSE184869（RNA‑seq，已提供 log2 TMM CPM）、GSE125989（Affy GPL571）。
  - 转移灶（脑/非脑）：GSE43837（脑），GSE14017/14018。
  - 协变量：中性粒浸润（MCPcounter）、肿瘤纯度（tidyestimate/ESTIMATEScore→TumorPurity），批次/平台分层；不跨平台混合。
  - 评分法：singscore 与 ssGSEA（GSVA 2.2 新 API，param 模式；对重复行名按“列最大值”聚合，避免行名冲突），ΔFM 统一按样本内 z 标准化后计算。
- scRNA（机制承载）：
  - GSE186344 多肿瘤脑转移单细胞，使用 10x 原始计数（.h5/.h5.gz、或 mtx 三件套；优先 h5）。
  - 中性粒严格门控（高 NEUT ≥ q，低 TNK/Mono ≤ 60% 分位），在中性粒内计算 ΔFM.cell 并聚类锁定“低ΔFM”亚群，输出样本层比例；
  - 肿瘤like（EPCAM/KRT*）细胞内计算 THBS1 中位量，与“低ΔFM 中性粒亚群比例”做样本层相关。

验证（Validation）
- ΔFM 配对差异（met−primary）：
  - GSE184869（RNA‑seq）：singscore Wilcoxon V=14, p=2.10e−4, FDR=3.22e−4，Cliff’s δ≈−0.70，n=20；ssGSEA 方向与显著性一致。
  - GSE125989（Affy）：V=91, p=0.252, FDR=0.298，Cliff’s δ≈0.25，n=16（不显著）。
  - 图：results/figures/Figure2_GSE184869_Paired_DFM.pdf（含 GSE125989）。
- ΔFM 与 THBS1 偏相关（控制 Neutrophil + TumorPurity）：
  - GSE184869：singscore ρ=−0.3004, p=0.00415, FDR≈0.00831, n=90；ssGSEA ρ=−0.2732, p=0.00936, FDR≈0.00936，方向一致。
  - GSE125989：ρ=0.2269, p=0.2109, FDR≈0.3149, n=32（不显著）。
  - 图：results/figures/Figure3_GSE184869_THBS1_Residuals.pdf（残差散点，ΔFM 与 THBS1 均残差化后相关）。
- F/M 相关性（Spearman）跨队列总览（5 队列 × 2 方法）：
  - GSE184869：singscore ρ=−0.2081, p=0.0491, FDR=0.098；ssGSEA ρ=−0.1216, p=0.2537。
  - GSE125989：ρ=0.4941, p=0.00449, FDR=0.01495；ssGSEA ρ=0.4025, p=0.02314。
  - GSE43837：ρ≈0.07–0.11，p>0.5；GSE14017：ρ≈0.23–0.26，p>0.18；GSE14018：ρ≈0.65–0.67，p≪0.001。
  - 图：results/figures/Figure4_FM_Corr_Overview.pdf；数据：results/tables/fm_correlation_overview.tsv。
- scRNA（GSE186344，Colab 运行，严格门控，高 NEUT≥0.70 且低 TNK/Mono≤0.60；肿瘤like 上皮 q=0.60）：
  - 低ΔFM 中性粒亚群比例（样本层）：如 GSM5645901≈0.0188，GSM5645905≈0.0435；另两样本≈0。
  - THBS1 关联：当前覆盖 n≈4，THBS1 中位值跨样本变异度不足，相关系数未定义（constant input），已在 results/figures/Figure5_scRNA_THBS1_Assoc.txt 记录。
  - 图：results/figures/Figure5_scRNA_LowDFM_Fraction.pdf；表：results/tables/gse186344_neutrophil_lowDFM_fraction.tsv 与 gse186344_subpop_thbs1_assoc.tsv。

反思（Reflection）
- scRNA 的“负结果”并未否定 Bulk 的强证据；反而提示“低ΔFM”更可能是“细胞外、由 NETs 破裂后残留物定义”的稳态微环境信号（NET‑M 标记端成分相对 NET‑F 功能端更持久）。标准 scRNA 仅捕捉细胞内转录，无法直接记录这一“细胞外/非转录”的持久信号。
- 这解释了：Bulk 层面 ΔFM 与 THBS1/ECM/外生长的耦合，以及单细胞层面难以在“活细胞转录”中稳定复现。

推论（Inference）与后续验证
- 机制模型（Figure6_Mechanism_Model.pdf）：NETs 解聚后，NET‑M（PADI4 相关）在组织微环境中相对 NET‑F（弹性蛋白酶复合体）更稳定，形成“低ΔFM 微环境”，与 THBS1/ECM/血管生成轴耦合，促外生长。
- 建议实验：
  - 病理/空间蛋白：IHC/IF 与空间蛋白组标记 NET‑M 残留、THBS1/TSP‑1 与 ECM/血管生成共定位。
  - 蛋白组/质谱：NETs 降解产物与抑制因子平衡（SERPINA1/B1）与 ΔFM 的一致性；
  - 多模态单细胞：CITE‑seq/空间蛋白捕捉“非转录”信号。

图表生成与路径
- 运行：`Rscript scripts/r/07_figures_story.R`
- 输出：`results/figures/`
  - Figure1_DFM_Schematic.pdf：ΔFM 概念示意（ASCII 标签，避免字体编码问题）。
  - Figure2_GSE184869_Paired_DFM.pdf（含 GSE125989）：配对 ΔFM 连线。
  - Figure3_GSE184869_THBS1_Residuals.pdf：ΔFM 与 THBS1 偏相关（残差化）散点。
  - Figure4_FM_Corr_Overview.pdf：F/M 相关性总览热图。
  - Figure5_scRNA_LowDFM_Fraction.pdf 与 Figure5_scRNA_THBS1_Assoc.txt：单细胞低ΔFM 亚群比例与 THBS1 关联结论。
  - Figure6_Mechanism_Model.pdf：机制模型示意（NET‑M 残留 > NET‑F）。

可复现（本地/Colab）
- 本地（严格门控已同步参数化）：
  - 中性粒严格门控与亚群：
    - `python3 scripts/py/sc_neutrophil_subsets.py --h5ad data/processed/scRNA/gse186344_colab.h5ad --net_f resources/modules/net_f_v1.tsv --net_m resources/modules/net_m_v2.tsv --outdir results/tables --neut_quantile 0.70`
  - THBS1 关联（可调上皮门控）：
    - `python3 scripts/py/sc_tumor_thbs1_assoc.py --h5ad data/processed/scRNA/gse186344_colab.h5ad --frac results/tables/gse186344_neutrophil_lowDFM_fraction.tsv --out results/tables/gse186344_subpop_thbs1_assoc.tsv --epi_quantile 0.60`
- Colab（推荐用于海量 GSM 下载与处理；支持 .h5.gz、aria2c 并行）：
  - `!python colab/gse186344_colab.py --gsm <12个GSM> --project_dir "/content/drive/MyDrive/NetsAnalysisProject" --neut_quantile 0.70 --epi_quantile 0.60`

原始技术指引与数据清单（保留如下）
NET‑F vs NET‑M 分析计划与交付（本地执行指引）

概览
- 目标：在多队列（bulk + scRNA）统一框架下构建并比较 NET‑F（ELANE/PRTN3/CTSG）与 NET‑M（以 PADI4 为锚）的活性评分，形成核心指标 ΔFM = z(score_F) − z(score_M)。
- 证据链：在转移灶层面验证 ΔFM 与上游/中游/下游通路与表型的共变关系，并在 scRNA 鉴定 PADI4^low & ELANE/PRTN3^high 亚群；进行鲁棒性与对照检验。
- 执行顺序：配对转移优先（GSE184869/GSE125989）→ 转移灶 bulk（GSE43837/GSE14017/GSE14018）→ 单细胞（GSE186344）。

本地已下载数据（根目录）
- 平台注释：`GPL96.annot.gz`, `GPL570.annot.gz`, `GPL571.annot.gz`, `GPL1352.annot.gz`
- 配对/转移 bulk：
  - `GSE184869_rna_seq_batch_corrected_log2_TMM_normalised_CPM_protein_coding_genes.xlsx`
  - `GSE125989_series_matrix.txt.gz`, `GSE125989_RAW.tar`
  - `GSE14017_series_matrix.txt.gz`, `GSE14017_RAW.tar`
  - `GSE14018_series_matrix.txt.gz`, `GSE14018_RAW.tar`
  - `GSE43837_series_matrix.txt.gz`, `GSE43837_RAW.tar`
- 参考队列：`TCGA-BRCA.star_fpkm-uq.tsv.gz`, `TcgaTargetGtex_rsem_gene_tpm.gz`, `TCGA-CDR-SupplementalTableS1.xlsx`, `brca_metabric.tar.gz`
- 单细胞：`GSE186344_RAW.tar`

项目目录规划（不移动现有原始文件）
- 原始数据（保持在根目录）：上述 GEO/TCGA/METABRIC 压缩或矩阵文件。
- 处理与结果：
  - `data/processed/`：各数据集标准化/整合后的表达矩阵与元数据（RDS/RData/TSV/H5AD）。
  - `resources/`：基因集合（GMT/TSV）、模块清单、对照基因集、配置 YAML。
  - `scripts/r/`：bulk 分析与统计建模（R）。
  - `scripts/py/`：scRNA 与图谱相关（Python/Scanpy）。
  - `workflow/`：Snakemake 或 targets/drake 流程文件。
  - `results/figures/`：F1–F8 图表（PDF/SVG）。
  - `results/tables/`：主结果表与补充表（CSV/TSV）。
  - `logs/`：运行日志与会话信息（sessionInfo/conda env）。

环境与依赖（建议）
- Conda 环境：
  - R：`limma`, `edgeR`, `sva`, `GSVA`, `singscore`, `estimate`, `MCPcounter`, `ppcor`, `metafor`, `psych`（相关比较），`ggplot2`, `ComplexHeatmap`, `data.table`, `readxl`。
  - Python：`scanpy`, `anndata`, `harmonypy`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `gseapy`, `pyAUCell`（或 R AUCell via reticulate/替代）。
- 版本锁定与记录：`env/environment.yml`（conda），`R/renv.lock`（可选）；运行时导出 `logs/sessionInfo_*.txt`。

模块与通路资源（resources/）
- `resources/modules/net_f_v1.tsv`：ELANE, PRTN3, CTSG（功能端种子）。
- `resources/modules/net_m_v1.tsv`：PADI4（标记端锚点）。
- `resources/modules/net_f_v2.tsv` 与 `net_m_v2.tsv`：扩展版（相关网络筛选，互斥）。
- `resources/pathways/msigdb_hallmark.gmt`：HALLMARK 集合（TNFA/NFKB, IL6/JAK/STAT3, ROS, G2M, E2F, ANGIOGENESIS, EMT）。
- `resources/pathways/reactome_go_p38_il1.gmt`：p38 MAPK 与 IL‑1 响应。
- `resources/controls/`：随机同大小集合、无关粒细胞基因集、抑制剂平衡对（ELANE:SERPINA1, PRTN3:SERPINB1）。

统一评分与核心指标
- Bulk 评分：并行计算 `singscore` 与 `ssGSEA(GSVA)`；对每个数据集内做 z 标准化。
- scRNA 评分：`AUCell` + `AddModuleScore`（Scanpy/Seurat 等价）并做深度校正；样本层聚合。
- 核心量：`ΔFM = z(score_F) − z(score_M)`；对比 `score_M` 的增益（ΔR²、ΔAIC、偏相关）。
- 协变量：中性粒浸润（MCP-counter/xCell/CIBERSORT 选1–2）、肿瘤纯度（ESTIMATE 或现成字段）、批次/平台（阵列：RMA+ComBat；RNA‑seq：TMM/TPM+voom；不跨平台合并）。

分阶段执行步骤（可复现流水线）
1) 数据准备与标准化（bulk）
   - 输入：根目录的 `GSE*_series_matrix.txt.gz`、`*_RAW.tar`、`GSE184869_*.xlsx`、`TCGA*`、`METABRIC`。
   - 处理：
     - 阵列：RMA（raw CEL 在 `*_RAW.tar`；若只用 series matrix 则直接 log2 标准化与批次标注）、按 GPL 注释统一基因符号、ComBat 去批。
     - RNA‑seq：TMM/TPM + voom；`GSE184869` 已给 log2 TMM CPM，可直接入模但保留批次字段。
     - 生成样本元数据：样本类型（原发/转移/脑转移/部位）、配对 ID、平台/批次、PAM50（TCGA/METABRIC）、是否肿瘤细胞比例字段。
   - 输出：
     - `data/processed/<dataset>/<dataset>.expr.tsv.gz`（gene × sample）。
     - `data/processed/<dataset>/<dataset>.pheno.tsv`（样本信息）。
     - `data/processed/<dataset>/<dataset>.platform.tsv`（平台/批次/注释）。

2) 协变量估计（bulk）
   - 输入：步骤1输出表达矩阵。
   - 处理：MCP‑counter/xCell/CIBERSORT 计算中性粒分数；ESTIMATE 计算 Stromal/Immune/Tumor purity。
   - 输出：`data/processed/<dataset>/<dataset>.covariates.tsv`（Neutrophil, Purity, Batch, PAM50 等）。

3) 模块扩展与基因集准备
   - 输入：种子模块、转移灶（优先）中性粒富集样本/集群的表达矩阵。
   - 处理：对种子基因做相关网络，筛选稳健共表达基因（阈值建议 |ρ|>0.3，跨样本/集群一致；互斥 NET‑F/NET‑M），形成 v2。
   - 输出：`resources/modules/net_f_v2.tsv`, `resources/modules/net_m_v2.tsv`；并生成版本说明 `resources/modules/README.tsv`。

4) 评分计算与 ΔFM（bulk）
   - 输入：步骤1表达矩阵、步骤3基因集。
   - 处理：singscore 与 ssGSEA 并行；z 标准化；生成 `score_F`, `score_M`, `ΔFM`；同时计算通路分数（MSigDB/Reactome）。
   - 输出：
     - `data/processed/<dataset>/<dataset>.scores.tsv`（样本 × 指标：score_F/score_M/ΔFM/通路）。
     - `results/tables/<dataset>_scores_summary.tsv`（描述统计）。

5) 转移灶主线检验（bulk）
   - 配对差异：同例转移 vs 原发的 ΔFM（配对检验 + Cliff’s δ）。
   - 分组差异：脑转移 vs 其他部位（Mann‑Whitney/线性模型）。
   - 分段相关：上游（CTSC/PRTN3/IL1B、NF‑κB/IL‑1）、中游（p38/ROS）、下游（THBS1↓、ECM/血管/增殖↑、休眠↓）。
   - 部分相关/多变量：控制中性粒浸润、纯度、PAM50、批次，比较 ΔFM 与 score_M 的独立效应与增益（Steiger/Meng 相关比较）。
   - 输出：
     - `results/tables/<dataset>_paired_tests.tsv`, `<dataset>_group_tests.tsv`。
     - `results/tables/<dataset>_segment_cor_partial.tsv`（ρ/β、FDR、95%CI）。
     - `results/tables/<dataset>_model_compare.tsv`（ΔR²/ΔAIC/偏相关）。

6) 解耦与对照检验（bulk）
   - 残差相关：`res_F = lm(score_F ~ Neutrophil + Purity + PAM50 + Batch)` 残差 vs `res_M` 残差；跨队列森林图。
   - 对照：随机同大小集、无关粒细胞集，基因置换 1,000 次经验 p 值；抑制剂平衡（ELANE:SERPINA1, PRTN3:SERPINB1）。
   - 输出：`results/tables/<dataset>_residual_corr.tsv`, `results/tables/<dataset>_null_perm.tsv`，以及元分析表 `results/tables/meta_*`。

7) 单细胞（GSE186344）
   - 预处理：QC（空滴、线粒体比例、双tscore）、归一化、HVG、批次整合（Harmony/CCA）、注释中性粒（S100A8/A9, CSF3R, FCGR3B, CXCR2, MPO）。
   - 模块评分：在中性粒细胞计算 NET‑F/NET‑M/ΔFM.cell（AUCell + module score），深度校正；Leiden/Louvain 聚类。
   - 亚群：筛选 PADI4^low & ELANE/PRTN3^high；样本层计算该亚群比例；关联外生长 proxy（同队列 bulk 或肿瘤细胞增殖/ECM/血管分数）。
   - 输出：
     - `data/processed/scRNA/gse186344.h5ad`（整合对象）。
     - `results/tables/gse186344_neutrophil_cluster_markers.tsv`，`gse186344_sample_subpop_fraction.tsv`。

8) 图表与最终交付
   - F1 流程图：`results/figures/F1_workflow.pdf`。
   - F2 配对 ΔFM：`results/figures/F2_paired_deltaFM.pdf`。
   - F3 脑转移 vs 其他部位：`results/figures/F3_site_deltaFM.pdf`。
   - F4 分段偏相关热图：`results/figures/F4_segment_partial_corr.pdf`。
   - F5 THBS1 vs ΔFM：`results/figures/F5_thbs1_deltaFM_regression.pdf`。
   - F6 残差相关森林图：`results/figures/F6_residual_corr_forest.pdf`。
   - F7 scRNA UMAP 与亚群：`results/figures/F7_scRNA_umap_subpop.pdf`；样本比例相关：`F7b_subpop_fraction_assoc.pdf`。
   - F8 鲁棒性与负对照：`results/figures/F8_robustness_controls.pdf`。
   - 附表：模块基因表、协变量定义、各队列 QC 指标：`results/tables/Supp_*`。

MVP 优先包（先行完成）
- GSE184869/GSE125989：验证转移灶 ΔFM↑（配对），与 THBS1↓ / 增殖/ECM↑ 偏相关（控浸润/纯度）。
- GSE186344：定义 PADI4^low & ELANE/PRTN3^high 亚群，样本比例与外生长 proxy 关联。
- 对应最少生成文件：
  - `data/processed/GSE184869/*.tsv`, `data/processed/GSE125989/*.tsv`。
  - `results/tables/GSE184869_*`, `GSE125989_*`；`results/figures/F2/F5`。
  - `data/processed/scRNA/gse186344.h5ad`；`results/tables/gse186344_sample_subpop_fraction.tsv`；`results/figures/F7*`。

代码组织与脚本入口（示例，不立即执行）
- R（scripts/r/）
  - `00_utils_io.R`：IO/注释/对齐基因符号。
  - `01_preprocess_bulk.R`：阵列 RMA/ComBat、RNA‑seq TMM/voom；导出 expr/pheno。
  - `02_covariates_scores.R`：MCP‑counter/ESTIMATE、singscore、GSVA、通路分数、ΔFM。
  - `03_stats_bulk.R`：配对/分组差异、分段偏相关、模型比较、残差相关、置换与对照。
  - `04_meta_analysis.R`：跨队列元分析与森林图。
  - `05_figures_bulk.R`：F2–F6–F8 图表。
- Python（scripts/py/）
  - `sc_preprocess_gse186344.py`：QC、整合、注释中性粒、保存 H5AD。
  - `sc_scores_clusters.py`：AUCell + module score、ΔFM.cell、聚类与亚群识别、样本比例计算。
  - `sc_figures.py`：UMAP、分数分布、样本比例相关图。
- 配置与流程
  - `resources/config.yml`：数据集路径、平台/批次字段、协变量方法、评分参数、阈值。
  - `workflow/Snakefile` 或 `workflow/_targets.R`：串联步骤 1–8；产出到 `results/`。

统计与阈值（统一标准）
- 双侧检验，FDR<0.05；报告效应量（Cliff’s δ/ρ/β）与 95% CI。
- 元分析：随机效应（DL 法），I²<50% 视为一致；异质性高使用方向一致 + 过半显著。
- 相关比较：Steiger/Meng；置换 1,000 次经验 p 值。

数据集要点与备注
- `GSE184869`：已批次校正的 log2 TMM CPM（蛋白编码），可直接评分；仍需提取配对与部位信息。
- `GSE125989/GSE14017/GSE14018/GSE43837`：若使用 `*_RAW.tar`，需 RMA；若先用 series matrix，可直接进入评分但保留平台/批次。
- `GSE186344`：scRNA RAW 包体积较大，建议先解压到临时目录，输出整合对象到 `data/processed/scRNA/gse186344.h5ad`。
- `TCGA-BRCA` 与 `METABRIC`：用于亚型/协变量与参照，不参与原发预测主线合并；不跨平台混合表达矩阵。

运行顺序（建议）
1. 先完成 MVP：GSE184869/GSE125989（bulk）→ GSE186344（scRNA）。
2. 扩展到 GSE43837/GSE14017/GSE14018（转移灶 bulk）。
3. 整体元分析与鲁棒性网格（评分方法 × 模块版本 × 协变量组合）。

产出清单对标（F1–F8）
- 详见“图表与最终交付”章节的文件路径；所有图表 PDF/SVG 同步导出，表格 CSV/TSV 同步导出。

注意与风险缓解
- 中性粒计数混杂：统一纳入 Neutrophil 分数与纯度；报告残差与偏相关。
- 平台异质性：不跨平台合并表达；队列内得效应量后做元分析。
- scRNA 中性粒易丢失：严格 QC、样本层聚合，并可引入外部同癌种/同部位 scRNA 方向验证。
- 标记/功能边界：保留 v1 种子与 v2 扩展两版，主文以 v1 种子为主，扩展为鲁棒性验证。

附：常用命令（示例，仅作为执行参考）
- 创建目录：`mkdir -p data/processed resources/modules resources/pathways resources/controls scripts/r scripts/py workflow results/figures results/tables logs`
- 生成评分（R）：`Rscript scripts/r/02_covariates_scores.R --dataset GSE184869 --config resources/config.yml`
- scRNA 预处理（Py）：`python scripts/py/sc_preprocess_gse186344.py --raw GSE186344_RAW.tar --out data/processed/scRNA/gse186344.h5ad`

**执行进展与思路更新（Bulk 收官 → scRNA 启动）**
- **GSVA 修复与统一口径**: 适配 GSVA>=2.2 新 API（`ssgseaParam` + `gsva(param)`），在评分与通路脚本内统一处理；对表达矩阵行名去重（按样本列取最大值）并实现 Ensembl→Symbol 映射与聚合，避免重复行名报错。
- **协变量（官方口径）**: 优先使用 `MCPcounter` 与 `tidyestimate` 计算协变量，生成 `Neutrophils` 与 `TumorPurity`；RNA‑seq 队列在 `tidyestimate` 仅输出 `ESTIMATEScore` 时，按 ESTIMATE 原公式推导 `TumorPurity`。见 `data/processed/<ds>/<ds>.covariates.tsv`。
- **模块与 ΔFM 稳定性**: 在 GSE125989 的中性粒富集样本中，以 PADI4 为锚做共表达扩展，生成 `resources/modules/net_m_v2.tsv`（互斥于 NET‑F）；统计脚本统一“现场重算 ΔFM”，避免常数化偏差。
- **MVP 关键结果（配对与偏相关）**:
  - `GSE184869`（RNA‑seq，脑转移配对）：ΔFM（met−primary）显著下降（singscore: V=14, p=2.10e−4, FDR=3.22e−4；ssGSEA 类似，Cliff’s δ≈−0.70，n=20）；偏相关 ΔFM ↔ THBS1（控 Neutrophil + TumorPurity）显著负相关（singscore: ρ=−0.3004, FDR≈0.0083；ssGSEA: ρ=−0.2732, FDR≈0.0094，n=90）。
  - `GSE125989`（Affy，配对）：配对差异不显著；偏相关与 THBS1 未达显著（FDR>0.2）。
  - 结果表均已标注来源：`neutrophil_source` 与 `purity_source`，见 `results/tables/*_paired_tests.tsv`, `*_thbs1_partial_cor.tsv`。
- **F/M 相关性（解耦证据）**: 新增总览表 `results/tables/fm_correlation_overview.tsv`（5 队列 × 2 方法）。`GSE184869` 中 `F/M` 相关性最低且为负（singscore: ρ≈−0.21, FDR≈0.10），而多数其他队列为正相关（如 `GSE14018`：ρ≈0.65–0.67, FDR≪0.01），支持“脑转移微环境特异性解耦”的核心假说。
- **通路资源**: 通过 `msigdbr` 缓存 `resources/pathways/msig_*.rds`，并导出筛选清单 `resources/pathways/selected_sets.tsv`（Reactome p38 与 GO:BP IL‑1 相关集合）。各队列通路分数见 `data/processed/<ds>/<ds>.pathways.tsv`。
- **关键脚本更新**:
  - 评分与 ΔFM：`scripts/r/03_scores_deltafm.R`（GSVA 新 API、去重、映射、并行 singscore/ssGSEA）。
  - 协变量：`scripts/r/02_covariates_estimate.R`（MCPcounter + tidyestimate，RNA‑seq 纯度推导，代理回退）。
  - 统计（MVP）：`scripts/r/04_bulk_stats_mvp.R`（配对检验、THBS1 偏相关、来源注释、现场重算 ΔFM）。
  - 通路评分：`scripts/r/05_pathway_scores.R`（ssGSEA + 缓存基因集）。
  - F/M 相关性总览：`scripts/r/06_fm_correlation_overview.R`。

**Step 3｜单细胞（GSE186344）执行计划与脚本**
- **预处理与细胞鉴定**（阶段 1）
  - 脚本：`scripts/py/sc_preprocess_gse186344.py`
  - 输入：`GSE186344_RAW.tar`（或解压根目录）
  - 处理：自动发现多个 10x 样本，合并；QC（min_genes≥200、线粒体<20%）、Normalize/log1p、HVG、PCA、Harmony(by sample, 若可用)、邻域/UMAP/Leiden。
  - 输出：`data/processed/scRNA/gse186344.h5ad`
- **中性粒与低ΔFM亚群**（阶段 2）
  - 脚本：`scripts/py/sc_neutrophil_subsets.py`
  - 动作：中性粒标记打分→门控中性粒→计算 NET‑F/NET‑M 分数与 ΔFM.cell→在中性粒内聚类并锁定“NET‑M高 & NET‑F低”的低ΔFM 亚群；输出样本层亚群比例与 markers。
  - 输出：
    - `results/tables/gse186344_neutrophil_lowDFM_fraction.tsv`
    - `results/tables/gse186344_lowDFM_markers.tsv`
- **样本层关联（关键验证）**（阶段 3）
  - 脚本：`scripts/py/sc_tumor_thbs1_assoc.py`
  - 动作：识别上皮/肿瘤like 细胞（EPCAM/KRT*）、计算样本层 THBS1 中位表达；与低ΔFM 中性粒亚群比例做 Spearman 相关。
  - 输出：`results/tables/gse186344_subpop_thbs1_assoc.tsv`
- **细胞通讯（可选机制）**（阶段 4）
  - 方案：后续可导出对象到 CellChat/NicheNet/LIANA，验证与肿瘤/内皮细胞的互作轴是否与外生长相关。

**Step 3｜一键运行命令（建议顺序）**
- 预处理：`python3 scripts/py/sc_preprocess_gse186344.py --raw GSE186344_RAW.tar --out data/processed/scRNA/gse186344.h5ad`
- 中性粒与亚群：`python3 scripts/py/sc_neutrophil_subsets.py --h5ad data/processed/scRNA/gse186344.h5ad --net_f resources/modules/net_f_v1.tsv --net_m resources/modules/net_m_v2.tsv --outdir results/tables`
- 样本层关联：`python3 scripts/py/sc_tumor_thbs1_assoc.py --h5ad data/processed/scRNA/gse186344.h5ad --frac results/tables/gse186344_neutrophil_lowDFM_fraction.tsv --out results/tables/gse186344_subpop_thbs1_assoc.tsv`

**Step 3｜环境建议（避免 NumPy ABI 冲突）**
- 推荐使用 conda 环境（已提供锁定文件）：
  - 创建：`conda env create -f env/environment-scrna.yml`
  - 激活：`conda activate nets-scrna`
  - 运行上述 Step 3 三个脚本
- 如果使用 pip，请确保 NumPy<2 与相容版本：
  - `pip install 'numpy<2' 'scanpy==1.9.8' 'anndata==0.9.2' 'pandas<2.2' 'numexpr<2.9' 'zarr<2.17' 'numcodecs<=0.12.1'`
