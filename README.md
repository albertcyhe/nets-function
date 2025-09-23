# NET 标记与功能解耦在脑转移中的证据链（含可复现记录）

本 README 以“研究思路 → 执行 → 结果”的顺序，重构并整合先前在 RNA 与 Proteomics/Metabolomics 层完成的全部工作与产出（含失败与修正）。

## 1. 研究动机与假设（标记≠功能，可否解耦？）

- 背景线索：既往研究提示 CTSC 诱导炎症 → NETosis → 释放 NETs → 裂解 TSP1（THBS1）并促进转移；多篇工作将 NETs“标记（NET‑M）”与“功能（NET‑F）”强绑定。
- 核心问题：标记并不等于功能，是否存在“标记与功能可解耦”的场景？在急性炎症期，文献多见“标记＜功能”。若在癌症转移的预防中标记与功能可分离，临床意义与用药分流将截然不同。
- 我们的可量化表述：ΔFM = z(F) − z(M)
  - NET‑F（功能端）：ELANE/PRTN3/CTSG（避免上游 ROS 生成组分混入）。
  - NET‑M（标记端）：以 PADI4 为锚，在中性粒富集样本中做共表达扩展（v2，与 F 互斥）。
- 研究假设：解耦为“脑转移微环境特异”，表现为“标记＞功能”（ΔFM 下降），并与 THBS1/ECM/血管生成等外生长轴耦合。

统计口径与统一阈值
- 检验：配对 Wilcoxon（Cliff’s δ 报效应量）、非配对 Mann‑Whitney/线性模型；相关采用 Spearman + BH FDR（α=0.05）。
- 偏相关/混合模型：控制 Neutrophil 分数与 TumorPurity，分层或加入 (1|cohort/dataset)。
- 报告：双侧检验，FDR<0.05；给出效应量与 95% CI；元分析随机效应（DL），I²<50% 为一致。

小结与下一步
- 结果要点：提出 ΔFM 作为“标记–功能失衡”的定量指标，明确定义 NET‑F/NET‑M 及统计口径，为后续 RNA/蛋白两层验证奠定统一框架。
- 临床意义：若证实“标记＞功能”的脑特异解耦，将直接影响预防性干预与随访分层（功能锚点需单独监测）。
- 下一步：在更多独立队列复核 ΔFM 稳定性，细化 NET‑M v2 模块与阴性对照集合，形成锁定版资源清单。

## 2. 转录组证据：肺不解耦，脑解耦（标记＞功能）

- 数据与设计（Bulk）
  - 配对转移：GSE184869（RNA‑seq，已提供 log2 TMM CPM）、GSE125989（Affy GPL571）。
  - 转移灶分层：GSE43837（脑），GSE14017/14018（非脑）。
  - 协变量：Neutrophils（MCPcounter）、TumorPurity（tidyestimate/ESTIMATEScore→TumorPurity）。
  - 评分：singscore 与 ssGSEA（GSVA≥2.2 新 API param 模式，行名去重聚合）。统一在样本内 z 标准化后计算 ΔFM。

- 主要结果（Bulk）
  - ΔFM 配对差异（met−primary）
    - GSE184869（RNA‑seq）：singscore Wilcoxon V=14, p=2.10e−4, FDR=3.22e−4，Cliff’s δ≈−0.70，n=20；ssGSEA 方向一致。
      - 图：results/figures/Figure2_GSE184869_Paired_DFM.pdf
      - 表：results/tables/GSE184869_paired_tests.tsv
    - GSE125989（Affy）：V=91, p=0.252, FDR=0.298，Cliff’s δ≈0.25，n=16（不显著）。
      - 表：results/tables/GSE125989_paired_tests.tsv
  - ΔFM 与 THBS1 偏相关（控 Neutrophil + TumorPurity）
    - GSE184869：singscore ρ=−0.3004, p=0.00415, FDR≈0.00831, n=90；ssGSEA ρ=−0.2732, p=0.00936, FDR≈0.00936。
      - 图：results/figures/Figure3_GSE184869_THBS1_Residuals.pdf
      - 表：results/tables/GSE184869_THBS1_residuals.tsv（偏相关备份：results/tables/GSE184869_thbs1_partial_cor.tsv）
    - GSE125989：ρ=0.2269, p=0.2109, FDR≈0.3149（不显著）。
  - F/M 相关性跨队列总览（5 队列 × 2 方法）
    - GSE184869：singscore ρ=−0.2081, p=0.0491, FDR=0.098；ssGSEA ρ=−0.1216, p=0.2537。
    - GSE125989：ρ=0.4941, p=0.00449, FDR=0.01495；ssGSEA ρ=0.4025, p=0.02314。
    - GSE43837：ρ≈0.07–0.11；GSE14017：ρ≈0.23–0.26；GSE14018：ρ≈0.65–0.67（p≪0.001）。
      - 图：results/figures/Figure4_FM_Corr_Overview.pdf
      - 表：results/tables/fm_correlation_overview.tsv

- 单细胞（机制承载）GSE186344
  - 严格门控（高 NEUT≥0.70 且低 TNK/Mono≤0.60）后，本次统计未找到稳定的中性粒细胞群（或门控结果为空/近 0），因此无法在中性粒内计算 ΔFM.cell 并锁定“低 ΔFM”亚群。
  - 解释：更符合“中性粒已裂解成 NETs、功能端被抑制”的场景；scRNA 仅捕捉细胞内转录，难以直接记录这种细胞外/降解残留信号。
    - 图：results/figures/Figure5_scRNA_LowDFM_Fraction.pdf（如为占位图请结合日志解读）
    - 记录：results/figures/Figure5_scRNA_THBS1_Assoc.txt（样本量与常量化导致相关不可估）
    - 表：
      - results/tables/gse186344_neutrophil_lowDFM_fraction.tsv（可能为空/近 0，占位记录）
      - results/tables/gse186344_sample_subpop_fraction.tsv（可能为空/近 0）
      - results/tables/gse186344_subpop_thbs1_assoc.tsv（相关不可估的记录表）

- 执行脚本（转录组）
  - R：scripts/r/02_covariates_estimate.R、scripts/r/03_scores_deltafm.R、scripts/r/04_bulk_stats_mvp.R、scripts/r/05_pathway_scores.R、scripts/r/06_fm_correlation_overview.R
  - Py（scRNA/辅助）：scripts/py/sc_preprocess_gse186344.py、scripts/py/sc_neutrophil_subsets.py、scripts/py/sc_tumor_thbs1_assoc.py、scripts/py/sc_figures.py、scripts/py/process_gse184869_excel.py

小结与下一步
- 结果要点：Bulk 层面在脑转移显示 ΔFM 显著下降并与 THBS1 负相关；肺/非脑未见同等强度的解耦。scRNA 在本次严格门控统计中未检出稳定的中性粒细胞群（或门控为空/近 0），因此“低 ΔFM 中性粒亚群比例”不可稳定估计；这与“中性粒已裂解为 NETs、功能端受抑”的解释一致。
- 临床意义：提示“低 ΔFM 微环境”更可能是细胞外/降解残留层面的稳态信号，更适合体液蛋白/肽与空间蛋白读出，而非依赖单细胞转录。
- 下一步：扩充脑/非脑转移队列（GSE43837/14017/14018 全量），完善偏相关/残差图；考虑空间蛋白/空间转录验证；在 scRNA 侧放宽门控或采用含粒细胞保存度更高的数据源做方向性复核。

## 3. 蛋白/肽层证据与机制：Serpin 抑制导致“标记＞功能”（脑特异）

- 挑战：多数公开队列中 NET‑F（ELANE/PRTN3/CTSG）蛋白/肽难以稳定检出。为避免“未检出=无功能”的误解，我们构建三类功能读出与两类桥接：
  - ISI（Inhibitor–Protease Stoichiometry）指数：ISI = log2(Σ 抑制子 / Σ 靶蛋白)；抑制子 core：SERPINA1/A3/SERPINB1/B6/B8（扩展：SLPI/ELAFIN/A2M）。
  - MNAR（Missingness‑as‑Signal）：以“是否检出 F/M”作因变量，Serpin_score 为自变量，(1|cohort) 随机项；M 作为阴性对照。
  - 底物足迹（NE/PR3 footprint）：半/非胰切口的 N‑端富集位点汇总成足迹指数。
  - THBS1 裂解指数：log2(neo‑N sum) − log2(tryptic N sum)（半胰重搜）。
  - RNA→蛋白桥接：Proteo‑ΔFM ~ Serpin_RNA * organ + 协变量；控制中性粒/纯度。

- 已执行与结果（严格证据链）
  - ISI（scripts/py/isi_model.py）
    - 输出：results/tables/isi_per_sample.tsv、results/tables/isi_models.tsv
    - 结论：PXD046330/PXD005719/PXD051579 中 F 目标和为 0，ISI 不可定义，符合“F 被抑制/难检出，而 M 持续”。
  - MNAR（scripts/py/mnar_detection_model.py）
    - 输出：results/tables/mnar_detection_table.tsv、results/tables/mnar_logit_results.tsv
    - 结论：detect_F=0、detect_M=1（几乎全样本）；Logit 不可估但模式明确。
  - 足迹（scripts/py/footprint_indices.py）
    - 输出：results/tables/footprints.tsv
    - 结论：指数整体为负；PXD005719 脑变体低于亲本（≈ −0.786 vs −0.633），与 THBS1 方向一致；PXD046330/PXD051579 足迹偏负。
  - 分层相关（results/tables/assoc_summary.tsv）
    - 足迹 vs Serpin：合并 ρ≈−0.563, q≈0.00182；非脑 ρ≈−0.560, q≈0.00130；脑层 n 小未稳。
    - THBS1 abundance vs Serpin：合并 ρ≈+0.583, q≈0.00218；非脑 ρ≈+0.611, q≈0.00130。
    - Proteo‑ΔFM vs Serpin：弱且不显著（ρ≈+0.225, p≈0.29），因 F 全缺信号弱化。
  - 混合效应（scripts/py/mixed_effects_a3.py）
    - 模型：endpoint ~ Serpin_score_core * organ + log_total_PSMs + (1 | dataset)
    - 输出：results/tables/mixed_effects_a3.tsv；诊断图：results/figures/mixed_effects_diagnostics.pdf
    - 结果：THBS1 abundance β_Serpin≈+0.028（未显著）；足迹模型脑层 n 小出现 singular（表内注明）。
  - 成长模型（非脑，scripts/py/mixed_effects_growth.py）
    - 输出：results/tables/mixed_effects_growth.tsv、results/tables/mixed_effects_messages.txt
    - 足迹 ~ Serpin：β≈+0.0027, p≈4e−20（指数为负，小正 β 表示抑制增强→足迹更负）。
    - THBS1 裂解 ~ Serpin：β≈−0.011（p≈0.41）；THBS1 log abundance ~ Serpin：β≈+0.0125（p≈0.30）。
  - SEM（探索性；scripts/py/prepare_sem_data.py + scripts/r/14_sem_model.R）
    - 输出：results/tables/sem_input.tsv、results/models/sem_results.json、results/figures/A3_sem_path.svg
    - 说明：n≈24，有不可逆信息矩阵与负方差警告，仅作概念示意。

- 跨层整合与图示
  - 样本对齐：results/tables/multiomics_alignment.tsv、results/tables/multiomics_pairwise_deltas.tsv（PXD046330↔GSE96860；PXD005719↔GSE12237）。
  - 特征汇总：results/tables/multilayer_sample_map.tsv、results/tables/multilayer_features.tsv（RNA ΔFM、Proteo‑ΔFM、Serpin、THBS1 裂解、ROS/MPO 占位）。
  - 机制总图：“双层三明治”（脑 vs 肺/非脑）：results/figures/F_double_sandwich.pdf。

- 执行脚本（蛋白/肽层）
  - Py：scripts/py/isi_model.py、scripts/py/mnar_detection_model.py、scripts/py/footprint_indices.py、scripts/py/thbs1_cleavage_index.py、scripts/py/proteomics_deltafm_score.py、scripts/py/serpin_scores.py、scripts/py/mixed_effects_a3.py、scripts/py/mixed_effects_growth.py、scripts/py/prepare_sem_data.py、scripts/py/module_scores.py、scripts/py/module_assoc.py
  - R：scripts/r/12_mixed_effects.R、scripts/r/14_sem_model.R

小结与下一步
- 结果要点：F 端在多个蛋白队列中呈“完全缺失/难检出”，而 M 端稳定；Serpin 分数与 NE/PR3 底物足迹负相关、与 THBS1 蛋白量正相关；混合效应的方向一致，SEM 概念模型与“Serpin 抑制→F 活性受抑→THBS1 裂解下降→ΔFM↓”相符。
- 临床意义：功能锚点（足迹/THBS1 裂解）在体液样本中具可行性，可与标记端联合用于风险识别；亦提示在脑微环境中优先针对“抑制–蛋白酶平衡”的调节策略。
- 下一步：系统化半胰重搜以量化 THBS1 裂解指数；实现肽水平零膨胀/层级模型；跨队列元分析“缺失即信号”；补充 SLPI/ELAFIN/A2M 的灵敏度分析。

## 4. 临床意义：早期液体生物标志物原型（l‑NFS）

- 目标：在 CSF/Plasma 中用“功能近端锚点（足迹/THBS1）+ 标记锚点（NET‑M/Serpin）”组合，预测脑转移高风险，指导预防用药与影像随访频率。

- 统一入口脚本：scripts/py/site_triage_liquid.py
  - --mode maxquant：读取 MaxQuant txt/ + 补充临床表，产出 run+patient 特征与单队列模型。
  - --mode csf：解析三层表头矩阵，仅保留 Area，构建 Serpin、NET‑F/M、THBS1/footprint 代理与 CSF 模型。
  - --mode plasma：跨队列融合（主 PXD032767；辅 PXD018301），加入 dataset_indicator，并输出稳健 Logit 与 (1|dataset) 混合效应。

- 队列与结果
  - CSF（MSV000089062，n=72，BrainMet=17）
    - 特征：data/processed/proteomics/MSV000089062/csf_patient_features.tsv
    - 模型：results/tables/msv000089062_csf_metrics.tsv（5 折 AUC≈0.70，整体 AUC≈0.66）
    - 诊断：results/tables/msv000089062_csf_{roc,pr,calibration,dca}.tsv；Spearman：..._spearman.tsv
    - 备注：源矩阵缺乏 THBS1/NE 裂解位点 → 稳健 Logit 奇异，日志中注明。
  - Plasma‑EV（PXD018301，n=512，BrainMet=12；三层表头自动抽取 Area）
    - 特征：data/processed/proteomics/PXD018301/csf_patient_features.tsv（含 sample_id/dataset）
    - 模型：results/tables/pxd018301_csf_metrics.tsv（5 折 AUC≈0.62，整体 AUC≈0.60）
    - 稳健 Logit（HC3 手工实现）：results/tables/pxd018301_csf_robust_logit.tsv（THBS1_log↓、footprint_index↓、THBS1_cleave_idx↑ 均 p<0.05）
  - 跨队列 Plasma（PXD032767×PXD018301，n=565，BrainMet=29）
    - 融合：results/tables/pxd032767_pxd018301_plasma_metrics.tsv（5 折 AUC≈0.79 / AP≈0.21；整体 AUC≈0.75）
    - 稳健 Logit：results/tables/..._plasma_robust_logit.tsv（THBS1_log↓、footprint_index↓、THBS1_cleave_idx↑ 显著；Proteo‑ΔFM 正向；dataset_indicator<0）
    - 混合效应（1|dataset）：results/tables/..._plasma_mixed_effects.tsv（方向与稳健 Logit 一致；示例：dataset_indicator β≈−0.65, p≈0.003）
  - 图：
    - ROC：results/figures/site_triage_roc.{pdf,png}
    - 效应量：results/figures/site_triage_effects.{pdf,png}

小结与下一步
- 结果要点：CSF 与 Plasma 原型在多数据集上达到可用的 AUC 水平（约 0.60–0.79 不等），稳健 Logit 与 (1|dataset) 混合效应一致指向“THBS1_log 下降、footprint 下降、THBS1 裂解上升”为风险特征；Proteo‑ΔFM 在跨队列融合中为正向贡献。
- 临床意义：为脑转移高风险早筛与随访频率分层提供可实施管线；可与影像学协同，降低不必要随访。
- 下一步：外部验证（独立血浆/CSF 队列）、阈值与决策曲线优化、与炎症/免疫共变项的特异性检验、前瞻性样本设计。

- 现场问题与修正（液体原型）
  - PXD032767 未公开非脑 RAW → 通过 Supplementary Table S1 映射患者分组；构建 run_manifest.tsv 并以 MaxQuant 直接计算足迹/THBS1。
    - 临床表：data/interim/proteomics/PXD032767/metadata/clinical_metadata.tsv
    - 运行表：data/interim/proteomics/PXD032767/metadata/run_manifest.tsv
    - 中间/特征：data/interim/proteomics/PXD032767/metadata/{footprint_maxquant.tsv, thbs1_cleavage.tsv}
  - MSV000089062 缺乏 THBS1/NE 裂解位点 → 近端锚点退化；保留为透明记录，主展示 ROC/效应量。
  - PXD018301 多层表头 → 修复为“仅保留 Area，以第 2 层样本名重命名列”。
  - statsmodels HC3：Logit 无现成交互 → 手动实现 HC3（杠杆 h 与残差构造 meat 矩阵）。
  - 混合效应：采用 BinomialBayesMixedGLM，报告 fe_mean/fe_sd；与稳健 Logit 方向一致。

## 5. 失败与调整（透明记录）

- 转录组：scRNA THBS1 相关因样本量与常量化不可估，且本次严格门控未检出稳定的中性粒细胞群（或门控为空/近 0），因此“低 ΔFM 中性粒亚群比例”无法稳定估计；更符合“中性粒已裂解成 NETs、功能端受抑”的解释。
- 蛋白/肽层：NET‑F 目标稀疏，ISI 与 Logit 在若干队列不可估；改以“缺失即信号 + 足迹/THBS1 裂解”双锚点承载功能信号。
- 混合效应：statsmodels==0.14.5 更新后仍有 singular 警告（低 n 合理），系数方向稳定；信息记录在 results/tables/mixed_effects_messages.txt。
- SEM：n≈24，模型识别度受限（不可逆信息矩阵、负方差），仅作机制示意。
- Metabolomics：Workbench（ST000745/ ST001104/ ST002921）仅 mwTab 元信息 + RAW，缺命名矩阵，红氧化比（GSH/GSSG）暂不可算；脚本先归档待恢复。

## 6. 数据、目录与环境（统一约定）

- 目录结构
  - 原始：data/raw/<dataset>/（保留 checksums.sha256）
  - 中间：data/interim/<modality>/<dataset>/（mzML、搜索工作区、临时 TSV）
  - 处理后：data/processed/<dataset>/（基因/蛋白矩阵、裂解指数、协变量）
  - 资源：resources/modules/、resources/manifests/、resources/annotation/（包括 gencode.v43.basic.annotation.gtf.gz）、resources/msfragger/
  - 结果：results/tables/、results/figures/

- 环境与工具
  - Proteomics 环境：env/proteomics.yml；激活：
    - export MAMBA_ROOT_PREFIX="$PWD/env/.mamba"
    - env/bin/micromamba env create -y -f env/proteomics.yml
    - eval "$(env/bin/micromamba shell hook -s bash)"
    - micromamba activate proteomics
    - export PATH="$PWD/env/bin:$PATH"；philosopher version（自检）
  - RNA bigWig：已集成 pybigwig 与 gffread；注释：resources/annotation/gencode/gencode.v43.basic.annotation.gtf.gz
  - R 包缓存（可选）：R_LIBS_USER=env/.R（msqrob2/lmerTest/plot 等）
  - 长作业：tmux/screen；日志到 logs/（例：logs/proteomics/pipeline_<timestamp>.log）
  - 脚本化环境：scripts/sh/setup_proteomics_env.sh（可一键拉起 proteomics 依赖）；RAW→mzML：scripts/sh/run_trfp_parallel.sh；端到端烟测：scripts/sh/msfragger_smoke.sh；GSE186344 补充抓取：scripts/sh/fetch_gse186344_suppl.sh。

- 数据获取（示例命令）
  - PRIDE：bash scripts/sh/fetch_pride_dataset.sh PXD011796 --categories RAW,RESULT --dest data/raw/PXD011796
  - PXD005719：scripts/sh/fetch_pride_dataset.sh PXD005719 --categories RAW,RESULT --dest data/raw/PXD005719
  - PXD046330：scripts/sh/fetch_pride_dataset.sh PXD046330 --categories RAW,RESULT --dest data/raw/PXD046330
  - PXD051579（暂停）：scripts/sh/fetch_pride_dataset.sh PXD051579 --categories RAW,RESULT --dest data/raw/PXD051579
  - GSE96860 bigWig：python scripts/py/process_gse_bigwig.py --bigwig-dir data/raw/transcriptomics/GSE96860 --gtf resources/annotation/gencode/gencode.v43.basic.annotation.gtf.gz --dataset GSE96860 --output-dir data/processed/GSE96860
  - GSE12237 芯片：python scripts/py/parse_series_matrix.py --series data/raw/transcriptomics/GSE12237/GSE12237_series_matrix.txt.gz --gpl data/raw/transcriptomics/GSE12237/GPL96.annot.gz --dataset GSE12237 --outdir data/processed/GSE12237
  - CSF 蛋白组（多来源）：scripts/sh/fetch_pride_dataset.sh <CSF_ACCESSION> --categories RESULT --match 'CSF' --dest data/raw/<CSF_ACCESSION>
  - CSF 代谢组：aria2c -c -x16 -s16 --max-tries=0 --retry-wait=30 -i resources/manifests/metabolomics_urls.txt -d data/raw/metabolomics

## 7. 可复现脚本入口与建议顺序

- Bulk（MVP）
  - Rscript scripts/r/02_covariates_estimate.R --dataset GSE184869 --config resources/config.yml
  - Rscript scripts/r/03_scores_deltafm.R --dataset GSE184869 --config resources/config.yml
  - Rscript scripts/r/04_bulk_stats_mvp.R --dataset GSE184869 --config resources/config.yml
  - Rscript scripts/r/05_pathway_scores.R --dataset GSE184869 --config resources/config.yml
  - Rscript scripts/r/06_fm_correlation_overview.R --config resources/config.yml

- scRNA（GSE186344）
  - python scripts/py/sc_preprocess_gse186344.py --raw GSE186344_RAW.tar --out data/processed/scRNA/gse186344.h5ad
  - python scripts/py/sc_neutrophil_subsets.py --h5ad data/processed/scRNA/gse186344.h5ad --net_f resources/modules/net_f_v1.tsv --net_m resources/modules/net_m_v2.tsv --outdir results/tables
  - python scripts/py/sc_tumor_thbs1_assoc.py --h5ad data/processed/scRNA/gse186344.h5ad --frac results/tables/gse186344_neutrophil_lowDFM_fraction.tsv --out results/tables/gse186344_subpop_thbs1_assoc.tsv

- 蛋白/肽层（A3）
  - python scripts/py/footprint_indices.py → results/tables/footprints.tsv
  - python scripts/py/thbs1_cleavage_index.py → results/tables/thbs1_cleave_idx.tsv
  - python scripts/py/serpin_scores.py → results/tables/serpin_scores.tsv
  - python scripts/py/proteomics_deltafm_score.py → data/processed/proteomics/<dataset>/proteo_deltafm.tsv
  - python scripts/py/isi_model.py → results/tables/isi_per_sample.tsv, isi_models.tsv
  - python scripts/py/mnar_detection_model.py → results/tables/mnar_detection_table.tsv, mnar_logit_results.tsv
  - python scripts/py/mixed_effects_a3.py → results/tables/mixed_effects_a3.tsv
  - python scripts/py/mixed_effects_growth.py → results/tables/mixed_effects_growth.tsv（诊断与消息：mixed_effects_messages.txt）
  - python scripts/py/prepare_sem_data.py；Rscript scripts/r/14_sem_model.R → sem_results.json, A3_sem_path.svg
  - python scripts/py/robust_assoc.py → 生成稳健相关/分层关联汇总（如 assoc_summary.tsv 等）
  - python scripts/py/neg_controls.py → results/tables/neg_controls.tsv, neg_controls_assoc.tsv（阴性/无关通路对照）
  - （可选流程）python scripts/py/proteomics_pipeline.py（端到端 orchestrate：RAW→mzML→搜索→汇总→指标；按需使用）

- 液体原型（CSF/Plasma）
  - CSF：python scripts/py/site_triage_liquid.py --mode csf --dataset MSV000089062 --input-xlsx data/interim/proteomics/MSV000089062/metadata/vdac161_suppl_supplementary_table_s1.xlsx
  - Plasma（单队列 MaxQuant）：python scripts/py/site_triage_liquid.py --mode maxquant --dataset PXD032767 --maxquant-dir data/interim/proteomics/PXD032767/txt --run-manifest data/interim/proteomics/PXD032767/metadata/run_manifest.tsv --clinical-xlsx data/interim/proteomics/PXD032767/metadata/vdac161_suppl_supplementary_table_s1.xlsx
  - Plasma（跨队列）：python scripts/py/site_triage_liquid.py --mode plasma --dataset PXD032767 --maxquant-dir ... --run-manifest ... --clinical-xlsx ... --aux-dataset PXD018301 --aux-features data/processed/proteomics/PXD018301/csf_patient_features.tsv
  - 图形产出：python scripts/py/plot_site_triage_figures.py → results/figures/site_triage_{roc,effects}.{pdf,png}
  - （资源构建）Rscript scripts/r/01_build_net_m_v2.R、Rscript scripts/r/export_selected_sets.R（模块/通路资源；如需重建）
  - （图像集）Rscript scripts/r/13_figures_a3.R（A3_forest/散点等发表图重绘）
  - （外展）python scripts/py/met500_fetch.py、python scripts/py/site_triage_met500.py（MET500 侧验证与资源获取）

## 8. 图与表（产出清单）

- Bulk 主结果
  - 图：results/figures/Figure2_GSE184869_Paired_DFM.pdf、results/figures/Figure3_GSE184869_THBS1_Residuals.pdf、results/figures/Figure4_FM_Corr_Overview.pdf、results/figures/Figure5_scRNA_LowDFM_Fraction.pdf
  - 附：results/figures/Figure5_scRNA_THBS1_Assoc.txt（样本层相关因常量化未估计的说明记录）
  - 表：results/tables/GSE184869_paired_tests.tsv、results/tables/GSE125989_paired_tests.tsv、results/tables/GSE184869_THBS1_residuals.tsv、results/tables/GSE184869_thbs1_partial_cor.tsv、results/tables/fm_correlation_overview.tsv、results/tables/gse186344_*.tsv

- 蛋白/肽层与整合
  - 图：results/figures/mixed_effects_diagnostics.pdf、results/figures/A3_sem_path.svg、results/figures/F_double_sandwich.pdf、results/figures/F_proteo_deltafm_boxplots.pdf、results/figures/F_proteo_deltafm_vs_transcriptome.pdf、results/figures/thbs1_spectrum_JIMT1_BR_R3.pdf、results/figures/A3_forest_assoc.pdf、results/figures/A3_scatter_serpin_vs_footprint.pdf、results/figures/A3_scatter_serpin_vs_ctrl_footprint.pdf、results/figures/A3_scatter_serpin_vs_deltaFM.pdf、results/figures/A3_scatter_serpin_vs_THBS1.pdf、results/figures/A3_scatter_serpin_vs_thbs1_cleave.pdf、results/figures/A3_scatter_serpin_vs_ecm_nc.pdf
  - 表：results/tables/footprints.tsv、results/tables/thbs1_cleave_idx.tsv、results/tables/serpin_scores.tsv、data/processed/proteomics/*/proteo_deltafm.tsv、results/tables/proteo_deltafm_summary.tsv、results/tables/isi_per_sample.tsv、results/tables/isi_models.tsv、results/tables/mnar_detection_table.tsv、results/tables/mnar_logit_results.tsv、results/tables/assoc_summary.tsv、results/tables/mixed_effects_a3.tsv、results/tables/mixed_effects_growth.tsv、results/tables/mixed_effects_messages.txt、results/tables/multiomics_alignment.tsv、results/tables/multiomics_pairwise_deltas.tsv、results/tables/multilayer_sample_map.tsv、results/tables/multilayer_features.tsv、results/tables/proteomics_transcript_sample_map.tsv、results/tables/proteomics_transcript_matched.tsv、results/tables/proteomics_transcript_unmatched.tsv、results/tables/multiomics_correlations.tsv、results/tables/proteo_serpin_assoc.tsv、results/tables/proteo_module_detection_summary.tsv

- 液体原型
  - 图：results/figures/site_triage_roc.{pdf,png}、results/figures/site_triage_effects.{pdf,png}、results/figures/A4_site_triage_roc_calib_dca.pdf
  - 表：results/tables/site_triage_metrics.tsv、results/tables/site_triage_predictions.tsv、results/tables/site_triage_coefficients.tsv、results/tables/site_triage_dca.tsv、results/tables/site_triage_dca_baseline.tsv、results/tables/site_triage_dca_serpin_extended.tsv、results/tables/site_triage_dca_null_inflammation.tsv、results/tables/site_triage_variant_metrics.tsv、results/tables/site_triage_variant_predictions.tsv、results/tables/site_triage_variant_coefficients.tsv、results/tables/msv000089062_csf_metrics.tsv、results/tables/msv000089062_csf_{roc,pr,calibration,dca}.tsv、results/tables/msv000089062_csf_{coefficients,scaler,spearman,predictions}.tsv、results/tables/pxd018301_csf_{metrics,roc,pr,calibration,coefficients,scaler,spearman,predictions}.tsv、results/tables/pxd018301_csf_{dca}.tsv、results/tables/pxd018301_csf_robust_logit.tsv、results/tables/pxd032767_liquid_{metrics,roc,pr,calibration,coefficients,scaler,predictions}.tsv、results/tables/pxd032767_liquid_{dca,spearman}.tsv、results/tables/pxd032767_liquid_{model,robust_logit}.tsv、results/tables/pzd032767_model_summary.tsv、results/tables/pxd032767_pxd018301_plasma_{metrics,roc,pr,calibration,predictions,coefficients,scaler,spearman}.tsv、results/tables/pxd032767_pxd018301_plasma_{robust_logit,mixed_effects}.tsv、results/tables/pxd032767_pxd018301_plasma_{dca}.tsv

- 其它图（示意/机制/总览）
  - 图：results/figures/Figure1_DFM_Schematic.pdf、results/figures/Figure6_Mechanism_Model.pdf

- 负/正对照与模块
  - 表：results/tables/neg_controls.tsv、results/tables/neg_controls_assoc.tsv、results/tables/module_scores.tsv、results/tables/module_assoc.tsv、results/tables/bulk_paired_tests.tsv、results/tables/GSE*_deltafm_scores.tsv

## 9. 下一步（两周可交付最小集）

- ISI（core/扩展）＋ 器官交互，导出森林图。
- Logit/零膨胀：F 检出 ~ Serpin（M 作阴性对照），汇总跨队列元分析。
- 底物足迹与 THBS1 裂解：指数与 Serpin/Proteo‑ΔFM 关联的稳健估计与可视化。
- 液体原型图：CSF/Plasma ROC 与效应量一页式整合；脑 vs 非脑 机制总图更新。

——

注：本 README 覆盖 `README_proteo_metab.md` 中的全部脚本、图表与失败记录。历史长文档仍保留为参考。
