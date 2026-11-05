# 数据分析方法（详细版，与最终图件一一对应）

> 本版以项目已有的分析记录与 `paper/Figure/` 的最终图件为准，按方法学写作的详细度展开：每一步写清输入、输出、函数、参数、阈值与图件归属。原始记录中未落盘的参数以 **`[待补]`** 标出，不臆造数值；未进入最终图的分析集中在最后一节。
---

## 1. 数据来源与样本

### 1.1 单细胞转录组与空间转录组

从 GEO 下载了 12 个 TNBC 单细胞样品（GSE176078）[1]，共 38,947 个细胞；从 GEO GSE210616 下载了 12 张 10x Genomics Visium 切片 [2]，来自 6 例 TNBC 患者（P9、P11–P15），每例 2 张，共 16,607 个 spot；切片编号为 092A、092B、093A–093D、094A–094D、095A、095B。

## 2. 单细胞数据预处理、聚类与注释

使用 Seurat（Seurat 4.4.0 / SeuratObject 4.1.4）[3] 对单细胞数据进行分析。过滤低质量细胞：`nFeature_RNA < 200` 或 `> 5,000`、线粒体比例 `> 20%` 的细胞被剔除；经 `NormalizeData`、`FindVariableFeatures`、`ScaleData`、`RunPCA`（50 个主成分）、harmony（1.2.3）批次校正与整合聚类 [4]，获得 **25 个 cluster**。

使用以下 marker 进行细胞类型注释：上皮（EPCAM、KRT8、KRT18、KRT19、KRT14、KRT5）、T 细胞（CD3D、CD3E、CD8A、CD4、FOXP3、CXCL13、GZMB、GZMK、TNFRSF4、LAG3、CCR7）、巨噬细胞（CD68、LYZ、CD163、CCL18、C1QA、SPP1、MKI67、FCN1）、成纤维细胞（COL1A1、DCN、LUM、MMP11、APOD、CFD）、B/浆细胞（MS4A1、CD79A、JCHAIN、MZB1）、内皮（PECAM1、VWF）与周细胞（RGS5、ACTA2）。随后对巨噬细胞、成纤维细胞与 T 细胞三个谱系内分别进行二次聚类，得到 7 个 TAM 状态（C1Q、CCL18、SPP1、FCN1、ISG15、MKI67、MMP9）、6 个 CAF 状态（MMP11、CD74、IGF1、CFD、APOD、STMN1）与 15 个 T 细胞状态。

## 3. 恶性细胞鉴定（inferCNV / SpatialInferCNV）

### 3.1 单细胞层面

以免疫细胞为参考，对上皮细胞做 inferCNV [5]：基因按染色体位置排序，以参考细胞平均表达为基线计算相对拷贝数；沿染色体滑窗平滑、HMM 分段并去噪，排除性染色体，得到 **CNV 负荷评分 `CNV_R_Score`** 与 **`Is_Malignant`（positive/negative）**。恶性 vs 正常上皮的 `CNV_R_Score` 用**双侧 Wilcoxon 秩和检验**（P < 2.2 × 10⁻¹⁶；Figure 1c）。

### 3.2 空间 inferCNV（SpatialInferCNV）与克隆判定

空间层面的 CNV 推断在**患者层面**完成：以上皮 spot 为观测组，将成纤维细胞、巨噬细胞、内皮细胞、T 细胞、B 细胞、周细胞与浆细胞七类非上皮 spot 合并为参考组；同一患者的两个切片合并后执行一次推断，使两片共享同一参考基线。使用 inferCNV 计算每个 spot 相对参考组的拷贝数，经滑窗平滑、**6 状态 HMM** 分段与去噪，得到基因 × spot 的 CNV 矩阵。

克隆判定采用"**树节点定义克隆**"策略，而非对 CNV 矩阵直接做无监督聚类：从 inferCNV 的去噪树状图重建有根树，在上皮细胞的 CNV 矩阵上做 PCA，对 k = 2–8 切树并以**轮廓系数**确定 k；每个簇取**最近公共祖先节点**的后代 spot 作为一个克隆，因而克隆在树上天然单系。随后按染色体臂汇总得到克隆级 consensus CNA，以参考组为外群、用**邻接法**构建克隆树。最终纳入 5 例患者、10 张切片、共 15 个克隆；克隆的 CNV score 定义为跨基因的平均绝对偏离（`mean |CNV − 1|`）。

### 3.3 三区判定

三区标签由 RCTD 反卷积得到的肿瘤上皮细胞比例（`tumor_prop`）按固定阈值划分：`tumor_prop > 0.50` 为肿瘤核心（Tumor Core），`0.20 ≤ tumor_prop ≤ 0.50` 为边界（Interface），`tumor_prop < 0.20` 为基质（Stroma）。`tumor_prop` 取RCTD反卷积中 9 大群水平的 Tumor Epithelial cells 比例；分区在 spot 水平进行，阈值全项目统一，不按切片分位数调整，也不做空间或邻域平滑。共得到 Core 9,515、Interface 6,687、Stroma 405 个 spot（合计 16,607）。

---

## 4. 空间反卷积（RCTD）

使用 RCTD（spacexr 2.2.1）[6]，以单细胞注释结果为参考对每个 Visium spot 进行两级反卷积：**主要细胞类型（9 类**：Tumor Epithelial、Normal Epithelial、T、Macrophages、Fibroblasts、Plasma、B、Pericytes、Endothelial**）** 用于定义肿瘤比例与三区；**细分亚群（35 类**：7 个 TAM 状态、6 个 CAF 状态、15 个 T 细胞状态等**）** 用于区分 TAM/CAF 的具体状态。运行采用 `RCTDmode = "full"`，为每个 spot 估计各类细胞的连续比例而非硬判定单一类型——三区划分依赖这一连续性——其余参数均为 spacexr 默认值。最终进入参考的基因为 2,083 个，输出为每个 spot 的两套比例矩阵（16,607 × 9 与 16,607 × 35）。

---

## 5. 空间域与三区一致性（BANKSY）

使用 BANKSY [7] 在每张切片内独立进行无监督空间域聚类，检验按肿瘤比例阈值划分的三区是否对应数据自身的空间域。输入为各切片的 SCT log-normalized 表达矩阵与 spot 阵列坐标，参数取 `k_geom = 18`、`λ = 0.2`、`npcs = 20`、Leiden `resolution = 0.8`；切片之间不合并、不做批次整合。以同一切片、同一表达矩阵的非空间 k-means 作为对照，簇数设为该切片 BANKSY 的域数。BANKSY 域与三区标签的一致性用调整兰德指数（ARI）逐切片及合并计算；另以域内平均肿瘤比例的域秩与 Spearman 相关检验域是否沿肿瘤比例梯度排列。

---

## 6. 空间可变基因（Squidpy Moran's I）

使用 Squidpy [8] 计算基因表达的空间自相关（Moran's I）并给出 FDR，用于筛选空间可变基因。候选基因为 SPP1、YBX1、C1QA、C1QB、C1QC 与 MMP11；另按三区汇总这 6 个基因的平均表达，比较其在不同区域的表达水平（Figure 2f、S2c）。

---

## 7. 区域组成与距核心梯度

**TAM 亚群的三区比例**：基于 RCTD 亚群比例矩阵，计算 7 个 TAM 亚群（FCN1、ISG15、C1Q、SPP1、CCL18、MKI67、MMP9）在"切片 × 区"内的平均比例并跨切片汇总；C1Q 与 SPP1 两个亚群以配对箱线图呈现，按切片配对比较它们相对 Tumor Core 的差异（配对 Wilcoxon 检验）。因部分切片不含 Stroma，Core 与 Interface 的样本数为 12、Stroma 为 8（Figure 2b、2c、S2a、S2b）。

**距肿瘤核心梯度**：对每个 spot 计算其到同切片内最近 Tumor Core spot 的欧氏距离，沿该距离考察 C1Q/(C1Q + SPP1) 比值及两个 TAM 亚群的占比。每张切片单独拟合平滑曲线，跨切片取均值并以 95% bootstrap 置信带表示离散度（Figure 2d）。距离单位为 Visium 阵列坐标，未换算为物理距离。

---

## 8. 邻域富集（Ro/e）

以"切片 × 区"为分析单元，检验指定细胞类型在空间邻域中的共现是否偏离随机期望。每个类型在该切片-区内取丰度前 25%（Q75）的 spot 定义为阳性；在每张切片内以 6 近邻构建邻域图，邻域期望由片内阳性标签随机置换得到。统计量为观测邻域阳性比例与置换均值的比值（Ro/e），每个单元置换 300 次，跨切片以 Stouffer z 合并并作 BH 校正。分析的类型为 C1QC TAM、SPP1 TAM、MMP11 CAF、CD8 与 CD4（另含 TAM 与 CAF 总群），结果以三区富集矩阵呈现，并计算 Stroma 与 Interface 之间的差值（Figure 3b–d）。该统计量反映的是阈值化阳性 spot 的观测/期望比，与第 9 节的共定位相关属不同度量。

---

## 9. 区域分层空间共定位（Squidpy）

在每个区内、按切片分别构建 7 近邻空间图（`squidpy.gr.spatial_neighbors`，`coord_type = "generic"`），计算 C1QC TAM、SPP1 TAM、MMP11 CAF、CD8 与 CD4 五类丰度与其邻域平均丰度之间的 Pearson 相关，构成对称的共定位矩阵。显著性由片内标签置换 1,000 次得到经验 P 值。结果呈现 Interface 区的共定位矩阵及其与 Stroma 的差值（Figure S3b、S3c）。由于 Stroma 仅 405 个 spot，部分切片-区不可估；各切片的空间图独立构建，组织大小不可跨片直接比较。

---

## 10. CD8 T 细胞功能评分

每张切片内，以 RCTD 反卷积得到的 CD8 T 细胞比例（`CD8_sum`）前 25% 的 spot 定义为 CD8 T 细胞 spot。在此基础上按基因集计算功能评分：细胞毒（GZMB、PRF1、NKG7、GNLY）与 TCR 活化（IL7R、CD69、NR4A1）；Interface 区的活化评分单独取 CD69、GZMB 与 PRF1 三个基因。评分方式为各基因在该片 CD8 spot 内逐基因 z-score 后取均值，再按"切片 × 生态位"汇总。生态位分为 Double+（C1QC 与 MMP11 均阳性）、C1QC only 与 MMP11 only，Other 不参与比较（Figure 3e、3f）。统计以切片为配对单位（n = 12），仅纳入三个生态位均含 ≥5 个 spot 的切片；三个生态位两两做配对 Wilcoxon 符号秩检验并作 BH 校正，图中标注的是未校正 P 值。

---

## 11. 细胞通讯

### 11.1 CellChat

使用 CellChat 2.1.2 [9] 与 CellChatDB.human，通讯概率以 `triMean` 估计（`nboot = 100`、`population.size = FALSE`，即按信号强度而非细胞数加权），并按 `min.cells = 10` 过滤。空间数据按"样本 × 分区"为单元独立建网，每个 spot 依反卷积比例归入单一细胞群（双阳性时归比例较高者），并启用空间距离限制（`distance.use = TRUE`，作用范围 500 个阵列单位），与后续空间流分析的邻域口径一致；单细胞数据与克隆水平分析按各细胞群直接建网。分析聚焦 C1QC TAM、MMP11 CAF、SPP1 TAM 与 CD8/CD4 T 细胞之间的通讯（Figure 5b、S5d、6a/6c）。

克隆水平的分析在两位患者各自的 **Interface 区**内独立进行，不合并、不跨患者比较。分组规则为：有克隆标签的 spot 优先归入 Clone1 或 Clone2；其余 spot 按 C1QC TAM 与 SPP1 TAM 比例的前 25%（片内 q75）定义为相应的 TAM 富集 spot，同一 spot 两者同时命中时归比例较高者，由此得到 Clone1、Clone2、C1Q TAM 与 SPP1 TAM 四组。该分析中过表达基因的识别使用标准 Wilcoxon 实现（环境无 presto 加速）。

需注意 CellChat 的 `net$weight` 与 `netP$prob` 会先将 p ≥ 0.05 的配受体对置零再求和，因此这些汇总值等于"效应量 × 显著条目数"，不等同于通讯强度；涉及克隆间比较时以逐配受体对的概率差与显著性分类为准。

### 11.2 NicheNet

使用 NicheNet v2 [10] 的先验（配体–靶基因矩阵与配体–受体网络）评估配体对受体细胞转录状态的调控潜力。受体基因集定义为：在受体阳性 spot 中表达与发送方丰度正相关（Pearson r ≥ 0.30、BH 校正 P < 0.05）、在受体阳性 spot 中富集（log 倍数变化 ≥ 0.25）且表达率 ≥ 10% 的基因，按相关性取前 200 个；配体则须在发送方 spot 中表达率 ≥ 10% 且平均表达 > 0.05。配体活性以 **AUPR** 衡量；由于 AUPR 的绝对值随受体基因集大小变化，仅用于配体之间的相对排序，不用于跨方向或跨区比较。CAF → C1Q TAM 与 C1Q TAM → CAF 两个方向分别独立分析（Figure 5a、S5c）。

由于空间 spot 为多细胞混合，受体基因集可能混入发送方自身的转录本（在界面区尤为明显），因此配体–靶基因面板在展示时剔除了胶原与间质结构基因，靶基因的解释以受体侧成分为准。

---

## 12. 肿瘤克隆与 TAM 生态位（Figure 4 / S4）

克隆的推断与判定方法见 **3.2**，本节不重复。最终纳入 5 例患者、10 张切片、共 15 个克隆（P9、P12、P13、P14、P15）；克隆编号带患者前缀，**不跨患者比较**。克隆的 CNV score 定义为跨基因的平均绝对偏离（`mean |CNV − 1|`），患者内共用纵轴（Figure 4b / S4b）。

### 12.1 克隆与 TAM 的共定位

在每张切片内分别评估克隆 spot 与 TAM 亚群 spot 的空间邻近程度，并列两种互为补充的统计量：

**半径邻域（观测边数）**：以片内 1-近邻距离中位数的 2 倍为邻域半径；TAM 阳性 spot 定义为该亚群比例 ≥ max(1%, 2 × 片内均值)。统计量为半径邻域内"克隆–TAM 阳性"spot 对的观测边数，以 TAM 阳性标签随机置换 1,000 次构建零分布，得到 z 值与经验 P 值。

**k 近邻 Ro/e**：以各 TAM 亚群比例前 25% 的 spot 为富集 spot，计算克隆 spot 及其 6 个最近邻落入富集 spot 的观测数与期望数之比（Ro/e）；以片内随机重贴克隆标签 1,000 次得到经验 P 值。

两种方法均在"切片 × 区 × 克隆 × TAM"层面统计，按 Tumor Core 与 Interface 分层（Stroma 因克隆 spot 过少不可估）；样本量不足的组合（区内 spot < 100、克隆 spot < 20 或 TAM 阳性 spot < 20）标为不可估。结果以 log2(Ro/e) 热图呈现，点的大小表示经验 P 值（Figure 4c、S4c–S4f）。P13 与 P15 各自归一、分别展示，不跨患者比较。

### 12.2 克隆间差异表达

差异表达在 P13 与 P15 两例患者中**分别进行，不合并、不跨患者比较**。输入为 spot 水平的 Visium 表达矩阵，原始 counts 经 LogNormalize（`scale.factor = 1e4`）标准化后，以 Seurat 的 `FindMarkers`（Wilcoxon 秩和检验）比较两个克隆，`ident.1 = Clone1`、`ident.2 = Clone2`，`logfc.threshold = 0` 以导出完整基因表而不预筛小效应基因。显著性以 BH 校正后的 `padj < 0.05` 判定，不额外设 log2FC 阈值；方向由 `avg_log2FC` 符号决定——`> 0` 表示 Clone1 相对 Clone2 高表达，`< 0` 反之（已用锚定基因校验方向）。除全片比较外，另在 Interface 区内对四张切片分别重复同一检验（Figure 4d / S4d）。

由于克隆与切片在部分患者中存在共线性，该检验反映的是两群 spot 的整体表达差异，不能单独归因于克隆内在的肿瘤细胞程序。

---

## 13. 生存分析（Figure S6）

- **队列与样本**：GSE21653（GEO，发现队列，DFS）与 GSE58812（GEO，验证队列，OS）；仅保留原发肿瘤与三阴性样本。GSE21653 以 ER、PR、ERBB2 均阴性定义三阴（87 例 → 83 例有 DFS 数据，25 个事件）；GSE58812 队列本身全为 TNBC（107 例，29 个 OS 事件）。
- **基因集**：从单细胞 Mac_sub / Fib_sub 差异表中按 `log2FC` 取 top20，排除技术/性别/线粒体基因（`MT-`、`RPL`、`RPS`、`HB[AB]`、`MALAT1`、`XIST`、`RPS4Y1`、`DDX3Y`）与 `IGKV`/`IGLV`；脚本 `05_survival/01_build_gene_lists.py`。三套基因集各 20 个基因（C1Q TAMs、SPP1 TAMs、MMP11 CAFs）。
- **打分**：队列内 ssGSEA（GSVA 1.50.0，`ssgsea.norm = TRUE`，α = 0.25）。基因名按 HGNC 新旧符号映射；缺失基因不插补。
- **分组**：中位数二分（`High = score > median`）；双高组 = C1Q 与 MMP11 均高于各自中位数。
- **统计**：KM + log-rank；Cox 模型估计 High vs Low 的 HR；BH 校正（族 = 队列 × 终点，4 个分析组）。
- **Figure S6**：GSE21653 DFS（发现）与 GSE58812 OS（验证），均为 ssGSEA + median，展示 C1Q、双高、SPP1、MMP11 四组的方向一致性；参考输出见 `reference_outputs/survival_v2/`。
- **说明**：代码中保留了 TCGA-BRCA 与 METABRIC 的探索性分支（maxstat 切点、z-score 敏感性、cox.zph），但本文 Figure S6 只展示上述两个队列。

---

## 14. 软件与版本

| 分析 | 软件/包 | 版本 |
|---|---|---|
| 单细胞 | Seurat / SeuratObject | 4.4.0 / 4.1.4 |
| 批次整合 | harmony | 1.2.3 |
| 空间反卷积 | spacexr（RCTD） | 2.2.1 |
| 恶性判定 / 克隆 | inferCNV | 1.22.0 |
| 空间域 | BANKSY | 0.99.13 |
| 空间自相关（Moran's I） | Squidpy | 1.2.2 |
| 区域共定位 | Squidpy | 1.2.2 |
| 邻域富集 | FNN + 自写置换 | FNN 1.1.4.1 |
| 细胞通讯 | CellChat | 2.1.2 |
| 配体–靶基因 | nichenetr | 2.2.1 |
| 差异表达 | Seurat（FindMarkers，Wilcoxon） | 4.4.0 |
| 生存分析 | survival / survminer / maxstat | 3.8.3 / 0.5.0 / 0.7.26 |
| 基因集打分 | GSVA | 1.50.0 |
| 绘图 | ggplot2 / ComplexHeatmap / circlize / patchwork / ggsci | 3.5.2 / 2.18.0 / 0.4.16 / — / — |

---

## 参考文献

1. Wu SZ, Al-Eryani G, Roden DL, et al. A single-cell and spatially resolved atlas of human breast cancers. *Nat Genet*. 2021;53(9):1334–1347. doi:10.1038/s41588-021-00911-1. GEO: GSE176078.
2. Bassiouni R, Idowu MO, Gibbs LD, et al. Spatial Transcriptomic Analysis of a Diverse Patient Cohort Reveals a Conserved Architecture in Triple-Negative Breast Cancer. *Cancer Res*. 2023;83(1):34–48. doi:10.1158/0008-5472.CAN-22-2682. PMID: 36283023; PMCID: PMC9812886. GEO: GSE210616.
3. Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data. *Cell*. 2021;184(13):3573–3587.e29. doi:10.1016/j.cell.2021.04.048.
4. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods*. 2019;16(12):1289–1296. doi:10.1038/s41592-019-0619-0.
5. Puram SV, Tirosh I, Parikh AS, et al. Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer. *Cell*. 2017;171(7):1611–1624.e24. doi:10.1016/j.cell.2017.10.044.
6. Cable DM, Murray E, Zou LS, et al. Robust decomposition of cell type mixtures in spatial transcriptomics. *Nat Biotechnol*. 2022;40(4):517–526. doi:10.1038/s41587-021-00830-w.
7. Singhal V, Chou N, Lee J, et al. BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis. *Nat Genet*. 2024;56(3):431–441. doi:10.1038/s41588-024-01664-3.
8. Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. *Nat Methods*. 2022;19(2):171–178. doi:10.1038/s41592-021-01358-2.
9. Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell-cell communication using CellChat. *Nat Commun*. 2021;12:1088. doi:10.1038/s41467-021-21246-9.
10. Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. *Nat Methods*. 2020;17(2):159–162. doi:10.1038/s41592-019-0667-5.
