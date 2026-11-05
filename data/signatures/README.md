# 从 Mac_sub / Fib_sub 差异列表提取 top20 基因集

## 输入

| 文件 | 行数 | 列 |
|---|---:|---|
| `Mac_sub_差异列表.tsv` | 204,560 | `genes, log2FC, cluster, pct, avgExpr` |
| `Fib_sub_差异列表.tsv` | 113,862 | 同上 |

- 每个 cluster 覆盖全部检测基因（Mac 20,456；Fib 18,977），**无 p 值/padj 列**，也未做显著性过滤。
- 目标细胞群：`C1Q TAMs`、`SPP1 TAMs`（Mac 表）；`MMP11 CAFs`（Fib 表）。

## 筛选规则

1. 取目标 cluster 的全部基因；
2. 排除技术/性别/线粒体基因：`^(MT-|RPL|RPS|HB[AB]|MALAT1$|XIST$|RPS4Y1$|DDX3Y$)`
3. 排除免疫球蛋白可变区基因：`^(IGKV|IGLV)`
4. 按 `log2FC` 降序（同值按 `pct` 降序）取 **top 20**（排除后重新取，保持每组 20 个）。

各群排除情况：

| 细胞群 | 检测基因 | 技术/性别/线粒体 | IGKV/IGLV | 剩余 |
|---|---:|---:|---:|---:|
| C1Q TAMs | 20,456 | 114 | 59 | 20,283 |
| SPP1 TAMs | 20,456 | 114 | 59 | 20,283 |
| MMP11 CAFs | 18,977 | 114 | 50 | 18,813 |

### 因 IGKV/IGLV 规则被移出 top20 的基因（`removed_IG_genes.tsv`）

| 细胞群 | 基因 | log2FC | pct |
|---|---|---:|---:|
| C1Q TAMs | IGLV4-69 | 1.2721 | 0.023 |
| C1Q TAMs | IGLV1-47 | 1.1715 | 0.125 |
| C1Q TAMs | IGKV3-20 | 1.1527 | 0.317 |

C1Q TAMs 因此顺延补入：`TCEB2`、`CTSC`、`LIPA`。SPP1 TAMs 与 MMP11 CAFs 无 IGKV/IGLV 进入原 top20，列表不变。

## 输出

| 文件 | 内容 |
|---|---|
| `gene_list_top20_long.tsv` | 长表：cell_type, rank, genes, log2FC, pct, avgExpr, source_file |
| `gene_list_top20_wide.csv` | 宽表：rank × 三个细胞群 |
| `gene_list_unique.txt` | 三群 top20 的并集（去重，59 个基因） |
| `gene_list_stats.csv` | 每群检测/排除/剩余基因数 |
| `removed_IG_genes.tsv` | 被 IGKV/IGLV 规则移除且原本会进 top20 的基因 |

## 注意事项

1. **无显著性列**：只能按 log2FC 排序，无法要求 `padj<0.05`。
2. **低阳性率基因仍存在**：SPP1 的 `FABP4`（pct=0.107）、`MMP1`（0.101）、`CCL7`（0.221）、`MARCO`（0.237），MMP11 的 `CST1`（0.120），C1Q 的 `PTGDS`（0.201）。若作为 bulk signature，建议加 `pct ≥ 0.25` 过滤后再取 top20。
3. **`MT1*`/`MT2A` 不是线粒体基因**：它们是金属硫蛋白，正则 `^MT-` 不会误排除；是否保留需按生物学判断。
4. **与旧 signature 不可比**：本列表的 DE 引擎/标准化与 `scripts/01_build_signatures.R`（presto wilcoxauc，log-normalized）不同，top 基因差异很大（C1Q 由 APOE/C1QB 变为 PTGDS/GPX1/SEPP1）。用于重做生存分析时应作为一套新 signature 重新打分，不能与旧结果混用。
5. **跨群重叠**：C1Q × SPP1 = 0；C1Q × MMP11 = 0；SPP1 × MMP11 = 1（`MMP1`）。
