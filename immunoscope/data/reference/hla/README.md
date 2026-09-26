# HLA 本地参考库

该目录用于承载 ImmunoScope 的本地 HLA 参考序列资源，供后续：

- HLA locus 预测
- top candidate allele 排序
- HLA 结构/序列注释

统一使用。

## 目录结构

- `raw/class_i/`
  - 原始输入库，按 locus 分文件保存全长蛋白序列
  - 当前已有：`A_prot.fasta`、`B_prot.fasta`、`C_prot.fasta`、`E_prot.fasta`

- `derived/`
  - 后续由预处理脚本生成的标准化库
  - 建议产物：
    - `class_i_extracellular.fasta`
    - `class_i_metadata.csv`
    - `class_i_index.json`

## 设计约定

1. `raw/` 保存原始来源，不直接作为在线注释的最终检索库
2. `derived/` 保存可直接用于比对和注释的标准化库
3. 第一阶段优先支持：
   - `HLA-A`
   - `HLA-B`
   - `HLA-C`
   - `HLA-E`
4. 第一阶段的注释目标优先级：
   - `best locus`
   - `top candidate alleles`
   - `identity / coverage / confidence`

## 后续实施建议

先从 `raw/class_i/*.fasta` 生成一个只保留胞外成熟区的派生库，再在代码里用 `BioPython PairwiseAligner` 做本地比对。

