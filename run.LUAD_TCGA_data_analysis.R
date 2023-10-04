### TCGA-LUAD数据分析 ###
#library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(cowplot)
library(dplyr)
library(maftools)
library(DESeq2)
library(clusterProfiler)

# set the work directory
setwd("/data/LUAD_PYGL/")
###################################################################################################
##### step1. data pre-processing #####
###################################################################################################
### 1.读取表达谱数据
TCGA_LUAD_Exp <- readRDS("./TCGA/TCGA-LUAD_expres.rds") ##注意文件路径要正确
TCGA_LUAD_Exp
#class: RangedSummarizedExperiment
#dim: 60660 600
#metadata(1): data_release
#assays(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand
#rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ...
#  ENSG00000288674.1 ENSG00000288675.1
#rowData names(10): source type ... hgnc_id havana_gene
#colnames(600): TCGA-78-7156-01A-11R-2039-07
#  TCGA-44-6774-01A-21R-1858-07 ... TCGA-50-8459-01A-11R-2326-07
#  TCGA-73-4675-01A-01R-1206-07
#colData names(87): barcode patient ... paper_Ploidy.ABSOLUTE.calls
#  paper_Purity.ABSOLUTE.calls

dim(TCGA_LUAD_Exp) # 60660个基因，600个样本
#[1] 60660   600

# 查看样本信息
head(colData(TCGA_LUAD_Exp))
#DataFrame with 6 rows and 87 columns
#                                            barcode      patient
#                                        <character>  <character>
#TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156-01A-11R.. TCGA-78-7156
#TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774-01A-21R.. TCGA-44-6774
#TCGA-91-6829-11A-01R-1858-07 TCGA-91-6829-11A-01R.. TCGA-91-6829

# 查看基因信息
head(rowData(TCGA_LUAD_Exp))
#DataFrame with 6 rows and 10 columns
#                     source     type     score     phase            gene_id
#                   <factor> <factor> <numeric> <integer>        <character>
#ENSG00000000003.15   HAVANA     gene        NA        NA ENSG00000000003.15
#ENSG00000000005.6    HAVANA     gene        NA        NA  ENSG00000000005.6
#ENSG00000000419.13   HAVANA     gene        NA        NA ENSG00000000419.13

# 使用assay函数提取表达数据
assays(TCGA_LUAD_Exp) #里面包含了6种表达矩阵
#List of length 6
#names(6): unstranded stranded_first ... fpkm_unstrand fpkm_uq_unstrand
names(assays(TCGA_LUAD_Exp))
#[1] "unstranded"       "stranded_first"   "stranded_second"  "tpm_unstrand"
#[5] "fpkm_unstrand"    "fpkm_uq_unstrand"

# 1为原始count数据("unstranded")，4为TPM表达值("tpm_unstrand")，5为FPKM表达值("fpkm_unstrand")
TCGA_LUAD_Exp_count <- SummarizedExperiment::assay(TCGA_LUAD_Exp,1)
TCGA_LUAD_Exp_fpkm <- SummarizedExperiment::assay(TCGA_LUAD_Exp,5)
dim(TCGA_LUAD_Exp_count)
#[1] 60660   600
dim(TCGA_LUAD_Exp_fpkm)
#[1] 60660   600
TCGA_LUAD_Exp_count[1:5,1:5]
#                   TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07
#ENSG00000000003.15                         2296                         1403
#ENSG00000000005.6                             2                            1
#ENSG00000000419.13                          951                          662
TCGA_LUAD_Exp_fpkm[1:5,1:5]
#                   TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07
#ENSG00000000003.15                      14.0995                      12.7765
#ENSG00000000005.6                        0.0377                       0.0280
#ENSG00000000419.13                      21.9471                      22.6558

# 提取gene信息
gene_id=data.frame(id = rowData(TCGA_LUAD_Exp)@listData[["gene_id"]], 
                   gene_name = rowData(TCGA_LUAD_Exp)@listData[["gene_name"]],
                   gene_type = rowData(TCGA_LUAD_Exp)@listData[["gene_type"]])
head(gene_id)
#                  id gene_name      gene_type
#1 ENSG00000000003.15    TSPAN6 protein_coding
#2  ENSG00000000005.6      TNMD protein_coding
#3 ENSG00000000419.13      DPM1 protein_coding

# 合并基因id和表达矩阵
counts <- cbind(gene_id, TCGA_LUAD_Exp_count)
fpkm <- cbind(gene_id, TCGA_LUAD_Exp_fpkm)
counts[1:5,1:5]
#                                   id gene_name      gene_type
#ENSG00000000003.15 ENSG00000000003.15    TSPAN6 protein_coding
#ENSG00000000005.6   ENSG00000000005.6      TNMD protein_coding
#ENSG00000000419.13 ENSG00000000419.13      DPM1 protein_coding
#ENSG00000000457.14 ENSG00000000457.14     SCYL3 protein_coding
#ENSG00000000460.17 ENSG00000000460.17  C1orf112 protein_coding
dim(counts)
#[1] 60660   603

# 提取蛋白编码基因(protein_coding gens)
counts_pcg <- counts[counts$gene_type == "protein_coding",-c(1,3)]
counts_pcg[1:5,1:5]
#                   gene_name TCGA-78-7156-01A-11R-2039-07
#ENSG00000000003.15    TSPAN6                         2296
#ENSG00000000005.6       TNMD                            2
#ENSG00000000419.13      DPM1                          951
#ENSG00000000457.14     SCYL3                         1196
#ENSG00000000460.17  C1orf112                          141
dim(counts_pcg)
#[1] 19962   601
fpkm_pcg <- fpkm[fpkm$gene_type == "protein_coding",-c(1,3)]
dim(fpkm_pcg)
#[1] 19962   601

# save RDS data
saveRDS(counts_pcg,"TCGA_LUAD_counts_pcg.rds")
saveRDS(fpkm_pcg,"TCGA_LUAD_fpkm_pcg.rds")

# 区分肿瘤样本和正常样本
table(colData(TCGA_LUAD_Exp)$shortLetterCode)
#NT  TP  TR 
#59 539   2 
#NT:Solid Tissue Normal
#TP:Primary Solid Tumor
#TR:Recurrent Solid Tumor

#2.提取临床信息
#可以直接从表达数据中提取metadata临床信息
TCGA_LUAD_clinData<-SummarizedExperiment::colData(TCGA_LUAD_Exp)
head(TCGA_LUAD_clinData)
#DataFrame with 6 rows and 87 columns
#                                            barcode      patient
#                                        <character>  <character>
#TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156-01A-11R.. TCGA-78-7156
#TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774-01A-21R.. TCGA-44-6774
#TCGA-91-6829-11A-01R-1858-07 TCGA-91-6829-11A-01R.. TCGA-91-6829

# 选择所需的meta信息
meta_col <- c("barcode","patient","sample","shortLetterCode","definition","days_to_death","days_to_last_follow_up","ajcc_pathologic_stage",
              "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m","years_smoked","gender","vital_status","age_at_index","progression_or_recurrence")
meta <- TCGA_LUAD_clinData[,meta_col]
head(meta)
# DataFrame with 6 rows and 16 columns
# barcode      patient           sample shortLetterCode          definition
# <character>  <character>      <character>     <character>         <character>
# TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156-01A-11R.. TCGA-78-7156 TCGA-78-7156-01A              TP Primary solid Tumor
# TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774-01A-21R.. TCGA-44-6774 TCGA-44-6774-01A              TP Primary solid Tumor
# TCGA-91-6829-11A-01R-1858-07 TCGA-91-6829-11A-01R.. TCGA-91-6829 TCGA-91-6829-11A              NT Solid Tissue Normal
# TCGA-69-A59K-01A-11R-A262-07 TCGA-69-A59K-01A-11R.. TCGA-69-A59K TCGA-69-A59K-01A              TP Primary solid Tumor

#肿瘤分组
meta$group <- ifelse(meta$shortLetterCode == "NT","Normal","Tumor")
#生存状态
meta$status <- ifelse(meta$vital_status == "Dead",1,0)
#是否吸烟
meta$smoking <- ifelse(meta$years_smoked >=1,"Yes","No")
#年龄分组
meta$age <- ifelse(meta$age_at_index > 70, ">70","<=70")
#TNM分型
meta$TNM <- paste0(meta$ajcc_pathologic_t,meta$ajcc_pathologic_n,meta$ajcc_pathologic_m)
#生存时间
meta$OS <- ifelse(meta$vital_status=='Alive',meta$days_to_last_follow_up,meta$days_to_death)
meta$OS_month <- round(meta$OS/30,2) #以month为单位，保留两位小数
meta$OS_year <- round(meta$OS/365,2) #以year为单位，保留两位小数
meta$stage <- NA
table(meta$ajcc_pathologic_stage)
# Stage I   Stage IA   Stage IB   Stage II  Stage IIA  Stage IIB Stage IIIA Stage IIIB   Stage IV 
#       5        152        170          1         55         83         85         12         28 
meta[meta$ajcc_pathologic_stage %in% c("Stage I","Stage IA","Stage IB"),"stage"] <- "Stage I"
meta[meta$ajcc_pathologic_stage %in% c("Stage II","Stage IIA","Stage IIB"),"stage"] <- "Stage II"
meta[meta$ajcc_pathologic_stage %in% c("Stage IIIA","Stage IIIB"),"stage"] <- "Stage III"
meta[meta$ajcc_pathologic_stage %in% c("Stage IV"),"stage"] <- "Stage IV"
table(meta$stage)
# Stage I  Stage II Stage III  Stage IV 
#     327       139        97        28
saveRDS(meta,"TCGA_LUAD_meta_data.rds")

#3.读取SNV突变数据
#下载的SNV_maf文件没有临床信息需要自己整理一下才能使用maftools
library(maftools)

# using maftools for data summary 
TCGA_LUAD_maf <- readRDS("./TCGA/TCGA-LUAD_SNV_Masked_Somatic_Mutation.rds")
TCGA_LUAD_maf[1:5,1:10]
#  X1 Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position
#1  1       ACOT7          11332     BI     GRCh38       chr1        6264665
#2  1     TNFRSF9           3604     BI     GRCh38       chr1        7937723
#3  1     DNAJC16          23341     BI     GRCh38       chr1       15566098
#4  1      AKR7A3          22977     BI     GRCh38       chr1       19282839
#5  1      EIF4G3           8672     BI     GRCh38       chr1       20980404
#  End_Position Strand Variant_Classification
#1      6264665      +      Missense_Mutation
#2      7937723      +      Missense_Mutation
#3     15566098      +                 Silent
#4     19282839      +      Missense_Mutation
#5     20980404      +      Nonsense_Mutation
maftools.input <- read.maf(TCGA_LUAD_maf)

#计算tmb值
tmb_table_wt_log = tmb(maf = maftools.input)
#查看tmb值
head(tmb_table_wt_log)
#            Tumor_Sample_Barcode total total_perMB total_perMB_log
# 1: TCGA-49-AARR-01A-11D-A410-08     1        0.02        -1.69897
# 2: TCGA-55-8513-01A-11D-2393-08     1        0.02        -1.69897
# 3: TCGA-17-Z019-01A-01W-0746-08     2        0.04        -1.39794
# 4: TCGA-17-Z054-01A-01W-0747-08     4        0.08        -1.09691
# 5: TCGA-L4-A4E6-01A-11D-A24D-08     4        0.08        -1.09691
# 6: TCGA-44-6148-01A-11D-1753-08     5        0.10        -1.00000

# Check summary
plotmafSummary(maf = maftools.input, 
               rmOutlier = TRUE, 
               addStat = 'median', 
               dashboard = TRUE)

library(ComplexHeatmap)
oncoplot(maf = maftools.input,top = 20, clinicalFeatures = "stage")
oncostrip(maf=maftools.input,genes=c("TP53","KRAS","EGFR","ALK"))

###################################################################################################
##### step2. Tumor vs. Normal DE analysis #####
###################################################################################################
library(DESeq2)

# load raw count data
########## tumor-vs-normal DESeq2差异表达分析 ##########
expDataCounts <- readRDS("TCGA_LUAD_counts_pcg.rds")
expDataCounts[1:5,1:5]
#                    gene_name TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07 TCGA-91-6829-11A-01R-1858-07
# ENSG00000000003.15    TSPAN6                         2296                         1403                          780
# ENSG00000000005.6       TNMD                            2                            1                            2
# ENSG00000000419.13      DPM1                          951                          662                          789
# ENSG00000000457.14     SCYL3                         1196                          256                          406
dim(expDataCounts)
#[1] 19962   601

table(duplicated(expDataCounts$gene_name))
#FALSE  TRUE
#19938    24

#remove duplicated genes
expDataCounts <- expDataCounts[!duplicated(expDataCounts$gene_name),]
rownames(expDataCounts) <- expDataCounts$gene_name
expDataCounts <- expDataCounts[,-1]
dim(expDataCounts)
#[1] 19938   600
expDataCounts[1:5,1:5]
#         TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07 TCGA-91-6829-11A-01R-1858-07 TCGA-69-A59K-01A-11R-A262-07
# TSPAN6                           2296                         1403                          780                         1338
# TNMD                                2                            1                            2                            0
# DPM1                              951                          662                          789                         2217
# SCYL3                            1196                          256                          406                         1101
# C1orf112                          141                          158                          146                          620

# load meta data
meta <- readRDS("TCGA_LUAD_meta_data.rds")
head(meta)
#DataFrame with 6 rows and 22 columns
#                                             barcode      patient           sample shortLetterCode          definition days_to_death
#                                         <character>  <character>      <character>     <character>         <character>     <integer>
# TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156-01A-11R.. TCGA-78-7156 TCGA-78-7156-01A              TP Primary solid Tumor           976
# TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774-01A-21R.. TCGA-44-6774 TCGA-44-6774-01A              TP Primary solid Tumor            NA
# TCGA-91-6829-11A-01R-1858-07 TCGA-91-6829-11A-01R.. TCGA-91-6829 TCGA-91-6829-11A              NT Solid Tissue Normal          1258
table(meta$group,meta$definition)
#        Primary solid Tumor Recurrent Solid Tumor Solid Tissue Normal
# Normal                   0                     0                  59
# Tumor                  539                     2                   0

expDataCounts <- expDataCounts[,rownames(meta)]
dim(expDataCounts)
#[1] 19938   600

colData <- as.data.frame(meta[,c("group","gender","stage")])
head(colData)
#                                    group      gender       stage
# TCGA-78-7156-01A-11R-2039-07       Tumor        male    Stage IV
# TCGA-44-6774-01A-21R-1858-07       Tumor      female   Stage III
# TCGA-91-6829-11A-01R-1858-07      Normal        male     Stage I

# 构建dds对象
dds <- DESeqDataSetFromMatrix(countData=round(expDataCounts, digits = 0),
                              colData=colData, 
                              design=~group)
dds
# class: DESeqDataSet 
#dim: 19938 600
#metadata(1): version
#assays(1): counts
#rownames(19938): TSPAN6 TNMD ... AL391628.1 AP006621.6
#rowData names(0):
#colnames(598): TCGA-78-7156-01A-11R-2039-07
#  TCGA-44-6774-01A-21R-1858-07 ... TCGA-50-8459-01A-11R-2326-07
#  TCGA-73-4675-01A-01R-1206-07
#colData names(3): group gender stage

# 数据过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 3  #过滤低表达基因，至少有3个样品都满足10个以上的reads数
dds <- dds[keep,]

# vst数据标准化
vsd <- vst(dds, blind=FALSE)
exprSet=assay(vsd)
exprSet[1:5,1:5]
#         TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07
#TSPAN6                      11.544649                    11.487493
#TNMD                         6.225937                     6.192885
#DPM1                        10.351699                    10.469047
#SCYL3                       10.655603                     9.271773
write.csv(exprSet,"TCGA-LUAD_count_vst_normalized.csv",quote=F)

# PCA降维
#plot by ggplot2
pcaData <- plotPCA(vsd, intgroup = c("group"), returnData = TRUE)
head(pcaData)
#                                   PC1        PC2 group group.1
#TCGA-78-7156-01A-11R-2039-07  11.97376  20.571341    TP      TP
#TCGA-44-6774-01A-21R-1858-07 -17.44269  -3.530650    TP      TP
#TCGA-91-6829-11A-01R-1858-07  32.78290   2.046207    NT      NT
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  #  coord_fixed() +
  ggtitle("PCA with VST data") + theme_bw()

# DESeq2差异分析
dds <- DESeq(dds)
resultsNames(dds)
#[1] "Intercept"      "group_Tumor_vs_Normal"
res <- results(dds,contrast=c("group","Tumor","Normal"))

# MA plot without shrink
plotMA(res, ylim = c(-5, 5))

# 结果保存
##筛选差异表达基因
#首先对表格排个序，先按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1 <- as.data.frame(res1)
head(res1)
#          baseMean log2FoldChange      lfcSE      stat        pvalue          padj  sig
# FAM83A  6243.1199       6.787874 0.23491570  28.89494 1.382284e-183 2.544784e-179   Up
# PYCR1   3599.8456       3.688964 0.13168909  28.01268 1.138578e-172 1.048061e-168   Up
# OTUD1   1102.3923      -2.109789 0.08083939 -26.09852 3.789375e-150 2.325413e-146 Down
# EPAS1  27000.2314      -2.714742 0.11040856 -24.58815 1.691479e-133 7.785033e-130 Down

#log2FC≥1 & padj<0.05 标识Up，代表显著上调的基因
#log2FC≤-1 & padj<0.05 标识Down，代表显著下调的基因
#其余标识 None，代表非差异的基因
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),'sig'] <- 'Up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05),'sig'] <- 'Down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.05),'sig'] <- 'None'
write.csv(res1, file = 'TCGA-LUAD_tumor_vs_normal.DESeq2.all.csv', quote = FALSE)

#输出差异基因总表
res1_sig <- subset(res1, sig %in% c('Up', 'Down'))
write.csv(res1_sig, file = 'TCGA-LUAD_tumor_vs_normal.DESeq2.DEs.csv', quote = FALSE)

###################################################################################################
##### step3. functional enrichment analysis #####
###################################################################################################
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

res1 <- read.csv(file = 'TCGA-LUAD_tumor_vs_normal.DESeq2.all.csv',row.names = 1)
head(res1)
#          baseMean log2FoldChange      lfcSE      stat        pvalue          padj  sig
# FAM83A  6243.1199       6.787874 0.23491570  28.89494 1.382284e-183 2.544784e-179   Up
# PYCR1   3599.8456       3.688964 0.13168909  28.01268 1.138578e-172 1.048061e-168   Up
# OTUD1   1102.3923      -2.109789 0.08083939 -26.09852 3.789375e-150 2.325413e-146 Down
# EPAS1  27000.2314      -2.714742 0.11040856 -24.58815 1.691479e-133 7.785033e-130 Down

### GSEA基因集富集分析
gene_df <- data.frame(gene=rownames(res1),logFC=res1$log2FoldChange)
head(gene_df)
#     gene     logFC
# 1 FAM83A  6.787874
# 2  PYCR1  3.688964
# 3  OTUD1 -2.109789
gene_id=bitr(gene_df$gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
head(gene_id)
#   SYMBOL ENTREZID
# 1 FAM83A    84985
# 2  PYCR1     5831
# 3  OTUD1   220213

genes <- merge(gene_id,gene_df,by.x="SYMBOL",by.y="gene")
head(genes)
#    SYMBOL ENTREZID       logFC
# 1    A1BG        1  0.36168132
# 2    A1CF    29974  3.64453703
# 3     A2M        2 -1.87464592

geneList <- genes$logFC
names(geneList) <- genes$SYMBOL
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)
#     TFF2   MAGEA3   MAGEA6     PDX1  MAGEA12      FGB 
# 9.600072 9.177694 9.157925 8.798172 8.728547 8.638270

### KEGG代谢通路富集
kegg <- read.gmt("c2.cp.kegg.v2023.1.Hs.symbols.gmt")
### HALLMARK代谢通路富集
hallmark <- read.gmt("h.all.v2023.1.Hs.symbols.gmt")

### 合并KEGG和HALLMARK的代谢通路
pathway_all <- rbind(kegg,hallmark)
head(pathway_all)

gsea_pathway <- GSEA(geneList,TERM2GENE = pathway_all, pvalueCutoff = 0.05, pAdjustMethod = "BH")
dim(gsea_pathway)
dotplot(gsea_pathway)

gsea_res <- as.data.frame(gsea_pathway)
head(gsea_res)
# 挑选感兴趣得通路展示
pathway <- c("KEGG_DNA_REPLICATION","KEGG_CELL_CYCLE","HALLMARK_IL2_STAT5_SIGNALING",
             "KEGG_STARCH_AND_SUCROSE_METABOLISM","KEGG_STEROID_HORMONE_BIOSYNTHESIS",
             "KEGG_DRUG_METABOLISM_OTHER_ENZYMES","KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY","KEGG_LYSOSOME",
             "KEGG_CHEMOKINE_SIGNALING_PATHWAY","KEGG_CELL_ADHESION_MOLECULES_CAMS",
             "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_G2M_CHECKPOINT","HALLMARK_INFLAMMATORY_RESPONSE",
             "HALLMARK_GLYCOLYSIS","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_MYC_TARGETS_V1"
             )
gsea_res_selected <- gsea_res[gsea_res$Description %in% pathway,]
write.csv(gsea_res_selected,"KEGG_HALLMARK_TvsN_GSEA_selected_16pathway.csv",quote = F,row.names = F)

###################################################################################################
##### step4. gene expression profiling #####
###################################################################################################
expDataFPKM <- readRDS("TCGA_LUAD_fpkm_pcg.rds")
table(duplicated(expDataFPKM$gene_name))
#FALSE  TRUE
#19938    24
expDataFPKM <- expDataFPKM[!duplicated(expDataFPKM$gene_name),]
rownames(expDataFPKM) <- expDataFPKM$gene_name
expDataFPKM <- expDataFPKM[,-1]
expDataFPKM[1:5,1:5]
#         TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07
#TSPAN6                        46.2924                      43.6109
#TNMD                           0.1239                       0.0955
#DPM1                          72.0583                      77.3323

# 读取metadata信息
meta <- readRDS("TCGA_LUAD_meta_data.rds")
head(meta)
#                                                  barcode      patient
#TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156
#TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774
#TCGA-91-6829-11A-01R-1858-07 TCGA-91-6829-11A-01R-1858-07 TCGA-91-6829

data2 <- t(expDataFPKM)
dim(data2)
#[1]   600 19938
data2 <- data2[rownames(meta),]

# 合并表达矩阵和metadata信息
data3 <- as.data.frame(cbind(meta,data2))
table(data3$group)
# Normal  Tumor 
#     59    541
saveRDS(data3,"TCGA_LUAD_mRNA_FPKM_metadata.rds")

# 定义一种主题，方便后面重复使用
theme_boxplot <- theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
                       axis.line=element_line(colour="black",size=0.25),
                       axis.title=element_text(size=15,face="plain",color="black"),
                       axis.text = element_text(size=15,face="plain",color="black"),
                       legend.position="none")

# PYGL expression in normal and tumor samples
my_comparisons <- list(c("Normal","Tumor"))

# boxplot
for(gene in c("GYS1","GYG1","GYG2","UGP2","PYGL","PYGB","PYGM","AGL")){
  pdf(paste0("Boxplot_wilcoxon_test_",gene,".pdf"),height=4.5,width=4.3)
  p <- ggplot(data3, aes(group, log2(get(gene)+1))) + 
    # geom_violin(aes(fill = group), width=0.6,trim = F) +
    geom_boxplot(aes(fill = group),width=0.5,notch = T,outlier.size = 0.5) +
    geom_jitter(size=0.5,position = position_jitter(0.2)) + 
    scale_fill_manual(values = c("#2DB2EB","#EB4232"))+ 
    theme_boxplot + stat_compare_means(comparisons = my_comparisons) + 
    ylab("Log2(FPKM+1)") + xlab("Sample") + ggtitle(paste0(gene," expression")) +
    theme(plot.title = element_text(size=15,hjust=0.5))
  print(p)
  dev.off()
}

table(data3$group,data3$stage)
#        Stage I Stage II Stage III Stage IV
# Normal      30       13        13        2
# Tumor      297      126        84       26

data3$Stage <- NA
data3[data3$group == "Normal","Stage"] <- "Normal"
data3[data3$group == "Tumor" & data3$ajcc_pathologic_stage %in% c("Stage I","Stage IA","Stage IB"),"Stage"] <- "Stage I"
data3[data3$group == "Tumor" & data3$ajcc_pathologic_stage %in% c("Stage II","Stage IIA","Stage IIB"),"Stage"] <- "Stage II"
data3[data3$group == "Tumor" & data3$ajcc_pathologic_stage %in% c("Stage IIIA","Stage IIIB"),"Stage"] <- "Stage III"
data3[data3$group == "Tumor" & data3$ajcc_pathologic_stage %in% c("Stage IV"),"Stage"] <- "Stage IV"
table(data3$Stage)
# Normal   Stage I  Stage II Stage III  Stage IV 
# 59       297       126        84        26

# PYGL表达与肿瘤分期的关系
for(gene in c("GYS1","GYG1","GYG2","UGP2","PYGL","PYGB","PYGM","AGL")){
  pdf(paste0("Boxplot_kw_test_normal_",gene,".pdf"),height=5,width=5)
  p <- ggplot(data3[data3$Stage %in% c("Normal","Stage I","Stage II","Stage III","Stage IV"),],
              aes(Stage, log2(get(gene)+1))) +
    #  geom_violin(aes(fill = group), width=0.6,trim = F) +
    geom_boxplot(aes(fill = Stage),width=0.5,notch = T,outlier.size = 0.5) +
    geom_jitter(size=0.5,position = position_jitter(0.2)) + scale_fill_npg() +
    theme_boxplot + 
    stat_compare_means(comparisons = list(c("Normal","Stage I"),c("Normal","Stage II"),c("Normal","Stage III"),c("Normal","Stage IV"))) +
    stat_compare_means(label.x = 3.5 ) + # Add global p-value
    ylab("Log2(FPKM+1)") + ggtitle(paste0(gene," expression")) +
    theme(plot.title = element_text(size=15,hjust=0.5))
  print(p)
  dev.off()
}

##### 提取肿瘤样本
data4 <- data3[data3$group=="Tumor",]
dim(data4)
#[1]   541 19960
table(data4$stage)
# Stage I  Stage II Stage III  Stage IV 
# 297       126        84        26

# 根据PYGL表达分组
data4$type <- ifelse(data4$PYGL > median(data4$PYGL),"PYGL_high","PYGL_low")
table(data4$type)
#PYGL_high  PYGL_low 
#270       271

# 添加边际图，可以使用ggExtra或psych
library(ggExtra)

p <- ggplot(data4, aes(log2(PYGL+1),log2(MKI67+1))) + geom_point(aes(color=type)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB"))+
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson") + theme(legend.position="bottom")
pdf("Correlation_PYGL-vs-MKI67_scatterplot.pdf",height=6.4,width=6)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(data4, aes(log2(PYGL+1),log2(MYC+1))) + geom_point(aes(color=type)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB"))+
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson") + theme(legend.position="bottom")
pdf("Correlation_PYGL-vs-MYC_scatterplot.pdf",height=6.4,width=6)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(data4, aes(log2(PYGL+1),log2(PDCD1+1))) + geom_point(aes(color=type)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB"))+
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson") + theme(legend.position="bottom")
pdf("Correlation_PYGL-vs-PDCD1_scatterplot.pdf",height=6.4,width=6)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(data4, aes(log2(PYGL+1),log2(CD274+1))) + geom_point(aes(color=type)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB"))+
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson") + theme(legend.position="bottom")
pdf("Correlation_PYGL-vs-PDL1_scatterplot.pdf",height=6.4,width=6)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(data4, aes(log2(PYGL+1),log2(CD276+1))) + geom_point(aes(color=type)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB"))+
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson") + theme(legend.position="bottom")
pdf("Correlation_PYGL-vs-CD276_scatterplot.pdf",height=6.4,width=6)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(data4, aes(log2(PYGL+1),log2(HAVCR2+1))) + geom_point(aes(color=type)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB"))+
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson") + theme(legend.position="bottom")
pdf("Correlation_PYGL-vs-HAVCR2_scatterplot.pdf",height=6.4,width=6)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()

###################################################################################################
### step5. survival analysis ########
###################################################################################################
# 加载包
library(survival)
library(survminer)

data3 <- readRDS("TCGA_LUAD_mRNA_FPKM_metadata.rds")
# 提取肿瘤样本
data4 <- data3[data3$group=="Tumor",]

### Starch通路基因集做生存分析
starch <- kegg[kegg$term =="KEGG_STARCH_AND_SUCROSE_METABOLISM",]$gene
data_starch <- data4[,colnames(data4) %in% starch]
#取基因集的均值
data_starch$set <- rowMeans(data_starch)
data_starch$group <- ifelse(data_starch$set > median(data_starch$set),"Starch_high","Starch_low")
table(data_starch$group)
#Starch_high  Starch_low 
#270         271
data_starch$OS_month <- data4$OS_month
data_starch$status <- data4$status

# 使用survfit()函数拟合KM生存曲线
fit <- survfit(Surv(OS_month, status) ~ group, data = data_starch)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
ggsurvplot(fit, data = data_starch,
           surv.median.line = "hv", # 添加中位数生存时间线
           # Change legends: title & labels
           legend = "right",
           legend.title = "Group", # 设置图例标题
           legend.labs = c("Starch high", "Starch low"), # 指定图例分组标签
           # Add p-value and tervals
           pval = TRUE, # 设置添加P值
           pval.method = TRUE, #设置添加P值计算方法
           conf.int = T, # 设置添加置信区间
           conf.int.style = "step",
           # Add risk table
           risk.table = TRUE, # 设置添加风险因子表
           tables.height = 0.2, # 设置风险表的高度
           tables.theme = theme_cleantable(), # 设置风险表的主题
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#EB4232","#2DB2EB"), # 设置颜色画板
           #palette = "lancet", 
           ggtheme = theme_bw() # Change ggplot2 theme
)

### Glycolysis通路基因集做生存分析
glycolysis <- hallmark[hallmark$term =="HALLMARK_GLYCOLYSIS",]$gene
data_glycolysis  <- data4[,colnames(data4) %in% glycolysis]
#取基因集的均值
data_glycolysis$set <- rowMeans(data_glycolysis)
data_glycolysis$group <- ifelse(data_glycolysis$set > median(data_glycolysis$set),"Glycolysis_high","Glycolysis_low")
table(data_glycolysis$group)
#Glycolysis_high  Glycolysis_low 
#270         271
data_glycolysis$OS_month <- data4$OS_month
data_glycolysis$status <- data4$status

# 使用survfit()函数拟合KM生存曲线
fit <- survfit(Surv(OS_month, status) ~ group, data = data_glycolysis)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
ggsurvplot(fit, data = data_glycolysis,
           surv.median.line = "hv", # 添加中位数生存时间线
           # Change legends: title & labels
           legend = "right",
           legend.title = "Group", # 设置图例标题
           legend.labs = c("Glycolysis high", "Glycolysis low"), # 指定图例分组标签
           # Add p-value and tervals
           pval = TRUE, # 设置添加P值
           pval.method = TRUE, #设置添加P值计算方法
           conf.int = T, # 设置添加置信区间
           conf.int.style = "step",
           # Add risk table
           risk.table = TRUE, # 设置添加风险因子表
           tables.height = 0.2, # 设置风险表的高度
           tables.theme = theme_cleantable(), # 设置风险表的主题
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#EB4232","#2DB2EB"), # 设置颜色画板
           #palette = "lancet", 
           ggtheme = theme_bw() # Change ggplot2 theme
)

data4$PYGL_group <- ifelse(data4$PYGL > median(data4$PYGL),"PYGL_high","PYGL_low")
table(data4$PYGL_group)
# PYGL_high  PYGL_low 
#     270       271

# 使用survfit()函数拟合KM生存曲线
fit <- survfit(Surv(OS_month, status) ~ PYGL_group, data = data4)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
ggsurvplot(fit, data = data4,
           surv.median.line = "hv", # 添加中位数生存时间线
           # Change legends: title & labels
           legend = "right",
           legend.title = "Group", # 设置图例标题
           legend.labs = c("PYGL high", "PYGL low"), # 指定图例分组标签
           # Add p-value and tervals
           pval = TRUE, # 设置添加P值
           pval.method = TRUE, #设置添加P值计算方法
           conf.int = T, # 设置添加置信区间
           conf.int.style = "step",
           # Add risk table
           risk.table = TRUE, # 设置添加风险因子表
           tables.height = 0.2, # 设置风险表的高度
           tables.theme = theme_cleantable(), # 设置风险表的主题
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#EB4232","#2DB2EB"), # 设置颜色画板
           #palette = "lancet", 
           ggtheme = theme_boxplot # Change ggplot2 theme
)

###################################################################################################
########## step6. 根据PYGL的表达分组进行差异分析 ##########
###################################################################################################
expDataCounts <- readRDS("TCGA_LUAD_counts_pcg.rds")
expDataCounts[1:5,1:5]
#                    gene_name TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07 TCGA-91-6829-11A-01R-1858-07
# ENSG00000000003.15    TSPAN6                         2296                         1403                          780
# ENSG00000000005.6       TNMD                            2                            1                            2
# ENSG00000000419.13      DPM1                          951                          662                          789
# ENSG00000000457.14     SCYL3                         1196                          256                          406
dim(expDataCounts)
#[1] 19962   601

table(duplicated(expDataCounts$gene_name))
#FALSE  TRUE
#19938    24
#remove duplicated genes
expDataCounts <- expDataCounts[!duplicated(expDataCounts$gene_name),]
rownames(expDataCounts) <- expDataCounts$gene_name
expDataCounts <- expDataCounts[,-1]
dim(expDataCounts)
#[1] 19938   600
expDataCounts[1:5,1:5]
#         TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07 TCGA-91-6829-11A-01R-1858-07 TCGA-69-A59K-01A-11R-A262-07
# TSPAN6                           2296                         1403                          780                         1338
# TNMD                                2                            1                            2                            0
# DPM1                              951                          662                          789                         2217
# SCYL3                            1196                          256                          406                         1101
# C1orf112                          141                          158                          146                          620

data4[1:5,1:5]
#                                                   barcode      patient           sample shortLetterCode          definition
# TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156 TCGA-78-7156-01A              TP Primary solid Tumor
# TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774 TCGA-44-6774-01A              TP Primary solid Tumor
# TCGA-69-A59K-01A-11R-A262-07 TCGA-69-A59K-01A-11R-A262-07 TCGA-69-A59K TCGA-69-A59K-01A              TP Primary solid Tumor
dim(data4)
#[1]   541 19969

expDataCounts <- expDataCounts[,rownames(data4)]
dim(expDataCounts)
# dim(expDataCounts)
# [1] 19938   541

colData <- data4[,c("gender","stage")]
table(data4$PYGL_group)
# PYGL_low PYGL_high 
# 271       270
colData$group <- data4$PYGL_group
head(colData)
#                              gender     stage     group
# TCGA-78-7156-01A-11R-2039-07   male  Stage IV PYGL_high
# TCGA-44-6774-01A-21R-1858-07 female Stage III PYGL_high
# TCGA-69-A59K-01A-11R-A262-07 female  Stage II PYGL_high

# 构建dds对象
dds <- DESeqDataSetFromMatrix(countData=round(expDataCounts, digits = 0),
                              colData=colData, 
                              design=~group)
dds
# class: DESeqDataSet 
# dim: 19938 541 
# metadata(1): version
# assays(1): counts
# rownames(55150): TSPAN6 TNMD ... LINC01144 AC007389.5
# rowData names(0):
#   colnames(594): TCGA-44-6777-11A-01R-1858-07 TCGA-49-6743-11A-01R-1858-07 ...
# TCGA-78-8655-01A-11R-2403-07 TCGA-38-4627-01A-01R-1206-07
# colData names(3): group gender stage

# 数据过滤
keep <- rowSums(counts(dds) >= 10) >= 3  #过滤低表达基因，至少有3个样品都满足10个以上的reads数
dds <- dds[keep,]

# 数据标准化
vsd <- vst(dds, blind=FALSE) #vst()函数效果和rlog()一样，且速度更快。
exprSet=assay(vsd)
exprSet[1:5,1:5]
#         TCGA-44-6777-11A-01R-1858-07 TCGA-49-6743-11A-01R-1858-07
# TSPAN6                      10.190136                    10.894689
# TNMD                         4.613970                     7.544313
# DPM1                        10.127907                    10.597370
# SCYL3                        9.166590                     9.290892
write.csv(exprSet,"TCGA-LUAD_tumor_count_vst_normalized.csv",quote=F)

# PCA降维
#plot by ggplot2
pcaData <- plotPCA(vsd, intgroup = c("group"), returnData = TRUE)
pcaData
#                                  PC1       PC2  group group.1                         name
# TCGA-44-6777-11A-01R-1858-07 38.85698 -5.120239 normal  normal TCGA-44-6777-11A-01R-1858-07
# TCGA-49-6743-11A-01R-1858-07 42.24485  6.620592 normal  normal TCGA-49-6743-11A-01R-1858-07
# TCGA-44-2662-11A-01R-1758-07 47.72076  6.561794 normal  normal TCGA-44-2662-11A-01R-1758-07
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  #  coord_fixed() +
  ggtitle("PCA with VST data") + theme_bw()

# 差异分析
dds <- DESeq(dds)
# using pre-existing size factors
# estimating dispersions
# found already estimated dispersions, replacing these
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 9072 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
resultsNames(dds)
#[1] "Intercept"             "group_PYGL_high_vs_PYGL_low"
res <- results(dds,contrast=c("group","PYGL_high","PYGL_low"))
summary(res)

# 结果保存
##筛选差异表达基因
#首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1 <- as.data.frame(res1)

#log2FC≥1 & padj<0.05 标识up，代表显著上调的基因
#log2FC≤-1 & padj<0.05 标识down，代表显著下调的基因
#其余标识 none，代表非差异的基因
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),'sig'] <- 'Up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05),'sig'] <- 'Down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.05),'sig'] <- 'None'
table(res1$sig)
#Down  None    Up 
#  551 17559   284
write.csv(res1, file = 'TCGA-LUAD_tumor_PYGL_high_vs_low.DESeq2.all.csv', quote = FALSE)

#输出选择的差异基因总表
res1_select <- subset(res1, sig %in% c('Up', 'Down'))
write.csv(res1_select, file = 'TCGA-LUAD_tumor_PYGL_high_vs_low.DESeq2.DEs.csv', quote = FALSE)

##### GSVA基因集变异分析 ######
library(GSVA)

data3 <- readRDS("TCGA_LUAD_mRNA_FPKM_metadata.rds")
# 提取肿瘤样本
data4 <- data3[data3$group=="Tumor",]
dim(data4)
#[1]   541 19960
table(data4$stage)
# Stage I  Stage II Stage III  Stage IV 
# 297       126        84        26
data4[1:5,1:5]
#                                                   barcode      patient           sample shortLetterCode          definition
# TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156-01A-11R-2039-07 TCGA-78-7156 TCGA-78-7156-01A              TP Primary solid Tumor
# TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774-01A-21R-1858-07 TCGA-44-6774 TCGA-44-6774-01A              TP Primary solid Tumor
# TCGA-69-A59K-01A-11R-A262-07 TCGA-69-A59K-01A-11R-A262-07 TCGA-69-A59K TCGA-69-A59K-01A              TP Primary solid Tumor

data_tumor <- t(data4[,-c(1:22)])
data_tumor[1:5,1:5]
#          TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07 TCGA-69-A59K-01A-11R-A262-07 TCGA-44-2665-01A-01R-A278-07
# TSPAN6                        14.0995                      12.7765                       6.3424                      22.3860
# TNMD                           0.0377                       0.0280                       0.0000                      21.5676
# DPM1                          21.9471                      22.6558                      39.4940                      46.6072
dim(data_tumor)
#[1] 19938   541

exp <- as.matrix(log2(data_tumor + 1))
exp[1:5,1:5]
#          TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07 TCGA-69-A59K-01A-11R-A262-07 TCGA-44-2665-01A-01R-A278-07
# TSPAN6                     3.91642887                   3.78413750                     2.876252                     4.547573
# TNMD                       0.05338942                   0.03984026                     0.000000                     4.496181
# DPM1                       4.52023994                   4.56412205                     5.339636                     5.573108
# SCYL3                      2.54599307                   1.34278228                     2.150365                     2.880842
# C1orf112                   0.72935699                   1.06571016                     1.692874                     1.632547

ferro <- read.gmt("GSVA/ferroptosis_genesets.gmt")
table(ferro$gene %in% rownames(exp))
#FALSE  TRUE 
#72   436 
ferro <- ferro[ferro$gene %in% rownames(exp),]
dim(ferro)
#[1] 436   2

# 提取基因集
geneset <- split(ferro$gene, ferro$term)

# 进行GSVA分析
es.max <- gsva(exp, geneset, 
               method="gsva",kcdf="Gaussian",
               parallel.sz=5)#线程数
# 写出GSVA结果到文件
write.table(es.max, file="GSVA/GSVA_ferroptosis_result.txt",sep="\t",quote=F,row.names = T)

data_df <- cbind(data4,t(es.max))

# 添加PYGL分组
data_df$PYGL_group <- ifelse(data_df$PYGL > median(data_df$PYGL),"PYGL_high","PYGL_low")
table(data_df$PYGL_group)
# PYGL_high  PYGL_low 
#     270       271
data_df$PYGL_group <- factor(data_df$PYGL_group,levels = c("PYGL_low","PYGL_high"))

# 添加铁死亡基因集分组
data_df$driver_group <- ifelse(data_df$ferroptosis_driver > median(data_df$ferroptosis_driver),"FerroDriver_high","FerroDriver_low")
table(data_df$driver_group)
# FerroDriver_high  FerroDriver_low 
# 270              271
data_df$suppressor_group <- ifelse(data_df$ferroptosis_suppressor > median(data_df$ferroptosis_suppressor),"FerroSuppressor_high","FerroSuppressor_low")
table(data_df$suppressor_group)
# FerroSuppressor_high  FerroSuppressor_low 
# 270                  271 

# PYGL表达与铁死亡的相关性
pdf("GSVA/Boxplot_wilcoxon_test_ferroptosis_driverr.pdf",height=4.5,width=4.5)
ggplot(data_df,aes(PYGL_group,ferroptosis_driver,fill=PYGL_group)) + 
  #  geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(aes(fill = PYGL_group),width=0.5,notch = T,outlier.size = 0.5) +
  geom_jitter(size=0.5,position = position_jitter(0.2)) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) + 
  theme_boxplot + stat_compare_means() + 
  theme(plot.title = element_text(size=15,hjust=0.5))
dev.off()
pdf("GSVA/Boxplot_wilcoxon_test_ferroptosis_suppressor.pdf",height=4.5,width=4.5)
ggplot(data_df,aes(PYGL_group,ferroptosis_suppressor,fill=PYGL_group)) + 
  #  geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(aes(fill = PYGL_group),width=0.5,notch = T,outlier.size = 0.5) +
  geom_jitter(size=0.5,position = position_jitter(0.2)) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) + 
  theme_boxplot + stat_compare_means() + 
  theme(plot.title = element_text(size=15,hjust=0.5))
dev.off()

data_df$driver_type <- "NA"
data_df[data_df$PYGL_group == "PYGL_high" & data_df$driver_group == "FerroDriver_high","driver_type"] <- "PYGL_high+FerroDriver_high"
data_df[data_df$PYGL_group == "PYGL_high" & data_df$driver_group == "FerroDriver_low","driver_type"] <- "PYGL_high+FerroDriver_low"
data_df[data_df$PYGL_group == "PYGL_low" & data_df$driver_group == "FerroDriver_high","driver_type"] <- "PYGL_low+FerroDriver_high"
data_df[data_df$PYGL_group == "PYGL_low" & data_df$driver_group == "FerroDriver_low","driver_type"] <- "PYGL_low+FerroDriver_low"
table(data_df$driver_type)
#PYGL_high+FerroDriver_high  PYGL_high+FerroDriver_low  PYGL_low+FerroDriver_high   PYGL_low+FerroDriver_low 
#152                        118                        118                        153
data_df$suppressor_type <- "NA"
data_df[data_df$PYGL_group == "PYGL_high" & data_df$suppressor_group == "FerroSuppressor_high","suppressor_type"] <- "PYGL_high+FerroSuppressor_high"
data_df[data_df$PYGL_group == "PYGL_high" & data_df$suppressor_group == "FerroSuppressor_low","suppressor_type"] <- "PYGL_high+FerroSuppressor_low"
data_df[data_df$PYGL_group == "PYGL_low" & data_df$suppressor_group == "FerroSuppressor_high","suppressor_type"] <- "PYGL_low+FerroSuppressor_high"
data_df[data_df$PYGL_group == "PYGL_low" & data_df$suppressor_group == "FerroSuppressor_low","suppressor_type"] <- "PYGL_low+FerroSuppressor_low"
table(data_df$suppressor_type)
# PYGL_high+FerroSuppressor_high  PYGL_high+FerroSuppressor_low  PYGL_low+FerroSuppressor_high   PYGL_low+FerroSuppressor_low 
# 154                            116                            116                            155

library(survival)
library(survminer)
# 使用survfit()函数拟合KM生存曲线
fit <- survfit(Surv(OS_month, status) ~ driver_type, data = data_df)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
pdf("GSVA/Survival_plot_ferroptosis_driver_type.pdf",height=6,width=7.5,onefile = F)
ggsurvplot(fit, data = data_df,
           surv.median.line = "hv", # 添加中位数生存时间线
           # Change legends: title & labels
           legend.title = "Group", # 设置图例标题
           #legend.labs = c("FerroDriver_high", "FerroDriver_low"), # 指定图例分组标签
           # Add p-value and tervals
           pval = TRUE, # 设置添加P值
           pval.method = TRUE, #设置添加P值计算方法
           conf.int = TRUE, # 设置添加置信区间
           conf.int.style = "step",
           # Add risk table
           risk.table = TRUE, # 设置添加风险因子表
           tables.height = 0.2, # 设置风险表的高度
           tables.theme = theme_cleantable(), # 设置风险表的主题
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           #palette = c("#EB4232","#2DB2EB","#E7B800", "#2E9FDF"), # 设置颜色画板
           palette = "nejm", 
           ggtheme = theme_bw() # Change ggplot2 theme
)
dev.off()

fit <- survfit(Surv(OS_month, status) ~ suppressor_type, data = data_df)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
pdf("GSVA/Survival_plot_ferroptosis_suppressor_type.pdf",height=6,width=8,onefile = F)
ggsurvplot(fit, data = data_df,
           surv.median.line = "hv", # 添加中位数生存时间线
           # Change legends: title & labels
           legend.title = "Group", # 设置图例标题
           #legend.labs = c("FerroSuppressor_high", "FerroSuppressor_low"), # 指定图例分组标签
           # Add p-value and tervals
           pval = TRUE, # 设置添加P值
           pval.method = TRUE, #设置添加P值计算方法
           conf.int = TRUE, # 设置添加置信区间
           conf.int.style = "step",
           # Add risk table
           risk.table = TRUE, # 设置添加风险因子表
           tables.height = 0.2, # 设置风险表的高度
           tables.theme = theme_cleantable(), # 设置风险表的主题
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           #palette = c("#EB4232","#2DB2EB","#E7B800", "#2E9FDF"), # 设置颜色画板
           palette = "nejm", 
           ggtheme = theme_bw() # Change ggplot2 theme
)
dev.off()

# 提取铁死亡Driver high组的样本
data_df2 <- data_df[data_df$driver_group == "FerroDriver_high",]
dim(data_df2)
#[1]   270 19967
table(data_df2$PYGL_group)
# PYGL_high  PYGL_low 
# 152       118

ggplot(data_df2,aes(PYGL_group,ferroptosis_driver,fill=PYGL_group)) + 
  #  geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(aes(fill = PYGL_group),width=0.5,notch = T,outlier.size = 0.5) +
  geom_jitter(size=0.5,position = position_jitter(0.2)) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) + 
  theme_boxplot + stat_compare_means() + 
  theme(plot.title = element_text(size=15,hjust=0.5))

ggplot(data_df2,aes(PYGL_group,ferroptosis_suppressor,fill=PYGL_group)) + 
  #  geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(aes(fill = PYGL_group),width=0.5,notch = T,outlier.size = 0.5) +
  geom_jitter(size=0.5,position = position_jitter(0.2)) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) + 
  theme_boxplot + stat_compare_means() + 
  theme(plot.title = element_text(size=15,hjust=0.5))

# 加载包
library(survival)
library(survminer)
# 使用survfit()函数拟合KM生存曲线
fit <- survfit(Surv(OS_month, status) ~ PYGL_group, data = data_df2)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
ggsurvplot(fit, data = data_df2,
           surv.median.line = "hv", # 添加中位数生存时间线
           # Change legends: title & labels
           legend.title = "Group", # 设置图例标题
           legend.labs = c("PYGL low", "PYGL high"), # 指定图例分组标签
           # Add p-value and tervals
           pval = TRUE, # 设置添加P值
           pval.method = TRUE, #设置添加P值计算方法
           conf.int = TRUE, # 设置添加置信区间
           conf.int.style = "step",
           # Add risk table
           risk.table = TRUE, # 设置添加风险因子表
           tables.height = 0.2, # 设置风险表的高度
           tables.theme = theme_cleantable(), # 设置风险表的主题
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#2DB2EB","#EB4232"), # 设置颜色画板
           #palette = "lancet", 
           ggtheme = theme_bw() # Change ggplot2 theme
)

###################################################################################################
##### step7. 免疫浸润分析 #####
###################################################################################################
library(preprocessCore)
library(immunedeconv)
library(tibble)

data3 <- readRDS("TCGA_LUAD_mRNA_FPKM_metadata.rds")
# 提取肿瘤样本
data4 <- data3[data3$group=="Tumor",]
data_tumor <- t(data4[,-c(1:22)])

### xcell deconvolute
result2 <- deconvolute(data_tumor, method = "xcell") 
result2 <- as.data.frame(result2)
rownames(result2) <- result2$cell_type
result2 <- as.data.frame(t(result2[,-1]))
write.csv(result2,"TCGA_LUAD_TME_celltype_xCell.csv",quote = F)

### 绘制箱线图
library(tibble)
library(tidyr)

# xCell results
dd1 <- result2 %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = 2:37,
               names_to = "CellType",
               values_to = "Composition")

plot.info <- dd1[,c(5,1,6)]

col <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 9,name = "Set1"))(40)

ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition") +
  theme_boxplot +
  theme(axis.text.x = element_text(
    angle = 60,
    hjust = 1,
    vjust = 1
  ),legend.position = "right") + scale_fill_manual(values = col) + ggpubr::stat_compare_means(label.x = 3)


###################################################################################################
##### step8. 免疫药物治疗预测分析 #####
###################################################################################################
### oncoPredict usage ###
library(oncoPredict)

#Set the seed for reproducibility. 
set.seed(12345)

dir='/data/LUAD_PYGL/drug/DataFiles/Training Data'

GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Expr[1:5,1:5]
#          COSMIC_906826 COSMIC_687983 COSMIC_910927 COSMIC_1240138 COSMIC_906792
# TSPAN6        7.632023      7.548671      8.712338       7.797142      7.074533
# TNMD          2.964585      2.777716      2.643508       2.817923      2.889677
# DPM1         10.379553     11.807341      9.880733       9.883471      9.773987
dim(GDSC2_Expr)
#[1] 17419   805
GDSC2_Res <- exp(GDSC2_Res)
GDSC2_Res[1:5,1:5]
#                Camptothecin_1003 Vinblastine_1004 Cisplatin_1005 Cytarabine_1006 Docetaxel_1007
# COSMIC_906826          0.3158373      0.208843106     1116.05899      18.5038719   0.0529577422
# COSMIC_687983          0.2827342      0.013664227       26.75839      16.2943594   0.0122399143
# COSMIC_910927          0.0295671      0.006684071       12.09379       0.3387418   0.0075478087
dim(GDSC2_Res)
#[1] 805 198

data_tumor[1:5,1:5]
#          TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07 TCGA-69-A59K-01A-11R-A262-07 TCGA-44-2665-01A-01R-A278-07
# TSPAN6                        14.0995                      12.7765                       6.3424                      22.3860
# TNMD                           0.0377                       0.0280                       0.0000                      21.5676
# DPM1                          21.9471                      22.6558                      39.4940                      46.6072

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = as.matrix(data_tumor),
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
testPtype[1:4, 1:4]
#                             V1 Camptothecin_1003 Vinblastine_1004 Cisplatin_1005
# 1 TCGA-55-7281-01A-11R-2039-07      2.913413e-04      0.012996676   4.752889e-02
# 2 TCGA-05-4433-01A-22R-1858-07      1.691975e+02               NA   8.181903e+09
# 3 TCGA-80-5611-01A-01R-1628-07      1.751992e-03      2.569609235   9.383037e-02

drugs <- c("Cisplatin", "Docetaxel", "Doxorubicin", "Gemcitabine")
drugs %in% gsub("_.*$","",colnames(testPtype))
#[1]  TRUE  TRUE FALSE  TRUE

# 添加PYGL分组
data4$PYGL_group <- ifelse(data4$PYGL > median(data4$PYGL),"PYGL_high","PYGL_low")
data4$PYGL_group <- factor(data4$PYGL_group,levels = c("PYGL_low","PYGL_high"))
group <- data.frame(id=rownames(data4),PYGL_group=data4$PYGL_group,PYGL=data4$PYGL)

testPtype_sub <- testPtype[group$id %in% testPtype$V1,]
testPtype_sub$group <- group[testPtype$V1 %in% group$id,"PYGL_group"]
testPtype_sub$PYGL <- group[testPtype$V1 %in% group$id,"PYGL"]
head(testPtype_sub)
#                             V1 Camptothecin_1003 Vinblastine_1004 Cisplatin_1005 Cytarabine_1006 Docetaxel_1007 Gefitinib_1010
# 1 TCGA-78-7156-01A-11R-2039-07        0.15055968       0.03752343       68.24723        4.885721    0.017764738       26.51723
# 2 TCGA-44-6774-01A-21R-1858-07        0.12517871       0.03299005       25.29981        8.763347    0.014328061       38.14545
# 3 TCGA-69-A59K-01A-11R-A262-07        0.08378431       0.04307948       35.51561       12.216783    0.028561700       41.06495
# 4 TCGA-44-2665-01A-01R-A278-07        0.08359785       0.01066778       13.23418        2.182081    0.008535937       29.43170

library(ggpubr)
library(ggExtra)

# 散点图
p <- ggplot(testPtype_sub,aes(x=log2(PYGL+1),y=log2(Cisplatin_1005+1))) + geom_point(aes(color=group)) + 
  scale_color_manual(values = c("#2DB2EB","#EB4232")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 15) + stat_cor(method = "pearson") + ylim(0,12)
pdf("Correlation_PYGL-vs-Cisplatin_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(testPtype_sub,aes(x=log2(PYGL+1),y=log2(Gemcitabine_1190+1))) + geom_point(aes(color=group)) +
  scale_color_manual(values = c("#2DB2EB","#EB4232")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 15) + stat_cor(method = "pearson") + ylim(0,3)
pdf("Correlation_PYGL-vs-Gemcitabine_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(testPtype_sub,aes(x=log2(PYGL+1),y=log2(Docetaxel_1007+1))) + geom_point(aes(color=group)) + 
  scale_color_manual(values = c("#2DB2EB","#EB4232")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 15) + stat_cor(method = "pearson") + ylim(0,0.05)
pdf("Correlation_PYGL-vs-Docetaxel_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(testPtype_sub,aes(x=log2(PYGL+1),y=log2(Gefitinib_1010+1))) + geom_point(aes(color=group)) + 
  scale_color_manual(values = c("#2DB2EB","#EB4232")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 15) + stat_cor(method = "pearson") + ylim(0,7)
pdf("Correlation_PYGL-vs-Gefitinib_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(testPtype_sub,aes(x=log2(PYGL+1),y=log2(Paclitaxel_1080+1))) + geom_point(aes(color=group)) + 
  scale_color_manual(values = c("#2DB2EB","#EB4232")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 15) + stat_cor(method = "pearson") + ylim(0,0.5)
pdf("Correlation_PYGL-vs-Paclitaxel_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(testPtype_sub,aes(x=log2(PYGL+1),y=log2(Erlotinib_1168+1))) + geom_point(aes(color=group)) + 
  scale_color_manual(values = c("#2DB2EB","#EB4232")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 15) + stat_cor(method = "pearson") + ylim(0,6)
pdf("Correlation_PYGL-vs-Erlotinib_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()

###################################################################################################
##### step9. ICB免疫治疗响应分析 #####
###################################################################################################
data <- read.table("./CTR_DB_immune drug/anti_PD1_NSCLC_datasets/GSE135222/GSE135222_PYGL_selected.txt",
                   header = T,sep = "\t",check.names = F)
head(data)
# PYGL_group                    type PFS PFS_event benefit        ID TIGIT  PYGL HAVCR2 CD274 PDCD1 HAVCR1 CTLA4 ICOS  TOX
# 1  PYGL_High     PYGL_High_Responder 412         0       Y  NSCLC378  5.52 41.20  21.64 16.70  2.30   2.02  3.70 3.53 1.22
# 2   PYGL_Low  PYGL_Low_non_Responder  29         1       N NSCLC1104  0.58 13.72   2.95  7.20  0.16   0.00  0.41 0.08 0.28
# 3   PYGL_Low      PYGL_Low_Responder 618         0       Y  NSCLC947  3.01 11.67  14.02  9.73  3.11   1.29  6.65 2.52 1.31
data$Responder <- ifelse(data$benefit == "Y","R","NR")
data_summary <- table(data$PYGL_group, data$Responder)
data_summary
#           NR  R
# PYGL_High  6  7
# PYGL_Low  13  1
chisq.test(data_summary)
# Pearson's Chi-squared test with Yates' continuity correction
# data:  data_summary
# X-squared = 4.9895, df = 1, p-value = 0.0255
data_summary <- reshape2::melt(data_summary)
data_summary
# Group Response Value
# 1 PYGL_High       NR     6
# 2  PYGL_Low       NR    13
# 3 PYGL_High        R     7
# 4  PYGL_Low        R     1
colnames(data_summary) <- c("Group","Response","Value")
data_summary$Response <- factor(data_summary$Response,levels = c("R","NR"))
data_summary$Group <- factor(data_summary$Group,levels = c("PYGL_Low","PYGL_High"))
data_summary$Percent <- c("46.15","92.86","53.85","7.14")
p <- ggplot(data_summary,aes(Group,Value,fill=Response)) + 
  geom_bar(stat = "identity",position = "fill") + theme_boxplot +
  scale_fill_manual(values = c("#FF1439","#40E0D0")) + theme(legend.position = "top") +
  xlab("X-squared = 4.9895, df = 1, p-value = 0.0255") + ylab("Relative percent (%)")
print(p)

p <- ggplot(data, aes(x = Responder, y = log1p(PDCD1), fill = Responder)) +
  geom_flat_violin(aes(fill = Responder),position = position_nudge(x = 0.1, y = 0), trim = TRUE, alpha = .8, colour = NA)+
  geom_point(aes(x = .65, y = log1p(PDCD1), colour = Responder),position = position_jitter(width = .1), size = 3, shape = 20)+
  geom_boxplot(aes(x = Responder, y = log1p(PDCD1), fill = Responder),outlier.shape = NA, alpha = .8, width = .1, colour = "black")+
  scale_colour_manual(values =  c("#40E0D0","#FF1439"))+
  scale_fill_manual(values = c("#40E0D0","#FF1439"))+
  theme_classic()+
  geom_signif(comparisons = list(c("R","NR")),map_signif_level = TRUE,test = "wilcox.test")+
  theme(axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text = element_text(size=14))+
  # ylim(0,1)+
  ylab('Expression')
p
p <- ggplot(data, aes(x = Responder, y = log1p(CD276), fill = Responder)) +
  geom_flat_violin(aes(fill = Responder),position = position_nudge(x = 0.1, y = 0), trim = TRUE, alpha = .8, colour = NA)+
  geom_point(aes(x = .65, y = log1p(CD276), colour = Responder),position = position_jitter(width = .1), size = 3, shape = 20)+
  geom_boxplot(aes(x = Responder, y = log1p(CD276), fill = Responder),outlier.shape = NA, alpha = .8, width = .1, colour = "black")+
  scale_colour_manual(values =  c("#40E0D0","#FF1439"))+
  scale_fill_manual(values = c("#40E0D0","#FF1439"))+
  theme_classic()+
  geom_signif(comparisons = list(c("R","NR")),map_signif_level = TRUE,test = "wilcox.test")+
  theme(axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text = element_text(size=14))+
  # ylim(0,1)+
  ylab('Expression')
p

##### ICB signature GSVA enrichment ######
markers <- clusterProfiler::read.gmt("ICB_signatures.gmt")
markers <- markers[markers$gene != "",]

geneSets <- split(markers$gene,markers$term)
names(geneSets)
# [1] "IFN-gamma"                                      "T cell-inflamed GEP(gene expression profiles) "
# [3] "IMPRES (immuno-predictive score)"               "Inflammatory.Sig"                              
# [5] "EMT.Sig"                                        "Blood.Sig"                                     
# [7] "ImmuneCells.Sig"                                "CIBERSORT.CD8"                                 
# [9] "TIS(T Cell Infiltration Score)"                 "CYT(immune cytolytic activity)"                
# [11] "AMP(Antigen Presenting Machinery)"              "IPS (Immunophenoscore)"     
gsva <- gsva(as.matrix(data_tumor), geneSets, method="gsva",kcdf="Gaussian",parallel.sz=5)#线程数
gsva[1:5,1:5]
# TCGA-78-7156-01A-11R-2039-07 TCGA-44-6774-01A-21R-1858-07 TCGA-69-A59K-01A-11R-A262-07
# IFN-gamma                                                        -0.2801620                  -0.08084813                   0.27570117
# T cell-inflamed GEP(gene expression profiles)                    -0.5252799                   0.09362897                   0.19158668
# IMPRES (immuno-predictive score)                                 -0.7403156                   0.50006056                  -0.03203424                                                      0.36125419                    0.4741209                   0.49612633

# PYGL分组
data4$PYGL_group <- ifelse(data4$PYGL > median(data4$PYGL),"PYGL_high","PYGL_low")
table(data4$PYGL_group)
# PYGL_high  PYGL_low 
# 270       271

gsva_res <- as.data.frame(t(gsva))
gsva_res$group <- data4$PYGL_group
head(gsva_res)
# IFN-gamma T cell-inflamed GEP(gene expression profiles)  IMPRES (immuno-predictive score) Inflammatory.Sig
# TCGA-78-7156-01A-11R-2039-07 -0.28016203                                    -0.52527993                      -0.74031563      -0.56433304
# TCGA-44-6774-01A-21R-1858-07 -0.08084813                                     0.09362897                       0.50006056      -0.17826471
# TCGA-69-A59K-01A-11R-A262-07  0.27570117                                     0.19158668                      -0.03203424       0.09090266
df <- reshape::melt(gsva_res)
head(df)
# group  variable       value
# 1 PYGL_high IFN-gamma -0.28016203
# 2  PYGL_low IFN-gamma -0.08084813
# 3 PYGL_high IFN-gamma  0.27570117
df$group <- factor(df$group,levels = c("PYGL_low","PYGL_high"))

ggplot(df, aes(group, value)) + 
  # geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(aes(fill = group),width=0.5,notch = T,outlier.size = 0.5) +
  geom_jitter(size=0.5,position = position_jitter(0.2)) + scale_fill_aaas() + 
  theme_boxplot + stat_compare_means(aes(group = group),
                                     method = "wilcox.test",
                                     label = "p.signif",
                                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                      symbols = c("***", "**", "*", "ns"))) +
  facet_wrap(.~variable) + 
  ylab("GSVA score") + xlab("Group") +
  theme(plot.title = element_text(size=12,hjust=0.5))

###################################################################################################
##### step9. 肿瘤突变负荷分析 #####
###################################################################################################
### TCIA数据库neoantigen分析
data <- read.table("TCIA-HeterogeneityData.tsv",header = T,sep = "\t",check.names = F)
data[1:5,1:10]
#        Patient Disease   AUC    KS fraction_clonal fraction_clonal_NEO fraction_sub_clonal fraction_sub_clonal_NEO num_clonal
# 1 TCGA-67-3773    LUAD 0.149 0.730           0.352               0.345               0.648                   0.655         43
# 2 TCGA-91-6848    LUAD 0.010 0.091           0.866               0.883               0.134                   0.117        505
# 3 TCGA-55-6986    LUAD 0.169 0.682           0.364               0.826               0.636                   0.174          8
# 4 TCGA-73-4677    LUAD 0.039 0.195           0.805               0.818               0.195                   0.182        215
# 5 TCGA-78-7540    LUAD 0.136 0.548           0.419               0.690               0.581                   0.310         13
# num_clonal_NEO
# 1             48
# 2             98
# 3             19
# 4            139
# 5             20

group <- data.frame(id=data4$barcode,PYGL_group=data4$PYGL_group,PYGL=data4$PYGL)
dim(group)
#[1] 541  23
head(group)
#                             id PYGL_group    PYGL
# 1 TCGA-78-7156-01A-11R-2039-07  PYGL_high 15.4226
# 2 TCGA-44-6774-01A-21R-1858-07   PYGL_low 10.2042
# 3 TCGA-69-A59K-01A-11R-A262-07  PYGL_high 34.8445
# 4 TCGA-44-2665-01A-01R-A278-07  PYGL_high 17.8359

group$barcode<-substr(group$id,1,12)
head(group)
#                             id PYGL_group      barcode    PYGL
# 1 TCGA-78-7156-01A-11R-2039-07  PYGL_high TCGA-78-7156 15.4226
# 2 TCGA-44-6774-01A-21R-1858-07   PYGL_low TCGA-44-6774 10.2042
# 3 TCGA-69-A59K-01A-11R-A262-07  PYGL_high TCGA-69-A59K 34.8445
NEO <- merge(data,group, by.x = "Patient",by.y = "barcode")

library(ggstatsplot)
library(ggplot2)
library(ggExtra)

#散点图
p <- ggplot(NEO,aes(x=log2(PYGL+1),y= log2(num_clonal_NEO+1))) + geom_point(aes(color=PYGL_group)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson")
pdf("Correlation_PYGL-vs-num_clonal_NEO_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(NEO,aes(x=log2(PYGL+1),y= log2(num_sub_clonal_NEO+1))) + geom_point(aes(color=PYGL_group)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson")
pdf("Correlation_PYGL-vs-num_sub_clonal_NEO_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(NEO,aes(x=log2(PYGL+1),y= log2(num_clonal+1))) + geom_point(aes(color=PYGL_group)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson")
pdf("Correlation_PYGL-vs-num_clonal_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()
p <- ggplot(NEO,aes(x=log2(PYGL+1),y= log2(num_sub_clonal+1))) + geom_point(aes(color=PYGL_group)) + 
  scale_color_manual(values = c("#EB4232","#2DB2EB")) +
  geom_smooth(method = "lm") + theme_bw(base_size = 14) + stat_cor(method = "pearson")
pdf("Correlation_PYGL-vs-num_sub_clonal_scatterplot.pdf",height=6,width=7.5)
ggMarginal(p,type ="densigram", margins = "both", fill="#A7B9D7", color="steelblue")
dev.off()

#箱式图
ggplot(NEO,aes(x=PYGL_group,y=num_clonal,fill=PYGL_group))+geom_boxplot()+
  theme_bw(base_size = 14)+ggpubr::stat_compare_means()
ggplot(NEO,aes(x=PYGL_group,y=num_clonal,fill=PYGL_group))+geom_boxplot()+
  theme_bw(base_size = 14)+ggpubr::stat_compare_means()+geom_jitter()+ylim(0,1000)

pdf("Boxplot_num_clonal_NEO_PYGL_group.pdf",height=6,width=6.8)
ggplot(NEO,aes(x=PYGL_group,y=log2(num_clonal_NEO+1),fill=PYGL_group))+
  geom_boxplot(width=0.5,notch = T)+theme_bw(base_size = 15)+
  scale_fill_manual(values = c("#2DB2EB","#EB4232"))+theme(panel.grid = element_blank()) + 
  ggpubr::stat_compare_means()+geom_jitter(width = 0.2)
dev.off()
pdf("Boxplot_num_sub_clonal_NEO_PYGL_group.pdf",height=6,width=6.8)
ggplot(NEO,aes(x=PYGL_group,y=log2(num_sub_clonal_NEO+1),fill=PYGL_group))+
  geom_boxplot(width=0.5,notch = T)+theme_bw(base_size = 15)+
  scale_fill_manual(values = c("#2DB2EB","#EB4232"))+theme(panel.grid = element_blank()) + 
  ggpubr::stat_compare_means()+geom_jitter(width = 0.2)
dev.off()
pdf("Boxplot_num_clonal_PYGL_group.pdf",height=6,width=6.8)
ggplot(NEO,aes(x=PYGL_group,y=log2(num_clonal+1),fill=PYGL_group))+
  geom_boxplot(width=0.5,notch = T)+theme_bw(base_size = 15)+
  scale_fill_manual(values = c("#2DB2EB","#EB4232"))+theme(panel.grid = element_blank()) + 
  ggpubr::stat_compare_means()+geom_jitter(width = 0.2)
dev.off()
pdf("Boxplot_num_sub_clonal_PYGL_group.pdf",height=6,width=6.8)
ggplot(NEO,aes(x=PYGL_group,y=log2(num_sub_clonal+1),fill=PYGL_group))+
  geom_boxplot(width=0.5,notch = T)+theme_bw(base_size = 18)+
  scale_fill_manual(values = c("#2DB2EB","#EB4232"))+theme(panel.grid = element_blank()) + 
  ggpubr::stat_compare_means()+geom_jitter(width = 0.2)
dev.off()

