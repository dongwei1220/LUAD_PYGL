### GEO-LUAD数据分析 ###
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

setwd("/data/LUAD_PYGL")

# 定义一种主题，方便后面重复使用
theme_boxplot <- theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
                       axis.line=element_line(colour="black",size=0.25),
                       axis.title=element_text(size=15,face="plain",color="black"),
                       axis.text = element_text(size=15,face="plain",color="black"),
                       legend.position="none")

###################################################################################################
##### cell paper LUAD data #####
###################################################################################################
### mRNA data
cell_mrna <- read.table("GEO/CELL/mRNA_FPKM.txt",header = T,sep = "\t",row.names = 1,check.names = F)
cell_mrna[1:5,1:5]
#            LUAD_17T   LUAD_18T    LUAD_19T    LUAD_20T    LUAD_22T
# A1BG    1.252622106 0.30719541 0.503480839 0.634031629 1.118910992
# A1CF    0.001838387         NA 0.003263509 0.008952105 0.006170243
# A2ML1   0.632022599 0.01360909 0.003225077          NA 0.004319388
dim(cell_mrna)
#[1] 16188   100

cell_meta <- read.table("GEO/CELL/clin_data.txt",header = T,sep = "\t",check.names = F)
head(cell_meta)
#   Sample ID Mutation counts TP53 mut (0:no; 1:yes) EGFR mut (0:no; 1:yes) KRAS mut (0:no; 1:yes) Gender (1:Male;2:Female) Age
# 1   LUAD_01              56                      1                      0                      0                        2  42
# 2   LUAD_02              67                      1                      1                      0                        1  54
# 3   LUAD_03              73                      0                      0                      0                        1  58
# 4   LUAD_04             776                      1                      0                      1                        1  70
dim(cell_meta)
#[1] 103  31

cell_mrna_df <- as.data.frame(t(cell_mrna))
cell_mrna_df[1:5,1:5]
cell_mrna_df$group <- gsub("LUAD_\\d+","",rownames(cell_mrna_df))
cell_mrna_df$group <- ifelse(cell_mrna_df$group == "T","Tumor","Normal")
table(cell_mrna_df$group)
# Normal  Tumor 
# 49     51 

genes <- c("GYS1","UGP2","GYG1","GYG2","PYGL","PYGB","PYGM","AGL")
cell_selected <- cell_mrna_df[,c(genes,"group")]
head(cell_selected)
#               GYS1     GYG1      GYG2      UGP2      PYGL     PYGB      PYGM       AGL group
# LUAD_17T 10.077363 7.929824 2.0891887  6.757389 14.107913 26.22127 0.2520733 1.5243676 Tumor
# LUAD_18T  7.223650 7.093139 0.5241234 11.741969  7.356141 30.10305 0.1168133 3.8900782 Tumor
# LUAD_19T  8.374919 6.074122 0.1807786  9.262701 23.937505 25.82397 0.1098050 0.8377212 Tumor

library(reshape2)
cell_selected <- melt(cell_selected,variable.name = "gene",value.name = "FPKM")
head(cell_selected)
#   group gene      FPKM
# 1 Tumor GYS1 10.077363
# 2 Tumor GYS1  7.223650
# 3 Tumor GYS1  8.374919

ggplot(cell_selected, aes(gene, log2(FPKM+1),fill=group)) + 
  geom_boxplot(width=1,notch = T,outlier.size = 0.3) +
  geom_jitter(size=0.6,position = position_jitter(0.5),alpha=0.8) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) + 
  theme_boxplot + stat_compare_means(aes(group = group),
                                     method = "wilcox.test",
                                     label = "p.signif",
                                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                      symbols = c("***", "**", "*", "ns"))) +
  ylab("Log2(FPKM+1)") + xlab("Gene") +
  theme(plot.title = element_text(size=15,hjust=0.5),
        axis.text.x  = element_text(size=15,face="italic",color="black"),
        legend.position="right")

# PYGL表达与肿瘤分期的关系
cell_mrna_df2 <- cell_mrna_df
table(cell_mrna_df2$group)
# Normal  Tumor 
# 49     51
cell_mrna_df2[cell_mrna_df2$group == "Normal","stage"] <- "Normal"
cell_mrna_df2[cell_mrna_df2$group == "Tumor","stage"] <- cell_meta[cell_meta$`Sample ID` %in% gsub("[T]","",rownames(cell_mrna_df2[cell_mrna_df2$group == "Tumor",])),]$`Tumor stage(A/B)`
table(cell_mrna_df2$stage)
# IA1    IA2    IA3     IB    IIA    IIB   IIIA   IIIB    IVA Normal 
# 2      1      6     13      2      7     14      5      1     49
cell_mrna_df2$Stage <- NA
cell_mrna_df2[cell_mrna_df2$group == "Normal","Stage"] <- "Normal"
cell_mrna_df2[cell_mrna_df2$group == "Tumor" & cell_mrna_df2$stage %in% c("IA1","IA2","IA3","IB"),"Stage"] <- "Stage I"
cell_mrna_df2[cell_mrna_df2$group == "Tumor" & cell_mrna_df2$stage %in% c("IIA","IIB"),"Stage"] <- "Stage II"
cell_mrna_df2[cell_mrna_df2$group == "Tumor" & cell_mrna_df2$stage %in% c("IIIA","IIIB"),"Stage"] <- "Stage III"
cell_mrna_df2[cell_mrna_df2$group == "Tumor" & cell_mrna_df2$stage %in% c("IVA"),"Stage"] <- "Stage IV"
table(cell_mrna_df2$Stage)
# Normal   Stage I  Stage II Stage III  Stage IV 
# 49        22         9        19         1 

for(gene in c("GYS1","GYG1","GYG2","UGP2","PYGL","PYGB","PYGM","AGL")){
  pdf(paste0("Boxplot_kw_test_normal_cell_",gene,".pdf"),height=5,width=5)
  p <- ggplot(cell_mrna_df2[cell_mrna_df2$Stage %in% c("Normal","Stage I","Stage II","Stage III"),],
              aes(Stage, log2(get(gene)+1))) +
    #  geom_violin(aes(fill = group), width=0.6,trim = F) +
    geom_boxplot(aes(fill = Stage),width=0.5,notch = T,outlier.size = 0.5) +
    geom_jitter(size=0.5,position = position_jitter(0.2)) + scale_fill_npg() +
    theme_boxplot + 
    stat_compare_means(comparisons = list(c("Normal","Stage I"),c("Normal","Stage II"),c("Normal","Stage III"))) +
    stat_compare_means(label.x = 3 ) + # Add global p-value
    ylab("Log2(FPKM+1)") + ggtitle(paste0(gene," expression")) +
    theme(plot.title = element_text(size=15,hjust=0.5))
  print(p)
  dev.off()
}

### protein data
cell_protein <- read.table("GEO/CELL/protein_intensity.txt",header = T,sep = "\t",check.names = F,row.names = 1)
cell_protein[1:5,1:5]
#        iBAQ LUAD_01T iBAQ LUAD_02T iBAQ LUAD_03T iBAQ LUAD_04T iBAQ LUAD_05T
# A1BG       140454759     327608842      64750505     143022749      73436305
# A1CF              NA            NA            NA            NA            NA
# A2M        222419458     265796590     244582962      98489218     132576705
dim(cell_protein)
#[1] 11056   206

cell_protein_df <- as.data.frame(t(cell_protein))
cell_protein_df[1:5,1:5]
cell_protein_df$group <- gsub("iBAQ LUAD_\\d+","",rownames(cell_protein_df))
cell_protein_df$group <- ifelse(cell_protein_df$group == "T","Tumor","Normal")
table(cell_protein_df$group)
# Normal  Tumor 
# 103     103 

genes <- c("GYS1","UGP2","GYG1","GYG2","PYGL","PYGB","PYGM","AGL")
cell_selected <- cell_protein_df[,c(genes,"group")]
head(cell_selected)
#                  GYS1     GYG1     GYG2      UGP2     PYGL      PYGB     PYGM       AGL group
# iBAQ LUAD_01T 6896093 40132235       NA  60787405 11910727  77761632 457957.5 8318903.7 Tumor
# iBAQ LUAD_02T 7897758 20022588       NA 129717596 19581026  41257395       NA 5210492.3 Tumor
# iBAQ LUAD_03T 8280315 14217435 840260.5  61064535 28319935  47229687 518928.6 8994832.9 Tumor
# iBAQ LUAD_04T 3308851 11691814       NA  11161743  1759434  33190798       NA  670604.4 Tumor

library(reshape2)
cell_selected <- melt(cell_selected,variable.name = "gene",value.name = "Intensity")
head(cell_selected)
#   group gene Intensity
# 1 Tumor GYS1   6896093
# 2 Tumor GYS1   7897758
# 3 Tumor GYS1   8280315

ggplot(cell_selected, aes(gene, log2(Intensity+1),fill=group)) + 
  #  geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(width=1,notch = T,outlier.size = 0.3) +
  geom_jitter(size=0.6,position = position_jitter(0.5),alpha=0.8) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) +
  theme_boxplot + stat_compare_means(aes(group = group),
                                     method = "wilcox.test",
                                     label = "p.signif",
                                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                      symbols = c("***", "**", "*", "ns"))) +
  ylab("Log2(Intensity+1)") + xlab("Gene") +
  theme(plot.title = element_text(size=15,hjust=0.5),
        axis.text.x  = element_text(size=15,face="italic",color="black"),
        legend.position="right")

### 生存预后分析
cell_mrna_tumor <- cell_mrna_df[cell_mrna_df$group == "Tumor",]
cell_mrna_tumor[1:5,1:5]
#               A1BG        A1CF       A2ML1    A4GALT      AAAS
# LUAD_17T 1.2526221 0.001838387 0.632022599 6.1695628 11.282293
# LUAD_18T 0.3071954          NA 0.013609092 7.1856388  7.060750
# LUAD_19T 0.5034808 0.003263509 0.003225077 1.3891861  5.970035
# LUAD_20T 0.6340316 0.008952105          NA 0.3995857  5.875602
# LUAD_22T 1.1189110 0.006170243 0.004319388 1.8500520  5.809641
dim(cell_mrna_tumor)
#[1]    51 16189
cell_mrna_tumor$sample <- gsub("T","",rownames(cell_mrna_tumor))

head(cell_meta)
#   Sample ID Mutation counts TP53 mut (0:no; 1:yes) EGFR mut (0:no; 1:yes) KRAS mut (0:no; 1:yes) Gender (1:Male;2:Female) Age
# 1   LUAD_01              56                      1                      0                      0                        2  42
# 2   LUAD_02              67                      1                      1                      0                        1  54
# 3   LUAD_03              73                      0                      0                      0                        1  58
# 4   LUAD_04             776                      1                      0                      1                        1  70
dim(cell_meta)
#[1] 103  31

cell_mrna_meta <- merge(cell_mrna_tumor,cell_meta,
                        by.x="sample",by.y="Sample ID")
cell_mrna_meta$status <- cell_mrna_meta$OS
colnames(cell_mrna_meta)[16218] <- "OS_day"
colnames(cell_mrna_meta)[16219] <- "OS_month"
dim(cell_mrna_meta)
#[1]    51 16221

# 加载包
library(survival)
library(survminer)

cell_mrna_meta$PYGL_group <- ifelse(cell_mrna_meta$PYGL > median(cell_mrna_meta$PYGL),"PYGL_high","PYGL_low")
table(cell_mrna_meta$PYGL_group)
# PYGL_high  PYGL_low 
#     25       26

# 使用survfit()函数拟合KM生存曲线
fit <- survfit(Surv(OS_month, status) ~ PYGL_group, data = cell_mrna_meta)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
ggsurvplot(fit, data = cell_mrna_meta,
           surv.median.line = "hv", # 添加中位数生存时间线
           # Change legends: title & labels
           legend.title = "Group", # 设置图例标题
           legend = "right",
           legend.labs = c("PYGL high", "PYGL low"), # 指定图例分组标签
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
           palette = c("#EB4232","#2DB2EB"), # 设置颜色画板
           #palette = "lancet", 
           ggtheme = theme_boxplot # Change ggplot2 theme
)

###################################################################################################
##### GSE37745 paper LUAD data #####
###################################################################################################
gse_mrna <- read.table("GEO/GSE37745/GSE37745_expression.txt",header = T,sep = "\t",row.names = 1,check.names = F)
gse_mrna[1:5,1:5]
#          GSM1019138 GSM1019139 GSM1019140 GSM1019141 GSM1019142
# A1BG       6.453640   6.056568   6.446127   7.109993   6.887309
# A1BG-AS1   5.368951   5.901213   5.678825   5.948762   5.778068
# A1CF       4.145956   4.210113   5.579812   4.230375   4.089493
# A2M        8.504073   9.276961   9.143159   8.496978   8.393224
# A2M-AS1    4.961126   5.442135   5.900450   4.720295   5.353107
dim(gse_mrna)
#[1] 23520   196

gse_meta <- read.table("GEO/GSE37745/GSE37745_clinical.txt",header = T,sep = "\t",check.names = F)
head(gse_meta)
#    accession                      title platform                       source     organism dead days_to_determined_death_status
# 1 GSM1019138  Patient 6, male, squamous   GPL570  cancer cells from patient 6 Homo sapiens   no                            3418
# 2 GSM1019139  Patient 11, female, adeno   GPL570 cancer cells from patient 11 Homo sapiens   no                            3186
# 3 GSM1019140    Patient 12, male, adeno   GPL570 cancer cells from patient 12 Homo sapiens  yes                             162
# 4 GSM1019141 Patient 13, male, squamous   GPL570 cancer cells from patient 13 Homo sapiens  yes                               6
dim(gse_meta)
#[1] 196  16
table(gse_meta$histology)
# adeno    large squamous 
# 106       24       66  
gse_meta <- gse_meta[gse_meta$histology == "adeno",]
dim(gse_meta)
#[1] 106  16
gse_mrna <- gse_mrna[,gse_meta$accession]
dim(gse_mrna)
#[1] 23520   106

gse_meta_df <- data.frame(sample=gse_meta$accession,age=gse_meta$age,sex=gse_meta$gender,
                          stage=paste0("Stage ",toupper(gse_meta$tumor_stage)),
                          status=gse_meta$dead,OS=gse_meta$days_to_determined_death_status)
gse_meta_df$OS_month <- round(gse_meta_df$OS/30,2) #以month为单位，保留两位小数
gse_meta_df$status <- ifelse(gse_meta_df$status == "yes",1,0)
head(gse_meta_df)
#       sample age    sex    stage status   OS OS_month
# 1 GSM1019139  66 female Stage 1A      0 3186   106.20
# 2 GSM1019140  77   male Stage 1A      1  162     5.40
# 3 GSM1019143  40 female Stage 1B      1  556    18.53
# 4 GSM1019144  68   male Stage 3A      1  304    10.13
gse_mrna_df <- as.data.frame(t(gse_mrna))
gse_mrna_df$sample <- rownames(gse_mrna_df)
gse_mrna_df[1:5,1:5]
#                A1BG A1BG-AS1     A1CF      A2M  A2M-AS1
# GSM1019139 6.056568 5.901213 4.210113 9.276961 5.442135
# GSM1019140 6.446127 5.678825 5.579812 9.143159 5.900450
# GSM1019143 7.551140 6.351549 4.328969 8.955866 5.254130
# GSM1019144 6.768441 5.861642 4.431411 8.803923 5.141022
# GSM1019145 7.312181 6.039282 4.312239 8.935653 5.909488

gse_mrna_meta <- merge(gse_mrna_df,gse_meta_df,
                       by.x="sample",by.y="sample")
gse_mrna_meta[1:5,1:5]

# 加载包
library(survival)
library(survminer)
gse_mrna_meta$PYGL_group <- ifelse(gse_mrna_meta$PYGL > mean(gse_mrna_meta$PYGL),"PYGL_high","PYGL_low")
table(gse_mrna_meta$PYGL_group)
# PYGL_high  PYGL_low 
#        53        53

# 使用survfit()函数拟合KM生存曲线
fit <- survfit(Surv(OS_month, status) ~ PYGL_group, data = gse_mrna_meta)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
ggsurvplot(fit, data = gse_mrna_meta,
           surv.median.line = "hv", # 添加中位数生存时间线
           # Change legends: title & labels
           legend.title = "Group", # 设置图例标题
           legend = "right",
           legend.labs = c("PYGL high", "PYGL low"), # 指定图例分组标签
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
           palette = c("#EB4232","#2DB2EB"), # 设置颜色画板
           #palette = "lancet", 
           ggtheme = theme_boxplot # Change ggplot2 theme
)

###################################################################################################
##### GSE42127 paper LUAD data #####
###################################################################################################
gse_mrna <- read.table("GEO/GSE42127/GSE42127_expression.txt",header = T,sep = "\t",row.names = 1,check.names = F)
gse_mrna[1:5,1:5]
#        GSM1032881 GSM1032882 GSM1032883 GSM1032884 GSM1032885
# 15E1.2      5.930      4.930      4.920       6.43      4.660
# 2PDE        7.870      6.710      7.500       7.42      7.800
# 76P         6.980      7.240      6.130       6.44      8.420
# 7A5         3.900      4.430      4.680       3.63      3.600
# A1BG        4.205      4.295      3.805       3.96      3.895
dim(gse_mrna)
#[1] 25438   176

gse_meta <- read.table("GEO/GSE42127/GSE42127_clinical.txt",header = T,sep = "\t",check.names = F)
head(gse_meta)
#    accession                          title platform                     source     organism had_adjuvant_chemo overall survival months
# 1 GSM1032881 human lung cancer sample 250-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                    33.6
# 2 GSM1032882 human lung cancer sample 259-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                    33.6
# 3 GSM1032883 human lung cancer sample 289-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                   120.0
# 4 GSM1032884 human lung cancer sample 296-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                   106.8
# 5 GSM1032885 human lung cancer sample 309-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                   132.0
dim(gse_meta)
#[1] 176  13
table(gse_meta$histology)
# Adenocarcinoma       Squamous 
# 133             43  
gse_meta <- gse_meta[gse_meta$histology == "Adenocarcinoma",]
dim(gse_meta)
#[1] 133  13
gse_mrna <- gse_mrna[,gse_meta$accession]
dim(gse_mrna)
#[1] 25438   133

gse_meta_df <- data.frame(sample=gse_meta$accession,age=gse_meta$`age at surgery`,sex=gse_meta$gender,
                          stage=paste0("Stage ",toupper(gse_meta$final.pat.stage)),
                          status=gse_meta$`survival status`,OS_month=gse_meta$`overall survival months`)
gse_meta_df$OS <- gse_meta_df$OS_month * 30
gse_meta_df$status <- ifelse(gse_meta_df$status == "D",1,0)
head(gse_meta_df)
#       sample  age sex      stage status OS_month
# 1 GSM1032881 63.8   M   Stage IB      1     33.6
# 2 GSM1032882 60.4   F   Stage IB      1     33.6
# 3 GSM1032884 67.6   M Stage IIIB      1    106.8
gse_mrna_df <- as.data.frame(t(gse_mrna))
gse_mrna_df$sample <- rownames(gse_mrna_df)
gse_mrna_df[1:5,1:5]
#             13CDNA73      15E1_2      182_FIP      3_HEXO     384D8_2
# GSM302996 0.05380707 -0.15709413  0.003104384 -0.04061761 -0.03437657
# GSM302997 0.13583142 -0.13950812 -0.083948146 -0.10526053 -0.06070256
# GSM302998 0.10346561 -0.15080859  0.002672228 -0.18103291 -0.07984835
# GSM302999 0.15128748 -0.18488011  0.035767355 -0.20748460  0.02976743

gse_mrna_meta <- merge(gse_mrna_df,gse_meta_df,
                       by.x="sample",by.y="sample")
gse_mrna_meta[1:5,1:5]

# 加载包
library(survival)
library(survminer)

gse_mrna_meta$PYGL_group <- ifelse(gse_mrna_meta$PYGL > median(gse_mrna_meta$PYGL),"PYGL_high","PYGL_low")
table(gse_mrna_meta$PYGL_group)
# PYGL_high  PYGL_low 
#        66        67

##### 生存分析 #####
# 使用survfit()函数拟合KM生存曲线
fit <- survfit(Surv(OS_month, status) ~ PYGL_group, data = gse_mrna_meta)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
ggsurvplot(fit, data = gse_mrna_meta,
           surv.median.line = "hv", # 添加中位数生存时间线
           # Change legends: title & labels
           legend.title = "Group", # 设置图例标题
           legend = "right",
           legend.labs = c("PYGL high", "PYGL low"), # 指定图例分组标签
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
           palette = c("#EB4232","#2DB2EB"), # 设置颜色画板
           #palette = "lancet", 
           ggtheme = theme_boxplot # Change ggplot2 theme
)

##### GSVA基因集变异分析 #####
library(GSVA)
###################################################################################################
##### GSE42127 paper LUAD data #####
###################################################################################################
gse_mrna <- read.table("GEO/GSE42127/GSE42127_expression.txt",header = T,sep = "\t",row.names = 1,check.names = F)
gse_mrna[1:5,1:5]
#        GSM1032881 GSM1032882 GSM1032883 GSM1032884 GSM1032885
# 15E1.2      5.930      4.930      4.920       6.43      4.660
# 2PDE        7.870      6.710      7.500       7.42      7.800
# 76P         6.980      7.240      6.130       6.44      8.420
# 7A5         3.900      4.430      4.680       3.63      3.600
# A1BG        4.205      4.295      3.805       3.96      3.895
dim(gse_mrna)
#[1] 25438   176

gse_meta <- read.table("GEO/GSE42127/GSE42127_clinical.txt",header = T,sep = "\t",check.names = F)
head(gse_meta)
#    accession                          title platform                     source     organism had_adjuvant_chemo overall survival months
# 1 GSM1032881 human lung cancer sample 250-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                    33.6
# 2 GSM1032882 human lung cancer sample 259-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                    33.6
# 3 GSM1032883 human lung cancer sample 289-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                   120.0
# 4 GSM1032884 human lung cancer sample 296-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                   106.8
# 5 GSM1032885 human lung cancer sample 309-T  GPL6884 Non-Small-Cell Lung Cancer Homo sapiens              FALSE                   132.0
dim(gse_meta)
#[1] 176  13
table(gse_meta$histology)
# Adenocarcinoma       Squamous 
# 133             43  
gse_meta <- gse_meta[gse_meta$histology == "Adenocarcinoma",]
dim(gse_meta)
#[1] 133  13
gse_mrna <- gse_mrna[,gse_meta$accession]
dim(gse_mrna)
#[1] 25438   133

gse_meta_df <- data.frame(sample=gse_meta$accession,age=gse_meta$`age at surgery`,sex=gse_meta$gender,
                          stage=paste0("Stage ",toupper(gse_meta$final.pat.stage)),
                          status=gse_meta$`survival status`,OS_month=gse_meta$`overall survival months`)
gse_meta_df$OS <- gse_meta_df$OS_month * 30
gse_meta_df$status <- ifelse(gse_meta_df$status == "D",1,0)
head(gse_meta_df)
#       sample  age sex      stage status OS_month
# 1 GSM1032881 63.8   M   Stage IB      1     33.6
# 2 GSM1032882 60.4   F   Stage IB      1     33.6
# 3 GSM1032884 67.6   M Stage IIIB      1    106.8
gse_mrna_df <- as.data.frame(t(gse_mrna))
gse_mrna_df$sample <- rownames(gse_mrna_df)
gse_mrna_df[1:5,1:5]
#       sample 15E1.2 2PDE  76P  7A5
# 1 GSM1032881   5.93 7.87 6.98 3.90
# 2 GSM1032882   4.93 6.71 7.24 4.43
# 3 GSM1032884   6.43 7.42 6.44 3.63
# 4 GSM1032885   4.66 7.80 8.42 3.60
# 5 GSM1032886   4.88 7.52 6.19 4.11

gse_mrna_meta <- merge(gse_mrna_df,gse_meta_df,
                       by.x="sample",by.y="sample")
gse_mrna_meta[1:5,1:5]
#       sample 15E1.2 2PDE  76P  7A5
# 1 GSM1032881   5.93 7.87 6.98 3.90
# 2 GSM1032882   4.93 6.71 7.24 4.43
# 3 GSM1032884   6.43 7.42 6.44 3.63
# 4 GSM1032885   4.66 7.80 8.42 3.60
# 5 GSM1032886   4.88 7.52 6.19 4.11

data_tumor <- t(gse_mrna_meta[,-c(1,25440:25445)])
colnames(data_tumor) <- gse_mrna_meta$sample
data_tumor[1:5,1:5]
#        GSM1032881 GSM1032882 GSM1032884 GSM1032885 GSM1032886
# 15E1.2      5.930      4.930       6.43      4.660      4.880
# 2PDE        7.870      6.710       7.42      7.800      7.520
# 76P         6.980      7.240       6.44      8.420      6.190
# 7A5         3.900      4.430       3.63      3.600      4.110
# A1BG        4.205      4.295       3.96      3.895      4.255
dim(data_tumor)
#[1] 25438   133

exp <- as.matrix(log2(data_tumor + 1))
exp[1:5,1:5]
#        GSM1032881 GSM1032882 GSM1032884 GSM1032885 GSM1032886
# 15E1.2   2.792855   2.568032   2.893362   2.500802   2.555816
# 2PDE     3.148934   2.946731   3.073820   3.137504   3.090853
# 76P      2.996389   3.042644   2.895303   3.235727   2.845992
# 7A5      2.292782   2.440952   2.211012   2.201634   2.353323
# A1BG     2.379898   2.404631   2.310340   2.291309   2.393691
ferro <- read.gmt("GSVA/ferroptosis_genesets.gmt")

table(ferro$gene %in% rownames(exp))
#FALSE  TRUE 
#72   436 
ferro <- ferro[ferro$gene %in% rownames(exp),]

# 提取基因集
geneset <- split(ferro$gene, ferro$term)    

# 进行GSVA分析
es.max <- gsva(exp, geneset, 
               method="gsva",kcdf="Gaussian",
               parallel.sz=5)#线程数
# 写出GSVA结果到文件
write.table(es.max, file="GSVA/GSVA_GSE42127_ferroptosis_result.txt",sep="\t",
            quote=F,row.names = T)

data_df <- cbind(gse_mrna_meta,t(es.max))
dim(data_df)
#[1]   541 19963
data_df[1:5,19961:19963]
# PSKH1 PSKH2    PSMA1
# GSM1032881  3.76  3.37 7.956667
# GSM1032882  4.15  2.79 7.610000
# GSM1032884  4.08  3.98 7.666667
# GSM1032885  4.08  4.04 7.226667

# 添加PYGL分组
data_df$PYGL_group <- ifelse(data_df$PYGL > median(data_df$PYGL),"PYGL_high","PYGL_low")
table(data_df$PYGL_group)
# PYGL_high  PYGL_low 
#     270       271
data_df$PYGL_group <- factor(data_df$PYGL_group,levels = c("PYGL_low","PYGL_high"))


data_df$driver_group <- ifelse(data_df$ferroptosis_driver > median(data_df$ferroptosis_driver),"FerroDriver_high","FerroDriver_low")
table(data_df$driver_group)
# FerroDriver_high  FerroDriver_low 
# 270              271
data_df$suppressor_group <- ifelse(data_df$ferroptosis_suppressor > median(data_df$ferroptosis_suppressor),"FerroSuppressor_high","FerroSuppressor_low")
table(data_df$suppressor_group)
# FerroSuppressor_high  FerroSuppressor_low 
# 270                  271 

# PYGL表达与铁死亡的相关性
pdf("GSVA/Boxplot_GSE42127_wilcoxon_test_ferroptosis_driverr.pdf",height=4.5,width=4.5)
ggplot(data_df,aes(PYGL_group,ferroptosis_driver,fill=PYGL_group)) + 
  #  geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(aes(fill = PYGL_group),width=0.5,notch = T,outlier.size = 0.5) +
  geom_jitter(size=0.5,position = position_jitter(0.2)) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) + 
  theme_boxplot + stat_compare_means() + 
  theme(plot.title = element_text(size=15,hjust=0.5))
dev.off()
pdf("GSVA/Boxplot_GSE42127_wilcoxon_test_ferroptosis_suppressor.pdf",height=4.5,width=4.5)
ggplot(data_df,aes(PYGL_group,ferroptosis_suppressor,fill=PYGL_group)) + 
  #  geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(aes(fill = PYGL_group),width=0.5,notch = T,outlier.size = 0.5) +
  geom_jitter(size=0.5,position = position_jitter(0.2)) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) + 
  theme_boxplot + stat_compare_means() + 
  theme(plot.title = element_text(size=15,hjust=0.5))
dev.off()

table(data_df$PYGL_group,data_df$driver_group)
#           FerroDriver_high FerroDriver_low
# PYGL_high              152             118
# PYGL_low               118             153
table(data_df$PYGL_group,data_df$suppressor_group)
#           FerroSuppressor_high FerroSuppressor_low
# PYGL_high                  154                 116
# PYGL_low                   116                 155

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

# 加载包
library(survival)
library(survminer)

fit <- survfit(Surv(OS_month, status) ~ driver_type, data = data_df)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
pdf("GSVA/Survival_plot_GSE42127_ferroptosis_driver_type.pdf",height=6,width=8.8,onefile = F)
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
           ggtheme = theme_boxplot # Change ggplot2 theme
)
dev.off()
fit <- survfit(Surv(OS_month, status) ~ suppressor_type, data = data_df)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
pdf("GSVA/Survival_plot_GSE42127_ferroptosis_suppressor_type.pdf",height=6,width=9.6,onefile = F)
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
           ggtheme = theme_boxplot # Change ggplot2 theme
)
dev.off()

# 提取铁死亡Driver high组的样本
data_df2 <- data_df[data_df$driver_group == "FerroDriver_high",]
dim(data_df2)
#[1]   270 19967
table(data_df2$PYGL_group)
# PYGL_high  PYGL_low 
# 152       118

pdf("GSVA/Boxplot_GSE42127_wilcoxon_test_FerroDriver_high_ferroptosis_driver.pdf",height=4.5,width=4.5)
ggplot(data_df2,aes(PYGL_group,ferroptosis_driver,fill=PYGL_group)) + 
  #  geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(aes(fill = PYGL_group),width=0.5,notch = T,outlier.size = 0.5) +
  geom_jitter(size=0.5,position = position_jitter(0.2)) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) + 
  theme_boxplot + stat_compare_means() + 
  theme(plot.title = element_text(size=15,hjust=0.5))
dev.off()
pdf("GSVA/Boxplot_GSE42127_wilcoxon_test_FerroDriver_high_ferroptosis_suppressor.pdf",height=4.5,width=4.5)
ggplot(data_df2,aes(PYGL_group,ferroptosis_suppressor,fill=PYGL_group)) + 
  #  geom_violin(aes(fill = group), width=0.6,trim = F) +
  geom_boxplot(aes(fill = PYGL_group),width=0.5,notch = T,outlier.size = 0.5) +
  geom_jitter(size=0.5,position = position_jitter(0.2)) + 
  scale_fill_manual(values = c("#2DB2EB","#EB4232")) + 
  theme_boxplot + stat_compare_means() + 
  theme(plot.title = element_text(size=15,hjust=0.5))
dev.off()

### 基因相关性分析
genes <- c('PYGL', 'MKI67', 'TOP2A', 'MCM5', 'CDC6', 'CCNB1', 'BID', 'BAX', 'BAD', 'BAK1', 'BCL2L11', 'BCL2')
data_cor <- data_df2[,genes]
head(data_cor)
pygl_cor <- cor(data_cor,method = "pearson")
head(pygl_cor)

# 绘制相关性热图
library(corrplot)
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue"))
pdf("Corrplot_correlation_FerroDriver_high_PYGL.pdf")
corrplot(pygl_cor,type = "lower",col = rev(col1(100)),diag = F)
dev.off()


# 加载包
library(survival)
library(survminer)

# 使用survfit()函数拟合KM生存曲线
fit <- survfit(Surv(OS_month, status) ~ PYGL_group, data = data_df2)
fit
# 使用ggsurvplot()函数绘制基础KM生存曲线
pdf("GSVA/Survival_plot_GSE42127_FerroDriver_high_PYGL.pdf",height=6,width=6,onefile = F)
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
           ggtheme = theme_boxplot # Change ggplot2 theme
)
dev.off()
