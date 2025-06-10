
###### load required packages ######
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)
library(tidyr)
library(paletteer)
library(reshape2)


##############################################################################################################################
# 13 paired Diagnosis and R/R dataset in this study
##############################################################################################################################

library(limma)
library(edgeR)

###### expression matrix and group info ######

# expression matrix and group info

expression_data <- read.table("RNAseq_count_analysis.norm.txt",header=T)
row.names(expression_data) = expression_data$Gene
expression_data <- expression_data[,-1]

group_data <- read.table("patient.info.txt",header=T,fill = TRUE,sep = "\t")



###### differential analysis, CSNK1E high vs low group ######

# add CSNK1E high vs low group
CSNK1E_data <- expression_data[rownames(expression_data) == "CSNK1E", ]

transposed_data <- as.data.frame(t(CSNK1E_data))
colnames(transposed_data) <- rownames(CSNK1E_data)
rownames(transposed_data) <- colnames(expression_data)

# merge into the table with patient table
transposed_data <- data.frame(Seq_ID = rownames(transposed_data), transposed_data)
group_data <- merge(group_data, transposed_data, by.x = "Seq_ID", by.y = "Seq_ID")

# CSNK1E groups, high vs low by median
group_data$CSNK1E_group <- ifelse(group_data$CSNK1E >= quantile(group_data$CSNK1E,0.5), "High", "Low")
group_data$CSNK1E_group = factor(group_data$CSNK1E_group, levels=c("Low", "High"))

samples <- colnames(expression_data)
group <- group_data_rchop$CSNK1E_group[match(samples, group_data_rchop$Seq_ID)]
group <- factor(group)

# design for limma
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# limma differential expression analysis
fit <- lmFit(expression_data_rchop, design)

contrast.matrix <- makeContrasts(High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# limma results
results_CK1E <- topTable(fit2, adjust = "fdr", number = Inf)
head(results_CK1E)


##### volcano plot

# cut-off for significant differential expression
results_CK1E$logQ <- -log10(results_CK1E$adj.P.Val)
results_CK1E$Significance <- "Not Significant"
results_CK1E$Significance[results_CK1E$logFC > log2(1.5) & results_CK1E$adj.P.Val < 0.01] <- "Upregulated"
results_CK1E$Significance[results_CK1E$logFC < -log2(1.5) & results_CK1E$adj.P.Val < 0.01] <- "Downregulated"

results_CK1E$Label <- ifelse(rownames(results_CK1E) %in% c("TNFRSF13B"), rownames(results_CK1E), NA)

p <- ggplot(results_CK1E, aes(x = logFC, y = logQ)) +
  geom_point(aes(color = Significance), alpha = 0.6, show.legend = FALSE) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw(base_size=18) +
  theme(axis.text=element_text(colour="black"), aspect.ratio = 1) +
  labs(title = "CSNK1E high vs low", x = "Log2 Fold Change", y = "-Log10 Adj.P")+ylim(0,15)+xlim(-2.2,2.2)

p + geom_text_repel(aes(label = Label), min.segment.length = 0.1, seed = 42, max.overlaps = 10, box.padding = 0.5, point.padding = 0.5, segment.color = 'black')



###### GSEA analysis for cluster signatures ######

library(GSVA)
library(GSEABase)

# gene sets
gset <- getGmt("gene_sets_bulk_RNAseq.gmt")
gset <- subsetGeneSets(gset, rownames(expression_data)) 

# gene sets enrichment analysis
ssgseaPar <- ssgseaParam(as.matrix(expression_data), gset)
gset_ssgsea <- gsva(ssgseaPar, verbose=FALSE)
gset_ssgsea = as.data.frame(gset_ssgsea)


# add GSEA results into patient table

transposed_data <- as.data.frame(t(gset_ssgsea))
colnames(transposed_data) <- rownames(gset_ssgsea)
rownames(transposed_data) <- colnames(gset_ssgsea)

# add gene expression into the patient info table
transposed_data <- data.frame(Seq_ID = rownames(transposed_data), transposed_data)
merged_data <- merge(group_data, transposed_data, by.x = "Seq_ID", by.y = "Seq_ID")



###### survival analysis ######


##### CSNK1E
merged_data$CSNK1E_group <- ifelse(merged_data$CSNK1E >= quantile(merged_data$CSNK1E,0.5), "High", "Low")
merged_data$CSNK1E_group = factor(merged_data$CSNK1E_group, levels=c("Low", "High"))

# OS
surv_object <- Surv(time = (merged_data$OS)/12, event = merged_data$OSS)
fit <- survfit(surv_object ~ CSNK1E_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by CSNK1E Expression",
           xlab = "Time (Months)", ylab = "OS (%)",
           legend.title = "CSNK1E",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)

## PFS
surv_object <- Surv(time = (merged_data$PFS)/12, event = merged_data$PFSS)
fit <- survfit(surv_object ~ CSNK1E_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by CSNK1E Expression",
           xlab = "Time (Months)", ylab = "PFS (%)",
           legend.title = "CSNK1E",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### TNFRSF13B

# adding group information
TNFRSF13B_data <- expression_data[rownames(expression_data) == "TNFRSF13B", ]
transposed_data <- as.data.frame(t(TNFRSF13B_data))
colnames(transposed_data) <- rownames(TNFRSF13B_data)
rownames(transposed_data) <- colnames(expression_data)

# merge into the table with patient table
transposed_data <- data.frame(Seq_ID = rownames(transposed_data), transposed_data)
group_data <- merge(group_data, transposed_data, by.x = "Seq_ID", by.y = "Seq_ID")

# TNFRSF13B group
merged_data$TNFRSF13B_group <- ifelse(merged_data$TNFRSF13B >= quantile(merged_data$TNFRSF13B,0.5), "High", "Low")
merged_data$TNFRSF13B_group = factor(merged_data$TNFRSF13B_group, levels=c("Low", "High"))


# OS
surv_object <- Surv(time = (merged_data$OS)/12, event = merged_data$OSS)
fit <- survfit(surv_object ~ TNFRSF13B_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by TNFRSF13B Expression",
           xlab = "Time (Months)", ylab = "OS (%)",
           legend.title = "TNFRSF13B",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


## PFS
surv_object <- Surv(time = (merged_data$PFS)/12, event = merged_data$PFSS)
fit <- survfit(surv_object ~ TNFRSF13B_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by TNFRSF13B Expression",
           xlab = "Time (Months)", ylab = "PFS (%)",
           legend.title = "TNFRSF13B",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### Cluster0
merged_data$scRNA_Cluster0_group <- ifelse(merged_data$scRNA_Cluster0 >= quantile(merged_data$scRNA_Cluster0,0.5), "High", "Low")
merged_data$scRNA_Cluster0_group = factor(merged_data$scRNA_Cluster0_group, levels=c("Low", "High"))

# OS
surv_object <- Surv(time = (merged_data$OS)/12, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster0_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster0",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "scRNA_Cluster0",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.2), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)

# PFS
surv_object <- Surv(time = (merged_data$PFS)/12, event = merged_data$PFSS)
fit <- survfit(surv_object ~ scRNA_Cluster0_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster0",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "scRNA_Cluster0",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster3
merged_data$scRNA_Cluster3_group <- ifelse(merged_data$scRNA_Cluster3 >= quantile(merged_data$scRNA_Cluster3,0.5), "High", "Low")
merged_data$scRNA_Cluster3_group = factor(merged_data$scRNA_Cluster3_group, levels=c("Low", "High"))

# OS
surv_object <- Surv(time = (merged_data$OS)/12, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster3_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster3",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "scRNA_Cluster3",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)

# PFS
surv_object <- Surv(time = (merged_data$PFS)/12, event = merged_data$PFSS)
fit <- survfit(surv_object ~ scRNA_Cluster3_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster3",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "scRNA_Cluster3",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster5
merged_data$scRNA_Cluster5_group <- ifelse(merged_data$scRNA_Cluster5 >= quantile(merged_data$scRNA_Cluster5,0.5), "High", "Low")
merged_data$scRNA_Cluster5_group = factor(merged_data$scRNA_Cluster5_group, levels=c("Low", "High"))

# OS
surv_object <- Surv(time = (merged_data$OS)/12, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster5_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster5",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "scRNA_Cluster5",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)

# PFS
surv_object <- Surv(time = (merged_data$PFS)/12, event = merged_data$PFSS)
fit <- survfit(surv_object ~ scRNA_Cluster5_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster5",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "scRNA_Cluster5",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster6
merged_data$scRNA_Cluster6_group <- ifelse(merged_data$scRNA_Cluster6 >= quantile(merged_data$scRNA_Cluster6,0.5), "High", "Low")
merged_data$scRNA_Cluster6_group = factor(merged_data$scRNA_Cluster6_group, levels=c("Low", "High"))

# OS
surv_object <- Surv(time = (merged_data$OS)/12, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster6_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster6",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "scRNA_Cluster6",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)

# PFS
surv_object <- Surv(time = (merged_data$PFS)/12, event = merged_data$PFSS)
fit <- survfit(surv_object ~ scRNA_Cluster6_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster6",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "scRNA_Cluster6",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster0+6, signature for cluster 0+6 as a combined subset
merged_data$scRNA_Cluster06_group <- ifelse(merged_data$scRNA_Cluster06 >= quantile(merged_data$scRNA_Cluster06,0.5), "High", "Low")
merged_data$scRNA_Cluster06_group = factor(merged_data$scRNA_Cluster06_group, levels=c("Low", "High"))

# OS

surv_object <- Surv(time = (merged_data$OS)/12, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster06_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster06",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "scRNA_Cluster06",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)

# PFS
surv_object <- Surv(time = (merged_data$PFS)/12, event = merged_data$PFSS)
fit <- survfit(surv_object ~ scRNA_Cluster06_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "K-M Curve by scRNA_Cluster06",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "scRNA_Cluster06",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)





##############################################################################################################################
# DLBCL cohort in Lacy et al. 2020 Blood
##############################################################################################################################


###### expression matrix ######

# reading the expression matrix from GEO, with acc# GSE181063
expression_data = read.table("GSE181063_series_matrix.nor.txt",header=T,sep = "\t")
rownames(expression_data) <- expression_data$ID_REF
expression_data <- expression_data[,-1]

# patient information in the series.matrix file
group_data <- read.table("GSE181063_series_matrix.sample.info.txt",header=T,fill = TRUE,sep = "\t")


###### survival analysis of CSNK1E and TNFRSF13B ######

# CSNK1E: ILMN_1808913, ILMN_1708858, ILMN_2415235
# TNFRSF13B: ILMN_1759075

# Extract gene expression
genelist = c("ILMN_1808913","ILMN_1708858","ILMN_2415235","ILMN_1759075","ILMN_1763398")
gene_data <- expression_data[rownames(expression_data) %in% genelist, ]

# merge with patient table
transposed_data <- as.data.frame(t(gene_data))
colnames(transposed_data) <- rownames(gene_data)
rownames(transposed_data) <- colnames(expression_data)

# use CSNK1E average for analysis
transposed_data$CSNK1E_avg = rowMeans(transposed_data[, c("ILMN_1808913", "ILMN_1708858","ILMN_2415235")], na.rm = TRUE)
transposed_data <- data.frame(ACC = rownames(transposed_data), transposed_data)

# merge tables, and only include CHOP/RCHOP
merged_data <- merge(group_data, transposed_data, by.x = "GSE_acc", by.y = "ACC")
merged_data = subset(merged_data, disease == "DLBCL" & FL %in% c("CHOP-R","CHOP"))


##### CSNK1E survival analysis
merged_data$CSNK1Eavg_group <- ifelse(merged_data$CSNK1E_avg >= quantile(merged_data$CSNK1E_avg,0.5), "High", "Low")
merged_data$CSNK1Eavg_group <- factor(merged_data$CSNK1Eavg_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$OS_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ CSNK1Eavg_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "KM Curve by CSNK1E Expression",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "CSNK1E",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)

##### TNFRSF13B survival analysis
merged_data$TNFRSF13B_group <- ifelse(merged_data$ILMN_1759075 >= quantile(merged_data$ILMN_1759075,0.5), "High", "Low")
merged_data$TNFRSF13B_group <- factor(merged_data$TNFRSF13B_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$OS_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ TNFRSF13B_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "KM Curve by TNFRSF13B Expression",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "TNFRSF13B",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)



###### annotation wit gene SYMBOL for GSEA ######

library("AnnotationDbi")
library("illuminaHumanv4.db")
library("org.Hs.eg.db")
columns(illuminaHumanv4.db)
library(tibble)

id_str <- rownames(expression_data)

annotable = AnnotationDbi::select(illuminaHumanv4.db,id_str,c("SYMBOL"), keytype="PROBEID")

# add gene symbol to the matrix
expression_data_df <- expression_data %>%  tibble::rownames_to_column(var = "PROBEID")
expression_data_df <- expression_data_df %>% left_join(annotable, by = "PROBEID")

df_final <- expression_data_df %>% filter(!is.na(SYMBOL))

# collapse probe sets to gene symbol with mean values
df_final <- df_final %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), 
   ~ if(n() > 1) mean(.x, na.rm = TRUE) else .x[1])) %>% ungroup()

df_final <- df_final %>% column_to_rownames(var = "SYMBOL")



###### gene set analysis ######

library(GSVA)
library(GSEABase)

# gene sets
gset <- getGmt("gene_sets_bulk_RNAseq.gmt")
gset <- subsetGeneSets(gset, rownames(df_final)) 

# gene sets enrichment analysis
ssgseaPar <- ssgseaParam(as.matrix(df_final), gset)
gset_ssgsea <- gsva(ssgseaPar, verbose=FALSE)
gset_ssgsea = as.data.frame(gset_ssgsea)

# merge GSEA results with patient table
transposed_data <- as.data.frame(t(gset_ssgsea))
transposed_data <- data.frame(ACC = rownames(transposed_data), transposed_data)
merged_data <- merge(group_data, transposed_data, by.x = "GSE_acc", by.y = "ACC")
merged_data = subset(merged_data, disease == "DLBCL" & FL %in% c("CHOP-R","CHOP"))



###### survival analysis ######

##### scRNA_Cluster0
merged_data$scRNA_Cluster0_group <- ifelse(merged_data$scRNA_Cluster0 >= quantile(merged_data$scRNA_Cluster0,0.5), "High", "Low")
merged_data$scRNA_Cluster0_group <- factor(merged_data$scRNA_Cluster0_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$OS_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster0_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_0",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_0",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster3
merged_data$scRNA_Cluster3_group <- ifelse(merged_data$scRNA_Cluster3 >= quantile(merged_data$scRNA_Cluster3,0.5), "High", "Low")
merged_data$scRNA_Cluster3_group <- factor(merged_data$scRNA_Cluster3_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$OS_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster3_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_3",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_3",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster5
merged_data$scRNA_Cluster5_group <- ifelse(merged_data$scRNA_Cluster5 >= quantile(merged_data$scRNA_Cluster5,0.5), "High", "Low")
merged_data$scRNA_Cluster5_group <- factor(merged_data$scRNA_Cluster5_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$OS_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster5_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_5",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_5",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster6
merged_data$scRNA_Cluster6_group <- ifelse(merged_data$scRNA_Cluster6 >= quantile(merged_data$scRNA_Cluster6,0.5), "High", "Low")
merged_data$scRNA_Cluster6_group <- factor(merged_data$scRNA_Cluster6_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$OS_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster6_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_6",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_6",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster0+6, signature for cluster 0+6 as a combined subset
merged_data$scRNA_Cluster06_group <- ifelse(merged_data$scRNA_Cluster06 >= quantile(merged_data$scRNA_Cluster06,0.5), "High", "Low")
merged_data$scRNA_Cluster06_group <- factor(merged_data$scRNA_Cluster06_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$OS_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster06_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_06",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_06",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)





##############################################################################################################################
# DLBCL cohort in Lenz et al. 2008 NEJM
##############################################################################################################################

###### expression matrix ######

# reading the expression matrix from GEO, with acc# GSE10846
expression_data = read.table("GSE10846_series_matrix.nor.txt",header=T,sep = "\t")
rownames(expression_data) <- expression_data$ID_REF
expression_data <- expression_data[,-1]

# patient information in the series.matrix file
group_data <- read.table("GSE10846_series_matrix.patient.group.txt",header=T,fill = TRUE,sep = "\t")



###### survival analysis of CSNK1E and TNFRSF13B ######

# CSNK1E: 202332_at, 222015_at
# TNFRSF13B: 207641_at

# Extract gene expression
genelist = c("202332_at","222015_at","207641_at")
gene_data <- expression_data[rownames(expression_data) %in% genelist, ]

# merge with patient table
transposed_data <- as.data.frame(t(gene_data))
colnames(transposed_data) <- rownames(gene_data)
rownames(transposed_data) <- colnames(expression_data)

# use CSNK1E average for analysis
transposed_data$CSNK1E_avg = rowMeans(transposed_data[, c("202332_at", "222015_at")], na.rm = TRUE)
transposed_data <- data.frame(ACC = rownames(transposed_data), transposed_data)

# merge tables, and only include CHOP/RCHOP
merged_data <- merge(group_data, transposed_data, by.x = "GSE_acc", by.y = "ACC")
merged_data = subset(merged_data, treatment %in% c("CHOP","R-CHOP"))
merged_data = subset(merged_data, OSS != "NA")

##### CSNK1E survival analysis
merged_data$CSNK1Eavg_group <- ifelse(merged_data$CSNK1E_avg >= quantile(merged_data$CSNK1E_avg,0.5), "High", "Low")
merged_data$CSNK1Eavg_group <- factor(merged_data$CSNK1Eavg_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$FU_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ CSNK1Eavg_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "KM Curve by CSNK1E Expression",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "CSNK1E",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400")
           )
print(p, newpage = FALSE)


##### TNFRSF13B survival analysis
merged_data$TNFRSF13B_group <- ifelse(merged_data$TNFRSF13B >= quantile(merged_data$TNFRSF13B,0.5), "High", "Low")
merged_data$TNFRSF13B_group <- factor(merged_data$TNFRSF13B_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$FU_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ TNFRSF13B_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "KM Curve by TNFRSF13B Expression",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "TNFRSF13B",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400")
           )
print(p, newpage = FALSE)



###### annotation wit gene SYMBOL for GSEA ######

library("AnnotationDbi")
library("hgu133plus2.db")
library("org.Hs.eg.db")                                          
columns(org.Hs.eg.db)
library(tibble)


id_str <- rownames(expression_data)
annotable = AnnotationDbi::select(hgu133plus2.db,id_str,c("SYMBOL"), keytype="PROBEID")

# add gene symbol to the matrix
expression_data_df <- expression_data %>%  tibble::rownames_to_column(var = "PROBEID")
expression_data_df <- expression_data_df %>% left_join(annotable, by = "PROBEID")

# remove rows without annotation results
df_filtered <- expression_data_df %>% filter(!is.na(SYMBOL))

# collapse probe sets to gene symbol with mean values
df_final <- df_filtered %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), 
   ~ if(n() > 1) mean(.x, na.rm = TRUE) else .x[1])) %>% ungroup()

df_final <- df_final %>% column_to_rownames(var = "SYMBOL")



###### gene set analysis ######

gset <- getGmt("gene_sets_bulk_RNAseq.gmt")
gset <- subsetGeneSets(gset, rownames(df_final)) 

ssgseaPar <- ssgseaParam(as.matrix(df_final), gset, )
gset_ssgsea <- gsva(ssgseaPar, verbose=FALSE)
gset_ssgsea = as.data.frame(gset_ssgsea)

# merge with table containing patient information
transposed_data <- as.data.frame(t(gset_ssgsea))
transposed_data <- data.frame(ACC = rownames(transposed_data), transposed_data)
merged_data <- merge(group_data, transposed_data, by.x = "GSE_acc", by.y = "ACC")
merged_data = subset(merged_data, OSS != "NA")


###### survival analysis ######


##### scRNA_Cluster0
merged_data$scRNA_Cluster0_group <- ifelse(merged_data$scRNA_Cluster0 >= quantile(merged_data$scRNA_Cluster0,0.5), "High", "Low")
merged_data$scRNA_Cluster0_group <- factor(merged_data$scRNA_Cluster0_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$FU_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster0_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_0",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_0",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster3
merged_data$scRNA_Cluster3_group <- ifelse(merged_data$scRNA_Cluster3 >= quantile(merged_data$scRNA_Cluster3,0.5), "High", "Low")
merged_data$scRNA_Cluster3_group <- factor(merged_data$scRNA_Cluster3_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$FU_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster3_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_3",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_3",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster5
merged_data$scRNA_Cluster5_group <- ifelse(merged_data$scRNA_Cluster5 >= quantile(merged_data$scRNA_Cluster5,0.5), "High", "Low")
merged_data$scRNA_Cluster5_group <- factor(merged_data$scRNA_Cluster5_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$FU_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster5_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_5",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_5",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster6
merged_data$scRNA_Cluster6_group <- ifelse(merged_data$scRNA_Cluster6 >= quantile(merged_data$scRNA_Cluster6,0.5), "High", "Low")
merged_data$scRNA_Cluster6_group <- factor(merged_data$scRNA_Cluster6_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$FU_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster6_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_6",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_6",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


##### scRNA_Cluster0+6, signature for cluster 0+6 as a combined subset
merged_data$scRNA_Cluster6_group <- ifelse(merged_data$scRNA_Cluster06 >= quantile(merged_data$scRNA_Cluster06,0.5), "High", "Low")
merged_data$scRNA_Cluster6_group <- factor(merged_data$scRNA_Cluster6_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$FU_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ scRNA_Cluster6_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Curve by Cluster_06",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Cluster_06",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


