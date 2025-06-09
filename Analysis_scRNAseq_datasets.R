
##############################################################################################################################
# Diagnosis and paired R/R dataset in this study
##############################################################################################################################

##############################################################################################################################
###### cellranger analysis with fastq ######
ref="./ref/cellranger/refdata-gex-GRCh38-2024-A"
cellranger="./bin/cellranger-7.1.0/cellranger"
RUN_id="run_id"

echo cellranger multi ${RUN_id} at `date`

${cellranger} multi --mempercore 6 --id=${RUN_id} --csv=./multi.config.csv

echo cellranger multi ${RUN_id} at `date`



##############################################################################################################################
###### load required R packages ######

library(ggplot2)
library(ggridges)
library(ggpubr)
library(ggthemes)
library(cowplot)
library(dplyr)
library(plyr)
library(reshape2)
library(dittoSeq)
library(RColorBrewer)
library(grDevices)
library(viridis)
library(scDataviz)
library(GSEABase)
library(GGally)

library(scran)
library(Seurat)
library(SeuratWrappers)
library(scCustomize)
library(umap)
library(harmony)
library(SingleR)
library(celldex)

library(UCell)
library(AUCell)
library(CellChat)
library(monocle3)
library(Nebulosa)

###### color codes for cell clusters ######

# all clusters
color.lib <- c("#E31A1C", "#55c2fc", "#A6761D", "#F1E404", "#33A02C", "#1F78B4", 
               "#FB9A99", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#00E5DF",
               "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
               "#F4A11D", "#8DC8ED", "#4C6CB0", "#8A1C1B", "#CBCC2B", "#EA644C",
               "#634795", "#005B1D", "#26418A", "#CB8A93", "#B2DF8A", "#E22826",
               "#A6CEE3", "#F4D31D", "#F4A11D", "#82C800", "#8B5900", "#858ED1",
               "#FF72E1", "#CB50B2", "#007D9B", "#26418A", "#8B495F", "#FF394B")

# B cell clusters with non-malignant B
color.lib.bcels <- c("#E31A1C", "#F1E404", "#1F78B4", "#FB9A99", "#FF7F00", "#CAB2D6", "#00E5DF", "#1B9E77", "#F4A11D")

# B cell clusters with without non-malignant B
color.lib.mbcels <- c("#E31A1C", "#F1E404", "#1F78B4", "#FB9A99", "#FF7F00", "#CAB2D6", "#1B9E77", "#F4A11D")



###### load cellranger output files ######

## SD patients
pt01_dx = Read10X(data.dir="./run_id/outs/per_sample_outs/PS_PT01_DX/count/sample_filtered_feature_bc_matrix")
pt01_rr = Read10X(data.dir="./run_id/outs/per_sample_outs/PS_PT01_RR/count/sample_filtered_feature_bc_matrix")

pt02_dx = Read10X(data.dir="./run_id/outs/per_sample_outs/PS_PT02_DX/count/sample_filtered_feature_bc_matrix")
pt02_rr = Read10X(data.dir="./run_id/outs/per_sample_outs/PS_PT02_RR/count/sample_filtered_feature_bc_matrix")

## CR patients
pt03_dx = Read10X(data.dir="./run_id/outs/per_sample_outs/PS_PT03_DX/count/sample_filtered_feature_bc_matrix")
pt03_rr = Read10X(data.dir="./run_id/outs/per_sample_outs/PS_PT03_RR/count/sample_filtered_feature_bc_matrix")

pt04_dx = Read10X(data.dir="./run_id/outs/per_sample_outs/PS_PT04_DX/count/sample_filtered_feature_bc_matrix")
pt04_rr = Read10X(data.dir="./run_id/outs/per_sample_outs/PS_PT04_RR/count/sample_filtered_feature_bc_matrix")



###### create seurat objects ######

# create a count_list file with all sample
count_list <- list(
  sce_pt01_dx = pt01_dx,
  sce_pt01_rr = pt01_rr,  
  sce_pt02_dx = pt02_dx,
  sce_pt02_rr = pt02_rr,
  sce_pt03_dx = pt03_dx,
  sce_pt03_rr = pt03_rr,
  sce_pt04_dx = pt04_dx,
  sce_pt04_rr = pt04_rr
)

# create a list file with seurat objects
seurat_objects <- list()

seurat_objects <- lapply(names(count_list), function(sample_name) {
    CreateSeuratObject(
        counts = count_list[[sample_name]],
        project = sample_name,
        min.cells = 5,
        min.features = 500
    )
})

names(seurat_objects) <- names(count_list)
print(names(seurat_objects))
original_names <- names(seurat_objects)



###### data QC and filter ######

# MT percentage
seurat_objects <- lapply(seurat_objects, function(object) {
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  return(object)
})

## QC plots before filtering
plots <- lapply(names(seurat_objects), function(sample_name){
  object <- seurat_objects[[sample_name]]
  plot <- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(plot)
  ggsave(filename = paste0(sample_name, "_vlnplot.png"), plot = plot, width = 10, height = 6)
  return(plot)
})


# filter
seurat_objects <- lapply(seurat_objects, function(object) {
  object = subset(object, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA < 30000 & nCount_RNA > 1000)
  return(object)
})


## QC plots after filtering
plots <- lapply(names(seurat_objects), function(sample_name){
  object <- seurat_objects[[sample_name]]
  plot <- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(plot)
  ggsave(filename = paste0(sample_name, "_filtered_vlnplot.png"), plot = plot, width = 10, height = 6)
  return(plot)
})


# add meta.data
seurat_objects <- lapply(seurat_objects, function(object) {
  sample_name = object@project.name
  
  name_parts <- strsplit(sample_name, "_")[[1]]
  
  # patient ID and disease stages
  patient_id <- name_parts[2]  
  disease_stage <- name_parts[3]
  sample_id <- paste(name_parts[2], name_parts[3], sep="_")

  object$patient <- patient_id
  object$stage <- disease_stage
  object$sample <- sample_id

  return(object)
})



###### merge samples ######

# normalization
seurat_objects <- lapply(seurat_objects, function(object) {
  object = NormalizeData(object)
  return(object)
})

# add sample identity to cell ID
seurat_objects <- lapply(names(seurat_objects), function(i) {
  RenameCells(seurat_objects[[i]], add.cell.id = i)
})
names(seurat_objects) <- original_names
smerge <- Reduce(function(x, y) merge(x, y), seurat_objects)

print(smerge)
rm(seurat_objects)

# add groups of front-line response
smerge$response = ifelse(smerge$patient %in% c("pt01","pt02"),"SD","CR")
smerge$response = factor(smerge$response, levels=c("SD","CR"))



###### UMAP clustering ######

# clustering, with harmony correction
smerge = FindVariableFeatures(smerge_orig, selection.method = "vst")
smerge = ScaleData(smerge, verbose = FALSE)
smerge = RunPCA(smerge, verbose = FALSE, npcs = 30)
ElbowPlot(object = smerge, ndims = 40)          # ElbowPlot to decide the dim number

smerge = RunHarmony(smerge, "patient")
smerge = FindNeighbors(smerge, reduction = "harmony", dims = 1:18)
smerge = FindClusters(smerge,resolution =0.5)
smerge = RunUMAP(smerge, reduction = "harmony", dims = 1:18)

# UMAP plots by cluster and patients
DimPlot(smerge, reduction = "umap", pt.size = 0.1, group.by="seurat_clusters", label = TRUE, label.size = 4) + theme_few()
DimPlot(smerge, reduction = "umap", pt.size = 0.1, group.by="patient", label = TRUE, label.size = 4) + theme_few()



###### BCR-seq results ######

# using shell command line
# create a file with B cell clone types for all samples, and merge B cell clone
find ./run_id/outs/per_sample_outs/PS_*/vdj_b/filtered_contig_annotations.csv > filelist.filtered_contig_annotations.txt
less filelist.filtered_contig_annotations.txt | perl -e 'print"barcode\tclonotype\n"; while(<>){chomp; ($pt,$stage) =($1,$2) if(/PS\_(.*?)\_(.*?)/); $tag="sce_"."$pt"."_$stage"; open I,"$_"; $he=<I>; while(<I>){chomp; next if(/^barcode/); @a=split/\,/,$_; $a[0]="$tag"."_"."$a[0]"; $hash{$a[0]}=$a[-3];} $n++; close O; close I;} foreach(sort keys %hash){print "$_\t$hash{$_}\n";}' > vdj_b.clonotypes.allsamples.txt

# back to R
# load clonotypes
clonotype = read.table("./vdj_b.clonotypes.allsamples.txt",header=T)
rownames(clonotype) = clonotype$barcode

# group clonotypes, with five major clonotypes, and add into seurat object
clonotype$clonotype_br = ifelse(clonotype$clonotype %in% c("clonotype1","clonotype2","clonotype3","clonotype4","clonotype5"),clonotype$clonotype,"others")
smerge = AddMetaData(object =smerge, clonotype)

# dimplot with clonotypes
clone_color = c("#7CA0D4","#A48AD3","#E995EB","#BADE86","#2B8AAE","#624894")
DimPlot(smerge, reduction = "umap", pt.size = 0.1, group.by="clonotype_br", label = TRUE, label.size = 4, cols = clone_color) + theme_classic()

# clonotypes in each cluster
table(smerge@meta.data$RNA_snn_res.0.5, smerge@meta.data$clonotype_br)



###### SingleR annotation ######

# v5 is not supported, needs to be converted to v3 before SingleR
smerge[["RNA"]] <- as(smerge[["RNA"]], Class="Assay")

# main annotation
ref <- fetchReference("blueprint_encode", "2024-02-26")
pred.BlueprintEncodeData <- SingleR(test = smerge@assays$RNA@data, ref = ref, labels = ref$label.main)
smerge@meta.data$CellType.BlueprintEncodeData <- pred.BlueprintEncodeData$labels

# cell types across seurat_clusters
DimPlot(smerge, reduction = "umap", pt.size = 0.5, group.by="CellType_BlueprintEncodeData", label = TRUE, label.size = 4) + theme_few()
table(smerge$seurat_clusters, smerge$CellType_BlueprintEncodeData)

# add cell types information based on SingleR results

# cell type and cluster ID
Idents(smerge) = smerge$RNA_snn_res.0.5
smerge <- RenameIdents(object = smerge,
'0'='B_0',
'1'='CD4T_1',
'2'='CD8T_2',
'3'='B_3',
'4'='MONO_4',
'5'='B_5',
'6'='B_6',
'7'='CD8T_7',
'8'='B_8',
'9'='B_9',
'10'='NK_10',
'11'='B_11',
'12'='B_12',
'13'='CD8T_13',
'14'='ENDO_14',
'15'='DC_15',
'16'='FIB_16',
'17'='MAC_17',
'18'='B_18',
'19'='pDC_19'
)
smerge$cell_type_ID = Idents(smerge)

# broader cell types
Idents(smerge) = smerge$RNA_snn_res.0.5
smerge <- RenameIdents(object = smerge,
'0'='Malignant_B',
'1'='T_CD8',
'2'='T_CD4',
'3'='Malignant_B',
'4'='Monocytes',
'5'='Malignant_B',
'6'='Malignant_B',
'7'='T_CD8',
'8'='Malignant_B',
'9'='Malignant_B',
'10'='NK',
'11'='Malignant_B',
'12'='NonMalignant_B',
'13'='T_CD8',
'14'='Endothelial',
'15'='DC',
'16'='Fibroblasts',
'17'='Macrophages',
'18'='Malignant_B',
'19'='pDC'
)
smerge$broad_cell_type = Idents(smerge)
Idents(smerge) = smerge$RNA_snn_res.0.5

# UMAP by broader cell types, all samples
DimPlot(smerge, reduction = "umap", pt.size = 0.5, group.by="broad_cell_type", label = TRUE, label.size = 4) + theme_few()

# UMAP by broader cell types, for each samples
plots <- lapply(c("pt01","pt02","pt03","pt04"), function(pt_name){
  plot <- DimPlot(subset(smerge, patient == pt_name), reduction = "umap", pt.size = 0.01, split.by = "stage", group.by="major_cell_name", label = F, label.size = 4)
  print(plot)
  ggsave(filename = paste0(pt_name, ".umap.by-cell-type.pdf"), plot = plot, width = 6, height = 3.5)
  return(plot)
})



# plots of numbers and percentages for different cell types across samples
cell_type_num = as.data.frame(table(smerge$major_cell_name, smerge$sample))
colnames(cell_type_num) = c("Cell_type","Sample","Cell_Number")

# cell number
ggplot(cell_type_num, aes(x = Sample, y = Cell_Number, fill = Cell_type)) +
    geom_bar(stat = "identity") +
    labs(title = "Numbers by cell types", y = "Number of cells", x="") +
    theme_classic(base_size=18) +
    theme(axis.text=element_text(color="black"),axis.text.x = element_text(size = 18, angle=45, hjust=1),axis.text.y = element_text(size = 18))
# cell percentages
ggplot(cell_type_num, aes(x=Sample, y=Cell_Number, fill=Cell_type)) +
    geom_col(position="fill") +
    labs(title="Percentage by cell types", x="", y="Percentage (%)") +
    scale_y_continuous(labels=scales::percent) +
    theme_classic(base_size=18) +
    theme(axis.text=element_text(color="black"),axis.text.x = element_text(size = 18, angle=45, hjust=1),axis.text.y = element_text(size = 18))



###### Cytotrace analysis ######
library(CytoTRACE2)
smerge <- cytotrace2(smerge, is_seurat = TRUE, slot_type = "counts", species = 'human')



###### geneset analysis ######

counts <- LayerData(smerge, assay = "RNA", layer = "counts")
expressed_gene = rownames(counts)
rm(counts)

gset_cur <- getGmt("gene_sets_scRNAseq.gmt")
gset_cur <- subsetGeneSets(gset_cur, expressed_gene) 
smerge <- AddModuleScore_UCell(smerge, features=geneIds(gset_cur), name=NULL)



###### monocle3 analysis ######

DefaultAssay(smerge) <- "RNA"

# convert seurat to single cell object for monocle analysis
cds <- as.cell_data_set(smerge)
cds <- cluster_cells(cds, reduction_method = "UMAP")

# trajectory graph
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 5)

# order B cells clusters based on CytoTRACE2 and HSPC signature score
cds <- order_cells(cds)

# pseudotime UMAP
plot_cells(cds, color_cells_by = "pseudotime",  group_cells_by = "cluster", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, label_roots = T, trajectory_graph_color = "grey80", trajectory_graph_segment_size = 1,cell_size = 0.1)

# add pseudotime into seurat meta.data
cds_pst = pseudotime(cds,reduction_method = "UMAP")
smerge = AddMetaData(object =smerge, cds_pst, col.name ="B_pstime")



###### plot_density with Nebulosa ######

plot_density(smerge, c("CytoTRACE2_Score"), reduction = "umap")
plot_density(smerge, "KEGG_WNT_SIGNALING", reduction = "umap")
plot_density(smerge, c("CSNK1E","CSNK1D","CSNK1A1","CSNK1G3"), reduction = "umap")
plot_density(smerge, c("TNFRSF13B","TNFSF13"), reduction = "umap")



###### extracting meta.data information, and statistical test ######

###### meta.data
# all cells
metadata = smerge@meta.data
# B cells
metadata_b = subset(metadata, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18"))

###### T test

# T.test of SD diagnosis vs RR
t.test(subset(metadata_b, sample %in% c("pt01_dx","pt02_dx"))$CytoTRACE2_Score, subset(metadata_b, sample %in% c("pt01_rr","pt02_rr"))$CytoTRACE2_Score)
# mean of x mean of y: 0.3994222 0.4898647, p-value < 2.2e-16

# T.test of CR diagnosis vs RR
t.test(subset(metadata_b, sample %in% c("pt03_dx","pt04_dx"))$CytoTRACE2_Score, subset(metadata_b, sample %in% c("pt03_rr","pt04_rr"))$CytoTRACE2_Score)
# mean of x mean of y: 0.2340939 0.3080788, p-value < 2.2e-16

# T.test of SD diagnosis vs CR diagnosis
t.test(subset(metadata_b, sample %in% c("pt01_dx","pt02_dx"))$CytoTRACE2_Score, subset(metadata_b, sample %in% c("pt03_dx","pt04_dx"))$CytoTRACE2_Score)
# mean of x mean of y: 0.3994222 0.2340939, p-value < 2.2e-16


# Refractory signature
t.test(subset(metadata, RNA_snn_res.0.5 %in% c("8","9"))$RCHOP_Refractory_Up, subset(metadata, RNA_snn_res.0.5 %in% c("0"))$RCHOP_Refractory_Up)
# p-value < 2.2e-16
t.test(subset(metadata, RNA_snn_res.0.5 %in% c("8","9"))$RCHOP_Refractory_Up, subset(metadata, RNA_snn_res.0.5 %in% c("3"))$RCHOP_Refractory_Up)
# p-value < 2.2e-16
t.test(subset(metadata, RNA_snn_res.0.5 %in% c("8","9"))$RCHOP_Refractory_Up, subset(metadata, RNA_snn_res.0.5 %in% c("5"))$RCHOP_Refractory_Up)
# p-value < 2.2e-16
t.test(subset(metadata, RNA_snn_res.0.5 %in% c("8","9"))$RCHOP_Refractory_Up, subset(metadata, RNA_snn_res.0.5 %in% c("6"))$RCHOP_Refractory_Up)
# p-value < 2.2e-16

# KEGG_WNT_SIGNALING
t.test(subset(metadata, RNA_snn_res.0.5 %in% c("8","9"))$KEGG_WNT_SIGNALING, subset(metadata, RNA_snn_res.0.5 %in% c("0"))$KEGG_WNT_SIGNALING)
# p-value < 2.2e-16
t.test(subset(metadata, RNA_snn_res.0.5 %in% c("8","9"))$KEGG_WNT_SIGNALING, subset(metadata, RNA_snn_res.0.5 %in% c("3"))$KEGG_WNT_SIGNALING)
# p-value < 2.2e-16
t.test(subset(metadata, RNA_snn_res.0.5 %in% c("8","9"))$KEGG_WNT_SIGNALING, subset(metadata, RNA_snn_res.0.5 %in% c("5"))$KEGG_WNT_SIGNALING)
# p-value < 2.2e-16
t.test(subset(metadata, RNA_snn_res.0.5 %in% c("8","9"))$KEGG_WNT_SIGNALING, subset(metadata, RNA_snn_res.0.5 %in% c("6"))$KEGG_WNT_SIGNALING)
# p-value < 2.2e-16


###### correlation test
cor.test(metadata_b$RCHOP_Refractory_Up, metadata_b$KEGG_WNT_SIGNALING)
# cor = 0.3372518; p-value < 2.2e-16

cor.test(metadata_b$HSPC, metadata_b$KEGG_WNT_SIGNALING)
# cor = 0.4520674; p-value < 2.2e-16


###### Kolmogorov-Smirnov test

ks.test(subset(metadata_b, sample == ("pt01_rr"))$B_pstime, subset(metadata_b, sample == ("pt01_dx"))$B_pstime)
# D = 0.20011, p-value < 2.2e-16
ks.test(subset(metadata_b, sample == ("pt02_rr"))$B_pstime, subset(metadata_b, sample == ("pt02_dx"))$B_pstime)
# D = 0.48027, p-value < 2.2e-16

ks.test(subset(metadata_b, sample == ("pt03_rr"))$B_pstime, subset(metadata_b, sample == ("pt03_dx"))$B_pstime)
# D = 0.3846, p-value < 2.2e-16
ks.test(subset(metadata_b, sample == ("pt04_rr"))$B_pstime, subset(metadata_b, sample == ("pt04_dx"))$B_pstime)
# D = 0.3599, p-value < 2.2e-16



###### contour plot comparing disease stages ######
library(gridExtra)

embed_umap = Embeddings(smerge, reduction = "umap")
metadata_emm <- merge(metadata, embed_umap, by = "row.names")

# SD patients
ggplot(subset(metadata_emm, patient %in% c("pt01","pt02")), aes(x = umap_1, y = umap_2)) + geom_density_2d_filled(bins=9,alpha = 1) + scale_fill_brewer() + facet_grid(. ~ stage)+ theme_classic(base_size=18)+ theme(legend.position = "none") + xlim(-8,15) + ylim(-18,15)
# CR patients
ggplot(subset(metadata_emm, patient %in% c("pt03","pt04")), aes(x = umap_1, y = umap_2)) + geom_density_2d_filled(bins=9,alpha = 1) + scale_fill_brewer() + facet_grid(. ~ stage)+ theme_classic(base_size=18)+ theme(legend.position = "none") + xlim(-8,15) + ylim(-18,15)



###### add gene expression to metadata ######
ft_gene = FetchData(smerge, vars = c("CSNK1E","TNFRSF13B","TNFSF13")) 
ft_gene$CSNK1E_pos = ifelse(ft_gene$CSNK1E>0, "Positive","Negative")
ft_gene$TNFRSF13B_pos = ifelse(ft_gene$TNFRSF13B>0, "Positive","Negative")
smerge = AddMetaData(object =smerge, ft_gene[,c("CSNK1E","TNFRSF13B","TNFSF13","CSNK1E_pos","TNFRSF13B_pos")])



###### distribution/density plot ######

# B cell clusters: ("0","3","5","6","8","9","11","12","18")

# Comparing CytoTRACE2_Score at RR vs diagnosis in SD and CR patients
VlnPlot(subset(smerge, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18")), features = c("CytoTRACE2_Score"), pt.size = 0, group.by = "response", split.by = "stage", ncol = 1,split.plot = TRUE) & scale_fill_manual(values = c("#00BFC4","#F8766D"))


# Comparing pseudotime at RR vs diagnosis in each patient
plots <- lapply(c("pt01","pt02","pt03","pt04"), function(pt_name){
  plot <- ggplot(subset(metadata, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18") & patient == pt_name), aes(x=B_pstime, col=stage,fill = stage)) + geom_density(linewidth=1,alpha = 0.1) + theme_few(base_size=18) + scale_color_manual(values = c("#00BFC4","#F8766D")) + scale_fill_manual(values = c("#00BFC4","#F8766D"))+ labs(title="",x="Pseudotime", y = "Density")+ theme(legend.position = c(0.8, 0.8))+ ylim(0,0.17)
  print(plot)
  ggsave(filename = paste0(pt_name, ".B_pstime-density.pdf"), plot = plot, width = 4.5, height = 2.6)
  return(plot)
})



###### Cell cycle score ######
smerge <- CellCycleScoring(smerge, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
FeaturePlot_scCustom(smerge, colors_use = viridis_magma_light_high, features = c("G2M.Score"), num_columns = 1, na_cutoff = 0.1)



###### findmarkers to extract cluster signature ######

# signatures for persistent cluster
marker <- lapply(c("0","3","5","6"), function(cluster_id){
  cluster_marker = FindMarkers(smerge, ident.1 = cluster_id)
  cluster_marker = subset(cluster_marker, avg_log2FC > 0.5 & p_val_adj < 0.01 & pct.1 >0.1)
  write.table(cluster_marker, file = paste0("cluster.", cluster_id, ".signature.txt"), sep="\t", quote=F)
})

#signatures for persistent cluster 0 and 6 as a combined subset
write.table(FindMarkers(smerge, ident.1 = c("0","6")), "cluster.0and6.signature.txt", sep="\t", quote=F)



###### findmarkers for differential expression analysis ######

# differential expression in RR vs diagnosis
marker_RRvsDN = FindMarkers(smerge, ident.1 = "rr", ident.2 = "dx", group.by = 'stage',subset.ident = c("0","3","5","6","8","9","11","12","18"))
write.table(marker_RRvsDN, "RR-vs-Diagnosis.signature.txt", sep="\t", quote=F)

# differential expression in persistent clusters (0,3,5, and 6) vs sensitive clusters (8 and 9)
marker <- lapply(c("0","3","5","6"), function(cluster_id){
  cluster_marker = FindMarkers(smerge, ident.1 = cluster_id, ident.2 = c("8","9"))
  write.table(cluster_marker, file = paste0("cluster.", cluster_id, ".vs-cluster89.signature.txt"), sep="\t", quote=F)
})



###### violin plots ######
library(reshape2)
library(introdataviz)


##### RR versus diagnosis in B cells
# HSPC and progenitor signatures
gsetlist = c("HSPC","CLP_1","CLP_2","preB")
df_long = melt(metadata_b[, c("stage",gsetlist)], id.vars = "stage", variable.name ="geneset", value.name = "ES")
ggplot(df_long, aes(x=stage, y=ES, fill=stage)) + geom_violin(alpha = 0.9) + geom_boxplot(width=0.2, fill="white",outliers = F) + 
    scale_fill_manual(values = c("#00BFC4","#F8766D")) + facet_wrap(~geneset, scales="free_y",ncol =4) +
    theme_classic(base_size=18) + labs(x="", y = "Enrichment score") + 
    theme(axis.text=element_text(color="black"),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))

# HSPC and progenitor signatures
gsetlist = c("KEGG_WNT_SIGNALING","DIAPAUSE_UP_BOROVIAK","SAUL_SEN_MAYO","DUY_CISG_UP")
df_long = melt(metadata_b[, c("stage",gsetlist)], id.vars = "stage", variable.name ="geneset", value.name = "ES")
ggplot(df_long, aes(x=stage, y=ES, fill=stage)) + geom_violin(alpha = 0.9) + geom_boxplot(width=0.2, fill="white",outliers = F) + 
    scale_fill_manual(values = c("#00BFC4","#F8766D")) + facet_wrap(~geneset, scales="free_y",ncol =4) +
    theme_classic(base_size=18) + labs(x="", y = "Enrichment score") + 
    theme(axis.text=element_text(color="black"),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))


##### CSNK1E
# CSNK1E expression across B cell clusters
VlnPlot(subset(smerge, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18")), features = c("CSNK1E"), pt.size = 0.01, ncol = 1) & scale_fill_manual(values = color.lib.bcels)

# comparison by stage and drug response group
VlnPlot(subset(smerge, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18")), features = c("CSNK1E"), pt.size = 0, group.by = "response", split.by = "stage", ncol = 1) & scale_fill_manual(values = c("#00BFC4","#F8766D")) & ylim(0,2.8)

# enrichment score of KEGG_WNT_SIGNALING in CSNK1E positive and negative cells
ggplot(metadata_b, aes(x=CSNK1E_pos, y=KEGG_WNT_SIGNALING, fill=CSNK1E_pos)) +
    geom_violin() + geom_boxplot(width=0.1, fill="white") + scale_fill_brewer(palette="Dark2") + 
    theme_classic(base_size=18) + labs(title="CSNK1E",x="", y = "WNT SIGNALING") +
    theme(legend.position='none', axis.text=element_text(color="black"),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))

# HSPC and progenitor signatures in CSNK1E positive versus negative
gsetlist = c("HSPC","CLP_1","CLP_2","preB")
df_long = melt(metadata_b[, c("CSNK1E_pos",gsetlist)], id.vars = "CSNK1E_pos", variable.name ="geneset", value.name = "ES")
ggplot(df_long, aes(x=CSNK1E_pos, y=ES, fill=CSNK1E_pos)) + geom_violin(alpha = 0.9) + geom_boxplot(width=0.2, fill="white",outliers = F) + 
    scale_fill_brewer(palette="Dark2") + facet_wrap(~geneset, scales="free_y",ncol =4) + 
    theme_classic(base_size=18) + labs(x="", y = "Enrichment score") + 
    theme(axis.text=element_text(color="black"),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))

# CytoTRACE2_Score in CSNK1E positive versus negative
ggplot(metadata_b, aes(x=CSNK1E_pos, y=CytoTRACE2_Score, fill=CSNK1E_pos)) + 
    geom_violin() + geom_boxplot(width=0.1, fill="white") + scale_fill_brewer(palette="Dark2") + 
    theme_classic(base_size=18) + labs(title="CSNK1E",x="", y = "CytoTRACE2_Score") +
    theme(legend.position='none', axis.text=element_text(color="black"),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))


##### TNFRSF13B

# TNFRSF13B expression across B cell clusters
VlnPlot(subset(smerge, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18")), features = c("TNFRSF13B"), pt.size = 0.01, ncol = 1) & scale_fill_manual(values = color.lib.bcels)

# CSNK1E expression in TNFRSF13B positive and negative cells
VlnPlot(subset(smerge, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18")), features = c("CSNK1E"), pt.size = 0.1, group.by = "TNFRSF13B_pos", ncol = 1) & scale_fill_brewer(palette="Dark2")

# enrichment score of KEGG_WNT_SIGNALING in TNFRSF13B positive and negative cells
ggplot(metadata_b, aes(x=TNFRSF13B_pos, y=KEGG_WNT_SIGNALING, fill=TNFRSF13B_pos)) +
    geom_violin() + geom_boxplot(width=0.1, fill="white") + scale_fill_brewer(palette="Dark2") + 
    theme_classic(base_size=18) + labs(title="TNFRSF13B",x="", y = "WNT SIGNALING") +
    theme(legend.position='none', axis.text=element_text(color="black"),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))

# CytoTRACE2_Score in TNFRSF13B positive versus negative
ggplot(metadata_b, aes(x=TNFRSF13B_pos, y=CytoTRACE2_Score, fill=TNFRSF13B_pos)) + 
    geom_violin() + geom_boxplot(width=0.1, fill="white") + scale_fill_brewer(palette="Dark2") + 
    theme_classic(base_size=18) + labs(title="TNFRSF13B",x="", y = "CytoTRACE2_Score") +
    theme(legend.position='none', axis.text=element_text(color="black"),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18))



###### boxplot of ES in B cell clusters ######
# WNT signaling
ggplot(subset(metadata, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18")), aes(x=KEGG_WNT_SIGNALING, y=RNA_snn_res.0.5, fill=RNA_snn_res.0.5)) + 
  geom_boxplot(outlier.size = 0.2,outliers = F) + scale_fill_manual(values =color.lib.bcels) + 
  theme_classic(base_size = 18) + theme(legend.position="none", axis.text=element_text(color="black")) + 
  labs(x = "WNT SIGNALING", y="Cluster ID")+ xlim(0.07,0.3)

# Refractory signature
ggplot(subset(metadata, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18")), aes(x=RCHOP_Refractory_Up, y=RNA_snn_res.0.5, fill=RNA_snn_res.0.5)) + 
  geom_boxplot(outlier.size = 0.2,outliers = F) + scale_fill_manual(values =color.lib.bcels) + 
  theme_classic(base_size = 18) + theme(legend.position="none", axis.text=element_text(color="black")) + 
  labs(x = "Refractory signature", y="Cluster ID") + xlim(0.01,0.06)

# Pseudotime
ggplot(subset(metadata, RNA_snn_res.0.5 %in% c("0","3","5","6","8","9","11","12","18")), aes(x=B_pstime, y=RNA_snn_res.0.5, fill=RNA_snn_res.0.5)) + 
  geom_boxplot(outlier.size = 0.2,outliers = F) + scale_fill_manual(values =color.lib.bcels) + 
  theme_classic(base_size = 18) + theme(legend.position="none", axis.text=element_text(color="black")) + 
  labs(x = "Pseudotime", y="Cluster ID") + xlim(0.01,0.06)



###### trendline plot ######
list = c("KEGG_WNT_SIGNALING","RCHOP_Refractory_Up","CSNK1E")

# subset B cells with Pseudotime results
metadata_fin = metadata_b[is.finite(metadata_b$B_pstime),]

# The initial cell number was too high for fitness, use half of the data
set.seed(123)
sample_data <- metadata_fin[sample(nrow(metadata_fin), nrow(metadata_fin)/2), ]  

# WNT SIGNALING with Pseudotime
ggplot(sample_data) + geom_smooth(aes(B_pstime,KEGG_WNT_SIGNALING), method="loess",linewidth=1,level=0.95, span = 0.7) + 
  theme_few(base_size=18) + scale_color_manual(values = c("#00BFC4","#F8766D")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    text = element_text(size = 18), axis.text=element_text(colour="black"), aspect.ratio = 0.75) + xlab("Pseudotime") + ylab("WNT SIGNALING")
# Refractory signature with Pseudotime
ggplot(sample_data) + geom_smooth(aes(B_pstime,RCHOP_Refractory_Up), method="loess",linewidth=1,level=0.95, span = 0.7) + 
  theme_few(base_size=18) + scale_color_manual(values = c("#00BFC4","#F8766D")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    text = element_text(size = 18), axis.text=element_text(colour="black"), aspect.ratio = 0.75) + xlab("Pseudotime") + ylab("Refractory signature")
# CSNK1E with Pseudotime
ggplot(sample_data) + geom_smooth(aes(B_pstime,CSNK1E), method="loess",linewidth=1,level=0.95, span = 0.7) + 
  theme_few(base_size=18) + scale_color_manual(values = c("#00BFC4","#F8766D")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    text = element_text(size = 18), axis.text=element_text(colour="black"), aspect.ratio = 0.75) + xlab("Pseudotime") + ylab("CSNK1E")



###### correlation plot ######

# correlation plots with multiple gene sets
pair = c("KEGG_WNT_SIGNALING","B_pstime","RCHOP_Refractory_Up","HSPC","CLP_1","CLP_2","preB","DIAPAUSE_UP_BOROVIAK","SAUL_SEN_MAYO","DUY_CISG_UP")

corrplot(cor(metadata_fin[,pair]), type="lower", tl.col="black", tl.srt=35,insig = "blank", addCoef.col = "black", method="ellipse", sig.level = 0.05,cl.ratio = 0.15, col=rev(brewer.pal(n=10, name="RdBu")))

# correlation plots with two gene sets, WNT_SIGNALING and Refractory Signature
ggplot(metadata_b), aes(x = RCHOP_Refractory_Up, y = KEGG_WNT_SIGNALING)) + 
  geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Viridis") + geom_smooth(method=lm) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 18), axis.text=element_text(colour="black"), aspect.ratio = 1) + 
  xlab("Refractory Signature") + ylab("WNT SIGNALING") + coord_fixed()

# correlation plots with two gene sets, WNT_SIGNALING and HSPC Signature
ggplot(metadata_b), aes(x = HSPC, y = KEGG_WNT_SIGNALING)) + 
  geom_bin2d() + paletteer::scale_fill_paletteer_c("grDevices::Viridis") + geom_smooth(method=lm) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", text = element_text(size = 18), axis.text=element_text(colour="black"), aspect.ratio = 1) + 
  xlab("HSPC Signature") + ylab("WNT SIGNALING") + coord_fixed()



###### cellChat APRIL-TNFRSF13B communication analysis ######

# cellChat differential analysis starts with separate analysis for each group

# subset seurat object with B cells and APRIL-expression cells
sce_dx = subset(smerge, cell_type_ID %in% c("B_0","B_3","B_5","B_6","B_8","B_9","B_11","B_18","MONO_4","MAC_17","pDC_19") & sample %in% c("pt01_dx","pt02_dx","pt03_dx","pt04_dx"))
sce_B_rr = subset(smerge, cell_type_ID %in% c("B_0","B_3","B_5","B_6","B_8","B_9","B_11","B_18","MONO_4","MAC_17","pDC_19") & sample %in% c("pt01_rr","pt02_rr","pt03_rr","pt04_rr"))

sce_dx$cell_type_ID = factor(sce_dx$cell_type_ID, levels=c("B_0","B_3","B_5","B_6","B_8","B_9","B_11","B_18","MONO_4","MAC_17","pDC_19"))
sce_rr$cell_type_ID = factor(sce_rr$cell_type_ID, levels=c("B_0","B_3","B_5","B_6","B_8","B_9","B_11","B_18","MONO_4","MAC_17","pDC_19"))

# create cellchat object
cellchat_dx <- createCellChat(object = sce_dx, group.by = "cell_type_ID", assay = "RNA")
cellchat_rr <- createCellChat(object = sce_rr, group.by = "cell_type_ID", assay = "RNA")

##### common parameters
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# use the Secreted Signaling of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
options(future.globals.maxSize= 2891289600)

##### cellChat analysis for diagnosis and RR samples separately

# samples at diagnosis
cellchat_dx@DB <- CellChatDB.use

cellchat_dx <- subsetData(cellchat_dx)
future::plan("multisession", workers = 4) # do parallel
cellchat_dx <- identifyOverExpressedGenes(cellchat_dx)
cellchat_dx <- identifyOverExpressedInteractions(cellchat_dx)
cellchat_dx <- computeCommunProb(cellchat_dx, type = "triMean")
cellchat_dx <- filterCommunication(cellchat_dx, min.cells = 20)
cellchat_dx <- computeCommunProbPathway(cellchat_dx)
cellchat_dx <- aggregateNet(cellchat_dx)

# significant interactin
cellchat_dx@netP$pathways

# samples at RR
cellchat_rr@DB <- CellChatDB.use

cellchat_rr <- subsetData(cellchat_rr)
future::plan("multisession", workers = 4) # do parallel
cellchat_rr <- identifyOverExpressedGenes(cellchat_rr)
cellchat_rr <- identifyOverExpressedInteractions(cellchat_rr)
cellchat_rr <- computeCommunProb(cellchat_rr, type = "triMean")
cellchat_rr <- filterCommunication(cellchat_rr, min.cells = 20)
cellchat_rr <- computeCommunProbPathway(cellchat_rr)
cellchat_rr <- aggregateNet(cellchat_rr)

# significant interactin
cellchat_rr@netP$pathways

##### Merge analysis

object.list <- list(Diag = cellchat_dx, RR = cellchat_rr)
cellchat_all <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

pathways.show <- c("APRIL")

# plot the heatmap comparing cell-cell interaction through APRIL
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))





##############################################################################################################################
# DLBCL dataset by Roider et al. Nat Cell Biol. 2020
##############################################################################################################################

###### seurat object and QC ######

# load cellranger files, and create seurat object
# cellranger files were downloaded from: https://heidata.uni-heidelberg.de/dataset.xhtml?persistentId=doi:10.11588/data/VRJUNV
bcl1=Read10X(data.dir="./DLBCL1")
bcl2=Read10X(data.dir="./DLBCL2")
bcl3=Read10X(data.dir="./DLBCL3")

sbl1 = CreateSeuratObject(counts = bcl1, project = "bcl1", min.cells = 5, min.features = 500)
sbl2 = CreateSeuratObject(counts = bcl2, project = "bcl2", min.cells = 5, min.features = 500)
sbl3 = CreateSeuratObject(counts = bcl3, project = "bcl3", min.cells = 5, min.features = 500)

# QC and filter
sbl1[["percent.mt"]] = PercentageFeatureSet(sbl1, pattern = "^MT-")
sbl2[["percent.mt"]] = PercentageFeatureSet(sbl2, pattern = "^MT-")
sbl3[["percent.mt"]] = PercentageFeatureSet(sbl3, pattern = "^MT-")

sbl1 = subset(sbl1, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA < 40000 & nCount_RNA > 2000)
sbl2 = subset(sbl2, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA < 40000 & nCount_RNA > 2000)
sbl3 = subset(sbl3, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA < 40000 & nCount_RNA > 2000)

# normalization
sbl1 = NormalizeData(sbl1)
sbl2 = NormalizeData(sbl2)
sbl3 = NormalizeData(sbl3)

# merge objects of different samples
sbl1$case = "DLBCL1"
sbl2$case = "DLBCL2"
sbl3$case = "DLBCL3"

smerge = merge(sbl1, y = c(sbl2, sbl3), add.cell.ids = c("s1", "s2", "s3"), project = "DLBCL")



###### UMAP clustering ######

smerge = FindVariableFeatures(smerge, selection.method = "vst")
smerge = ScaleData(smerge, verbose = FALSE)
smerge = RunPCA(smerge, npcs = 30, verbose = FALSE)
ElbowPlot(object = smerge)

length(VariableFeatures(smerge))
smerge = RunHarmony(smerge, "case")
smerge = FindNeighbors(smerge, reduction = "harmony", dims = 1:20)
smerge = FindClusters(smerge)
smerge = RunTSNE(smerge, reduction = "harmony", dims = 1:15)



###### SingleR annotation ######
# v5 is not supported, needs to be converted to v3 before SingleR
smerge[["RNA"]] <- as(smerge[["RNA"]], Class="Assay")

# main annotation
ref <- fetchReference("blueprint_encode", "2024-02-26")
pred.BlueprintEncodeData <- SingleR(test = smerge@assays$RNA@data, ref = ref, labels = ref$label.main)
smerge@meta.data$CellType_BlueprintEncode <- pred.BlueprintEncodeData$labels

# dimplot
DimPlot(smerge, reduction = "umap", group.by = c("case", "seurat_clusters","CellType_BlueprintEncode"),label = T,label.size = 5)



###### geneset analysis ######

counts <- LayerData(smerge, assay = "RNA", layer = "counts")
expressed_gene = rownames(counts)
rm(counts)

gset_cur <- getGmt("gene_sets_scRNAseq.gmt")
gset_cur <- subsetGeneSets(gset_cur, expressed_gene) 
smerge <- AddModuleScore_UCell(smerge, features=geneIds(gset_cur), name=NULL)



###### plot_density with Nebulosa ######

plot_density(smerge, c("CSNK1E"), reduction = "umap")
plot_density(smerge, c("HSPC_2","CLP_1","CLP_2","prePB"), reduction = "umap")



###### monocle3 analysis ######

DefaultAssay(smerge) <- "RNA"

# convert seurat to single cell object for monocle analysis
cds <- as.cell_data_set(smerge)
cds <- cluster_cells(cds, reduction_method = "UMAP")

# trajectory graph
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 5)

# order B cells clusters based on CytoTRACE2 and HSPC signature score
cds <- order_cells(cds)

# pseudotime UMAP
plot_cells(cds, color_cells_by = "pseudotime",  group_cells_by = "cluster", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, label_roots = T, trajectory_graph_color = "grey80", trajectory_graph_segment_size = 1,cell_size = 0.1)

# add pseudotime into seurat meta.data
cds_pst = pseudotime(cds,reduction_method = "UMAP")
smerge = AddMetaData(object =smerge, cds_pst, col.name ="B_pstime")



###### violin plots ######
library(reshape2)
library(introdataviz)

metadata = smerge@meta.data
metadata$CSNK1E_pos = ifelse(metadata$CSNK1E>0, "Positive","Negative")
smerge = AddMetaData(object =smerge, metadata[,c("CSNK1E_pos")])


checklist = c("HSPC_2","CLP_1","CLP_2","preB")
df_long = melt(metadata[, c("CSNK1E_pos",checklist)], id.vars = "CSNK1E_pos", variable.name ="geneset", value.name = "ES")

ggplot(df_long, aes(x = geneset, y = ES, fill = CSNK1E_pos)) +
    introdataviz::geom_split_violin(alpha = .6, trim = FALSE) +
    geom_boxplot(width = .2, alpha = .8, fatten = NULL, show.legend = FALSE, , outlier.colour = NA) +
    stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, position = position_dodge(.175)) +
    scale_fill_brewer(palette="Dark2") +
    theme_classic() + theme(axis.text = element_text(colour="black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.ticks = element_line(colour = "black",linewidth=0.5), axis.ticks.length = unit(.2, "cm"), axis.line = element_line(colour = "black",linewidth=0.5), panel.grid.major.x = element_blank())


###### trend line plots ######

ggplot(subset(metadata, CellType_BlueprintEncode %in% c("B-cells"))) + geom_smooth(aes(B_pstime,CSNK1E), method="loess",linewidth=2,level=0.99) + theme_classic(base_size=18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 18), axis.text=element_text(colour="black"), aspect.ratio = 1) + xlab("Pseudotime") + ylab("CSNK1E Expression")

