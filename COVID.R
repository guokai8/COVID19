library(Seurat)
C141<-Read10X_h5("GSM4339769_C141_filtered_feature_bc_matrix.h5")
C142<-Read10X_h5("GSM4339770_C142_filtered_feature_bc_matrix.h5")
C143<-Read10X_h5("GSM4339771_C143_filtered_feature_bc_matrix.h5")
C144<-Read10X_h5("GSM4339772_C144_filtered_feature_bc_matrix.h5")
C145<-Read10X_h5("GSM4339773_C145_filtered_feature_bc_matrix.h5")
C146<-Read10X_h5("GSM4339774_C146_filtered_feature_bc_matrix.h5")
C148<-Read10X_h5("GSM4475051_C148_filtered_feature_bc_matrix.h5")
C149<-Read10X_h5("GSM4475052_C149_filtered_feature_bc_matrix.h5")
C152<-Read10X_h5("GSM4475053_C152_filtered_feature_bc_matrix.h5")
######
C141<-CreateSeuratObject(counts = C141, project = "C141",min.cells = 3, min.features = 200)
C142<-CreateSeuratObject(counts = C142, project = "C142",min.cells = 3, min.features = 200)
C143<-CreateSeuratObject(counts = C143, project = "C143",min.cells = 3, min.features = 200)
C144<-CreateSeuratObject(counts = C144, project = "C144",min.cells = 3, min.features = 200)
C145<-CreateSeuratObject(counts = C145, project = "C145",min.cells = 3, min.features = 200)
C146<-CreateSeuratObject(counts = C146, project = "C146",min.cells = 3, min.features = 200)
C148<-CreateSeuratObject(counts = C148, project = "C148",min.cells = 3, min.features = 200)
C149<-CreateSeuratObject(counts = C149, project = "C149",min.cells = 3, min.features = 200)
C152<-CreateSeuratObject(counts = C152, project = "C152",min.cells = 3, min.features = 200)
############################
C141$group<-"mild"
C142$group<-"mild"
C143$group<-"severe"
C144$group<-"mild"
C145$group<-"severe"
C146$group<-"severe"
C148$group<-"healthy"
C149$group<-"healthy"
C152$group<-"healthy"
#################################
C141[["percent.mt"]]<-PercentageFeatureSet(C141,pattern = "^MT")
p141<-VlnPlot(C141, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C142[["percent.mt"]]<-PercentageFeatureSet(C142,pattern = "^MT")
p142<-VlnPlot(C142, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C143[["percent.mt"]]<-PercentageFeatureSet(C143,pattern = "^MT")
p143<-VlnPlot(C143, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C144[["percent.mt"]]<-PercentageFeatureSet(C144,pattern = "^MT")
p144<-VlnPlot(C144, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C145[["percent.mt"]]<-PercentageFeatureSet(C145,pattern = "^MT")
p145<-VlnPlot(C145, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C146[["percent.mt"]]<-PercentageFeatureSet(C146,pattern = "^MT")
p146<-VlnPlot(C146, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C148[["percent.mt"]]<-PercentageFeatureSet(C148,pattern = "^MT")
p148<-VlnPlot(C148, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C149[["percent.mt"]]<-PercentageFeatureSet(C149,pattern = "^MT")
p149<-VlnPlot(C149, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C152[["percent.mt"]]<-PercentageFeatureSet(C152,pattern = "^MT")
p152<-VlnPlot(C152, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
####################
C141 <- subset(C141, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C142 <- subset(C142, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C143 <- subset(C143, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C144 <- subset(C144, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C145 <- subset(C145, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C146 <- subset(C146, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C148 <- subset(C148, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C149 <- subset(C149, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
C152 <- subset(C152, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
#############################################################################
C141 <- NormalizeData(C141)
C142 <- NormalizeData(C142)
C143 <- NormalizeData(C143)
C144 <- NormalizeData(C144)
C145 <- NormalizeData(C145)
C146 <- NormalizeData(C146)
C148 <- NormalizeData(C148)
C149 <- NormalizeData(C149)
C152 <- NormalizeData(C152)
##########################
C141 <- FindVariableFeatures(C141, selection.method = "vst", nfeatures = 2000)
C142 <- FindVariableFeatures(C142, selection.method = "vst", nfeatures = 2000)
C143 <- FindVariableFeatures(C143, selection.method = "vst", nfeatures = 2000)
C144 <- FindVariableFeatures(C144, selection.method = "vst", nfeatures = 2000)
C145 <- FindVariableFeatures(C145, selection.method = "vst", nfeatures = 2000)
C146 <- FindVariableFeatures(C146, selection.method = "vst", nfeatures = 2000)
C148 <- FindVariableFeatures(C148, selection.method = "vst", nfeatures = 2000)
C149 <- FindVariableFeatures(C149, selection.method = "vst", nfeatures = 2000)
C152 <- FindVariableFeatures(C152, selection.method = "vst", nfeatures = 2000)
#####################################
sample.anchors <- FindIntegrationAnchors(object.list = list(C141,C142,C143,C144,C145,C146,C148,C149,C152), dims = 1:50)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:50)
#########################################
DefaultAssay(sample.combined) <- "integrated"
all.genes <- rownames(sample.combined)
sample.combined <- ScaleData(sample.combined, verbose = TRUE,features = all.genes)
sample.combined <- RunPCA(sample.combined, npcs = 50, features = VariableFeatures(object = sample.combined))
sample.combined <- JackStraw(sample.combined, num.replicate = 50)
sample.combined <- ScoreJackStraw(sample.combined, dims = 1:20)
VizDimLoadings(sample.combined, dims = 1:2, reduction = "pca")
DimPlot(sample.combined, reduction = "pca")
######################################################
sample.combined <- RunUMAP(sample.combined, reduction = "pca", dims = 1:50)
sample.combined <- FindNeighbors(sample.combined, reduction = "pca", dims = 1:50)
pp1 <- DimPlot(sample.combined, reduction = "umap", group.by = "group")
pp2 <- DimPlot(sample.combined, reduction = "umap", label = TRUE)
sample.combined <- FindClusters(sample.combined, resolution = 0.85)
DimPlot(sample.combined, reduction = "umap", split.by = "group",label=T)
########################################################
#sample.combined.bk<-sample.combined
#############################################################
DefaultAssay(sample.combined) <- "RNA"
sample.markers <- FindAllMarkers(sample.combined, min.pct = 0.25, logfc.threshold = 0.25)
FeaturePlot(sample.combined, features = c("AGER","SFTPC","SCGB3A2","TPPP3","KRT5","CD68","FCN1","CD1C","TPSB2","CD3D","IGHG4","MS4A1","VWF","DCN"), min.cutoff = "q9")
######
clusterid<-c("Macrophage","Macrophage","Macrophage","Macrophage","CD4+ T","Macrophage","CD8+ T",
             "Macrophage","Macrophage","Macrophage","NK","Macrophage","Epithelial",
             "Macrophage","Macrophage","NK","Epithelial","B","NK","Macrophage",
             "Macrophage","NK","Macrophage","B","Macrophage","NK")
names(clusterid) <- levels(sample.combined)
sample.combined <- RenameIdents(sample.combined, clusterid)
DimPlot(sample.combined, reduction = "umap", split.by = "group",label=T,repel = T)
####################
VlnPlot(sample.combined, features = c("CD68","IL7R","CD4","CD8A","TPPP3","MS4A1"),pt.size = 0)
dev.print(pdf,file="cluster_expression_cd4.pdf")
#############
VlnPlot(sample.combined, features = c("CD68","IL7R","CD4","CD8A","TPPP3","MS4A1"),pt.size = 0,split.by = "group")
VlnPlot(sample.combined, features = c("APOC1"),pt.size = 0,split.by = "group")
####################MAST part
sample<-sample.combined
Idents(sample)<-"condition"
mildvshea_macrophage<-FindMarkers(sample, ident.1 = "mild_Macrophage", ident.2 = "healthy_Macrophage", test.use = "MAST",verbose = FALSE)
severevshea_macrophage<-FindMarkers(sample, ident.1 = "severe_Macrophage", ident.2 = "healthy_Macrophage", test.use = "MAST",verbose = FALSE)
severevsmild_macrophage<-FindMarkers(sample, ident.1 = "severe_Macrophage", ident.2 = "mild_Macrophage", test.use = "MAST",verbose = FALSE)
#
mildvshea_NK<-FindMarkers(sample, ident.1 = "mild_NK", ident.2 = "healthy_NK", test.use = "MAST",verbose = F)
severevshea_NK<-FindMarkers(sample, ident.1 = "severe_NK", ident.2 = "healthy_NK", test.use = "MAST",verbose = F)
severevsmild_NK<-FindMarkers(sample, ident.1 = "severe_NK", ident.2 = "mild_NK", test.use = "MAST",verbose = F)
#
mildvshea_ccr<-FindMarkers(sample, ident.1 = "mild_CCR+", ident.2 = "healthy_CCR+", test.use = "MAST",verbose = F)
severevshea_ccr<-FindMarkers(sample, ident.1 = "severe_CCR+", ident.2 = "healthy_CCR+", test.use = "MAST",verbose = F)
severevsmild_ccr<-FindMarkers(sample, ident.1 = "severe_CCR+", ident.2 = "mild_CCR+", test.use = "MAST",verbose = F)
##
mildvshea_cd8<-FindMarkers(sample, ident.1 = "mild_CD8+ T", ident.2 = "healthy_CD8+ T", test.use = "MAST",verbose = F)
severevshea_cd8<-FindMarkers(sample, ident.1 = "severe_CD8+ T", ident.2 = "healthy_CD8+ T", test.use = "MAST",verbose = F)
severevsmild_cd8<-FindMarkers(sample, ident.1 = "severe_CD8+ T", ident.2 = "mild_CD8+ T", test.use = "MAST",verbose = F)
###
mildvshea_b<-FindMarkers(sample, ident.1 = "mild_B", ident.2 = "healthy_B", test.use = "MAST",verbose = F)
severevshea_b<-FindMarkers(sample, ident.1 = "severe_B", ident.2 = "healthy_B", test.use = "MAST",verbose = F)
severevsmild_b<-FindMarkers(sample, ident.1 = "severe_B", ident.2 = "mild_B", test.use = "MAST",verbose = F)
###
mildvshea_epi<-FindMarkers(sample, ident.1 = "mild_Epithelial", ident.2 = "healthy_Epithelial", test.use = "MAST",verbose = F)
severevshea_epi<-FindMarkers(sample, ident.1 = "severe_Epithelial", ident.2 = "healthy_Epithelial", test.use = "MAST",verbose = F)
severevsmild_epi<-FindMarkers(sample, ident.1 = "severe_Epithelial", ident.2 = "mild_Epithelial", test.use = "MAST",verbose = F)
######################################################
write.csv(mildvshea_macrophage,file="mildvshea_macrophage.csv")
write.csv(mildvshea_b,file="mildvshea_b.csv")
write.csv(mildvshea_ccr,file="mildvshea_ccr.csv")
write.csv(mildvshea_cd8,file="mildvshea_cd8.csv")
write.csv(mildvshea_epi,file="mildvshea_epi.csv")
write.csv(mildvshea_NK,file="mildvshea_NK.csv")
################################################
write.csv(severevsmild_macrophage,file="severevsmild_macrophage.csv")
write.csv(severevsmild_b,file="severevsmild_b.csv")
write.csv(severevsmild_ccr,file="severevsmild_ccr.csv")
write.csv(severevsmild_cd8,file="severevsmild_cd8.csv")
write.csv(severevsmild_epi,file="severevsmild_epi.csv")
write.csv(severevsmild_NK,file="severevsmild_NK.csv")
####################################################

