library(tidyverse)
library(Seurat)
library(Matrix)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)
library(ggthemes)
library(future)
library(tidyverse)

# plan("multiprocess", workers=4)
plan("sequential")

setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_Hnf4")

date="20221212"

cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value
cmap <- c(cmap, "firebrick3", "navyblue")

#####
#read gene expression data
dat <- Read10X( data.dir="data/gex" )
seurat_table <- CreateSeuratObject(dat$`Gene Expression`, min.cells = 3, project="MEClo_Hnf4", min.features=0)

#integrate hashes
hashtag_table <- as.matrix(dat$`Antibody Capture`)
hashtag_table <- hashtag_table[c("HT1","HT2","HT3","HT5","HT6","HT8","HT9"),]
length(intersect(colnames(seurat_table), colnames(hashtag_table))) / length(union(colnames(seurat_table), colnames(hashtag_table)))

seurat_table[["Hash"]] <- CreateAssayObject( counts=hashtag_table )

rm(dat)

seurat_table <- NormalizeData(seurat_table, assay = "Hash", normalization.method = "CLR")
seurat_table <- ScaleData(seurat_table, assay = "Hash")

#####
#remove doublets and assign hash identities

pdf( paste0("figures/qc_hash_",date,".pdf" ), height=5, width=8 )
for (i in 1:7) {
  x <- ggplot( as.data.frame((seurat_table@assays$Hash[i,])) ) +
    geom_histogram( aes( x=seurat_table@assays$Hash[i,] ),
                    fill=cmap[i], color="black" ) +
    geom_vline( xintercept=1 ) +
    geom_vline( xintercept=2 ) +
    geom_vline( xintercept=3 ) +
    geom_vline( xintercept=4 ) +
    xlim(0,6) +
    theme_classic() +
    ggtitle( rownames(seurat_table@assays$Hash[i,]) ) +
    xlab( rownames(seurat_table@assays$Hash[i,]) )
  print(x)
}
dev.off()

png(paste0("figures/hash_pairs_plot_",date,".png"), height=1000, width=1000 )
pairs( t(seurat_table@assays$Hash[,]), cex=0.2, pch=20 )
dev.off()

seurat_table <- HTODemux(seurat_table, assay = "Hash", positive.quantile = 0.99)
table(seurat_table$Hash_classification.global)
seurat_table <- seurat_table[,seurat_table$hash.ID %in% c("Doublet","Negative") == F]

ident <- vector(mode="logical", length=ncol(seurat_table))
for( i in 1:ncol(seurat_table) ) {
  if( seurat_table$hash.ID[i] == "HT1" ) {
    ident[i] <- "wt_r1"
  } else if ( seurat_table$hash.ID[i] == "HT2" ) {
    ident[i] <- "wt_r2"
  } else if ( seurat_table$hash.ID[i] == "HT3" ) {
    ident[i] <- "wt_r3"
  } else if ( seurat_table$hash.ID[i] == "HT5" ) {
    ident[i] <- "hnf4_r1"
  } else if ( seurat_table$hash.ID[i] == "HT6" ) {
    ident[i] <- "hnf4_r2"
  } else if ( seurat_table$hash.ID[i] == "HT8" ) {
    ident[i] <- "hnf4_r3"
  } else if ( seurat_table$hash.ID[i] == "HT9" ) {
    ident[i] <- "hnf4_r4"
  } else {
    ident[i] <- "na"
  }
}

names(ident) <- colnames(seurat_table)
seurat_table@meta.data$hash.ident <- factor(ident)

ident <- vector(mode="logical", length=ncol(seurat_table))
for( i in 1:ncol(seurat_table) ) {
  if( seurat_table$hash.ID[i] == "HT1" ) {
    ident[i] <- "wt"
  } else if ( seurat_table$hash.ID[i] == "HT2" ) {
    ident[i] <- "wt"
  } else if ( seurat_table$hash.ID[i] == "HT3" ) {
    ident[i] <- "wt"
  } else if ( seurat_table$hash.ID[i] == "HT5" ) {
    ident[i] <- "hnf4"
  } else if ( seurat_table$hash.ID[i] == "HT6" ) {
    ident[i] <- "hnf4"
  } else if ( seurat_table$hash.ID[i] == "HT8" ) {
    ident[i] <- "hnf4"
  } else if ( seurat_table$hash.ID[i] == "HT9" ) {
    ident[i] <- "hnf4"
  } else {
    ident[i] <- "na"
  }
}

names(ident) <- colnames(seurat_table)
seurat_table@meta.data$geno.ident <- factor(ident, levels=c("wt","hnf4"))

#####
#filter data

seurat_table[["percent.mt"]] <- PercentageFeatureSet(seurat_table, pattern = "^mt-")

x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  geom_hline( yintercept=1250, color="red" ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="steelblue" ) +
  theme_bw() +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="pink" ) +
  theme_bw() +
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=7.5, color="red" ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="orchid4" ) +
  theme_bw() +
  xlab("")

pdf( paste0("figures/qc_violin_prefilter_",date,".pdf" ), height=3, width=5 )
plot_grid( x,y,z, ncol=3 )
dev.off()

# VlnPlot(seurat_table, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_table <- subset(seurat_table, subset = percent.mt < 7.5)
seurat_table <- subset(seurat_table, subset = nFeature_RNA > 1250)

x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  geom_hline( yintercept=100, color="red" ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="steelblue" ) +
  theme_bw() +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="pink" ) +
  theme_bw() +
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=10, color="red" ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="orchid4" ) +
  theme_bw() +
  xlab("")

pdf( paste0("figures/qc_violin_postfilter_",date,".pdf" ), height=3, width=5 )
plot_grid( x,y,z, ncol=3 )
dev.off()


#####
#normalize data
seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)

nfeatures = 2000
seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)

top10 <- head(VariableFeatures(seurat_table), 10)

plot1 <- VariableFeaturePlot(seurat_table)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf( paste0("figures/variable_feature_postfilter_",date,".pdf" ), height=4, width=6 )
plot1
plot2
dev.off()

#####
#run dimensionality reduction
all.genes <- rownames(seurat_table)
seurat_table <- ScaleData(seurat_table, features = VariableFeatures(object = seurat_table))
seurat_table <- RunPCA(seurat_table, features = VariableFeatures(object = seurat_table))

sink( "text_outputs/top_PC_genes_1-30_postfilter.txt", append=F )
print(seurat_table[["pca"]], dims = 1:30, nfeatures = 5)
sink()

pdf( paste0("figures/PC_dim_loading_postfilter_",date,".pdf" ), height=40, width=12 )
VizDimLoadings(seurat_table, dims = 1:30, reduction = "pca")
dev.off()

DimPlot(seurat_table, reduction = "pca", group.by="geno.ident")

pdf( paste0("figures/PC_dim_heatmap_postfilter_",date,".pdf" ), height=40, width=12 )
DimHeatmap(seurat_table, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

seurat_table <- JackStraw(seurat_table, num.replicate = 100, dims=50)
seurat_table <- ScoreJackStraw(seurat_table, dims = 1:50)

pdf( paste0("figures/PC_jackstraw_",date,".pdf" ), height=5, width=12 )
JackStrawPlot(seurat_table, dims = 1:50)
dev.off()

pdf( paste0("figures/PC_elbow_postfilter_",date,".pdf" ), height=5, width=6 )
ElbowPlot(seurat_table, ndims=50)
dev.off()

#####
# cluster and umap
seurat_table <- FindNeighbors(seurat_table, dims = 1:30)
seurat_table <- FindClusters(seurat_table, resolution = 2.2)
seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 123)

pdf( paste0("figures/umap_postfilter_",date,".pdf" ), height=8, width=12 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, pt.size=1)
dev.off()

#####
#analyze cluster composition
cluster.markers <- FindMarkers(seurat_table, ident.1=c(16,20), min.pct=0.01, logfc.threshold=0.25, only.pos=T, test.use="wilcox")
head(cluster.markers, n = 20)

#####
#remove contaminating cells
keep <- seurat_table@active.ident!=22 #t
seurat_table <- seurat_table[,keep]

#####
#combine homologous clusters
keep <- seurat_table@active.ident==1 | seurat_table@active.ident==2 | seurat_table@active.ident==12 | seurat_table@active.ident==18 | 
  seurat_table@active.ident==4 | seurat_table@active.ident==19 | seurat_table@active.ident==7 | seurat_table@active.ident==0 |
  seurat_table@active.ident==5 | seurat_table@active.ident==25
seurat_table@active.ident[keep] <- 1

keep <- seurat_table@active.ident==24 | seurat_table@active.ident==26
seurat_table@active.ident[keep] <- 24

keep <- seurat_table@active.ident==20 | seurat_table@active.ident==16
seurat_table@active.ident[keep] <- 16

keep <- seurat_table@active.ident==14 | seurat_table@active.ident==6
seurat_table@active.ident[keep] <- 6

keep <- seurat_table@active.ident==8 | seurat_table@active.ident==15
seurat_table@active.ident[keep] <- 8

keep <- seurat_table@active.ident==9 | seurat_table@active.ident==17
seurat_table@active.ident[keep] <- 9


#####
#plot per hash umap
keep1 <- seurat_table@meta.data$hash.ident %in% c( "wt_r1","wt_r2","wt_r3" )
keep2 <- seurat_table@meta.data$hash.ident %in% c( "hnf4_r1","hnf4_r2", "hnf4_r3", "hnf4_r4" )

w <- DimPlot(seurat_table, group.by="geno.ident", pt.size=0.5, order=F, cols=cmap )
x <- DimPlot(seurat_table[,keep1], group.by="hash.ident", pt.size=0.5, order=F, cols=cmap[c(3,4,5,6)] )
y <- DimPlot(seurat_table[,keep2], group.by="hash.ident", pt.size=0.5, order=F, cols=cmap )

# pdf( paste0("figures/umap_perGenotype_",date,".pdf" ), height=5, width=20 )
plot_grid(w,x,y, ncol=3)
# dev.off()

w <- DimPlot(seurat_table, group.by="geno.ident", pt.size=0.5, order=F, cols=cmap )
x <- DimPlot(seurat_table[,keep1], group.by="geno.ident", pt.size=0.5, order=F, cols="gray" )
y <- DimPlot(seurat_table[,keep2], group.by="geno.ident", pt.size=0.5, order=F, cols="gray" )

# pdf( paste0("figures/umap_perGenotype_monochrome_",date,".pdf" ), height=5, width=20 )
plot_grid(w,x,y, ncol=3)
# dev.off()

w <- DimPlot(seurat_table, pt.size=0.5, order=F, label=T, repel=T, cols=cmap ) + theme(legend.position="none") + ggtitle("Merged")
x <- DimPlot(seurat_table[,keep1], pt.size=0.5, order=F, cols=cmap, ) + theme(legend.position="none") + ggtitle("Hnf4_ctrl")
y <- DimPlot(seurat_table[,keep2], pt.size=0.5, order=F, cols=cmap ) + theme(legend.position="none") + ggtitle("Hnf4_dTEC")

pdf( paste0("figures/umap_perGenotype_",date,".pdf" ), height=6, width=17 )
plot_grid(w,x,y, ncol=3)
dev.off()

#####
#calculate fractional representation of clusters in wt vs hnf4

source("functions/calculate_cluster_fraction.R")

cluster_frac_list <- calculate_cluster_fraction(n_hash=7, 
                                                names_hash=unique(seurat_table@meta.data$hash.ident),
                                                n_cluster=length(unique(seurat_table@active.ident)),
                                                names_cluster=unique(seurat_table@active.ident),
                                                n_geno=2,
                                                names_geno=c("wt","hnf4"),
                                                hash_path=seurat_table@meta.data$hash.ident,
                                                cluster_path=seurat_table@active.ident)

dat <- cluster_frac_list[[1]][cluster_frac_list[[1]]$cluster!="Tuft",]

pdf( paste0("figures/cluster_abundance_boxplot_",date,".pdf"), height=4, width=12 )
ggplot( dat, aes( x=cluster, y=frac, fill=geno ) ) +
  geom_boxplot( width=0.4, position=position_dodge( width = 0.75 ), lwd=0.3, outlier.shape = NA ) +
  geom_jitter( aes(), position=position_dodge( width = 0.75 ) ) +
  # stat_summary(geom = "pointrange", aes(color=genotype), shape=3,
  #              fun.data = mean_se, position=position_dodge( width = 0.75 ) ) +
  scale_fill_manual( values=as.vector(cmap[2:3]) ) +
  scale_color_manual( values=cmap[2:3] ) +
  geom_vline(xintercept = seq(1.5,13.5, by=1), lwd=0.25, lty=3 ) +
  theme_bw() +
  theme( axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1), legend.position=c(0.9,0.8) ) +
  xlab("") +
  ylab( "Fraction of total Pdpn-CD104- MEClo")
dev.off()

dat <- cluster_frac_list[[2]][cluster_frac_list[[2]]$cluster%in%c("Tuft","TA","Aire-expressing","Transitional","Aire-adjacent","Immature")==F,]
dat$cluster <- factor(dat$cluster, levels=c( "Entero-hepato",
                                             "Microfold",
                                             "Neuroendocrine",
                                             "Keratinocyte",
                                             "Secretory",
                                             "Ionocyte",
                                             "Ciliated",
                                             "Tuft",
                                             "Muscle",
                                             "Klk-rich",
                                             "TA",
                                             'Transitional',
                                             "Immature",
                                             "Aire-expressing",
                                             "Aire-adjacent"),
                      ordered=T )

pdf( paste0("figures/cluster_abundance_barplot_",date,".pdf"), height=6, width=4 )
ggplot(dat, aes(fill=cluster, y=frac, x=genotype, color="black")) + 
  geom_bar(position="fill", stat="identity", width=0.65) +
  scale_fill_manual( values=cmap ) +
  scale_color_manual( values="black") +
  xlab("Cluster") +
  ylab("Fraction of mimetic cells") +
  theme_few() +
  labs(color="") +
  guides(color="none") +
  theme(legend.position="right")
dev.off()

#####
#save work
# saveRDS(seurat_table, file = "Rdata/meclo_hnf4_seurat_table_v2_20221212.rds")

library(tidyverse)
library(Seurat)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)
library(ggthemes)
library(BuenColors)
library(ggrepel)

setwd( "C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_Hnf4" )
seurat_table <- readRDS( "Rdata/meclo_hnf4_seurat_table_v2_20221212.rds" )

date="20221212"

cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value
names(cmap) <- unique(seurat_table@active.ident)
cmap2 <- ggthemes_data[["stata"]]$colors$schemes$s2color$value
cmap <- cmap2
cmap[14] <- "darkgray"

seurat_table@active.ident <- factor(seurat_table@active.ident, 
                                    levels=c( "Entero-hepato",
                                              "Microfold",
                                              "Neuroendocrine",
                                              "Keratinocyte",
                                              "Secretory",
                                              "Ionocyte",
                                              "Ciliated",
                                              "Tuft",
                                              "Muscle",
                                              "Klk-rich",
                                              "TA",
                                              'Transitional',
                                              "Immature",
                                              "Aire-expressing",
                                              "Aire-adjacent"),
                                    ordered=T )

#####
#viz data

feature <- "Atp6v1b2"

a <- DimPlot(seurat_table, reduction="umap", label=T, repel=T, cols=cmap ) + theme(legend.position="none")
x <- FeaturePlot(seurat_table, features = c(feature), pt.size=1, order=T ) + theme(legend.position="none")
y <- VlnPlot( seurat_table, features=c(feature), cols=cmap ) + theme(legend.position="none")

# pdf( paste0("figures/umap_",feature, "_",date,".pdf"), height=4, width=5 )
# plot_grid(a,x,y, ncol=3)
x
# dev.off()

# FeaturePlot(seurat_table, features = c(feature), pt.size=1, order=T ) + theme(legend.position="none")

#####
#analyze cluster composition
cluster.markers <- FindMarkers(seurat_table, ident.1="14", min.pct=0.1, logfc.threshold=0.5, only.pos=T, test.use="wilcox")
head(cluster.markers, n = 50)

#####
#analyze cluster composition high thru

pdf(paste0("figures/wt_vs_hnf4_wilcox_volcano_cluster_lab_", date ,".pdf"), height=4, width=6)

for ( cluster_feature in unique(seurat_table@active.ident) ) {

  idx1 <- (seurat_table@meta.data$geno.ident %in% c("wt") & seurat_table@active.ident %in% cluster_feature )
  idx2 <- (seurat_table@meta.data$geno.ident %in% c("hnf4") & seurat_table@active.ident %in% cluster_feature )
  
  set.seed(12345)
  
  cluster.markers <- FindMarkers(seurat_table@assays$RNA, 
                                 cells.1=colnames(seurat_table)[idx1], 
                                 cells.2=colnames(seurat_table)[idx2], 
                                 only.pos=F, logfc.threshold=0, min.pct=0, test.use="wilcox")
  
  idx <- abs(cluster.markers$avg_log2FC) > 1 & -log10(cluster.markers$p_val_adj) > 1.3
  label_idx <- ifelse( abs(cluster.markers$avg_log2FC) > 1 & -log10(cluster.markers$p_val_adj) > 1.3, rownames(cluster.markers), "" )
  color_idx <- ifelse( label_idx!="","1","0")
    
  # pdf( paste0("figures/wt_vs_hnf4_microfold_volcano_",date,".pdf"), height=4, width=6 )
  p <- ggplot( cluster.markers, aes(x=avg_log2FC, y=-log10(p_val_adj) ) ) +
    geom_point( aes(color=color_idx) ) +
    geom_point( data=cluster.markers[idx,], 
                aes(x=cluster.markers[idx,"avg_log2FC"], 
                    y=-log10(cluster.markers[idx,"p_val_adj"])), 
                color="#3F007D" ) +
    geom_text_repel(label=label_idx ) +
    theme_few() +
    xlab(paste0("log2 fold change, wt vs hnf4,", cluster_feature)) +
    ylab("-log10 adj p-value") +
    scale_color_manual( values=c("lightgray","#3F007D") ) +
    theme_bw() +
    scale_x_continuous(limits=c(-5,5), oob=squish) +
    scale_y_continuous(limits=c(0,7), oob=squish) +
    theme( legend.position="none" ) +
    ggtitle( paste0("wt vs hnf4, ", cluster_feature) )
  # dev.off()
    
  print(p)
  
  idx <- abs(cluster.markers$avg_log2FC) > 1 & -log10(cluster.markers$p_val) > 1.3
  label_idx <- ifelse( abs(cluster.markers$avg_log2FC) > 1 & -log10(cluster.markers$p_val) > 1.3, rownames(cluster.markers), "" )
  color_idx <- ifelse( label_idx!="","1","0")
  
  # pdf( paste0("figures/wt_vs_hnf4_microfold_volcano_",date,".pdf"), height=4, width=6 )
  q <- ggplot( cluster.markers, aes(x=avg_log2FC, y=-log10(p_val) ) ) +
    geom_point( aes(color=color_idx) ) +
    geom_point( data=cluster.markers[idx,], 
                aes(x=cluster.markers[idx,"avg_log2FC"], 
                    y=-log10(cluster.markers[idx,"p_val"])), 
                color="#3F007D" ) +
    geom_text_repel(label=label_idx ) +
    theme_few() +
    xlab(paste0("log2 fold change, wt vs hnf4,", cluster_feature)) +
    ylab("-log10 nominal p-value") +
    scale_color_manual( values=c("lightgray","#3F007D") ) +
    theme_bw() +
    scale_x_continuous(limits=c(-5,5), oob=squish) +
    scale_y_continuous(limits=c(0,10), oob=squish) +
    theme( legend.position="none" ) +
    ggtitle( paste0("wt vs hnf4, ", cluster_feature) )
  # dev.off()
  
  print(q)

}

dev.off()

#####
#edger approach
library(edgeR)

pdf(paste0("figures/wt_vs_hnf4_edgeR_volcano_cluster_lab_", date ,".pdf"), height=4, width=6)

cluster_list <- unique(seurat_table@active.ident)

for (cluster_feature in cluster_list) {
  
  if ( 
    (sum(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="wt_r1" & seurat_table@active.ident==cluster_feature])!=0) +
    (sum(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="wt_r2" & seurat_table@active.ident==cluster_feature])!=0) +
    (sum(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="wt_r3" & seurat_table@active.ident==cluster_feature])!=0) +
    (sum(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="hnf4_r1" & seurat_table@active.ident==cluster_feature])!=0) +
    (sum(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="hnf4_r2" & seurat_table@active.ident==cluster_feature])!=0) +
    (sum(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="hnf4_r3" & seurat_table@active.ident==cluster_feature])!=0) +
    (sum(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="hnf4_r4" & seurat_table@active.ident==cluster_feature])!=0) == 7 ) {

    dat <- cbind( rowSums(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="wt_r1" & seurat_table@active.ident==cluster_feature]),
                  rowSums(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="wt_r2" & seurat_table@active.ident==cluster_feature]),
                  rowSums(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="wt_r3" & seurat_table@active.ident==cluster_feature]),
                  rowSums(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="hnf4_r1" & seurat_table@active.ident==cluster_feature]),
                  rowSums(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="hnf4_r2" & seurat_table@active.ident==cluster_feature]),
                  rowSums(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="hnf4_r3" & seurat_table@active.ident==cluster_feature]),
                  rowSums(seurat_table@assays$RNA@counts[,seurat_table$hash.ident=="hnf4_r4" & seurat_table@active.ident==cluster_feature]) )
    colnames(dat) <- c("wt_r1","wt_r2","wt_r3","hnf4_r1","hnf4_r2","hnf4_r3","hnf4_r4")
    
    sample_names <- colnames(dat)
    
    i = "hnf4"
    j = "wt"
    
    idx <- substr( colnames(dat), 1, nchar(colnames(dat)) - 3 ) == i
    ctrl_idx <- substr( colnames(dat), 1, nchar(colnames(dat)) - 3 ) == j
    
    group <- factor(c(rep(1,times=sum(idx)),rep(2,times=sum(ctrl_idx))))
    y <- DGEList( counts=cbind( dat[,idx], dat[,ctrl_idx]), 
                  group=group, 
                  genes=rownames(dat) )
    
    keep <- filterByExpr(y, min.count=5)
    y <- y[keep, , keep.lib.sizes=FALSE]
    
    y <- calcNormFactors(y)
    y$samples
    
    design <- model.matrix(~group)
    y <- estimateDisp(y, design)
    
    fit <- glmQLFit(y,design)
    
    qlf <- glmQLFTest(fit,coef=2)
    
    qlf$table <- cbind( qlf$table, p.adjust( qlf$table$PValue, method="BH" ) )
    colnames(qlf$table)[5] = "q_val"
  
    features <- rownames(topTags(qlf, n=25))[qlf$table[rownames(topTags(qlf, n=25)),"q_val"] < 0.05]
    label_idx <- ifelse( rownames(qlf$table) %in% features & qlf$table$q_val < 0.05, rownames(qlf$table), "" )
    color_idx <- ifelse( rownames(qlf$table) %in% features & qlf$table$q_val < 0.05, "1", "0" )
    
    # pdf(paste0("figures/",j,"_vs_", i,"_volcano_lab_", date ,".pdf"), height=3.25, width=4.5)
    p <- ggplot( qlf$table, aes( x=qlf$table$logFC, y=-log10(qlf$table$q_val) ) ) +
      geom_point( aes(color=color_idx) ) +
      geom_point( data=qlf$table[features,], 
                  aes(x=qlf$table[features,"logFC"], 
                      y=-log10(qlf$table[features,"q_val"])), 
                  color="#3F007D" ) +
      geom_text_repel( label=label_idx, color="black", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
      xlab(paste0("log2 fold change, ", j, " vs ", i,", ", cluster_feature)) +
      ylab("-log10 adj p-value") +
      scale_color_manual( values=c("lightgray","#3F007D") ) +
      theme_bw() +
      scale_x_continuous(limits=c(-5,5), oob=squish) +
      scale_y_continuous(limits=c(0,7), oob=squish) +
      theme( legend.position="none" ) +
      ggtitle( paste0(j, " vs ", i,", ", cluster_feature) )
    # dev.off()
    
    print(p)
    
    de_genes <- c("Muc13","Trf","Far2","Gstm1","Apoa4","Ccl6","Alox15","Cobl","Pglyrp1",
                  "Ifitm1","Apoc3","Aoah","Cdhr5","Nostrin","Gabrp","Sbsn","Vil1","Btnl4")
    
    features <- rownames(qlf$table[abs(qlf$table$logFC) > 1 & qlf$table$q_val < 0.05,])
    label_idx <- ifelse( rownames(qlf$table) %in% de_genes, rownames(qlf$table), "" )
    color_idx <- ifelse( rownames(qlf$table) %in% features & qlf$table$logFC > 1, "up", "neutral")
    color_idx <- ifelse( rownames(qlf$table) %in% features & qlf$table$logFC < -1, "down", color_idx)
    
    pdf(paste0("figures/",i,"_vs_", j,"_volcano_lab_", date ,".pdf"), height=3.25, width=4.5)
    q <- ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
      geom_point( aes(color=color_idx) ) +
      geom_point( data=qlf$table[features,], 
                  aes(x=qlf$table[features,"logFC"], 
                      y=-log10(qlf$table[features,"PValue"]),
                      color=color_idx[features]) ) +
      geom_text_repel( label=label_idx, color="black", segment.alpha=0.75, segment.color="black", segment.size=0.25, max.overlaps=100 ) +
      xlab(paste0("log2 fold change, ", i, " vs ", j,", ", cluster_feature)) +
      ylab("-log10 nominal p-value") +
      scale_color_manual( values=c("firebrick3","lightgray","navyblue") ) +
      theme_bw() +
      scale_x_continuous(limits=c(-6,6), oob=squish) +
      scale_y_continuous(limits=c(0,8), oob=squish) +
      theme( legend.position="none" ) +
      ggtitle( paste0(i, " vs ", j,", ", cluster_feature) )
    print(q)
    dev.off()
    
    print(q)
    
  } else { next }

}

dev.off()

#####
#build cluster tree
seurat_table <- BuildClusterTree(seurat_table, dims=T)
PlotClusterTree(seurat_table)

cluster.markers <- FindMarkers(seurat_table, ident.1="10", min.pct=0.1, logfc.threshold=0.5, only.pos=T, test.use="wilcox")
head(cluster.markers, n = 20)

#####
#rename clusters
seurat_table <- RenameIdents( object=seurat_table,
                              "16"="TA",
                              "10"="Transitional",
                              "13"="Immature",
                              "6"="Aire-expressing",
                              "11"="Aire-adjacent",
                              "22"="Klk-rich",
                              "9"="Neuroendocrine",
                              "1"="Tuft",
                              "8"="Keratinocyte",
                              "3"="Secretory",
                              "28"="Ionocyte",
                              "27"="Entero-hepato",
                              "21"="Microfold",
                              "23"="Muscle",
                              "24"="Ciliated")

seurat_table@active.ident <- factor(seurat_table@active.ident, 
                                    levels=c( "TA","Transitional","Immature","Aire-expressing","Aire-adjacent","Tuft",
                                              "Entero-hepato","Microfold","Secretory","Ionocyte","Klk-rich",
                                              "Neuroendocrine","Ciliated","Keratinocyte","Muscle" ) )

#####
#subset post-Aire and reselect variable genes, reduce dimensionality and viz
keep <- seurat_table@active.ident%in%c("Entero-hepato","Microfold","Secretory")
subset_table <- seurat_table[,keep]

keep <- subset_table@active.ident%in%c(7)==F
subset_table <- subset_table[,keep]

subset_table <- NormalizeData(subset_table, normalization.method = "LogNormalize", scale.factor = 10000)
nfeatures = 1000
subset_table <- FindVariableFeatures(subset_table, selection.method = "vst", nfeatures = nfeatures)
all.genes <- rownames(subset_table)
subset_table <- ScaleData(subset_table, features = all.genes)
subset_table <- RunPCA(subset_table, features = VariableFeatures(object = subset_table))
subset_table <- FindNeighbors(subset_table, dims = 1:10)
subset_table <- FindClusters(subset_table, resolution = 1)
subset_table <- RunUMAP(subset_table, dims = 1:10, seed.use = 123)

subset_table$og.ident <- subset_table@active.ident
subset_table <- FindClusters(subset_table, resolution = 1)

subset_table@active.ident <- subset_table$og.ident

DimPlot(subset_table, reduction="umap", label=T, pt.size=2 )
FeaturePlot(subset_table, features=c("Dntt"), order=T, pt.size=2 )
DimPlot(subset_table, reduction="umap", label=T, pt.size=2, split.by="geno.ident" )

#####
#boxplot of cell freq
dat <- c(sum(seurat_table@active.ident=="Entero-hepato" & seurat_table$hash.ident=="wt_r1") /
    sum(seurat_table$hash.ident=="wt_r1"),
  sum(seurat_table@active.ident=="Entero-hepato" & seurat_table$hash.ident=="wt_r2") /
    sum(seurat_table$hash.ident=="wt_r2"),
  sum(seurat_table@active.ident=="Entero-hepato" & seurat_table$hash.ident=="wt_r3") /
    sum(seurat_table$hash.ident=="wt_r3"),
  sum(seurat_table@active.ident=="Entero-hepato" & seurat_table$hash.ident=="hnf4_r1") /
    sum(seurat_table$hash.ident=="hnf4_r1"),
  sum(seurat_table@active.ident=="Entero-hepato" & seurat_table$hash.ident=="hnf4_r2") /
    sum(seurat_table$hash.ident=="hnf4_r2"),
  sum(seurat_table@active.ident=="Entero-hepato" & seurat_table$hash.ident=="hnf4_r3") /
    sum(seurat_table$hash.ident=="hnf4_r3"),
  sum(seurat_table@active.ident=="Entero-hepato" & seurat_table$hash.ident=="hnf4_r4") /
    sum(seurat_table$hash.ident=="hnf4_r4"))

dat <- data.frame(geno=factor(c(rep("wt",times=3),rep("hnf4",times=4)), levels=c("wt","hnf4")),frac=dat)

pdf(paste0("figures/enterohepato_boxplot_hnf4_",date,".pdf"), height=4, width=3)
ggplot(dat, aes(x=geno, y=frac, fill=geno))+
  geom_boxplot(width=0.5) +
  geom_point() +
  theme_bw() +
  scale_fill_manual(values=c("lightgray","lightgray")) +
  theme(legend.position="none", panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Fraction of Pdpn-CD104- mTEClo") +
  ylim(0,0.021)
dev.off()

t.test(dat$frac~dat$geno, alternative="greater")

dat <- c(sum(seurat_table@active.ident=="Microfold" & seurat_table$hash.ident=="wt_r1") /
           sum(seurat_table$hash.ident=="wt_r1"),
         sum(seurat_table@active.ident=="Microfold" & seurat_table$hash.ident=="wt_r2") /
           sum(seurat_table$hash.ident=="wt_r2"),
         sum(seurat_table@active.ident=="Microfold" & seurat_table$hash.ident=="wt_r3") /
           sum(seurat_table$hash.ident=="wt_r3"),
         sum(seurat_table@active.ident=="Microfold" & seurat_table$hash.ident=="hnf4_r1") /
           sum(seurat_table$hash.ident=="hnf4_r1"),
         sum(seurat_table@active.ident=="Microfold" & seurat_table$hash.ident=="hnf4_r2") /
           sum(seurat_table$hash.ident=="hnf4_r2"),
         sum(seurat_table@active.ident=="Microfold" & seurat_table$hash.ident=="hnf4_r3") /
           sum(seurat_table$hash.ident=="hnf4_r3"),
         sum(seurat_table@active.ident=="Microfold" & seurat_table$hash.ident=="hnf4_r4") /
           sum(seurat_table$hash.ident=="hnf4_r4"))

dat <- data.frame(geno=factor(c(rep("wt",times=3),rep("hnf4",times=4)), levels=c("wt","hnf4")),frac=dat)

pdf(paste0("figures/microfold_boxplot_hnf4_",date,".pdf"), height=4, width=3)
ggplot(dat, aes(x=geno, y=frac, fill=geno))+
  geom_boxplot(width=0.5) +
  geom_point() +
  theme_bw() +
  scale_fill_manual(values=c("lightgray","lightgray")) +
  theme(legend.position="none", panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Fraction of Pdpn-CD104- mTEClo") +
  ylim(0,0.03)
dev.off()

t.test(dat$frac~dat$geno, alternative="two.sided")

pdf(paste0("figures/allClusters_boxplot_hnf4_",date,".pdf"), height=4, width=3)
for (i in levels(seurat_table@active.ident)) {
  dat <- c(sum(seurat_table@active.ident==i & seurat_table$hash.ident=="wt_r1") /
             sum(seurat_table$hash.ident=="wt_r1"),
           sum(seurat_table@active.ident==i & seurat_table$hash.ident=="wt_r2") /
             sum(seurat_table$hash.ident=="wt_r2"),
           sum(seurat_table@active.ident==i & seurat_table$hash.ident=="wt_r3") /
             sum(seurat_table$hash.ident=="wt_r3"),
           sum(seurat_table@active.ident==i & seurat_table$hash.ident=="hnf4_r1") /
             sum(seurat_table$hash.ident=="hnf4_r1"),
           sum(seurat_table@active.ident==i & seurat_table$hash.ident=="hnf4_r2") /
             sum(seurat_table$hash.ident=="hnf4_r2"),
           sum(seurat_table@active.ident==i & seurat_table$hash.ident=="hnf4_r3") /
             sum(seurat_table$hash.ident=="hnf4_r3"),
           sum(seurat_table@active.ident==i & seurat_table$hash.ident=="hnf4_r4") /
             sum(seurat_table$hash.ident=="hnf4_r4"))
  
  dat <- data.frame(geno=factor(c(rep("wt",times=3),rep("hnf4",times=4)), levels=c("wt","hnf4")),frac=dat)
  
  p <- ggplot(dat, aes(x=geno, y=frac, fill=geno))+
    geom_boxplot(width=0.5) +
    geom_point() +
    theme_bw() +
    scale_fill_manual(values=c("lightgray","lightgray")) +
    theme(legend.position="none", panel.grid.minor = element_blank()) +
    xlab("") +
    ylab("Fraction of Pdpn-CD104- mTEClo") +
    ggtitle(i) +
    ylim(0,max(dat$frac)+0.05*max(dat$frac))
  
  print(p)
}
dev.off()

#####
#heatmap

dat <- as.matrix(seurat_table[,seurat_table@active.ident%in%"Microfold"]@assays$RNA[,])
sig <- read.delim("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures/20210924_Mcell_sig_logfc1_pct10_qval0.01.txt", header=F, sep="\t")

sig <- sig$V1[sig$V1 %in% rownames(dat)]

dat <- dat[sig,order(seurat_table$geno.ident[seurat_table@active.ident%in%"Microfold"])]

cor_mat <- cor(t(dat))
dist_mat <- as.dist( (1-cor_mat)/2 )
hclust_mat <- hclust(dist_mat, method="ward.D2")

kmeans_mat <- kmeans(cor_mat, centers=2, algorithm="Hartigan-Wong")

dat2 <- log2( (dat+0.01)/(rowMeans(dat[,1:76])+0.01) )
kmeans_mat <- kmeans(dat2, centers=2, algorithm="Hartigan-Wong")

anno_cell <- data.frame("genotype"=seurat_table$geno.ident[seurat_table@active.ident%in%"Microfold"][order(seurat_table$geno.ident[seurat_table@active.ident%in%"Microfold"])])
anno_gene <- data.frame("kmeans"=factor(kmeans_mat$cluster))
breaksList <- seq(-2,0,by=0.1)
anno_cols <- list(genotype = c(wt="navajowhite3",hnf4="firebrick"),
                  kmeans = c("1"="#1B9E77","2"="#D95F02"))

pdf(paste0("figures/mcell_sig_hnf4_heatmap_monocolor_",date,".pdf"), height=10, width=12)
pheatmap(dat2[order(kmeans_mat$cluster,decreasing=T),],
         cluster_rows=F,
         cluster_cols=F,
         annotation_col=anno_cell,
         annotation_row=anno_gene,
         breaks=breaksList,
         # color=colorRampPalette((viridis(n=5, option="magma")))(length(breaksList)),
         color=colorRampPalette((brewer.pal(n=9, name="Purples")))(length(breaksList)),
         border_col=NA,
         show_colnames=F,
         annotation_colors=anno_cols)
dev.off()
