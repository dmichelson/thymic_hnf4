#####
#start
library(tidyverse)
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(pcaMethods)
library(ggrepel)
library(scales)
library(ggthemes)

setwd("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2022-06-23_hnf4_rnaseq/Batch130_Dan_Mouse")
date = "20220909"

count_table <- read.delim( "Genes_count_table.tsv", sep="\t", header=T, row.names=1 )
count_table <- count_table[,1:6]
count_table2 <- read.delim("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2020-08-21_mcell_meclo_rnaseq/Batch096_Michelson_Mouse/Genes_count_table.tsv", sep="\t", header=T, row.names=1)
count_table2 <- count_table2[,10:12]
count_table3 <- read.delim("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2021-04-02_foxj1_thy_lung_rnaseq/Batch110_Dan_Michelson_Mouse/Genes_count_table.tsv", sep="\t", header=T, row.names=1)
count_table3 <- count_table3[,10:12]

count_table_merge <- cbind(count_table, count_table2, count_table3)

sample_names <- colnames(count_table_merge)
anno_names <- vector(length=12)
anno_names <- substr(sample_names,1,nchar(sample_names)-13)
anno_names[7:9] <- substr(sample_names[7:9],1,nchar(sample_names[7:9])-6)
anno_names[10:12] <- substr(sample_names[10:12],1,nchar(sample_names[10:12])-7)

group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4))
y <- DGEList( counts=count_table_merge, 
              group=group, 
              genes=rownames(count_table) )

keep <- filterByExpr(y, min.count=5)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

design <- model.matrix(~group)
y <- estimateDisp(y, design)

norm_table <- cpm(y, normalized.lib.sizes = T, log=F)

#####
##QC

#remove genes with normalized expression < threshold
idx <- rowSums( norm_table > 5 ) < 2
norm_table <- norm_table[!idx,]


#####
##initial clustering

#look at sample clustering
sample_cor_matrix <- cor( norm_table, method="pearson" )
hclust_samples <- hclust( as.dist( 1-sample_cor_matrix ), method="ward.D2" )

scaled_norm_table <- log2( norm_table+1 )
scaled_norm_table <- scaled_norm_table - rowMeans(scaled_norm_table)

annotation_names = as.data.frame(sample_names)
rownames( annotation_names ) = rownames( sample_cor_matrix )

breaksList <- seq(0.7,1,by=0.01)

# pdf( paste0("figures/samples_distance_heatmap_unclustered_", date, ".pdf"), height=4, width=6)
pheatmap( sample_cor_matrix, 
          cluster_rows=F,
          cluster_cols=F,
          show_annotation_row=T,
          labels_row = as.vector(sample_names), 
          show_rownames=T, 
          show_colnames=F, 
          main="Samples distance matrix", 
          annotation_legend = F, 
          scale="none",
          breaks=breaksList, 
          color=colorRampPalette((brewer.pal(name="Reds",n=11)))(length(breaksList)),
          border_color = NA
)
# dev.off()


#####
##pca analysis

#run PCA
pca_object <- pca( t(scaled_norm_table), method="svd", nPcs=5, scale="uv", center=T )

#examine elbow plot
pdf( paste0("figures/elbow_plot_", date, ".pdf" ), width=3.5, height=4 )
plot( pca_object@R2, xlab="PC#", ylab="R^2", main="Elbow plot" )
lines( pca_object@R2 )
dev.off()

#examine PC scores
factor_idx <- factor( c(1,1,1,2,2,2,3,3,3,4,4,4) )
cmap <- ggthemes_data$tableau$`color-palettes`$regular$`Classic 10`$value

pdf( paste0("figures/PCA_plot_", date, ".pdf" ), height=3.5, width=4 )
for (i in 1:4) {
  temp_plot <- ggplot( as.data.frame(pca_object@scores), aes( x=pca_object@scores[,i],y=pca_object@scores[,i+1] ) ) +
    geom_point( aes( color=factor_idx ), size=4 ) +
    geom_text_repel( label=anno_names ) +
    theme_classic() +
    theme( legend.position="none" ) +
    xlab( paste0("PC",i, " (",pca_object@R2[i]*100,"%)") ) +
    ylab( paste0("PC",i+1, " (",pca_object@R2[i+1]*100,"%)") ) +
    scale_color_manual(values=cmap ) +
    theme( panel.border = element_rect(colour = "black", fill=NA, size=1) )
  print( temp_plot )
}
dev.off()

pdf( paste0("figures/PCA_plot_bespoke_", date, ".pdf" ), height=6, width=8 )
i=1
temp_plot <- ggplot( as.data.frame(pca_object@scores), aes( x=pca_object@scores[,i],y=pca_object@scores[,i+1] ) ) +
  geom_point( aes( color=factor_idx ), size=4 ) +
  geom_text_repel( label=anno_names ) +
  theme_classic() +
  theme( legend.position="none" ) +
  xlab( paste0("PC",i, " (",pca_object@R2[i]*100,"%)") ) +
  ylab( paste0("PC",i+1, " (",pca_object@R2[i+1]*100,"%)") ) +
  scale_color_manual(values=cmap ) +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1) )
print( temp_plot )
dev.off()

#find top loadings for PCs
pca_genes <- matrix( ncol=5, nrow=50 )
colnames( pca_genes ) = c("PC1","PC2","PC3","PC4","PC5")
for (i in 1:5) {
  temp <- sort( pca_object@loadings[, paste0("PC",i) ], decreasing=T )
  pca_genes[,i] = names( temp[1:50] )
}

#write top loadings to file
write_delim( as.data.frame( pca_genes ), paste0( "de_analysis/top_pca_genes_per_component_", date, ".txt" ), delim="\t", col_names=T )

#find bottom loadings for PCs
pca_genes <- matrix( ncol=5, nrow=50 )
colnames( pca_genes ) = c("PC1","PC2","PC3","PC4","PC5")
for (i in 1:5) {
  temp <- sort( pca_object@loadings[, paste0("PC",i) ], decreasing=F )
  pca_genes[,i] = names( temp[1:50] )
}

#write bottom loadings to file
write_delim( as.data.frame( pca_genes ), paste0( "de_analysis/bottom_pca_genes_per_component_", date, ".txt" ), delim="\t", col_names=T )

#####
##de analysis

# run edgeR pipeline
i = "EnteroHepato"
j = "Enterocyte"

idx <- substr( colnames(count_table_merge), 1, nchar(colnames(count_table_merge)) - 13 ) == i
# ctrl_idx <- substr( colnames(count_table_merge), 1, nchar(colnames(count_table_merge)) - 13 ) == j

# idx <- substr( colnames(count_table_merge), 1, 6 ) == i
ctrl_idx <- substr( colnames(count_table_merge), 1, 10 ) == j

group <- factor(c(rep(1,times=sum(idx)),rep(2,times=sum(ctrl_idx))))
y <- DGEList( counts=cbind( count_table_merge[,idx], count_table_merge[,ctrl_idx]), 
              group=group, 
              genes=rownames(count_table_merge) )

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

qlf.thy <- qlf

#####
##de analysis

# run edgeR pipeline
i = "Enterocyte"
j = "AT2"

count_table <- count_table_merge

idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 6 ) == i
ctrl_idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 7 ) == j

# idx <- substr( colnames(count_table), 1, 6 ) == i
# ctrl_idx <- substr( colnames(count_table), 1, 6 ) == j

group <- factor(c(rep(1,times=sum(idx)),rep(2,times=sum(ctrl_idx))))
y <- DGEList( counts=cbind( count_table[,idx], count_table[,ctrl_idx]), 
              group=group, 
              genes=rownames(count_table) )

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

qlf.gut <- qlf

#####
#fc fc
keep <- rownames(qlf.thy) %in% rownames(qlf.gut)
qlf.thy.shared <- qlf.thy[keep,]
keep <- rownames(qlf.gut) %in% rownames(qlf.thy)
qlf.gut.shared <- qlf.gut[keep,]

library(Hmisc)
rcorr(-qlf.thy.shared$table[,"logFC"], -qlf.gut.shared$table[,"logFC"])

pdf(paste0("figures/gut_vs_thymus_fc-fc_unlabeled_", date ,".pdf"), height=4, width=4.75)
ggplot( qlf.thy.shared$table, aes( x=-qlf.thy.shared$table$logFC, y=-qlf.gut.shared$table$logFC ) ) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  geom_abline(slope=1,intercept=0) +
  geom_point() +
  theme_bw() +
  theme( legend.position="none" ) +
  # geom_text_repel( label=label_idx, color="red", segment.alpha=0.75, segment.color="gray", segment.size=0.25  ) +
  scale_color_manual( values=c("black") ) +
  xlab( "log2 fold change, Entero-hepato vs mTEC" ) +
  ylab( "log2 fold change, Gut.Enterocyte vs Lung.AT2") +
  scale_x_continuous(limits=c(-20,20), oob=squish) +
  scale_y_continuous(limits=c(-20,20), oob=squish) +
  annotate(geom="text", label="r=0.36", x=15, y=-15) +
  ggtitle( "Gut vs Thymus" )
dev.off()

keep <- c("Apoa4","Apoc3","Reg3g","Muc13","Hnf4a","Hnf4g","Slc51b","Reg3g","Reg3b","Vil1","Guca2a","Guca2b","Ces2e","Ces2f")
label_idx <- ifelse( rownames(qlf.thy.shared$table) %in% keep, rownames(qlf.thy.shared$table), "" )
color_idx <- ifelse( rownames(qlf.thy.shared$table) %in% keep, "1", "0" )

pdf(paste0("figures/gut_vs_thymus_fc-fc_labeled_", date ,".pdf"), height=4, width=4.75)
ggplot( qlf.thy.shared$table, aes( x=-qlf.thy.shared$table$logFC, y=-qlf.gut.shared$table$logFC ) ) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  geom_abline(slope=1,intercept=0) +
  geom_point( color="lightgray") +
  geom_point( data=qlf.thy.shared$table[label_idx!="",], 
              aes( x=-qlf.thy.shared$table$logFC[label_idx!=""], y=-qlf.gut.shared$table$logFC[label_idx!=""] ), color="#3F007D" ) +
  theme_bw() +
  theme( legend.position="none" ) +
  geom_text_repel( label=label_idx, color="black", segment.alpha=0.75, segment.color="black", segment.size=0.25  ) +
  scale_color_manual( values=c("black") ) +
  xlab( "log2 fold change, Entero-hepato vs mTEC" ) +
  ylab( "log2 fold change, Gut.Enterocyte vs Lung.AT2") +
  scale_x_continuous(limits=c(-20,20), oob=squish) +
  scale_y_continuous(limits=c(-20,20), oob=squish) +
  ggtitle( "Gut vs Thymus" )
dev.off()

keep <- c("Chga","Chgb","Clca3b","Gip")
label_idx <- ifelse( rownames(qlf.thy.shared$table) %in% keep, rownames(qlf.thy.shared$table), "" )
color_idx <- ifelse( rownames(qlf.thy.shared$table) %in% keep, "1", "0" )

pdf(paste0("figures/gut_vs_thymus_fc-fc_labeled_extra_", date ,".pdf"), height=4, width=4.75)
ggplot( qlf.thy.shared$table, aes( x=-qlf.thy.shared$table$logFC, y=-qlf.gut.shared$table$logFC ) ) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  geom_abline(slope=1,intercept=0) +
  geom_point( color="lightgray") +
  geom_point( data=qlf.thy.shared$table[label_idx!="",], 
              aes( x=-qlf.thy.shared$table$logFC[label_idx!=""], y=-qlf.gut.shared$table$logFC[label_idx!=""] ), color="#3F007D" ) +
  theme_bw() +
  theme( legend.position="none" ) +
  geom_text_repel( label=label_idx, color="black", segment.alpha=0.75, segment.color="black", segment.size=0.25  ) +
  scale_color_manual( values=c("black") ) +
  xlab( "log2 fold change, Entero-hepato vs mTEC" ) +
  ylab( "log2 fold change, Gut.Enterocyte vs Lung.AT2") +
  scale_x_continuous(limits=c(-20,20), oob=squish) +
  scale_y_continuous(limits=c(-20,20), oob=squish) +
  ggtitle( "Gut vs Thymus" )
dev.off()


# write_delim( cbind( rownames(qlf$table),qlf$table ), paste0("de_analysis/",i,"_vs_", j,"_DE_table_", date ,".txt"), col_names=T, delim="\t" )

#plot exp-exp

label_idx <- ifelse( rownames(norm_table) %in% rownames(topTags(qlf, n=50)) &
                       !grepl("Rik", rownames(norm_table)), rownames(norm_table), "" )
color_idx <- ifelse( rownames(norm_table) %in% rownames(topTags(qlf, n=50)), "1", "0" )

# keep <- rowMeans(log2(norm_table[,idx])) > 0.95*rowMeans(log2(norm_table[,ctrl_idx])) + 3
# label_idx <- ifelse( keep==T, rownames(norm_table), "" )
# color_idx <- ifelse( keep==T, "1", "0" )

pdf(paste0("figures/",i,"_vs_", j,"_exp-exp_unlab", date ,".pdf"), height=8, width=9)
ggplot( norm_table, aes( x=rowMeans(log2(norm_table[,ctrl_idx])), y=rowMeans(log2(norm_table[,idx])) ) ) +
  geom_point() +
  # geom_text_repel( label=label_idx, color="red", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab( paste0( "log2 expression, ",j ) ) +
  ylab( paste0( "log2 expression, ",i ) ) +
  scale_color_manual( values=c("black","red") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  xlim(0,16) +
  ylim(0,16) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()

pdf(paste0("figures/",i,"_vs_", j,"_exp-exp_lab", date ,".pdf"), height=8, width=9)
ggplot( norm_table, aes( x=rowMeans(log2(norm_table[,ctrl_idx])), y=rowMeans(log2(norm_table[,idx])) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_text_repel( label=label_idx, color="red", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab( paste0( "log2 expression, ",j ) ) +
  ylab( paste0( "log2 expression, ",i ) ) +
  scale_color_manual( values=c("black","red") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  xlim(0,16) +
  ylim(0,16) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()

#plot volcano
label_idx <- ifelse( rownames(qlf$table) %in% rownames(topTags(qlf, n=50)), rownames(qlf$table), "" )
color_idx <- ifelse( rownames(qlf$table) %in% rownames(topTags(qlf, n=50)), "1", "0" )


pdf(paste0("figures/",i,"_vs_", j,"_volcano_unlab_", date ,".pdf"), height=6, width=8)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( ) +
  # geom_point( data=qlf$table[features,], aes(x=qlf$table[features,"logFC"], y=-log10(qlf$table[features,"PValue"])), color="red" ) +
  # geom_text_repel( label=label_idx, color="red", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("black","red") ) +
  theme_bw() +
  # scale_x_continuous(limits=c(-9,9), oob=squish) +
  scale_y_continuous(limits=c(0,10), oob=squish) +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()

pdf(paste0("figures/",i,"_vs_", j,"_volcano_lab_", date ,".pdf"), height=6, width=8)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( aes(color=color_idx) ) +
  # geom_point( data=qlf$table[features,], aes(x=qlf$table[features,"logFC"], y=-log10(qlf$table[features,"PValue"])), color="red" ) +
  geom_text_repel( label=label_idx, color="red", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("black","red") ) +
  theme_bw() +
  # scale_x_continuous(limits=c(-9,9), oob=squish) +
  scale_y_continuous(limits=c(0,10), oob=squish) +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()

fc_cutoff = 1
qval_cutoff = 0.05

de_genes <- abs( qlf$table$logFC ) > fc_cutoff & qlf$table$q_val < qval_cutoff
down_sig <- qlf$table$logFC > fc_cutoff & qlf$table$q_val < qval_cutoff
up_sig <- qlf$table$logFC < -fc_cutoff & qlf$table$q_val < qval_cutoff

write_delim( as.data.frame( rownames(qlf$table)[up_sig] ),
             paste0("de_analysis/de_genes/",i,"_vs_",j,"_UP_de_genes_logfc",fc_cutoff,"_qval",qval_cutoff,".txt"),
             col_names=F,
             delim="\t" )

write_delim( as.data.frame( rownames(qlf$table)[down_sig] ),
             paste0("de_analysis/de_genes/",i,"_vs_",j,"_DOWN_de_genes_logfc",fc_cutoff,"_qval",qval_cutoff,".txt"),
             col_names=F,
             delim="\t" )

#plots with de_gene overlay
# label_idx <- ifelse( rownames(qlf$table) %in% rownames(topTags(qlf, n=50)), "" )

color_idx <- vector( length=nrow(norm_table) )
for (k in 1:nrow(norm_table)) {
  if (rownames(norm_table)[k] %in% rownames(qlf$table)[up_sig]){
    color_idx[k] <- "up"
  } else if (rownames(norm_table)[k] %in% rownames(qlf$table)[down_sig]) {
    color_idx[k] <- "down"
  } else {
    color_idx[k] <- "neutral"
  }
}

pdf(paste0("figures/",i,"_vs_", j,"_exp-exp_up_down_sig_logfc",fc_cutoff,"_qval",qval_cutoff,"_", date ,".pdf"), height=5, width=5)
ggplot( norm_table, aes( x=rowMeans(log2(norm_table[,ctrl_idx])), y=rowMeans(log2(norm_table[,idx])) ) ) +
  geom_point( aes(color=color_idx) ) +
  xlab( paste0( "log2 expression, ",j ) ) +
  ylab( paste0( "log2 expression, ",i ) ) +
  scale_color_manual( values=c("red","black","blue") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  xlim(0,16) +
  ylim(0,16) +
  ggtitle( paste0(i, " vs ", j) ) +
  annotate( geom="text", label=sum(up_sig), x=2, y=15, color="blue" ) +
  annotate( geom="text", label=sum(down_sig), x=15, y=2, color="red" )
dev.off()

color_idx <- vector( length=nrow(qlf$table) )
for (k in 1:nrow(qlf$table)) {
  if (rownames(qlf$table)[k] %in% rownames(qlf$table)[up_sig]){
    color_idx[k] <- "up"
  } else if (rownames(qlf$table)[k] %in% rownames(qlf$table)[down_sig]) {
    color_idx[k] <- "down"
  } else {
    color_idx[k] <- "neutral"
  }
}

pdf(paste0("figures/",i,"_vs_", j,"_volcano_up_down_sig_logfc",fc_cutoff,"_qval",qval_cutoff,"_", date ,".pdf"), height=4, width=6)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( aes(color=color_idx) ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("blue","black","red") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) ) +
  scale_x_continuous(limits=c(-9,9), oob=squish) +
  scale_y_continuous(limits=c(0,10), oob=squish) +
  annotate( geom="text", label=sum(up_sig), x=6, y=7, color="red" ) +
  annotate( geom="text", label=sum(down_sig), x=-6, y=7, color="blue" )
dev.off()

#analyze signatures
path_name="C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures"
idx <- list.files(path=path_name)[grepl(".txt", list.files(path=path_name))]

# qlf <- qlf.gut
# i = "4aHet.4gHet.DN.mTEClo"
# j = "4aKO.4gKO.DN.mTEClo"

pdf( paste0("figures/",i,"vs",j,"_allSigs_",date,".pdf"), height=5, width=6 )
for (k in idx ) {
  sig <- read.delim(paste0(path_name,"/",k),header=T, sep="\t")
  sig <- sig[,1]
  sig <- sig[sig%in%rownames(qlf)]
  sig_name <- substr(k, 10, nchar(k)-30)
  
  color_idx <- ifelse( rownames(qlf) %in% sig, "1", "0" )
  
  p <- ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
    geom_hline( yintercept = 0 ) +
    geom_vline( xintercept = 0 ) +
    geom_point( aes( color=color_idx ) ) +
    geom_point( data=as.data.frame(qlf$table[sig,]), 
                aes( x=-qlf$table[sig,"logFC"], y=-log10(qlf$table[sig,"PValue"]) ), color="purple" ) +
    theme_bw() +
    xlab(paste0("log2FC, ",i, " vs ", j)) +
    ylab("-log10 PValue") +
    theme(legend.position="none") +
    scale_color_manual( values=c("gray","purple") ) +
    ggtitle(sig_name) +
    annotate( geom="text", label=sum(color_idx=="1"&qlf$table$logFC<0), x=3, y=1, color="purple" ) +
    annotate( geom="text", label=sum(color_idx=="1"&qlf$table$logFC>0), x=-3, y=1, color="purple" )
  print(p)
}
dev.off()
