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

setwd("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2022-06-23_hnf4_rnaseq/Batch129_Dan_Mouse")
date = "20220823"

norm_table <- read.delim( "Batch129_Mouse_Dan_1_adjust.gct", sep="\t", header=T, skip=2, row.names=1 )
norm_table <- norm_table[,2:29]
count_table <- read.delim( "Genes_count_table.tsv", sep="\t", header=T, row.names=1 )
count_table <- count_table[,1:28]

sample_names <- colnames(norm_table)
anno_names <- as.data.frame(substr( colnames(count_table[,1:28]), 2, nchar(colnames(count_table[,1:28])) - 8 ))

#####
##QC

#look at counts data
pdf( paste0("figures/counts_per_sample_",date,".pdf"), height=3, width=5 )
hist(log10(colSums(count_table)), breaks=10, xlim=range(0,10))
abline( v=3.3010, col="red" )
dev.off()

#look at normalized data
pdf( paste0("figures/genes_per_sample_",date,".pdf"), height=3, width=5 )
hist(colSums( norm_table > 1 ), breaks=20, xlim=range(0,40000))
dev.off()

pdf( paste0("figures/detection_per_gene_",date,".pdf"), height=3, width=5 )
hist(rowSums( norm_table > 20 ), breaks=20)
abline(v=2, col="red")
dev.off()

#remove genes with normalized expression < threshold
idx <- rowSums( norm_table > 20 ) < 2
norm_table <- norm_table[!idx,]

# calculate basic stats for samples and genes
means <- rowMeans(norm_table)
vars <- apply(norm_table,1,var)
cvs <- apply(norm_table,1,sd)/rowMeans(norm_table)

# plot means vs cvs 
pdf( paste0("figures/cv_vs_mean_",date,".pdf"), height=5, width=5 )
plot(log(means),log(cvs))
dev.off()

# pairs plot
png( paste0("figures/pairs_mtechi_",date,".png"), height=2000, width=2000 )
pairs(log2(norm_table[,1:14]))
dev.off()

png( paste0("figures/pairs_mteclo_",date,".png"), height=2000, width=2000 )
pairs(log2(norm_table[,15:28]))
dev.off()


#####
##initial clustering

#look at sample clustering
sample_cor_matrix <- cor( norm_table, method="pearson" )
hclust_samples <- hclust( as.dist( 1-sample_cor_matrix ), method="ward.D2" )

scaled_norm_table <- log2( norm_table )
scaled_norm_table <- scaled_norm_table - rowMeans(scaled_norm_table)

annotation_names = as.data.frame(sample_names)
rownames( annotation_names ) = rownames( sample_cor_matrix )

breaksList <- seq(0.9,1,by=0.01)

pdf( paste0("figures/samples_distance_heatmap_unclustered_", date, ".pdf"), height=5, width=6)
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
dev.off()

pdf( paste0("figures/samples_distance_heatmap_clustered_", date, ".pdf"), height=10, width=12)
pheatmap( sample_cor_matrix, 
          cluster_rows=hclust_samples,
          cluster_cols=hclust_samples,
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
dev.off()

#####
##pca analysis

#run PCA
pca_object <- pca( t(norm_table), method="svd", nPcs=5, scale="uv", center=T )

#examine elbow plot
pdf( paste0("figures/elbow_plot_", date, ".pdf" ), width=3.5, height=4 )
plot( pca_object@R2, xlab="PC#", ylab="R^2", main="Elbow plot" )
lines( pca_object@R2 )
dev.off()

#examine PC scores
factor_idx <- factor( c(1,1,1,1,2,2,2,3,3,3,4,4,4,4,
                        5,5,5,5,6,6,6,7,7,7,8,8,8,8) )
cmap <- ggthemes_data$tableau$`color-palettes`$regular$`Classic 10`$value

pdf( paste0("figures/PCA_plot_", date, ".pdf" ), height=3.5, width=4 )
for (i in 1:4) {
  temp_plot <- ggplot( as.data.frame(pca_object@scores), aes( x=pca_object@scores[,i],y=pca_object@scores[,i+1] ) ) +
    geom_point( aes( color=factor_idx ), size=4 ) +
    geom_text_repel( label=anno_names[,1] ) +
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
  geom_text_repel( label=anno_names[,1] ) +
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
j = "4aHet.4gHet.mTEChi"
i = "4aKO.4gKO.mTEChi"

idx <- substr( colnames(count_table), 2, nchar(colnames(count_table)) - 8 ) == i
ctrl_idx <- substr( colnames(count_table), 2, nchar(colnames(count_table)) - 8 ) == j

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

qlf.mechi <- qlf

# write_delim( cbind( rownames(qlf$table),qlf$table ), paste0("de_analysis/",i,"_vs_", j,"_DE_table_", date ,".txt"), col_names=T, delim="\t" )

# run edgeR pipeline
j = "4aHet.4gHet.DNmTEClo"
i = "4aKO.4gKO.DNmTEClo"

idx <- substr( colnames(count_table), 2, nchar(colnames(count_table)) - 8 ) == i
ctrl_idx <- substr( colnames(count_table), 2, nchar(colnames(count_table)) - 8 ) == j

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

qlf.meclo <- qlf

# write_delim( cbind( rownames(qlf$table),qlf$table ), paste0("de_analysis/",i,"_vs_", j,"_DE_table_", date ,".txt"), col_names=T, delim="\t" )

# run edgeR pipeline
j = "4aHet.4gHet.mTEChi"
i = "4aHet.4gKO.mTEChi"

idx <- substr( colnames(count_table), 2, nchar(colnames(count_table)) - 8 ) == i
ctrl_idx <- substr( colnames(count_table), 2, nchar(colnames(count_table)) - 8 ) == j

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

qlf.4g <- qlf

topTags(qlf.other, n=20)

write_delim( cbind( rownames(qlf$table),qlf$table ), paste0("de_analysis/",i,"_vs_", j,"_DE_table_", date ,".txt"), col_names=T, delim="\t" )

#plot exp-exp

qlf <- qlf.meclo
i = "4aKO.4gKO.DNmTEClo"
j = "4aHet.4gHet.DNmTEClo"

idx <- substr( colnames(count_table), 2, nchar(colnames(count_table)) - 8 ) == i
ctrl_idx <- substr( colnames(count_table), 2, nchar(colnames(count_table)) - 8 ) == j

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
  scale_x_continuous(limits=c(-6,6), oob=squish) +
  scale_y_continuous(limits=c(0,8), oob=squish) +
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
  scale_x_continuous(limits=c(-6,6), oob=squish) +
  scale_y_continuous(limits=c(0,8), oob=squish) +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()


features <- c("Vil1","Lypd8","Reg1","Reg3b","Reg3g","Ttr","Lyz1","Muc13","Apoa4","Apoc3","Fabp1")
label_idx <- ifelse( rownames(qlf$table) %in% features, rownames(qlf$table), "" )
color_idx <- ifelse( rownames(qlf$table) %in% features, "1", "0" )

pdf(paste0("figures/",i,"_vs_", j,"_volcano_lab_", date ,".pdf"), height=3.25, width=4.5)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_point( data=qlf$table[features,], aes(x=-qlf$table[features,"logFC"], y=-log10(qlf$table[features,"PValue"])), color="#3F007D" ) +
  geom_text_repel( label=label_idx, color="black", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("lightgray","#3F007D") ) +
  theme_bw() +
  scale_x_continuous(limits=c(-6,6), oob=squish) +
  scale_y_continuous(limits=c(0,8), oob=squish) +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()


fc_cutoff = 1
qval_cutoff = 0.05

de_genes <- abs( qlf$table$logFC ) > fc_cutoff & qlf$table$q_val < qval_cutoff
down_sig <- qlf$table$logFC > fc_cutoff & qlf$table$q_val < qval_cutoff
up_sig <- qlf$table$logFC < -fc_cutoff & qlf$table$q_val < qval_cutoff
# 
# write_delim( as.data.frame( rownames(qlf$table)[up_sig] ),
#              paste0("de_analysis/de_genes/",i,"_vs_",j,"_UP_de_genes_logfc",fc_cutoff,"_pval",qval_cutoff,".txt"),
#              col_names=F,
#              delim="\t" )
# 
# write_delim( as.data.frame( rownames(qlf$table)[down_sig] ),
#              paste0("de_analysis/de_genes/",i,"_vs_",j,"_DOWN_de_genes_logfc",fc_cutoff,"_pval",qval_cutoff,".txt"),
#              col_names=F,
#              delim="\t" )

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

pdf(paste0("figures/",i,"_vs_", j,"_volcano_up_down_sig_logfc",fc_cutoff,"_qval",qval_cutoff,"_", date ,".pdf"), height=3.25, width=4.5)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( aes(color=color_idx) ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("navyblue","lightgray","firebrick3") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) ) +
  scale_x_continuous(limits=c(-6,6), oob=squish) +
  scale_y_continuous(limits=c(0,8), oob=squish) +
  annotate( geom="text", label=sum(up_sig), x=6, y=7, color="firebrick3" ) +
  annotate( geom="text", label=sum(down_sig), x=-6, y=7, color="navyblue" )
dev.off()

#plot signatures

path_name="C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures"
idx <- list.files(path=path_name)[grepl(".txt", list.files(path=path_name))]

# qlf <- qlf.meclo
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

#####
#make mean gut sig expression

meanExpTable <- colMeans(log2(norm_table[signature,15:28]))
factor_idx <- factor( c("4aH\n4gH","4aH\n4gH","4aH\n4gH","4aH\n4gH",
                        "4aH\n4gK","4aH\n4gK","4aH\n4gK",
                        "4aK\n4gH","4aK\n4gH","4aK\n4gH",
                        "4aK\n4gK","4aK\n4gK","4aK\n4gK","4aK\n4gK"), 
                      levels=c("4aH\n4gH","4aH\n4gK","4aK\n4gH","4aK\n4gK"), ordered=T )

pdf("figures/20221212_enterohepato_sig_boxplot.pdf", height=5, width=3.5)
ggplot(as.data.frame(meanExpTable), aes(x=factor_idx,y=meanExpTable, fill=factor_idx)) +
  geom_boxplot(width=0.5, alpha=0.5) +
  geom_point() +
  theme_bw() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none" ) +
  xlab("Genotype") +
  ylab("Entero-hepato signature\nexpression (log2)") +
  scale_fill_manual(values=c("darkgray","steelblue","firebrick3","navyblue"))
dev.off()

#####
#gprofiler
gprofiler <- read.delim( "de_analysis/gProfiler_mmusculus_12-12-2022_9-57-46 AM__intersections.csv", header=T, sep="," )

head(gprofiler)
gprofiler <- gprofiler[gprofiler[,1]=="GO:BP",]
gprofiler$term_name <- factor(gprofiler$term_name, levels=rev(c(gprofiler$term_name)))

pdf("figures/4aK4gK_down_DESeq2_pathway_full_20221212.pdf", height=12, width=9)
ggplot(gprofiler, aes(x=negative_log10_of_adjusted_p_value, y=term_name)) +
  geom_col(color="black", fill="gray") +
  theme_few()
dev.off()
pdf("figures/4aK4gK_down_DESeq2_pathway_select_20221212.pdf", height=2, width=4.5)
ggplot(gprofiler[c(3,13,15,16,20),], aes(x=negative_log10_of_adjusted_p_value, y=term_name)) +
  geom_col(color="black", fill="gray") +
  theme_few()
dev.off()
