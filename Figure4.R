#Import dependencies
library(tidyverse)
library(rhdf5)
library(Seurat)
library(cowplot)
library(ggpubr)

#Set input/output directories
input_dir <- "data/"
output_dir <- "results/Figure4/"

#Figure 4A
##Reconstruct H & W matrices
genes <- read.delim("results/Subclustering_B/NMF/sa_input_genes_tumor.tsv", sep="\t")[,1]
H <- h5read(file = "results/Subclustering_B/NMF/NMF_Tumor_Poisson_HL1_WL1_n30/nmf_output.h5", name = "H")
H_df <- data.frame(t(H$block0_values))
rownames(H_df) <- gsub("^X", "", H$axis1)
rownames(H_df) <- gsub("\\.", "-", rownames(H_df))
colnames(H_df) <- H$block0_items
H_df$max_id <- data.frame(t(H$block1_values))[,1]
W <- h5read(file = "results/Subclustering_B/NMF/NMF_Tumor_Poisson_HL1_WL1_n30/nmf_output.h5", name = "W")
W_df <- data.frame(t(W$block0_values))
rownames(W_df) <- genes
colnames(W_df) <- W$block0_items
W_df$max_id <- data.frame(t(W$block1_values))[,1]
write.csv(H_df, paste0(output_dir, "H.csv"), row.names = T)
write.csv(W_df, paste0(output_dir, "W.csv"), row.names = T)

##Obtain markers
#UE = W*H.sum()
#FE = (UE.T/UE.T.sum()).T
#R=W*FE
nsigs <- length(grep("^S", colnames(W_df)))
UE <- data.frame(t(t(W_df[,1:nsigs])*colSums(H_df[,1:nsigs])))
FE <- data.frame(UE/rowSums(UE))
R <- W_df[,1:nsigs]*FE
write.csv(R, paste0(output_dir, "Markers.csv"), row.names = T)

##Select top 5 markers
R$Assignment <- apply(R, 1, function(x) {colnames(R)[which.max(x)]})
sigs <- colnames(W_df)[grep("^S", colnames(W_df))]
R <- R[-grep("^AC.*\\.|^AL.*\\.|RPS*|RPL*|^MT*", rownames(R)),]
for(i in 1:nsigs){
  tmp <- R[R$Assignment == sigs[i],]
  tmp <- tmp[order(tmp[,i], decreasing = T),]
  if(i==1){
    markers <- rownames(tmp)[1:5]
    out <- data.frame("Signature" = rep(sigs[i], length(markers)), "Gene" = markers)
  } else {
    markers <- rownames(tmp)[1:5]
    out <- rbind(out, data.frame("Signature" = rep(sigs[i], length(markers)), "Gene" = markers))
  }
}
write.csv(out, paste0(output_dir, "Subtype_markers_for_validation.csv"), row.names = F)

##Visualize signature activity across tumors in UMAP embedding
##Read in WM tumor object
wm <- readRDS("results/Subclustering_B/WMonly/integrated.rds")
##Extract UMAP coordinates and merge with H matrix
coord <- data.frame(wm[["umap"]]@cell.embeddings)
meta <- wm@meta.data
meta <- transform(merge(meta, coord, by="row.names"), row.names=Row.names, Row.names=NULL)
H_df <- read.csv(paste0(output_dir, "H.csv"), row.names = 1)
nsigs <- length(grep("^S", colnames(H_df)))
H_df$Cell_ID  <- rownames(H_df)
meta$Cell_ID <-meta$NewCellID
h_clust <- transform(merge(H_df, meta, by ="Cell_ID"), row.names=Cell_ID, Cell_ID=NULL)
write.csv(h_clust, paste0(output_dir, "Merged_Data.csv"), row.names = T)
nrow(h_clust)
##48,875
h_clust_long <- h_clust[, c("UMAP_1", "UMAP_2", "NewCellID", "CaTissueID", c(paste0("S", 1:nsigs)))]
h_clust_long <- gather(h_clust_long, "Signature", "Expression", -UMAP_1, -UMAP_2, -NewCellID, -CaTissueID)
h_clust_long$Signature <- factor(gsub("S", "GEX-", h_clust_long$Signature), levels = paste0("GEX-", 1:nsigs))
png(paste0(output_dir, "NMF_Tumor_Poisson_HL1_WL1_n30_GEXSignature_UMAP.png"), res = 300, units="in", width = 12, height = 10)
ggplot(h_clust_long) + 
  geom_point(aes(UMAP_1, UMAP_2, color=log10(Expression+1)), show.legend = T) + 
  scale_color_gradient(low="lightyellow2", high="midnightblue", name=expression("Log"[10]*"(GEX+1)")) +
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 18, face = "italic"), 
        legend.position = "bottom", 
        legend.text = element_text(size =16), 
        legend.title = element_text(size = 16)) +
  facet_wrap(.~Signature, ncol = 4)
dev.off()

##Plot signature dictionary
out <- read.csv(paste0(output_dir, "Subtype_markers_for_validation.csv"))
markers <- out$Gene
sigs <- paste0("S", 1:nsigs)
X <- GetAssayData(wm[markers,], slot = "data")
marker_table <- data.frame(t(data.frame(X)))
rownames(marker_table) <- gsub("^X","",rownames(marker_table))
rownames(marker_table) <- gsub("\\.","-",rownames(marker_table))
marker_table$max_id <- H_df$max_id[match(rownames(marker_table), H_df$Cell)]
colnames(marker_table) <- gsub("\\.","-",colnames(marker_table))
marker_means <- aggregate(scale(marker_table[, 1:(ncol(marker_table)-1)]), list(marker_table$max_id), mean)
colnames(marker_means)[1] <- "Signature"
marker_means$Signature <- paste0("S",marker_means$Signature)
marker_means_long <- gather(marker_means, "Gene","Exp",-Signature)
marker_means_long$Signature <- factor(marker_means_long$Signature, levels=sigs)
marker_means_long$Gene <- factor(marker_means_long$Gene, levels=markers)
marker_means_long <- marker_means_long[as.character(marker_means_long$Signature) %in% sigs,]
gene_num <- rep(5, nsigs)
for(i in seq_along(gene_num)){
  if(i==1){
    out <- c(gene_num[i])
  } else {
    out <- c(out, out[i-1]+gene_num[i])
  }
}
segment_data_top <- data.frame("x"=c(0,out[-length(out)])+0.5, "x_end"=out+0.5, "y" = c(1:nsigs)+0.5, "y_end" = c(1:nsigs)+0.5)
segment_data_bottom <- data.frame("x"=c(0,out[-length(out)])+0.5, "x_end"=out+0.5, "y" = c(0:(nsigs-1))+0.5, "y_end" = c(0:(nsigs-1))+0.5)
segment_data_left <- data.frame("x"=c(0,out[-length(out)])+0.5, "x_end"=c(0,out[-length(out)])+0.5, "y" = c(0:(nsigs-1))+0.5, "y_end" = c(1:nsigs)+0.5)
segment_data_right <- data.frame("x"=out+0.5, "x_end"=out+0.5,  "y" = c(0:(nsigs-1))+0.5, "y_end" = c(1:nsigs)+0.5)
png(paste0(output_dir, "NMF_Tumor_Poisson_HL1_WL1_n30_SignatureDictionary.png"), res = 300, units="in", width = 7, height = 5)
ggplot(marker_means_long[(!is.na(marker_means_long$Gene)),]) +
  geom_tile(aes(Gene,Signature,fill=Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Exp),0,max(marker_means_long$Exp))))+ 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill=NA, color="black"), 
        axis.text.x = element_text(angle = 65, hjust=1, color="black", size=12, face="italic"), 
        axis.text.y=element_text(size=12, color="black"), 
        legend.text = element_text(size=16), 
        legend.key.size = unit(1, "cm"), 
        legend.title = element_text(size=18), 
        legend.position = "top") + 
  xlab("") + 
  labs(fill="Mean Z-score") + 
  ylab("")+ 
  scale_y_discrete(labels=gsub("S","GEX-",marker_means_long$Signature)) +
  geom_segment(data=segment_data_top, aes(x=x, xend=x_end, y=y, yend=y_end), linetype="dashed", size=0.8) + 
  geom_segment(data=segment_data_bottom, aes(x=x, xend=x_end, y=y, yend=y_end), linetype="dashed", size=0.8) +
  geom_segment(data=segment_data_left, aes(x=x, xend=x_end, y=y, yend=y_end), linetype="dashed", size=0.8) +
  geom_segment(data=segment_data_right, aes(x=x, xend=x_end, y=y, yend=y_end), linetype="dashed", size=0.8)
dev.off()

#Figure 4B
##Plot heatmap of signature markers in GSE171739
##Import GSE171739 data and signature markers
data <- read.csv(paste0(input_dir, "GSE171739_allgenes_bcells.csv"))
markers <- read.csv(paste0(output_dir, "Subtype_markers_for_validation.csv"))
##Remove signature S4, which corresponds to monocyte-derived doublets and/or ambient, and S1 which is a general activation signature
markers <- markers[!markers$Signature %in% c("S4", "S1"),]
##ACTG1 == ACTB in this dataset, as the same probe corresponds to both genes
data$ACTG1 <- data$ACTB
##Subset data for markers
present <- as.character(markers$Gene)[markers$Gene %in% colnames(data)]
tmp <- data[, c("Sample", "Stage", present)]
##Scale gene expression
tmp_scaled <- cbind(tmp[,c("Sample", "Stage")], scale(tmp[, 3:ncol(tmp)]))
##Perform hierarchical clustering
hc <- tmp_scaled
rownames(hc) <- hc$Sample
hc$Stage <- NULL
res <- hclust(dist(hc, method = "euclidean"), method = "complete")
tmp$subtypes <- cutree(res, k=length(unique(markers$Signature)))
sample_order <- as.character(tmp$Sample)[order(tmp$subtypes, decreasing = F)]
##Plot heatmap
tmp_long <- gather(tmp_scaled, "Gene", "Expression", -Sample, -Stage)
markers <- markers[markers$Gene %in% present,]
tmp_long$Gene <- factor(tmp_long$Gene, levels = markers$Gene)
tmp_long$Sample <- factor(tmp_long$Sample, levels = sample_order)
segment.data.vertical <- data.frame("x" = seq(5, 5*5, 5) + 0.5, "x_end" = seq(5, 5*5, 5) +0.5, "y" = rep(0, 5), "y_end" = rep(nrow(tmp) + 1, 5))
segment.data.horizontal <- data.frame("x" = c(0, 5, 15) + 0.5, "x_end" = c(15, 30, 30) +0.5, "y" = c(10, 32, 41) + 0.5, "y_end" = c(10, 32, 41) + 0.5)
limits <- data.frame(table(tmp_long$Class)/30)
tmp_long$Dummy <- 1
tmp_long$Subtype <- tmp$subtypes[match(tmp_long$Sample, tmp$Sample)]
tmp_long$Class <- factor(ifelse(tmp_long$Subtype == "1", "ACTG1/S100A4", ifelse(tmp_long$Subtype == "2", "DUSP22/CD9", ifelse(tmp_long$Subtype %in% c("3", "4"), "CXCR4/AHNAK/BCL7A", "CXCR4/AHNAK"))), levels = c("ACTG1/S100A4", "DUSP22/CD9", "CXCR4/AHNAK/BCL7A", "CXCR4/AHNAK"))
a <- ggplot(tmp_long) +
  geom_tile(aes(Gene, Sample, fill = Expression)) +
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(tmp_long$Expression), 0, max(tmp_long$Expression)))) +
  geom_segment(data = segment.data.vertical, aes(x = x, xend = x_end, y=y, yend = y_end), linetype = "dashed", size = 0.8) +
  geom_segment(data = segment.data.horizontal, aes(x = x, xend = x_end, y=y, yend = y_end), linetype = "dashed", size = 0.8) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle=65, hjust = 1, color = "black", size = 12, face = "italic"), axis.title = element_text(size = 14)) +
  ylab("") +
  xlab("") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0))
b <- ggplot(tmp_long) +
  geom_tile(aes(Dummy, Sample, fill = factor(Subtype))) + 
  scale_fill_brewer(palette = "Set3", name = "Cluster") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
c <- ggplot(tmp_long) +
  geom_tile(aes(Dummy, Sample, fill = factor(Stage))) + 
  scale_fill_brewer(palette = "Paired", name = "Stage") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
d <- ggplot(tmp_long) +
  geom_tile(aes(Dummy, Sample, fill = factor(Class))) + 
  scale_fill_manual(values = c("steelblue", "lightyellow", "orange", "tomato2"), name = "Subtype") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
b_leg <- get_legend(b)
b <- b + theme(legend.position = "None")
a_leg <- get_legend(a)
a <- a + theme(legend.position = "None")
c_leg <- get_legend(c)
c <- c + theme(legend.position = "None")
d_leg <- get_legend(d)
d <- d + theme(legend.position = "None")
legs <- plot_grid(a_leg, b_leg, c_leg, NULL, d_leg, ncol = 1, align = "v", rel_heights = c(1, 1, 1, -0.1, 1))
plots <- plot_grid(d, b, c, a, align = "h", ncol = 4, rel_widths = c(0.05, 0.05, 0.05, 1))
png(paste0(output_dir, "NMF_Tumor_Poisson_HL1_WL1_n30_GSE171739_SubtypeHeatmap.png"), res = 300, units="in", width = 10, height = 8)
plot_grid(plots, legs, ncol = 2, align = "h", rel_widths = c(0.5, 0.1))
dev.off()

##Write GSE171739 classification
meta_class <- tmp_long[!duplicated(tmp_long$Sample), c("Sample", "Class")]
write.csv(meta_class, paste0(output_dir, "GSE171739_meta_class.csv"), row.names = F)

#Figure 4C
##Compare signature expression levels between IgM MGUS and WM
tmp$`GEX-2` <- rowMeans(tmp[, as.character(markers$Gene[markers$Signature == "S2"])])
tmp$`GEX-3` <- rowMeans(tmp[, as.character(markers$Gene[markers$Signature == "S3"])])
tmp$`GEX-5` <- rowMeans(tmp[, as.character(markers$Gene[markers$Signature == "S5"])])
tmp$`GEX-6` <- rowMeans(tmp[, as.character(markers$Gene[markers$Signature == "S6"])])
tmp$`GEX-7` <- rowMeans(tmp[, as.character(markers$Gene[markers$Signature == "S7"])])
tmp$`GEX-8` <- rowMeans(tmp[, as.character(markers$Gene[markers$Signature == "S8"])])
sig_prog <- tmp[,c(paste0("GEX-", c(2,3,5:8)), "Sample", "Stage")]
sig_prog_long <- gather(sig_prog, "Signature", "Activity", -Sample, -Stage)
sig_prog_long$Signature <- factor(sig_prog_long$Signature, levels = paste0("GEX-", c(2,3,5:8)))
write.csv(sig_prog_long, paste0(output_dir, "SourceData_Fig4C.csv"), row.names = F)
signatures <- paste0("GEX-", c(2,3,5:8))
sig_max <- aggregate(sig_prog_long$Activity, list(sig_prog_long$Signature), max)
stats <- compare_means(Activity ~ Stage, group.by = "Signature", data = sig_prog_long, p.adjust.method = "BH")
stats$y.pos <- sig_max$x[match(stats$Signature, sig_max$Group.1)] + 0.2
write.csv(stats, paste0(output_dir, "GSE171739_Stats.csv"), row.names = F)
png(paste0(output_dir, "GSE171739_Progression.png"), res = 300, units="in", width = 7, height = 5)
ggplot(sig_prog_long) +
  geom_boxplot(aes(Stage, Activity, fill = Stage), show.legend = F, alpha = 0.5, outlier.size = -1) +
  geom_point(aes(Stage, Activity, fill = Stage), show.legend = F, shape = 21, position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("steelblue", "tomato2")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(color = "black", size = 12), 
        axis.text.x = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 14, face = "italic")) +
  ylab("Signature Activity") +
  xlab("") +
  scale_x_discrete(labels = c("IgM MGUS", "WM")) +
  geom_bracket(aes(xmin = group1, xmax = group2, label = paste0("q=", signif(p.adj, 1)), y.position = y.pos), data = stats[stats$p.adj < 0.1,], size = 0.4, label.size = 4, vjust = 1.3) +
  facet_wrap(.~Signature, scales = "free_y", ncol = 3)
dev.off()

#Figure 4D
##Compare signature activity for signatures GEX-2, GEX-5, and GEX-7 by genotype within WM samples
wm <- read.csv("results/WM_BulkRNASeq/BulkRNASeq_Data.csv", check.names = F)
table(wm$Genotype)
##CXCR4-mut MYD88-mut 
##20        32
wm_long <- gather(wm, "Signature", "Score", -SampleID, -Batch, -Genotype, -PriorRX, -Adenopathy, -SM, -BM, -IgM)
wm_long$Genotype <- factor(wm_long$Genotype, levels = c("MYD88-mut", "CXCR4-mut"))
stats <- compare_means(Score ~ Genotype, group.by = "Signature", data = wm_long, p.adjust.method = "BH")
write.csv(stats, paste0(output_dir, "SignatureActivityByGenotype_Stats.csv"), row.names = F)
sig_max <- aggregate(wm_long$Score, list(wm_long$Signature), max)
stats$y.pos <- sig_max$x[match(stats$Signature, sig_max$Group.1)] + 0.2
png(paste0(output_dir, "SignatureActivityByGenotype.png"), res = 300, units="in", width = 6, height = 5)
ggplot(wm_long) +
  geom_boxplot(aes(Genotype, Score, fill = Genotype), show.legend = F, alpha = 0.5, outlier.size = -1) +
  geom_violin(aes(Genotype, Score, fill = Genotype), show.legend = F, alpha = 0.5) +
  geom_point(aes(Genotype, Score, fill = Genotype), show.legend = F, shape = 21, position = position_jitter(width = 0.2), size = 4) +
  scale_fill_manual(values = c("steelblue", "orange")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(color = "black", size = 12), 
        axis.text.x = element_text(size = 14, angle = 65, hjust = 1), 
        axis.title = element_text(size = 14), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 14, face = "italic")) +
  ylab("Signature Activity") +
  xlab("") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = paste0("q=", signif(p.adj, 2)), y.position = y.pos*1.1), data = stats[stats$p.adj < 0.05,], label.size = 5, vjust = 1.8) +
  facet_wrap(.~Signature, scales = "free_y", ncol = 3)
dev.off()

#Figure 4E
##Compare signature activity for signatures GEX-2, GEX-5, and GEX-7 by LAD within WM samples
wm <- read.csv("results/WM_BulkRNASeq/BulkRNASeq_Data.csv", check.names = F)
table(wm$Adenopathy)
##N  Y 
##25 27
wm_long <- gather(wm, "Signature", "Score", -SampleID, -Batch, -Genotype, -PriorRX, -Adenopathy, -SM, -BM, -IgM)
wm_long$LAD <- factor(ifelse(wm_long$Adenopathy %in% "Y", "LAD+", "LAD-"), levels = c("LAD-", "LAD+"))
stats <- compare_means(Score ~ LAD, group.by = "Signature", data = wm_long, p.adjust.method = "BH")
write.csv(stats, paste0(output_dir, "SignatureActivityByLAD_Stats.csv"), row.names = F)
sig_max <- aggregate(wm_long$Score, list(wm_long$Signature), max)
stats$y.pos <- sig_max$x[match(stats$Signature, sig_max$Group.1)] + 0.2
png(paste0(output_dir, "SignatureActivityByLAD.png"), res = 300, units="in", width = 6, height = 5)
ggplot(wm_long) +
  geom_boxplot(aes(LAD, Score, fill = LAD), show.legend = F, alpha = 0.5, outlier.size = -1) +
  geom_violin(aes(LAD, Score, fill = LAD), show.legend = F, alpha = 0.5) +
  geom_point(aes(LAD, Score, fill = LAD), show.legend = F, shape = 21, position = position_jitter(width = 0.2), size = 4) +
  scale_fill_manual(values = c("pink", "red")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(color = "black", size = 12), 
        axis.text.x = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 14, face = "italic")) +
  ylab("Signature Activity") +
  xlab("") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = paste0("q=", signif(p.adj, 2)), y.position = y.pos*1.1), data = stats[stats$p.adj < 0.05,], label.size = 5, vjust = 1.8) +
  facet_wrap(.~Signature, scales = "free_y", ncol = 3)
dev.off()

#Figure 4F
##Compare signature activity for signatures GEX-2, GEX-5, and GEX-7 with BM infiltration within WM samples
wm <- read.csv("results/WM_BulkRNASeq/BulkRNASeq_Data.csv", check.names = F)
stats <- data.frame("Signature" = c("GEX-2", "GEX-5", "GEX-7"), 
                    "r" = c(cor.test(wm$BM, wm$`GEX-2`)$estimate, cor.test(wm$BM, wm$`GEX-5`)$estimate, cor.test(wm$BM, wm$`GEX-7`)$estimate),
                    "p" = c(cor.test(wm$BM, wm$`GEX-2`)$p.value, cor.test(wm$BM, wm$`GEX-5`)$p.value, cor.test(wm$BM, wm$`GEX-7`)$p.value))
stats$padj <- p.adjust(stats$p, method = "BH")
write.csv(stats, paste0(output_dir, "BMBySignatureActivity_Stats.csv"), row.names = F)
stats$BM <- c(50, 50, 50)
stats$Score <- c(1750, 400, 700)
wm_long <- gather(wm[, c("GEX-2", "GEX-5", "GEX-7", "SampleID", "BM", "IgM")], "Signature", "Score", -SampleID, -BM, -IgM)
png(paste0(output_dir, "BMBySignatureActivity.png"), res = 300, units="in", width = 6, height = 5)
ggplot(wm_long) +
  geom_point(aes(BM, Score, fill = Signature), show.legend = F, shape = 21, position = position_jitter(width = 0.2), size = 4) +
  geom_smooth(aes(BM, Score, fill = Signature), color = "red", fill = "red", method = "lm") +
  scale_fill_manual(values = c("seagreen2", "lightblue", "lemonchiffon1")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(size = 14), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 14, face = "italic")) +
  ylab("Signature Activity") +
  xlab("(%) BM Infiltration") +
  geom_text(data = stats, aes(x = BM, y = Score, label = paste0("q=", signif(padj, 2))), size = 5, fontface = "italic") +
  facet_wrap(.~Signature, scales = "free_y", ncol = 3)
dev.off()

#Figure 4G
##Identify genes that are recurrently dysregulated across WM tumors
full <- read.csv("results/Subclustering_B/DE_vsMBC/Compiled_AllTumors_Wilcoxon.csv")
wm_full <- full[full$Tumor_Type %in% "WM",]
wm_genes_up <- wm_full$Gene[wm_full$padj < 0.05 & wm_full$Log2FC > 0.5]
wm_genes_up <- wm_genes_up[-grep("RPS*|RPL*|^MT*", wm_genes_up)]
wm_gene_up_freq <- data.frame(table(wm_genes_up))
write.csv(wm_gene_up_freq, paste0(output_dir, "WM_gene_up_freq.csv"), row.names = F)

wm_genes_down <- wm_full$Gene[wm_full$padj < 0.05 & wm_full$Log2FC < -0.5]
wm_genes_down <- wm_genes_down[-grep("RPS*|RPL*|^MT*", wm_genes_down)]
wm_gene_down_freq <- data.frame(table(wm_genes_down))
write.csv(wm_gene_down_freq, paste0(output_dir, "WM_gene_down_freq.csv"), row.names = F)

##Plot recurrently upregulated genes in WM
wm_gene_up_freq <- read.csv(paste0(output_dir, "WM_gene_up_freq.csv"))
wm_gene_up_freq$wm_genes_up <- gsub("\\.", "-", wm_gene_up_freq$wm_genes_up)
wm_gene_up_freq$Gene <- factor(wm_gene_up_freq$wm_genes_up, levels = wm_gene_up_freq$wm_genes_up[order(wm_gene_up_freq$Freq, decreasing = T)])


png(paste0(output_dir, "WM_genes_up.png"), res = 300, units="in", width = 15, height = 5)
ggplot(wm_gene_up_freq[wm_gene_up_freq$Freq > 0.3*max(wm_gene_up_freq$Freq),]) +
  geom_bar(aes(Gene, Freq), fill = "steelblue", stat = "identity", width = 0.2) +
  geom_point(aes(Gene, Freq), color = "orange", size =3) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA,  color = "black"), 
        axis.text.y = element_text(color = "black", size = 10), 
        axis.text.x = element_text(color = "black", size = 12, angle = 65, hjust = 1, face = "italic")) +
  ylab("# of tumors") +
  xlab("")
dev.off()

#Figure 4H
##Plot GSE171739 DE heatmap annotated by class
data <- read.csv(paste0(input_dir, "genes_with_adjpval_v2.csv"))
colnames(data) <- gsub("\\.", "-", colnames(data))
data <- cbind(data[,c(1:2)], scale(data[,c(3:ncol(data))]))
data$Stage <- ifelse(data$Stage %in% "MGUS", "IgM MGUS", "WM")
data$Stage <- factor(data$Stage, levels = c("IgM MGUS", "WM"))
##Get gene order
for(i in 3:ncol(data)){
  if(i == 3){
    meandiff <- data.frame("gene" = colnames(data)[i], "meandiff" = mean(data[data$Stage == "WM",i], na.rm = T) - mean(data[data$Stage == "IgM MGUS",i], na.rm = T))
  } else {
    meandiff <- rbind(meandiff, data.frame("gene" = colnames(data)[i], "meandiff" = mean(data[data$Stage == "WM",i], na.rm = T) - mean(data[data$Stage == "IgM MGUS",i], na.rm = T)))
  }
}
genes <- colnames(data[3:ncol(data)])
gene_order <- as.character(meandiff$gene)[order(meandiff$meandiff, decreasing = T)]
data_long <- gather(data, "Gene", "Expression", -Sample, -Stage)
data_long$Stage <- factor(data_long$Stage, levels = c("IgM MGUS", "WM"))
data_long$Sample <- factor(data_long$Sample, levels = data$Sample[order(data$Stage)])
data_long$Dummy <- 1
data_long$Gene <- factor(data_long$Gene, levels = gene_order)
data_long$Expression <- as.numeric(as.character(data_long$Expression))
a <- ggplot(data_long) +
  geom_tile(aes(Sample, Dummy, fill = Stage)) +
  scale_fill_manual(values = c("steelblue", "orange")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  scale_y_continuous(expand = c(0,0))
a_leg <- get_legend(a)
a <- a + theme(legend.position = "NULL")
b <- ggplot(data_long) + 
  geom_tile(aes(Sample, Gene, fill = Expression)) +
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(data_long$Expression),0,max(data_long$Expression)))) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text.y = element_text(color = "black", size = 12, face = "italic"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size = 14))
b_leg <- get_legend(b)
b <- b + theme(legend.position = "NULL")
meta_class <- read.csv(paste0(output_dir, "GSE171739_meta_class.csv"))
data_long$Class <- meta_class$Class[match(data_long$Sample, meta_class$Sample)]
c <- ggplot(data_long) +
  geom_tile(aes(Sample, Dummy, fill = Class)) +
  scale_fill_manual(values = c("steelblue", "lightyellow", "orange", "tomato2")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  scale_y_continuous(expand = c(0,0))
c_leg <- get_legend(c)
c <- c + theme(legend.position = "NULL")
legs <- plot_grid(a_leg, NULL, b_leg, NULL, c_leg, ncol = 1, align = "v", rel_heights = c(0.5, -0.1, 0.5, -0.1, 0.5))
p <- plot_grid(a, c, b, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.05, 0.05, 1))
png(paste0(output_dir, "GSE171739_DE_Heatmap.png"), res = 300, units="in", width = 10, height = 8)
plot_grid(p, legs, ncol = 2, align = "h", rel_widths = c(0.5, 0.1))
dev.off()
