#Import dependencies
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(ggpubr)
library(grid)

#Set output directory
input_dir <- "data/"
output_dir <- "results/SuppFigure1/"

#Supp. Figure 1A
##Myeloid cells
##Read in annotated myeloid object
renamed_sub <- readRDS("results/Figure1/Myeloid_Object.rds")
##Plot Myeloid marker heatmap
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
plotdf$Cell_ID <- rownames(plotdf)
markers <- c("LCN2", "PGLYRP1", "LTF", "NCF2", "CAMP", "FCGR3B", "CSF3R", "CMTM2", "LUCAT1", "MPO", "AZU1", "PRTN3", "ELANE", "RETN", "RNASE2", "CD14", "CD36", "MS4A6A", "LYZ", "FCN1", "SELL", "VCAN", "S100A8", "S100A9", "S100A10", "S100A12", "IL1B", "CXCL8", "CXCL2", "PLAUR", "ISG15", "IFI6", "IFI44L", "MX1", "STAT1", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CD74", "FCGR3A", "MS4A7", "CSF1R", "CDKN1C", "RHOC", "FCER1A", "CLEC10A", "CD1C", "CLEC9A", "C1orf54", "IDO1")
orig_cell_types <- c('LCN2+ Neut', 'FCGR3B+ Neut', 'MPO+ Neut', 'RETN+ CD14+ Mono', 'S100A+ CD14+ Mono', 'IL1B+ CD14+ Mono', 'IFN+ CD14+ Mono', 'HLA-DR++ CD14+ Mono', 'CD16+ Mono', 'cDC2', 'cDC1')
marker_counts <- data.frame(t(data.frame(GetAssayData(renamed_sub[markers,], slot = "data"))))
rownames(marker_counts) <- gsub("^X", "", rownames(marker_counts))
rownames(marker_counts) <- gsub("\\.", "-", rownames(marker_counts))
marker_counts <- data.frame(scale(marker_counts))
marker_counts$idents <- plotdf$Cell_Type[match(rownames(marker_counts), rownames(plotdf))]
marker_means <- aggregate(marker_counts[, 1:(ncol(marker_counts)-1)], list(marker_counts$idents), mean)
rownames(marker_means) <- marker_means$Group.1
marker_means$Group.1  <- NULL
marker_means_t <- data.frame(t(marker_means))
colnames(marker_means_t) <- rownames(marker_means)
rownames(marker_means_t) <- gsub("\\.","-", rownames(marker_means_t))
marker_means_t$Marker <- rownames(marker_means_t)
marker_means_long <- gather(marker_means_t, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels = orig_cell_types)
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)
png(paste0(output_dir, "Myeloid_Annotation_Heatmap.png"), res = 300, units="in", width = 11, height = 5)
ggplot(marker_means_long) +
  geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"), values = scales::rescale(c(min(marker_means_long$Mean_Exp), 0, max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 65, hjust=1, size=12, face="italic", color="black"), 
        axis.text.y = element_text(size=12, color="black")) + 
  xlab("") + 
  ylab("") + 
  labs(fill="Mean Exp")
dev.off()

##T cells
##Read in annotated T cell object
renamed_sub <- readRDS("results/Figure1/T_Object.rds")
##Plot T marker heatmap
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
plotdf$Cell_ID <- rownames(plotdf)
markers <- c("CD3D", "CD4", "LEF1", "TCF7", "SELL", "CCR7", "IL7R", "HNRNPLL", "FAS", "SESN3", "GATA3", "CCR4", "KRT1", "RORA", "CCR6", "KLRB1", "TNFRSF4", "IL2RA", "FOXP3", "CTLA4", "RTKN2", "IKZF2", "HLA-DRA", "HLA-DRB1", "CD74", "CXCR3", "GZMA", "GZMK", "LYAR", "CCL5", "JUN", "JUNB", "FOS", "FOSB", "DUSP1", "TSC22D3", "CD69", "MX1", "ISG15", "IFI6", "IFI44L", "STAT1", "IFIT1", "IFIT2", "IFIT3", "SLC4A10", "TRAV1-2", "TRBV20-1", "CD8A", "CD8B", "LINC02446", "NELL2", "CMC1", "XCL1", "XCL2", "GZMB", "GZMH", "PRF1", "GNLY", "FGFBP2", "FCGR3A", "TRGV9", "TRDV2", "TRDC", "TRGC1", "KLRD1", "TRDV1", "KLRC2", "KLRC3", "KIR3DL2")
orig_cell_types <- c('CD4+ TN', 'CD4+ TCM', 'Th17', 'Treg', 'Th1', 'aCD4+ TCM', 'IFN+ CD4+ TCM', 'MAIT', 'CD8+ TN', 'CD8+ TCM', 'aCD8+ TCM', 'GZMK+ CD8+ TEM', 'GZMB+ CD8+ TEM', 'IFN+ CD8+ TEM', 'Tgd', 'KIR+ CD8+ TEM')
marker_counts <- data.frame(t(data.frame(GetAssayData(renamed_sub[markers,], slot = "data"))))
rownames(marker_counts) <- gsub("^X", "", rownames(marker_counts))
rownames(marker_counts) <- gsub("\\.", "-", rownames(marker_counts))
marker_counts <- data.frame(scale(marker_counts))
marker_counts$idents <- plotdf$Cell_Type[match(rownames(marker_counts), rownames(plotdf))]
marker_means <- aggregate(marker_counts[, 1:(ncol(marker_counts)-1)], list(marker_counts$idents), mean)
rownames(marker_means) <- marker_means$Group.1
marker_means$Group.1  <- NULL
marker_means_t <- data.frame(t(marker_means))
colnames(marker_means_t) <- rownames(marker_means)
rownames(marker_means_t) <- gsub("\\.","-", rownames(marker_means_t))
marker_means_t$Marker <- rownames(marker_means_t)
marker_means_long <- gather(marker_means_t, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels = orig_cell_types)
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)
png(paste0(output_dir,"/T_Annotation_Heatmap.png"), res = 300, units="in",width = 13, height = 5)
ggplot(marker_means_long) +
  geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1, size=12, face="italic", color="black"), axis.text.y=element_text(size=12, color="black")) + 
  xlab("") + 
  ylab("") + 
  labs(fill="Mean Exp") 
dev.off()

##NK cells
##Read in annotated NK cell object
renamed_sub <- readRDS("results/Figure1/NK_Object.rds")
##Plot NK marker heatmap
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
plotdf$Cell_ID <- rownames(plotdf)
markers <- c("NCAM1", "SELL", "KLRC1", "GZMK", "CD160", "CCL3", "CCL4", "GZMB", "FGFBP2", "PRF1", "GNLY", "NKG7", "S100B", "GZMH", "ITGB7", "PTGDS", "JUN", "FOS", "FOSB", "DUSP1", "TNFAIP3", "ISG15", "IFIT1", "IFIT2", "IFIT3")
orig_cell_types <- c('CD56br NK', 'CD160+ NK', 'CD56dim NK', 'S100B+ CD56dim NK', 'aNK', 'IFN+ NK')
marker_counts <- data.frame(t(data.frame(GetAssayData(renamed_sub[markers,], slot = "data"))))
rownames(marker_counts) <- gsub("^X", "", rownames(marker_counts))
rownames(marker_counts) <- gsub("\\.", "-", rownames(marker_counts))
marker_counts <- data.frame(scale(marker_counts))
marker_counts$idents <- plotdf$Cell_Type[match(rownames(marker_counts), rownames(plotdf))]
marker_means <- aggregate(marker_counts[, 1:(ncol(marker_counts)-1)], list(marker_counts$idents), mean)
rownames(marker_means) <- marker_means$Group.1
marker_means$Group.1  <- NULL
marker_means_t <- data.frame(t(marker_means))
colnames(marker_means_t) <- rownames(marker_means)
rownames(marker_means_t) <- gsub("\\.","-", rownames(marker_means_t))
marker_means_t$Marker <- rownames(marker_means_t)
marker_means_long <- gather(marker_means_t, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels = orig_cell_types)
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)
png(paste0(output_dir,"/NK_Annotation_Heatmap.png"), res = 300, units="in",width = 8, height = 3)
ggplot(marker_means_long) +
  geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1, size=12, face="italic", color="black"), axis.text.y=element_text(size=12, color="black")) + 
  xlab("") + 
  ylab("") + 
  labs(fill="Mean Exp") 
dev.off()

##Progenitor cells
##Read in annotated Progenitor cell object
renamed_sub <- readRDS("results/Figure1/Progenitor_Object.rds")
##Remove pDC, C1Q+ Macrophage, and Stroma
renamed_sub <- subset(renamed_sub, subset = Cell_Type %in% c('HSC', 'Pro-E', 'Baso-E', 'Poly-E', 'Ortho-E', 'MK', 'GMP', 'MP', 'DP', 'BMP'))
##Plot Progenitor marker heatmap
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
plotdf$Cell_ID <- rownames(plotdf)
markers <- c("CD34", "NME2", "NPM1", "MIF", "STMN1", "APOC1", "STAT5", "ENG", "PKLR", "GATA1", "KLF1", "TFRC", "STOM", "HMGB1", "HMGB2", "ITGA4","GYPA", "GYPB", "ALAS2", "SNCA", "FOXO3", "ARL4A", "TMCC2", "EIF5", "XPO7", "IFIT1B", "PPBP", "PF4", "ITGA2B", "GP9", "NRGN", "CAVIN2", "MPO", "AZU1", "PRTN3", "ELANE", "LYZ", "S100A8", "LGALS1", "FCER1A", "CLC", "PRG2", "HDC", "TPSAB1", "TPSB2", "TPSD1", "CPA3")
orig_cell_types <- c('HSC', 'Pro-E', 'Baso-E', 'Poly-E', 'Ortho-E', 'MK', 'GMP', 'MP', 'DP', 'BMP')
marker_counts <- data.frame(t(data.frame(GetAssayData(renamed_sub[markers,], slot = "data"))))
rownames(marker_counts) <- gsub("^X", "", rownames(marker_counts))
rownames(marker_counts) <- gsub("\\.", "-", rownames(marker_counts))
marker_counts <- data.frame(scale(marker_counts))
marker_counts$idents <- plotdf$Cell_Type[match(rownames(marker_counts), rownames(plotdf))]
marker_means <- aggregate(marker_counts[, 1:(ncol(marker_counts)-1)], list(marker_counts$idents), mean)
rownames(marker_means) <- marker_means$Group.1
marker_means$Group.1  <- NULL
marker_means_t <- data.frame(t(marker_means))
colnames(marker_means_t) <- rownames(marker_means)
rownames(marker_means_t) <- gsub("\\.","-", rownames(marker_means_t))
marker_means_t$Marker <- rownames(marker_means_t)
marker_means_long <- gather(marker_means_t, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels = orig_cell_types)
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)
png(paste0(output_dir, "Progenitor_Annotation_Heatmap.png"), res = 300, units="in", width = 10, height = 4)
ggplot(marker_means_long) +
  geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1, size=12, face="italic", color="black"), axis.text.y=element_text(size=12, color="black")) + 
  xlab("") + 
  ylab("") + 
  labs(fill="Mean Exp") 
dev.off()

#Supp. Figure 1B
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
##Compare immune cell proportions between AWM and SMM
threshold <- 100
tmp <- read.csv(paste0(input_dir, "Full_count_table.csv"), row.names = 1, check.names = F)
##Remove stroma cells
tmp <- tmp[, !colnames(tmp) %in% "Stroma"]
##Aggregate cell counts by CaTissue ID
tmp$CaTissueID <- metadata$CaTissueID[match(rownames(tmp), metadata$Demux_Sample_ID_Final)]
df <- transform(aggregate(tmp[,!colnames(tmp) %in% c("CaTissueID")], list(tmp$CaTissueID), sum, na.rm=T), row.names=Group.1, Group.1=NULL)
colnames(df) <- colnames(tmp[,!colnames(tmp) %in% c("CaTissueID")])
##Remove samples with < 100 cells across subpopulations
df <- df[rowSums(df, na.rm = T)>= threshold,]
##Compute proportions
df <- df/rowSums(df, na.rm = T)
##Compare proportions between AWM and SMM using Wilcoxon tests
df$Status <- metadata$Disease[match(rownames(df), metadata$CaTissueID)]
df <- df[df$Status %in% c("AWM", "SMM"),]
nrow(df)
##49
table(df$Status)
##AWM SMM 
##25  24  
df_long <- gather(df, "Cell_Type", "Proportion", -Status)
out <- compare_means(Proportion ~ Status, group.by="Cell_Type", data=df_long, p.adjust.method = "BH")
write.csv(out, paste0(output_dir, "AWMvsSMM_Proportion_Stats.csv"), row.names = F)

##Plot overall volcano plot
volmat <- matrix(nrow=length(unique(df_long$Cell_Type)), ncol=6)
colnames(volmat) <- c("Cell_Type","Wilcoxon_p", "FDR", "mean_AWM", "mean_SMM", "LFC")
volmat <- data.frame(volmat)
volmat$Cell_Type <- unique(df_long$Cell_Type)
volmat$Wilcoxon_p <- out$p[match(volmat$Cell_Type, out$Cell_Type)]
volmat$FDR <- out$p.adj[match(volmat$Cell_Type, out$Cell_Type)]
for(ind in 1:length(unique(df_long$Cell_Type))){
  cl <- unique(df_long$Cell_Type)[ind]
  volmat[volmat$Cell_Type == cl,"mean_AWM"] <-   mean(as.numeric(as.character(df_long$Proportion[df_long$Cell_Type==cl & df_long$Status %in% "AWM"])), na.rm=T)
  volmat[volmat$Cell_Type == cl,"mean_SMM"] <-   mean(as.numeric(as.character(df_long$Proportion[df_long$Cell_Type==cl & df_long$Status %in% "SMM"])), na.rm=T)
}
volmat$LFC <- log2(volmat$mean_AWM/volmat$mean_SMM)
write.csv(volmat, paste0(output_dir, "AWMvsSMM_Proportion_Data.csv"), row.names = F)
crange <- t(matrix(c("#FC5A5A","#EAEAEA","#FF9933","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)
volmat$shape <- factor(ifelse(volmat$FDR < 0.1 & volmat$Wilcoxon_p < 0.05, 1, 0), levels = c("0", "1"))

png(paste0(output_dir, "AWMvsSMM_Proportion_Volcano.png"), res = 300, units="in",width = 7, height = 5)
ggplot(volmat,aes(x=LFC, y=-log10(Wilcoxon_p)), size=4) +
  annotation_custom(g, xmin = -(ceiling(max(c(abs(min(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)), abs(max(volmat$LFC[is.finite(volmat$LFC)], na.rm = T))))) + 1), xmax= ceiling(max(c(abs(min(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)), abs(max(volmat$LFC[is.finite(volmat$LFC)], na.rm = T))))) + 1, ymin=-.5, ymax= max(-log10(volmat$Wilcoxon_p)) + 0.5) + 
  xlim(-ceiling(max(c(abs(min(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)), abs(max(volmat$LFC[is.finite(volmat$LFC)], na.rm = T))))), ceiling(max(c(abs(min(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)), abs(max(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)))))) +
  geom_point(aes(shape=shape, alpha=shape), show.legend = T) + 
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color="black")+
  xlab(expression("Log"[2]*"fold-change of the mean")) +
  ylab(expression("-Log"[10]*"(p-value)")) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(color="black", size=12), 
        axis.title = element_text(size=14), 
        legend.position = c(0.92, 0.8)) + 
  geom_text_repel(data = volmat[volmat$shape == 1,], aes(x=LFC, y=-log10(Wilcoxon_p), label=Cell_Type), size = 5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)+
  scale_shape_manual(values=c(16,8), labels=c("NS", "<0.1"), name = "q-value") +
  scale_alpha_manual(values = c(0.3, 1), guide = "none") +
  annotate("text", x = ceiling(max(c(abs(min(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)), abs(max(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)))))*2/3, y = 0, label = paste0("AWM, n=", sum(df$Status %in% "AWM")), size = 8) +
  annotate("text", x = -ceiling(max(c(abs(min(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)), abs(max(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)))))*3/4, y = 0, label = paste0("SMM, n=", sum(df$Status %in% "SMM")), size = 8)
dev.off() 

#Supp. Figure 1C
##Plot T cell UMAP colored by the level of an IFN signature
integrated <- readRDS("results/Figure1/T_Object.rds")
integrated <- AddModuleScore(integrated, features = list(c("ISG15", "ISG20", "IFI6", "IFI44L", "MX1", "STAT1", "MX2", "OAS1", "OASL", "IFIT1", "IFIT2", "IFIT3", "IFI27", "EIF2AK2")), name = c("IFN_sig"))
integrated$IFN_sig <- integrated$IFN_sig1
plotdf <- FetchData(integrated, vars= c("Cell","UMAP_1", "UMAP_2", "Demux_Sample_ID", "Cell_Type", "Disease", "IFN_sig"))
plotdf <- plotdf[plotdf$Disease %in% c("NBM", "AWM", "SMM"),]
plotdf$Disease <- as.character(plotdf$Disease)
plotdf$Disease[plotdf$Disease %in% "NBM"] <- "HD"
plotdf$Disease <- factor(plotdf$Disease, levels = c("HD", "SMM", "AWM"))
write.csv(plotdf, paste0(output_dir, "IFNsig_T_Data.csv"), row.names = T)
png(paste0(output_dir, "IFNsig_T_UMAP.png"), res = 300, units="in",width = 8, height = 5)
ggplot(data = plotdf) + 
  geom_point(aes(UMAP_1, UMAP_2, fill=IFN_sig), show.legend = T, shape = 21, stroke = 0) +
  xlab("") +
  ylab("") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_text(hjust=0.05, vjust=0.5, size=16), axis.title.y = element_text(hjust=0.1, vjust=0.5, size=16), strip.background = element_blank(), strip.text = element_text(size=18, face="italic")) + 
  scale_fill_gradient(low="white", high="blue") +
  coord_cartesian(clip="off") +
  facet_grid(.~Disease)
dev.off()

#Supp. Figure 1D
##Compare IFN sig activity within IFN+ T cell subpopulations across diagnoses
plotdf <- read.csv(paste0(output_dir, "SourceData_SuppFig1D.csv"))
plotdf$Disease <- factor(plotdf$Disease, levels = c("HD", "SMM", "AWM"))
stats <- compare_means(IFN_sig ~ Disease, data = plotdf[grep("IFN", plotdf$Cell_Type),], p.adjust.method = "BH", ref.group = "AWM")
write.csv(stats, paste0(output_dir, "IFNsig_IFN_T_Stats.csv"), row.names = F)
stats$p.format <- ifelse(stats$p.adj < 2e-16, "q<2e-16", paste0("q=", signif(stats$p.adj, 2)))
df <- data.frame(table(plotdf$Disease))
df$x <- 1:nrow(df)

png(paste0(output_dir, "IFNsig_IFN_T_Boxplots.png"), res = 300, units="in",width = 5, height = 5)
ggplot(plotdf[grep("IFN", plotdf$Cell_Type),]) +
  geom_boxplot(aes(Disease, IFN_sig, fill = Disease), alpha = 0.5, show.legend = F) +
  geom_violin(aes(Disease, IFN_sig, fill = Disease), alpha = 0.5, show.legend = F) +
  scale_fill_manual(values = c("steelblue", "tomato2", "orange")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.title = element_text(size = 14)) +
  xlab("") +
  ylab("IFN sig activity in IFN+ T cells") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = p.format), y.position = 2.4, step.increase = 0.1, vjust = 1.5, label.size = 5, data = stats[stats$p.adj < 0.05,]) +
  scale_y_continuous(limits = c(min(plotdf$IFN_sig)-0.2, max(plotdf$IFN_sig) + 0.6)) +
  geom_text(data = df, aes(x = x, y = min(plotdf$IFN_sig)-0.15, label = paste0("n=", Freq)), size = 5)
dev.off()

#Supp. Figure 1E
##UMAP of Myeloid cells colored for CEBPD expression
integrated <- readRDS("results/Figure1/Myeloid_Object.rds")
plotdf <- FetchData(integrated, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type", "Disease"))
plotdf <- plotdf[plotdf$Disease %in% c("NBM", "AWM", "SMM"),]
plotdf$Disease <- as.character(plotdf$Disease)
plotdf$Disease[plotdf$Disease %in% "NBM"] <- "HD"
features <- c("CEBPD")
X <- data.frame(t(data.frame(GetAssayData(integrated[features,], slot="data"))))
rownames(X) <- gsub("^X", "", rownames(X))
rownames(X) <- gsub("\\.", "-", rownames(X))
tmp <- transform(merge(plotdf, X, by="row.names"), row.names=Row.names, Row.names=NULL)
tmp$Disease <- factor(tmp$Disease, levels = c("HD", "SMM", "AWM"))

png(paste0(output_dir, "CEBPD_Myeloid_UMAP.png"), res = 300, units="in", width = 8, height = 5)
ggplot() + 
  geom_point(data=tmp, aes(UMAP_1, UMAP_2, color=CEBPD), show.legend = T) +
  xlab("") +
  ylab("") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_text(hjust=0.05, vjust=0.5, size=16), axis.title.y = element_text(hjust=0.1, vjust=0.5, size=16), strip.background = element_blank(), strip.text = element_text(size=18, face="italic")) + 
  scale_color_gradient(low="lightblue", high="red") +
  coord_cartesian(clip="off") +
  facet_grid(.~Disease)
dev.off()

#Supp. Figure 1F
##UMAP of Myeloid cells colored for S100A8 expression
integrated <- readRDS("results/Figure1/Myeloid_Object.rds")
plotdf <- FetchData(integrated, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type", "Disease"))
plotdf <- plotdf[plotdf$Disease %in% c("NBM", "AWM", "SMM"),]
plotdf$Disease <- as.character(plotdf$Disease)
plotdf$Disease[plotdf$Disease %in% "NBM"] <- "HD"
features <- c("S100A8")
X <- data.frame(t(data.frame(GetAssayData(integrated[features,], slot="data"))))
rownames(X) <- gsub("^X", "", rownames(X))
rownames(X) <- gsub("\\.", "-", rownames(X))
tmp <- transform(merge(plotdf, X, by="row.names"), row.names=Row.names, Row.names=NULL)
tmp$Disease <- factor(tmp$Disease, levels = c("HD", "SMM", "AWM"))

png(paste0(output_dir, "S100A8_Myeloid_UMAP.png"), res = 300, units="in", width = 8, height = 5)
ggplot() + 
  geom_point(data=tmp, aes(UMAP_1, UMAP_2, color=S100A8), show.legend = T) +
  xlab("") +
  ylab("") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_text(hjust=0.05, vjust=0.5, size=16), axis.title.y = element_text(hjust=0.1, vjust=0.5, size=16), strip.background = element_blank(), strip.text = element_text(size=18, face="italic")) + 
  scale_color_gradient(low="lightblue", high="red") +
  coord_cartesian(clip="off") +
  facet_grid(.~Disease)
dev.off()

#Supp. Figure 1G
##UMAP of T cells colored for IL2RG expression
integrated <- readRDS("results/Figure1/T_Object.rds")
plotdf <- FetchData(integrated, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type", "Disease"))
plotdf <- plotdf[plotdf$Disease %in% c("NBM", "AWM", "SMM"),]
plotdf$Disease <- as.character(plotdf$Disease)
plotdf$Disease[plotdf$Disease %in% "NBM"] <- "HD"
features <- c("IL2RG")
X <- data.frame(t(data.frame(GetAssayData(integrated[features,], slot="data"))))
rownames(X) <- gsub("^X", "", rownames(X))
rownames(X) <- gsub("\\.", "-", rownames(X))
tmp <- transform(merge(plotdf, X, by="row.names"), row.names=Row.names, Row.names=NULL)
tmp$Disease <- factor(tmp$Disease, levels = c("HD", "SMM", "AWM"))

png(paste0(output_dir, "IL2RG_T_UMAP.png"), res = 300, units="in", width = 8, height = 5)
ggplot() + 
  geom_point(data=tmp, aes(UMAP_1, UMAP_2, color=IL2RG), show.legend = T) +
  xlab("") +
  ylab("") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_text(hjust=0.05, vjust=0.5, size=16), axis.title.y = element_text(hjust=0.1, vjust=0.5, size=16), strip.background = element_blank(), strip.text = element_text(size=18, face="italic")) + 
  scale_color_gradient(low="lightblue", high="red") +
  coord_cartesian(clip="off") +
  facet_grid(.~Disease)
dev.off()

#Supp. Figure 1H
##Plot variable importance of 17-feature SVM model
golub_ranks_long <- read.csv("results/DiseaseClassifier/VariableImportance_SVM_golub17f_Data.csv")
golub_ranks <- spread(golub_ranks_long[, !colnames(golub_ranks_long) %in% "sign"], Class, Score)
golub_ranks_long$Features <- factor(golub_ranks_long$Features, levels = golub_ranks$Features[order(golub_ranks$AWM, decreasing = T)])
golub_ranks_long$sign <- factor(golub_ranks_long$sign, levels = c("0", "1"))
png(paste0(output_dir, "VariableImportance_SVM_golub17f.png"), res = 300, units="in", width = 7, height = 4)
ggplot(golub_ranks_long) + 
  geom_bar(aes(x=Score, y=Features, fill = sign), stat="identity", color=NA, width = 0.3, show.legend = F) + 
  geom_point(aes(x=Score, y=Features), color="orange", size=2) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="black"), axis.text = element_text(size=12, color="black"), axis.title = element_text(size=14), strip.background = element_blank(), strip.text = element_text(size = 14, face = "italic")) +
  ylab("") + 
  xlab("Importance") +
  scale_fill_manual(values = c("steelblue", "tomato2")) +
  facet_wrap(.~Class)
dev.off()

#Supp. Figure 1I
##Identify misclassified cases by 17-feature SVM classifier and compare BM infiltration
cv_test <- read.csv("results/DiseaseClassifier/5fold_CV.csv", row.names = 1, check.names = F)
cv_test$status <- ifelse(cv_test$class != cv_test$predictions, "Incorrect", "Correct")
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
cv_test$BM <- metadata$BM_Infiltration[match(rownames(cv_test), metadata$CaTissueID)]
cv_test$Disease <- metadata$Disease[match(rownames(cv_test), metadata$CaTissueID)]
awm <- cv_test[cv_test$Disease %in% "AWM",]
awm$BM <- as.numeric(as.character(awm$BM))
awm$CaTissueID <- rownames(awm)
write.csv(awm[, c("CaTissueID", "status", "BM", "Disease")], paste0(output_dir, "MisclassificationData.csv"), row.names = F)
stats <- compare_means(BM ~ status, data = awm, p.adjust.method = "BH")
png(paste0(output_dir, "Misclassification_BM.png"), res = 300, units="in", width = 5, height = 5)
ggplot(awm) + 
  geom_boxplot(aes(status, BM*100, fill = status), show.legend = F, outlier.size = -1, alpha = 0.5) + 
  geom_violin(aes(status, BM*100, fill = status), show.legend = F, alpha = 0.5) + 
  geom_point(aes(status, BM*100, fill = status), shape = 21, size = 4, show.legend = F, position = position_jitter(width = 0.2), alpha = 0.8) +
  scale_fill_manual(values = c("seagreen3", "red")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.title = element_text(size = 14)) +
  ylab("(%) BM infiltration") +
  xlab("") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = paste0("p=", p.format)), y.position = 95, step.increase = 0.1, vjust = 1.5, label.size = 5, data = stats)
dev.off()

#How many patients with AWM out of 25 have < 10% infiltration?
sum(awm$BM < 0.1)
##8
prop.test(sum(awm$BM < 0.1), nrow(awm), correct = T)
#32%, 95% CI: 16-54

#Are misclassified cases associated with sample processing type?
cv_test <- read.csv(paste0(output_dir, "5fold_CV.csv"), row.names = 1, check.names = F)
cv_test$status <- ifelse(cv_test$class != cv_test$predictions, "Incorrect", "Correct")
cv_test$Method <- metadata$Sample_Processing_Method[match(rownames(cv_test), metadata$CaTissueID)]
cv_test$Disease <- metadata$Disease[match(rownames(cv_test), metadata$CaTissueID)]
awm <- cv_test[cv_test$Disease %in% "AWM",]
table(awm$status, awm$Method)
##          Ficoll RBCLB
##Correct       14     8
##Incorrect      1     2
fisher.test(table(awm$status, awm$Method))
##p-value = 0.5435