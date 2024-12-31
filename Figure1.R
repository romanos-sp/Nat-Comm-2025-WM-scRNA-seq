#Import dependencies
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(grid)
library(cowplot)
library(gridExtra)

#Set output directory
input_dir <- "data/"
output_dir <- "results/Figure1/"

#Set seed
set.seed(123)

#Figure 1A
##Myeloid cell UMAP
##Import RDS file
integrated <- readRDS("results/Subclustering_MonoDC2/integrated.rds")
##Create annotation vector, rename clusters and remove un-annotated clusters
annotation <- c('1' = 'S100A+ CD14+ Mono', '3' = 'RETN+ CD14+ Mono', '4' = 'HLA-DR++ CD14+ Mono', '15' = 'IFN+ CD14+ Mono', '16' = 'IL1B+ CD14+ Mono', '2' = 'CD16+ Mono', '13' = 'CD16+ Mono', '21' = 'cDC1', '6' = 'cDC2', '10' = 'LCN2+ Neut', '17' = 'FCGR3B+ Neut', '18' = 'MPO+ Neut')
renamed <- RenameIdents(integrated, annotation)
renamed_sub <- subset(renamed, idents=unname(annotation))
dim(renamed_sub)
##[1] 36601 49047
renamed_sub$Cell_Type <- Idents(renamed_sub)
##Write annotated object
saveRDS(renamed_sub, paste0(output_dir, "Myeloid_Object.rds"))
##Plot UMAP
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
plotdf$Cell_ID <- rownames(plotdf)
centroids <- aggregate(plotdf[, c("UMAP_1","UMAP_2")], list(factor(plotdf$Cell_Type)), median)
colnames(centroids)[1] <- "Cluster"
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
final_cols <- col_vector[1:length(unique(plotdf$Cell_Type))]
centroids$UMAP_2[centroids$Cluster=="IFN+ CD14+ Mono"] <- centroids$UMAP_2[centroids$Cluster=="IFN+ CD14+ Mono"]+0.3
centroids$UMAP_1[centroids$Cluster=="IFN+ CD14+ Mono"] <- centroids$UMAP_1[centroids$Cluster=="IFN+ CD14+ Mono"]-2
png(paste0(output_dir, "Myeloid_UMAP.png"), res = 300, units="in", width = 8, height = 5)
ggplot(plotdf,aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=factor(Cell_Type)), show.legend = F) +
  labs(color="Cluster") + 
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_blank()) + 
  annotate("text", x=centroids$UMAP_1, y=centroids$UMAP_2, label= centroids$Cluster, size=6) +
  scale_color_manual(values = final_cols) +
  coord_cartesian(clip="off")
dev.off()

##T cell UMAP
##Import RDS file
integrated <- readRDS("results/Subclustering_T/integrated.rds")
##Create annotation vector, rename clusters and remove un-annotated clusters
annotation <- c('4' = 'CD8+ TN', '0' = 'CD4+ TN', '6' = 'aCD4+ TCM', '14' = 'IFN+ CD8+ TEM', '17' = 'Tgd', '13' = 'MAIT', '12' = 'Th1', '7' = 'Th17', '1' = 'CD4+ TCM', '8' = 'Treg', '3' = 'GZMK+ CD8+ TEM', '10' = 'IFN+ CD4+ TCM', '2' = 'GZMB+ CD8+ TEM', '9' = 'CD8+ TCM', '15' = 'aCD8+ TCM', '16' = 'KIR+ CD8+ TEM')
renamed <- RenameIdents(integrated, annotation)
renamed_sub <- subset(renamed, idents=unname(annotation))
dim(renamed_sub)
##[1] 36601 95754
renamed_sub$Cell_Type <- Idents(renamed_sub)
##Write annotated object
saveRDS(renamed_sub, paste0(output_dir, "T_Object.rds"))
##Plot UMAP
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
plotdf$Cell_ID <- rownames(plotdf)
centroids <- aggregate(plotdf[, c("UMAP_1","UMAP_2")], list(factor(plotdf$Cell_Type)), median)
colnames(centroids)[1] <- "Cluster"
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
final_cols <- col_vector[1:length(unique(plotdf$Cell_Type))]
png(paste0(output_dir, "T_UMAP.png"), res = 300, units="in", width = 8, height = 5)
ggplot(plotdf,aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=factor(Cell_Type)), show.legend = F) +
  labs(color="Cluster") + 
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_blank()) + 
  annotate("text", x=centroids$UMAP_1, y=centroids$UMAP_2, label= centroids$Cluster, size=6) +
  scale_color_manual(values = final_cols) +
  coord_cartesian(clip="off")
dev.off()

##NK cell UMAP
##Import RDS file
integrated <- readRDS("results/Subclustering_NK2/integrated.rds")
##Create annotation vector, rename clusters and remove un-annotated clusters
annotation <- c('9' = 'CD56br NK', '6' = 'CD56br NK', '13' = 'CD160+ NK', '3' = 'aNK', '10' = 'aNK', '11' = 'IFN+ NK', '0' = 'S100B+ CD56dim NK', '1' = 'CD56dim NK', '2' = 'CD56dim NK', '4' = 'CD56dim NK')
renamed <- RenameIdents(integrated, annotation)
renamed_sub <- subset(renamed, idents=unname(annotation))
dim(renamed_sub)
##[1] 36601 24219
renamed_sub$Cell_Type <- Idents(renamed_sub)
##Write annotated object
saveRDS(renamed_sub, paste0(output_dir, "NK_Object.rds"))
##Plot UMAP
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
plotdf$Cell_ID <- rownames(plotdf)
centroids <- aggregate(plotdf[, c("UMAP_1","UMAP_2")], list(factor(plotdf$Cell_Type)), median)
colnames(centroids)[1] <- "Cluster"
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
final_cols <- col_vector[1:length(unique(plotdf$Cell_Type))]
centroids$UMAP_1[centroids$Cluster=="CD56int NK"] <- centroids$UMAP_1[centroids$Cluster=="CD56int NK"]+1
png(paste0(output_dir, "NK_UMAP.png"), res = 300, units="in", width = 8, height = 5)
ggplot(plotdf,aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=factor(Cell_Type)), show.legend = F) +
  labs(color="Cluster") + 
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_blank()) + 
  annotate("text", x=centroids$UMAP_1, y=centroids$UMAP_2, label= centroids$Cluster, size=6) +
  scale_color_manual(values = final_cols) +
  coord_cartesian(clip="off")
dev.off()

##Progenitor cell UMAP
##Import RDS file
integrated <- readRDS("results/Subclustering_Stem2/integrated.rds")
##Create annotation vector, rename clusters and remove un-annotated clusters
annotation <- c('3' = 'HSC', '24' = 'HSC', '16' = 'GMP', '0' = 'MP', '19' = 'DP', '14' = 'Pro-E', '2' = 'Baso-E', '11' = 'Baso-E', '5' = 'Poly-E', '8' = 'Poly-E', '9' = 'Ortho-E', '25' = 'BMP', '10' = 'MK', '7' = 'pDC', '27' = 'Stroma', '20' = 'C1Q+ Macrophage')
renamed <- RenameIdents(integrated, annotation)
renamed_sub <- subset(renamed, idents=unname(annotation))
dim(renamed_sub)
##[1] 36601 40999
renamed_sub$Cell_Type <- Idents(renamed_sub)
##Write annotated object
saveRDS(renamed_sub, paste0(output_dir, "Progenitor_Object.rds"))
##Remove pDC, C1Q+ Macrophage, and Stroma
renamed_sub <- subset(renamed_sub, subset = Cell_Type %in% c('HSC', 'Pro-E', 'Baso-E', 'Poly-E', 'Ortho-E', 'MK', 'GMP', 'MP', 'DP', 'BMP'))
##Plot UMAP
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
plotdf$Cell_ID <- rownames(plotdf)
centroids <- aggregate(plotdf[, c("UMAP_1","UMAP_2")], list(factor(plotdf$Cell_Type)), median)
colnames(centroids)[1] <- "Cluster"
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
final_cols <- col_vector[1:length(unique(plotdf$Cell_Type))]
png(paste0(output_dir, "Progenitor_UMAP.png"), res = 300, units="in", width = 8, height = 5)
ggplot(plotdf,aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=factor(Cell_Type)), show.legend = F) +
  labs(color="Cluster") + 
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_blank()) + 
  annotate("text", x=centroids$UMAP_1, y=centroids$UMAP_2, label= centroids$Cluster, size=6) +
  scale_color_manual(values = final_cols) +
  coord_cartesian(clip="off")
dev.off()

#Figure 1B
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
##Compare immune cell proportions between AWM and NBM
threshold <- 100
tmp <- read.csv(paste0(input_dir, "Full_count_table.csv"), row.names = 1, check.names = F)
##Remove stroma cells
tmp <- tmp[, !colnames(tmp) %in% "Stroma"]
##Aggregate cell counts by CaTissue ID
tmp$CaTissueID <- metadata$CaTissueID[match(rownames(tmp), metadata$Demux_Sample_ID_Final)]
df <- transform(aggregate(tmp[,!colnames(tmp) %in% c("CaTissueID")], list(tmp$CaTissueID), sum, na.rm=T), row.names=Group.1, Group.1=NULL)
colnames(df) <- colnames(tmp[,!colnames(tmp) %in% c("CaTissueID")])
##Remove samples with < 100 cells across immune subpopulations
df <- df[rowSums(df, na.rm = T)>= threshold,]
nrow(df)
##71
##Compute proportions
df <- df/rowSums(df, na.rm = T)
##Compare proportions between AWM and NBM using Wilcoxon tests
df$Status <- metadata$Disease[match(rownames(df), metadata$CaTissueID)]
df <- df[df$Status %in% c("AWM", "NBM"),]
nrow(df)
##43
table(df$Status)
##AWM NBM 
##25  18
df_long <- gather(df, "Cell_Type", "Proportion", -Status)
out <- compare_means(Proportion ~ Status, group.by="Cell_Type", data=df_long, p.adjust.method = "BH")
write.csv(out, paste0(output_dir, "AWMvsNBM_Proportion_Stats.csv"), row.names = F)
##Plot volcano plot
volmat <- matrix(nrow=length(unique(df_long$Cell_Type)), ncol=6)
colnames(volmat) <- c("Cell_Type", "Wilcoxon_p", "FDR", "mean_AWM", "mean_NBM", "LFC")
volmat <- data.frame(volmat)
volmat$Cell_Type <- unique(df_long$Cell_Type)
volmat$Wilcoxon_p <- out$p[match(volmat$Cell_Type, out$Cell_Type)]
volmat$FDR <- out$p.adj[match(volmat$Cell_Type, out$Cell_Type)]
for(ind in 1:length(unique(df_long$Cell_Type))){
  cl <- unique(df_long$Cell_Type)[ind]
  volmat[volmat$Cell_Type == cl, "mean_AWM"] <-   mean(as.numeric(as.character(df_long$Proportion[df_long$Cell_Type==cl & df_long$Status %in% "AWM"])), na.rm=T)
  volmat[volmat$Cell_Type == cl, "mean_NBM"] <-   mean(as.numeric(as.character(df_long$Proportion[df_long$Cell_Type==cl & df_long$Status %in% "NBM"])), na.rm=T)
}
volmat$LFC <- log2(volmat$mean_AWM/volmat$mean_NBM)
write.csv(volmat, paste0(output_dir, "AWMvsNBM_Proportion_Data.csv"), row.names = F)
crange <- t(matrix(c("#009BF4","#EAEAEA","#FF9933","#EAEAEA","#EAEAEA","#EAEAEA"), ncol=2))
g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"), interpolate = TRUE)
volmat$shape <- factor(ifelse(volmat$FDR < 0.1 & volmat$Wilcoxon_p < 0.05, 1, 0), levels = c("0", "1"))

png(paste0(output_dir, "AWMvsNBM_Proportion_Volcano.png"), res = 300, units="in",width = 7, height = 5)
ggplot(volmat,aes(x=LFC, y=-log10(Wilcoxon_p)), size=4) +
  annotation_custom(g, xmin = -(ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T))))) + 1), xmax= ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T))))) + 1, ymin=-.5, ymax= max(-log10(volmat$Wilcoxon_p)) + 0.5) + 
  xlim(-ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T))))), ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T)))))) +
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
        legend.position = c(0.9,0.5)) + 
  geom_text_repel(data = volmat[volmat$shape == 1,], aes(x=LFC, y=-log10(Wilcoxon_p), label=Cell_Type), size = 5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)+
  scale_shape_manual(values=c(16,8), labels=c("NS", "<0.1"), name = "q-value") +
  scale_alpha_manual(values = c(0.3, 1), guide = "none") +
  annotate("text", x = ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T)))))*3/4, y = 0, label = paste0("AWM, n=", sum(df$Status %in% "AWM")), size = 8) +
  annotate("text", x = -ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T)))))*3/4, y = 0, label = paste0("HD, n=", sum(df$Status %in% "NBM")), size = 8)
dev.off() 

#Figure 1C
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
##Compare immune cell proportions between SMM and NBM
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
##Compare proportions between SMM and NBM using Wilcoxon tests
df$Status <- metadata$Disease[match(rownames(df), metadata$CaTissueID)]
df <- df[df$Status %in% c("SMM", "NBM"),]
nrow(df)
##42
table(df$Status)
##NBM SMM 
##18  24
df_long <- gather(df, "Cell_Type", "Proportion", -Status)
out <- compare_means(Proportion ~ Status, group.by="Cell_Type", data=df_long, p.adjust.method = "BH")
write.csv(out, paste0(output_dir, "SMMvsNBM_Proportion_Stats.csv"), row.names = F)

##Plot overall volcano plot
volmat <- matrix(nrow=length(unique(df_long$Cell_Type)), ncol=6)
colnames(volmat) <- c("Cell_Type","Wilcoxon_p", "FDR", "mean_SMM", "mean_NBM", "LFC")
volmat <- data.frame(volmat)
volmat$Cell_Type <- unique(df_long$Cell_Type)
volmat$Wilcoxon_p <- out$p[match(volmat$Cell_Type, out$Cell_Type)]
volmat$FDR <- out$p.adj[match(volmat$Cell_Type, out$Cell_Type)]
for(ind in 1:length(unique(df_long$Cell_Type))){
  cl <- unique(df_long$Cell_Type)[ind]
  volmat[volmat$Cell_Type == cl,"mean_SMM"] <-   mean(as.numeric(as.character(df_long$Proportion[df_long$Cell_Type==cl & df_long$Status %in% "SMM"])), na.rm=T)
  volmat[volmat$Cell_Type == cl,"mean_NBM"] <-   mean(as.numeric(as.character(df_long$Proportion[df_long$Cell_Type==cl & df_long$Status %in% "NBM"])), na.rm=T)
}
volmat$LFC <- log2(volmat$mean_SMM/volmat$mean_NBM)
write.csv(volmat, paste0(output_dir, "SMMvsNBM_Proportion_Data.csv"), row.names = F)
crange <- t(matrix(c("#009BF4","#EAEAEA","#FC5A5A","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)
volmat$shape <- factor(ifelse(volmat$FDR < 0.1 & volmat$Wilcoxon_p < 0.05, 1, 0), levels = c("0", "1"))

png(paste0(output_dir, "SMMvsNBM_Proportion_Volcano.png"), res = 300, units="in",width = 7, height = 5)
ggplot(volmat,aes(x=LFC, y=-log10(Wilcoxon_p)), size=4) +
  annotation_custom(g, xmin = -(ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T))))) + 1), xmax= ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T))))) + 1, ymin=-.5, ymax= max(-log10(volmat$Wilcoxon_p)) + 0.5) + 
  xlim(-ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T))))), ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T)))))) +
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
  annotate("text", x = ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T)))))*3/4, y = 0, label = paste0("SMM, n=", sum(df$Status %in% "SMM")), size = 8) +
  annotate("text", x = -ceiling(max(c(abs(min(volmat$LFC, na.rm = T)), abs(max(volmat$LFC, na.rm = T)))))*3/4, y = 0, label = paste0("HD, n=", sum(df$Status %in% "NBM")), size = 8)
dev.off() 

#Figure 1D
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
##Compare immune cell proportions between IgM MGUS and NBM
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
##Compare proportions between IgM MGUS and NBM using Wilcoxon tests
df$Status <- metadata$Stage[match(rownames(df), metadata$CaTissueID)]
df <- df[df$Status %in% c("IgM MGUS", "NBM"),]
nrow(df)
##24
table(df$Status)
##IgM MGUS      NBM 
##6       18 
df_long <- gather(df, "Cell_Type", "Proportion", -Status)
out <- compare_means(Proportion ~ Status, group.by="Cell_Type", data=df_long, p.adjust.method = "BH")
write.csv(out, paste0(output_dir, "IgMMGUSvsNBM_Proportion_Stats.csv"), row.names = F)

##Plot overall volcano plot
volmat <- matrix(nrow=length(unique(df_long$Cell_Type)), ncol=6)
colnames(volmat) <- c("Cell_Type","Wilcoxon_p", "FDR", "mean_IgMMGUS", "mean_NBM", "LFC")
volmat <- data.frame(volmat)
volmat$Cell_Type <- unique(df_long$Cell_Type)
volmat$Wilcoxon_p <- out$p[match(volmat$Cell_Type, out$Cell_Type)]
volmat$FDR <- out$p.adj[match(volmat$Cell_Type, out$Cell_Type)]
for(ind in 1:length(unique(df_long$Cell_Type))){
  cl <- unique(df_long$Cell_Type)[ind]
  volmat[volmat$Cell_Type == cl,"mean_IgMMGUS"] <-   mean(as.numeric(as.character(df_long$Proportion[df_long$Cell_Type==cl & df_long$Status %in% "IgM MGUS"])), na.rm=T)
  volmat[volmat$Cell_Type == cl,"mean_NBM"] <-   mean(as.numeric(as.character(df_long$Proportion[df_long$Cell_Type==cl & df_long$Status %in% "NBM"])), na.rm=T)
}
volmat$LFC <- log2(volmat$mean_IgMMGUS/volmat$mean_NBM)
write.csv(volmat, paste0(output_dir, "IgMMGUSvsNBM_Proportion_Data.csv"), row.names = F)
crange <- t(matrix(c("#009BF4","#EAEAEA","#FFFACD","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)
volmat$shape <- factor(ifelse(volmat$FDR < 0.1 & volmat$Wilcoxon_p < 0.05, 1, 0), levels = c("0", "1"))

png(paste0(output_dir, "IgMMGUSvsNBM_Proportion_Volcano.png"), res = 300, units="in",width = 7, height = 5)
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
        legend.position = c(0.92, 0.5)) + 
  geom_text_repel(data = volmat[volmat$shape == 1,], aes(x=LFC, y=-log10(Wilcoxon_p), label=Cell_Type), size = 5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)+
  scale_shape_manual(values=c(16,8), labels=c("NS", "<0.1"), name = "q-value") +
  scale_alpha_manual(values = c(0.3, 1), guide = "none") +
  annotate("text", x = ceiling(max(c(abs(min(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)), abs(max(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)))))*2/3, y = 0, label = paste0("IgM MGUS, n=", sum(df$Status %in% "IgM MGUS")), size = 8) +
  annotate("text", x = -ceiling(max(c(abs(min(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)), abs(max(volmat$LFC[is.finite(volmat$LFC)], na.rm = T)))))*3/4, y = 0, label = paste0("HD, n=", sum(df$Status %in% "NBM")), size = 8)
dev.off() 

#Figure 1E
##Plot heatmap of proportions for cell types that change significantly in AWM across HD-IgM MGUS-SWM-WM individuals
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
##Compare immune cell proportions between IgM MGUS and NBM
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
##Prepare input matrix for plotting
df$Status <- metadata$Stage[match(rownames(df), metadata$CaTissueID)]
df <- df[df$Status %in% c("NBM", "IgM MGUS", "SWM", "WM"),]
df$Status <- factor(df$Status, levels = c("NBM", "IgM MGUS", "SWM", "WM"))
nrow(df)
##46
table(df$Status)
##NBM IgM MGUS      SWM       WM 
##18        6       19        3  
df$Sample <- rownames(df)
df_long <- gather(df, "Cell_Type", "Proportion", -Status, -Sample)
df_long$Sample <- factor(df_long$Sample, levels = df$Sample[order(df$Status)])
##Sort cell types by their LFC between AWM and NBM
volmat <- read.csv(paste0(output_dir, "AWMvsNBM_Proportion_Data.csv"))
sig <- volmat$Cell_Type[volmat$FDR < 0.1 & volmat$Wilcoxon_p < 0.05]
##Remove progenitor populations
sig <- sig[!sig %in% c("HSC", "GMP", "MP", "DP", "Pro-E", "Baso-E", "Poly-E", "Ortho-E", "BMP")]
length(sig)
##21
volmat <- volmat[volmat$Cell_Type %in% sig,]
df_long <- df_long[df_long$Cell_Type %in% sig,]
df_long$Cell_Type <- factor(df_long$Cell_Type, levels=volmat$Cell_Type[order(volmat$LFC, decreasing = F)])
df_long$Dummy <- 1
df_long$LFC <- volmat$LFC[match(df_long$Cell_Type, volmat$Cell_Type)]
df_long$Status <- as.character(df_long$Status)
df_long$Status[df_long$Status %in% "NBM"] <- "HD"
df_long$Status <- factor(df_long$Status, levels = c("HD", "IgM MGUS", "SWM", "WM"))
deid <- read.csv(paste0(input_dir, "Complete_SampleTable_DeIdentified.csv"))
df_long$DeIdentified_SampleID <- deid$SampleDeIdentified[match(df_long$Sample, deid$CaTissueID)]
sum(is.na(df_long$DeIdentified_SampleID))
##0
write.csv(df_long[!colnames(df_long) %in% c("Sample", "Dummy")], paste0(output_dir, "HD_MGUS_SWM_WM_Proportion_Data.csv"), row.names = F)

a <- ggplot(df_long) +
  geom_tile(aes(x=Sample, y=Cell_Type, fill=Proportion*100)) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(color="black", fill=NA), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size=12, color="black"), 
        axis.title = element_blank(), 
        legend.position = "bottom") +
  scale_fill_gradient(low="white", high="red", name = "(%) Proportion")
b <- ggplot(df_long) +
  geom_tile(aes(Sample, Dummy, fill=Status)) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(color="black", fill=NA), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank()) +
  scale_fill_manual(values=c("steelblue", "lemonchiffon1", "orange", "tomato2"), name = "")
b_leg <- get_legend(b)
b <- b + 
  theme(legend.position = "none")
c <- ggplot(df_long) +
  geom_tile(aes(Cell_Type, x=Dummy, fill=LFC)) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(color="black", fill=NA), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank()) +
  scale_fill_gradientn(colors = c("darkblue","white","orangered3"),values=scales::rescale(c(min(df_long$LFC),0,max(df_long$LFC))))
a_leg <- get_legend(a)
a <- a + 
  theme(legend.position = "none")
a <- a + 
  theme(axis.text.y = element_blank())
c_leg <- get_legend(c)
c <- c + 
  theme(legend.position = "none")
a_ticks <- ggplot(df_long) +
  geom_point(aes(x=NA, y=Cell_Type), size = -1) +
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size=12, color="black"), 
        axis.title = element_blank())
blankPlot <- ggplot() +
  geom_blank(aes(1,1)) +
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

png(paste0(output_dir, "HD_MGUS_SWM_WM_Proportion_Heatmap.png"), res = 300, units="in", width = 10, height = 5)
grid.arrange(blankPlot, 
             b, 
             blankPlot, 
             blankPlot, 
             blankPlot, 
             a_ticks, 
             a, 
             c, 
             c_leg, 
             b_leg, 
             blankPlot, 
             a_leg, 
             blankPlot, 
             blankPlot, 
             blankPlot, 
             nrow = 3, 
             ncol = 5, 
             widths = c(20, 60, 4, 10, 10), 
             heights = c(5, 60, 6))
dev.off()

#Figure 1F
##Compare Treg proportions between IgM MGUS and SWM 
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
##Compare immune cell proportions between IgM MGUS and SWM
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
##Compare proportions between IgM MGUS and NBM using Wilcoxon tests
df$Status <- metadata$Stage[match(rownames(df), metadata$CaTissueID)]
df <- df[df$Status %in% c("IgM MGUS", "SWM"),]
nrow(df)
##25
table(df$Status)
##IgM MGUS      SWM 
##6       19  
sum(is.na(df$Treg))
##0
df$Sample <- rownames(df)
df_long <- gather(df, "Cell_Type", "Proportion", -Status, -Sample)
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
df_long$DeIdentified_SampleID <- deid$SampleDeIdentified[match(df_long$Sample, deid$CaTissueID)]
sum(is.na(df_long$DeIdentified_SampleID))
##0
plotdf <- df_long[df_long$Cell_Type %in% c("Treg"),]
write.csv(plotdf[!colnames(plotdf) %in% "Sample"], paste0(output_dir, "IgMMGUSvsSWM_Proportion_Data.csv"), row.names = F)
stats <- compare_means(Proportion ~ Status, group.by = "Cell_Type", data = plotdf, p.adjust.method = "BH")
write.csv(stats, paste0(output_dir, "IgMMGUSvsSWM_Proportion_WilcoxonStats.csv"), row.names = F)
stats <- compare_means(Proportion ~ Status, group.by = "Cell_Type", data = plotdf, p.adjust.method = "BH", method = "t.test")
write.csv(stats, paste0(output_dir, "IgMMGUSvsSWM_Proportion_tStats.csv"), row.names = F)
maxes <- aggregate(plotdf$Proportion, list(plotdf$Cell_Type), max, na.rm = T)
stats$y.position <- maxes$x[match(stats$Cell_Type, maxes$Group.1)]

png(paste0(output_dir, "Treg_SWMvsIgM MGUS.png"), res = 300, units="in", width = 5, height = 5)
ggplot(plotdf) +
  geom_boxplot(aes(Status, Proportion*100, fill = Status), outlier.size = -1, alpha = 0.5, show.legend = F) +
  geom_violin(aes(Status, Proportion*100, fill = Status), alpha = 0.5, show.legend = F) +
  geom_point(aes(Status, Proportion*100, fill = Status), alpha = 0.8, shape = 21, size = 5, position = position_jitter(width = 0.2), show.legend = F) +
  scale_fill_manual(values = c("lemonchiffon1", "orange")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14)) +
  xlab("") +
  ylab("(%) Proportion of Tregs") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = paste0("p=", signif(p, 2)), y.position = y.position*100 + 0.5), data = stats, size = 0.4, label.size = 5, vjust = 1.5) +
  scale_y_continuous(limits = c(0, 6.8), breaks = c(0, 2, 4, 6))
dev.off()

#Figure 1G
##Plot boxplots for activated and IFN-stimulated populations
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
##Compare immune cell proportions between IgM MGUS and SWM
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
##Compare proportions between IgM MGUS and NBM using Wilcoxon tests
df$Status <- metadata$Stage[match(rownames(df), metadata$CaTissueID)]
df <- df[df$Status %in% c("NBM", "SMM", "IgM MGUS", "SWM"),]
nrow(df)
##67
df$Status <- as.character(df$Status)
df$Status[df$Status %in% "NBM"] <- "HD"
df$Status <- factor(df$Status, levels = c("HD", "SMM", "IgM MGUS", "SWM"))
table(df$Status)
##HD      SMM IgM MGUS      SWM 
##18       24        6       19 
df$Sample <- rownames(df)
df_long <- gather(df, "Cell_Type", "Proportion", -Status, -Sample)
##Isolate populations of interest
plotdf <- df_long[df_long$Cell_Type %in% c("aNK", "aCD4+ TCM", "aCD8+ TCM", "IFN+ CD4+ TCM", "IFN+ CD8+ TEM", "IFN+ NK"),]
plotdf$Cell_Type <- factor(plotdf$Cell_Type, levels = c("aCD4+ TCM", "aCD8+ TCM", "aNK", "IFN+ CD4+ TCM", "IFN+ CD8+ TEM", "IFN+ NK"))
write.csv(plotdf, paste0(output_dir, "Activated_And_IFN_Proportions_Data_Identifiable.csv"), row.names = F)
deid <- read.csv(paste0(input_dir, "Complete_SampleTable_DeIdentified.csv"))
plotdf$DeIdentified_SampleID <- deid$SampleDeIdentified[match(plotdf$Sample, deid$CaTissueID)]
sum(is.na(plotdf$DeIdentified_SampleID))
##0
write.csv(plotdf[, !colnames(plotdf) %in% "Sample"], paste0(output_dir, "Activated_And_IFN_Proportions_Data.csv"), row.names = F)
stats <- compare_means(Proportion ~ Status, group.by = "Cell_Type", data = plotdf, p.adjust.method = "BH")
maxes <- aggregate(plotdf$Proportion*100, list(plotdf$Cell_Type), max, na.rm = T)
stats$max <- maxes$x[match(stats$Cell_Type, maxes$Group.1)] + 0.3
stats$ypos <- NA
stats <- stats[stats$p.adj < 0.1,]
ct <- unique(stats$Cell_Type)
for(i in seq_along(ct)){
  tmp <- stats[stats$Cell_Type == ct[i],]
  for(x in 1:nrow(tmp)){
    tmp$ypos[x] <- tmp$max[x] + x*1.5
  }
  if(i==1){
    out <- tmp
  } else {
    out <- rbind(out, tmp)
  }
}
stats <- out
write.csv(stats, paste0(output_dir, "Activated_And_IFN_Proportions_Stats.csv"), row.names = F)
png(paste0(output_dir, "Activated_And_IFN_Proportions.png"), res = 300, units="in", width = 7, height = 7)
ggplot(plotdf, aes(Status, Proportion*100)) + 
  geom_boxplot(aes(fill=Status), show.legend = F, outlier.size = -1, alpha = 0.5, position = position_dodge()) +
  geom_violin(aes(fill=Status), show.legend = F, alpha = 0.5, position = position_dodge()) +
  geom_point(position=position_jitter(width = 0.2), size=2) +
  scale_fill_manual(values=c("steelblue", "tomato2", "lemonchiffon1", "orange")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="black"), axis.text.x = element_text(size=14, color="black", angle=65, hjust =1), axis.text.y = element_text(color="black", size=12), axis.title = element_text(size=14), strip.background = element_blank(), strip.text = element_text(size=14, face = "italic")) +
  xlab("") + 
  ylab("(%) Proportion") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = paste0("q=", signif(p.adj, 2)), y.position = ypos), label.size=4,  vjust = 1.8, data = stats[stats$p.adj < 0.1,]) +
  facet_wrap(.~Cell_Type, scales = "free_y")
dev.off()

#Figure 1H
##Scatterplot of LFC in myeloid cells from AWM vs NBM and AWM vs SMM
full <- read.csv("results/Subclustering_MonoDC2/Compiled_All.csv")
nbm <- grep("Log2FC_AWMvsNBM", colnames(full))
smm <- grep("Log2FC_AWMvsSMM", colnames(full))
full$Ave_LFC_AWMvsNBM <- rowMeans(as.matrix(full[, nbm]), na.rm=T)
full$Ave_LFC_AWMvsSMM <- rowMeans(as.matrix(full[, smm]), na.rm=T)
full$Diff_Ave <- full$Ave_LFC_AWMvsNBM-full$Ave_LFC_AWMvsSMM
tmp <- full
rem <- grep("RPL*|RPS*|\\.", tmp$Gene)
tmp <- tmp[-rem,]
labs <- tmp$Gene[order(tmp$Diff_Ave, decreasing = T)][1:30]
labs <- c(labs, tmp$Gene[order(tmp$Diff_Ave, decreasing = F)][1:30])
tmp$lab <- ifelse(tmp$Gene %in% c(labs, "MNDA", "FGL2", "PTPRC", "CXCL8"), as.character(tmp$Gene), NA)
write.csv(tmp[,c("Gene", "Ave_LFC_AWMvsNBM", "Ave_LFC_AWMvsSMM", "Diff_Ave")], paste0(output_dir, "DE_Diff_Compiled_Myeloid_Data.csv"), row.names = F)
tmp$lab2 <- ifelse(tmp$Gene %in% c("MNDA"), "bolditalic('MNDA')", 
                   ifelse(tmp$Gene %in% c("CEBPD"), "bolditalic('CEBPD')", 
                          ifelse(tmp$Gene %in% c("S100A8"), "bolditalic('S100A8')", 
                                 ifelse(tmp$Gene %in% c("S100A9"), "bolditalic('S100A9')", 
                                        paste0("italic('", as.character(tmp$lab), "')")))))
tmp$lab2[is.na(tmp$lab)] <- NA
png(paste0(output_dir, "DE_Diff_Compiled_Myeloid.png"), res = 300, units="in", width = 7, height = 5)
ggplot(tmp[,c("Gene", "Ave_LFC_AWMvsNBM", "Ave_LFC_AWMvsSMM", "Diff_Ave", "lab", "lab2")]) + 
  geom_point(aes(Ave_LFC_AWMvsNBM, Ave_LFC_AWMvsSMM, color=Diff_Ave)) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", alpha=0.5) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill=NA, color="black"), 
        axis.text = element_text(color = "black", size=12), 
        axis.title = element_text(size=14), 
        legend.text = element_text(size=10)) +
  geom_text_repel(aes(Ave_LFC_AWMvsNBM, Ave_LFC_AWMvsSMM, label=lab2), max.overlaps = Inf, segment.alpha=0.2, parse = T) +
  scale_color_gradientn(colors = c("steelblue","lavender","orangered3"),values=scales::rescale(c(min(tmp$Diff_Ave),0,max(tmp$Diff_Ave)))) +
  ylab(expression("Mean Log"[2]*"FC: AWM vs SMM")) +
  xlab(expression("Mean Log"[2]*"FC: AWM vs HD")) +
  labs(color="Difference") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  annotate("text", x = -0.5, y = 0.9, label = "Myeloid cells", size = 5, fontface = "italic")
dev.off() 

###Average adjusted p-values for selected genes
smmnbm <- grep("padj_SMMvsNBM", colnames(tmp))
tmp$Ave_Padj_SMMvsNBM <- rowMeans(as.matrix(tmp[, smmnbm]), na.rm=T)
tmp$Ave_Padj_SMMvsNBM[tmp$Gene == "CEBPD"]
##8.983144e-07
tmp$Ave_Padj_SMMvsNBM[tmp$Gene == "S100A8"]
##0.012
tmp$Ave_Padj_SMMvsNBM[tmp$Gene == "S100A9"]
##0.022
awmnbm <- grep("padj_AWMvsNBM", colnames(tmp))
tmp$Ave_Padj_AWMvsNBM <- rowMeans(as.matrix(tmp[, awmnbm]), na.rm=T)
tmp$Ave_Padj_AWMvsNBM[tmp$Gene == "MNDA"]
##5.395995e-10

#Figure 1I
###UMAP of Myeloid cells colored for MNDA expression
integrated <- readRDS(paste0(output_dir, "Myeloid_Object.rds"))
plotdf <- FetchData(integrated, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type", "Disease"))
plotdf <- plotdf[plotdf$Disease %in% c("NBM", "AWM", "SMM"),]
plotdf$Disease <- as.character(plotdf$Disease)
plotdf$Disease[plotdf$Disease %in% "NBM"] <- "HD"
features <- c("MNDA")
X <- data.frame(t(data.frame(GetAssayData(integrated[features,], slot="data"))))
rownames(X) <- gsub("^X", "", rownames(X))
rownames(X) <- gsub("\\.", "-", rownames(X))
tmp <- transform(merge(plotdf, X, by="row.names"), row.names=Row.names, Row.names=NULL)
tmp$Disease <- factor(tmp$Disease, levels = c("HD", "SMM", "AWM"))
table(tmp$Disease)
##HD   SMM   AWM 
##15957  5634 23521

png(paste0(output_dir, "MNDA_Myeloid_UMAP.png"), res = 300, units="in", width = 8, height = 5)
ggplot() + 
  geom_point(data=tmp, aes(UMAP_1, UMAP_2, color=MNDA), show.legend = T) +
  xlab("") +
  ylab("") + 
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.x = element_text(hjust=0.05, vjust=0.5, size=16), 
        axis.title.y = element_text(hjust=0.1, vjust=0.5, size=16), 
        strip.background = element_blank(), 
        strip.text = element_text(size=18, face="italic")) + 
  scale_color_gradient(low="lightblue", high="red") +
  coord_cartesian(clip="off") +
  facet_grid(.~Disease)
dev.off()

#Figure 1J
##Scatterplot of LFC in T cells from AWM vs NBM and AWM vs SMM
full <- read.csv("results/Subclustering_T/Compiled_All.csv")
nbm <- grep("Log2FC_AWMvsNBM", colnames(full))
smm <- grep("Log2FC_AWMvsSMM", colnames(full))
full$Ave_LFC_AWMvsNBM <- rowMeans(as.matrix(full[, nbm]), na.rm=T)
full$Ave_LFC_AWMvsSMM <- rowMeans(as.matrix(full[, smm]), na.rm=T)
full$Diff_Ave <- full$Ave_LFC_AWMvsNBM-full$Ave_LFC_AWMvsSMM
tmp <- full
rem <- grep("RPL*|RPS*|\\.", tmp$Gene)
tmp <- tmp[-rem,]
labs <- tmp$Gene[order(tmp$Diff_Ave, decreasing = T)][1:30]
labs <- c(labs, tmp$Gene[order(tmp$Diff_Ave, decreasing = F)][1:30])
write.csv(tmp[,c("Gene", "Ave_LFC_AWMvsNBM", "Ave_LFC_AWMvsSMM", "Diff_Ave")], paste0(output_dir, "DE_Diff_Compiled_T_Data.csv"), row.names = F)
tmp$lab <- ifelse(tmp$Gene %in% c(labs, "LCK", "CD2", "CD69"), as.character(tmp$Gene), NA)
tmp$lab[tmp$lab %in% c("TRAC", "TRBC1", "TRBC2", "PPDPF", "TRAF3IP3", "HNRNPA2B1", "DPSMB10", "RIPOR2", "PFN1")] <- NA
tmp$lab2 <- ifelse(tmp$Gene %in% c("PTPRC"), "bolditalic('PTPRC')", 
                   ifelse(tmp$Gene %in% c("CD2"), "bolditalic('CD2')", 
                          ifelse(tmp$Gene %in% c("LCK"), "bolditalic('LCK')", 
                                 ifelse(tmp$Gene %in% c("IL2RG"), "bolditalic('IL2RG')", 
                                        ifelse(tmp$Gene %in% c("IL7R"), "bolditalic('IL7R')", 
                                               ifelse(tmp$Gene %in% c("GIMAP4"), "bolditalic('GIMAP4')", 
                                                      ifelse(tmp$Gene %in% c("RAC2"), "bolditalic('RAC2')", 
                                                             ifelse(tmp$Gene %in% c("CXCR4"), "bolditalic('CXCR4')", 
                                                                    ifelse(tmp$Gene %in% c("CD69"), "bolditalic('CD69')", 
                                                                           ifelse(tmp$Gene %in% c("JUN"), "bolditalic('JUN')", 
                                                                                  ifelse(tmp$Gene %in% c("FOS"), "bolditalic('FOS')", 
                                                                                         ifelse(tmp$Gene %in% c("TNFAIP3"), "bolditalic('TNFAIP3')", 
                                                                                                paste0("italic('", as.character(tmp$lab), "')")))))))))))))

tmp$lab2[is.na(tmp$lab)] <- NA
png(paste0(output_dir, "DE_Diff_Compiled_T.png"), res = 300, units="in", width = 7, height = 5)
ggplot(tmp[,c("Gene", "Ave_LFC_AWMvsNBM", "Ave_LFC_AWMvsSMM", "Diff_Ave", "lab", "lab2")]) + 
  geom_point(aes(Ave_LFC_AWMvsNBM, Ave_LFC_AWMvsSMM, color=Diff_Ave)) +
  geom_abline(slope=1, intercept = 0, linetype="dashed", alpha=0.5) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill=NA, color="black"), 
        axis.text = element_text(color = "black", size=12), 
        axis.title = element_text(size=14), 
        legend.text = element_text(size=10)) +
  geom_text_repel(aes(Ave_LFC_AWMvsNBM, Ave_LFC_AWMvsSMM, label=lab2), max.overlaps = Inf, segment.alpha=0.2, parse = T) +
  scale_color_gradientn(colors = c("steelblue","lavender","orangered3"),values=scales::rescale(c(min(tmp$Diff_Ave),0,max(tmp$Diff_Ave)))) +
  ylab(expression("Mean Log"[2]*"FC: AWM vs SMM")) +
  xlab(expression("Mean Log"[2]*"FC: AWM vs HD")) +
  labs(color="Difference") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  annotate("text", x = -1, y = 0.65, label = "T cells", size = 5, fontface = "italic")
dev.off()

#Figure 1K
##Plot CERI ELISPOT ratios
ceri <- read.csv(paste0(input_dir, "CERI_ELISPOT_AWM.csv"))
ceri$Disease <- factor(ceri$Disease, levels = c("HD", "SMM", "AWM"))
table(ceri$Disease)
##HD SMM AWM 
##10  10   3
stats <- compare_means(Ratio ~ Disease, data = ceri, p.adjust.method = "BH")
write.csv(stats, paste0(output_dir, "CERI_ELISPOT_Stats.csv"), row.names = F)
png(paste0(output_dir, "CERI_ELISPOT.png"), res = 300, units="in", width = 5, height = 5)
ggplot(ceri) + 
  geom_boxplot(aes(Disease, Ratio, fill = Disease), alpha = 0.5, outlier.size = -1, show.legend = F) +
  geom_violin(aes(Disease, Ratio, fill = Disease), alpha = 0.5, show.legend = F) + 
  geom_point(aes(Disease, Ratio, fill = Disease), shape = 21, show.legend = F, size = 4, position = position_jitterdodge()) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 14)) +
  scale_fill_manual(values = c("steelblue", "tomato2", "orange")) +
  scale_y_log10() + 
  ylab("Ratio: CERI/DMSO") +
  xlab("") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = paste0("q=", signif(p.adj, 2))), y.position = 1.2, step.increase = 0.1, vjust = 1.5, label.size = 5, data = stats)
dev.off()

#Figure 1L
##Perform PCA on immune cell proportions and make density plot
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
##Read in full count table
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
##Compute principle components
pca <- data.frame(prcomp(df)$x)
write.csv(pca, paste0(output_dir, "PCA_ImmuneProportions.csv"), row.names = T)
##Plot density plot
pca <- read.csv(paste0(output_dir, "PCA_ImmuneProportions.csv"), row.names = 1)
pca$Disease <- metadata$Disease[match(rownames(pca), metadata$CaTissueID)]
pca$Disease <- as.character(pca$Disease)
pca$Disease[pca$Disease=="NBM"] <- "HD"
pca <- pca[pca$Disease %in% c("HD", "SMM", "AWM"),]
pca$Disease <- factor(pca$Disease, levels= c("HD", "SMM", "AWM"))
table(pca$Disease)
##HD SMM AWM 
##18  24  25
deid <- read.csv(paste0(input_dir, "Complete_SampleTable_DeIdentified.csv"))
pca$Deidentified_SampleID <- deid$SampleDeIdentified[match(rownames(pca), deid$CaTissueID)]
sum(is.na(pca$Deidentified_SampleID))
##0
write.csv(pca[, c("Deidentified_SampleID", "Disease", "PC1", "PC2")], paste0(output_dir, "PCA_ImmuneProportions_Data.csv"), row.names = F)
png(paste0(output_dir, "PCA_ImmuneProportions_DensityPlot.png"), res = 300, units="in", width = 6, height = 5)
ggplot(pca, aes(PC1, PC2, group = Disease))  + 
  geom_point(aes(fill = Disease), shape = 21, size = 4, alpha = 0.9) +
  stat_density2d(aes(alpha=..level.., fill=Disease), geom="polygon") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill=NA, color="black"), 
        axis.text = element_text(size=12, color="black"), 
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12), 
        legend.position = c(0.2, 0.85)) + 
  scale_fill_manual(values=c("steelblue", "tomato2", "orange"), labels = c("HD, n=18", "SMM, n=24", "AWM, n=25"), name = "") +
  guides(alpha = "none")
dev.off()

#Figure 1M
##Visualize confusion matrix for 5-fold cross-validation of disease classifier
cv_test <- read.csv("results/DiseaseClassifier/5fold_CV.csv", row.names = 1, check.names = F)
conf <- data.frame(table(cv_test$predictions, cv_test$class))
colnames(conf) <- c("Prediction", "Reference", "Count")
perc <- data.frame(prop.table(table(cv_test$predictions, cv_test$class), margin = 2))
colnames(perc) <- c("Prediction","Reference", "Agr")
anno_t <- data.frame(table(cv_test$predictions, cv_test$class))
colnames(anno_t) <- c("Prediction","Reference","Count")
perc$Prediction <- factor(perc$Prediction, levels = c("0", "2", "1"))
perc$Reference <- factor(perc$Reference, levels = c("0", "2", "1"))
anno_t$Prediction <- factor(anno_t$Prediction, levels = c("0", "2", "1"))
anno_t$Reference <- factor(anno_t$Reference, levels = c("0", "2", "1"))
png(paste0(output_dir, "ConfusionMatrix_SVM_golub17f_5foldCV.png"), res = 300, units="in", width = 6, height = 5)
p <- ggplot(perc) + 
  geom_tile(aes(Reference,Prediction,fill=Agr*100)) + 
  scale_fill_gradient(low="white",high="red") + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="black"), axis.title = element_text(size=16), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=10)) + 
  labs(fill="% of Ref") + 
  scale_x_discrete(label=c("HD", "SMM", "AWM")) + 
  scale_y_discrete(label=c("HD", "SMM", "AWM"))
p + 
  geom_text(aes(Reference, Prediction, label = Count), data=anno_t, size=10)
dev.off()

cv_test <- read.csv("results/DiseaseClassifier/5fold_CV.csv", row.names = 1, check.names = F)
deid <- read.csv(paste0(input_dir, "Complete_SampleTable_DeIdentified.csv"))
cv_test$Deidentified_SampleID <- deid$SampleDeIdentified[match(rownames(cv_test), deid$CaTissueID)]
sum(is.na(cv_test$Deidentified_SampleID))
##0
write.csv(cv_test[, c("Deidentified_SampleID", "class", "fold", "predictions")], paste0(output_dir, "5fold_CV_Data.csv"), row.names = F)

##Sensitivity/Specificity for cancer detection
prop.test(sum(anno_t$Count[anno_t$Prediction!="0" & anno_t$Reference!="0"]), sum(anno_t$Count[anno_t$Reference!="0"]), correct = T)
##49/49, 100% (95% CI: 91-100)
prop.test(sum(anno_t$Count[anno_t$Prediction=="0" & anno_t$Reference=="0"]), sum(anno_t$Count[anno_t$Reference=="0"]), correct = T)
##15/17, 88% (95% CI: 62-98)

##Sensitivity/Specificity for SMM diagnosis
prop.test(sum(anno_t$Count[anno_t$Prediction=="2" & anno_t$Reference=="2"]), sum(anno_t$Count[anno_t$Reference=="2"]), correct = T)
##20/24, 83% (95% CI: 62-95)
prop.test(sum(anno_t$Count[anno_t$Prediction!="2" & anno_t$Reference!="2"]), sum(anno_t$Count[anno_t$Reference!="2"]), correct = T)
##37/42, 88% (95% CI: 74-96)

##Sensitivity/Specificity for AWM diagnosis
prop.test(sum(anno_t$Count[anno_t$Prediction=="1" & anno_t$Reference=="1"]), sum(anno_t$Count[anno_t$Reference=="1"]), correct = T)
##22/25, 88% (95% CI: 68-97)
prop.test(sum(anno_t$Count[anno_t$Prediction!="1" & anno_t$Reference!="1"]), sum(anno_t$Count[anno_t$Reference!="1"]), correct = T)
##37/41, 90% (95% CI: 76-97)