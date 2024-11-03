#Import dependencies
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(Seurat)
library(ggpubr)

#Set input/output directories
input_dir <- "data/"
output_dir <- "results/SuppFigure3/"

#Supp. Figure 3A
##Write clonotype info on all tumors
##Read in tumor object
sub <- readRDS("results/Subclustering_B/tumor_integrated.rds")
sub_meta <- sub@meta.data
##Read in de-identification matrix
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
sum(sub_meta$CaTissueID %in% deid$CaTissueID)/nrow(sub_meta)
##1
sub_meta$Deidentified_SampleID <- deid$PatientDeIdentified[match(sub_meta$CaTissueID, deid$CaTissueID)]
sub_meta$Tumor <- ifelse(sub_meta$Annotation == "Tumor", paste0(sub_meta$Tumor_Type2, "_", sub_meta$Deidentified_SampleID), as.character(sub_meta$Annotation))
sub_meta$GeneLevelClonotype <- paste0(sub_meta$IgH_gene_ct, ";", sub_meta$IgLC_gene_ct)
clon_freq <- data.frame(table(sub_meta$GeneLevelClonotype, sub_meta$Tumor))
clon_freq <- clon_freq[!clon_freq$Var1 %in% c("NA;NA"),]
clon_freq <- clon_freq[clon_freq$Freq>0,]
for(i in 1:length(unique(clon_freq$Var2))){
  p <- as.character(unique(clon_freq$Var2))[i]
  tmp <- clon_freq[clon_freq$Var2==p,]
  if(length(grep("NA", tmp$Var1))>0){
    tmp <- tmp[-grep("NA", tmp$Var1),]
  }
  if(i==1){
    out <- tmp[which.max(tmp$Freq),]
  } else {
    out <- rbind(out, tmp[which.max(tmp$Freq),])
  }
}
out$Freq <- NULL
colnames(out) <- c("Clonotype", "Tumor")
write.csv(out, paste0(output_dir, "Tumor_Clonotypes.csv"), row.names = F)

##Plot clonotype concordance
out <- read.csv(paste0(output_dir, "Tumor_Clonotypes.csv"))
out$HC <- unlist(lapply(out$Clonotype, function(x) {str_split(x, "\\;")[[1]][1]}))
out$LC <- unlist(lapply(out$Clonotype, function(x) {str_split(x, "\\;")[[1]][2]}))
out$heavy_v_gene <- unlist(lapply(out$HC, function(x) {str_split(x, "\\.")[[1]][1]}))
out$heavy_j_gene <- unlist(lapply(out$HC, function(x) {str_split(x, "\\.")[[1]][2]}))
out$heavy_c_gene <- unlist(lapply(out$HC, function(x) {str_split(x, "\\.")[[1]][3]}))
out$light_v_gene <- unlist(lapply(out$LC, function(x) {str_split(x, "\\.")[[1]][1]}))
out$light_j_gene <- unlist(lapply(out$LC, function(x) {str_split(x, "\\.")[[1]][2]}))
out$light_c_gene <- unlist(lapply(out$LC, function(x) {str_split(x, "\\.")[[1]][3]}))
#Visualize in heatmap
data <- out
data$sample <- data$Tumor
ig_genes <- unique(as.character(c(data$heavy_v_gene, data$heavy_j_gene, data$heavy_c_gene, data$light_v_gene, data$light_j_gene, data$light_c_gene)))
out <- data.frame("Gene" = as.character(ig_genes))
for(i in 1:nrow(data)){
  tmp <- as.vector(t(data[i, c("heavy_v_gene", "heavy_j_gene", "heavy_c_gene", "light_v_gene", "light_j_gene", "light_c_gene")]))
  sample_id <- as.character(data$sample)[i]
  out[, sample_id] <- ifelse(out$Gene %in% tmp, 1, 0)
}
out_long <- gather(out, "Sample", "Count", -Gene)
out_long$clone <- out_long$Sample
data$lightchain <- "Kappa"
data$lightchain[grep("IGLC", data$light_c_gene)] <- "Lambda"
data$lightchain <- factor(data$lightchain, levels = c("Kappa", "Lambda"))
data <- data[order(data$sample, decreasing = T),]
data <- data[order(data$lightchain, decreasing = T),]
data$Type <- "WM"
data$Type[grep("CLL", data$Tumor)] <- "CLL"
data$Type[grep("UCL", data$Tumor)] <- "UCL"
out_long$Type <- data$Type[match(out_long$Sample, data$Tumor)]
data$sample <- factor(data$sample, levels = data$sample)
out_long$Sample <- factor(out_long$Sample, levels = data$sample)
out_long$dummy <- 1
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
final_cols <- c(col_vector, col_vector[1:7])
a <- ggplot(out_long) +
  geom_tile(aes(x=Sample, y=Gene, fill= factor(Count)), color = "black") +
  scale_fill_manual(values = c("white", "red"), name = "", guide = "none") + 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill=NA, color="black"), 
        axis.title = element_text(size=14), 
        axis.text.y = element_text(size=12, color="black", face = "italic"), 
        legend.text = element_text(size=10), 
        axis.text.x = element_text(angle=65, hjust = 1, size = 14)) +
  ylab("") +
  xlab("")
b <- ggplot(out_long) +
  geom_tile(aes(x=Sample, y=dummy, fill=Type), color = "black", show.legend = T) + 
  scale_fill_manual(values = final_cols) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill=NA, color="black"), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()) +
  scale_y_discrete(expand = c(0,0))
b_leg <- get_legend(b)
b <- b + theme(legend.position = "none")
plots <- plot_grid(b, a, align = "v", axis = "lr", ncol = 1, rel_heights = c(0.05, 1))
png(paste0(output_dir, "ClonotypeConcordance.png"), res = 300, units="in", width = 10, height = 10)
plot_grid(plots, b_leg, ncol = 2, rel_widths = c(1, 0.1))
dev.off()

#Supp. Figure 3B
##Plot tumor cell numbers per sample
##Read in tumor object
sub <- readRDS("results/Subclustering_B/tumor_integrated.rds")
##Keep only WM
wm <- subset(sub, Tumor_Type %in% "WM")
plotdf <- FetchData(wm, vars= c("Cell", "UMAP_1", "UMAP_2", "Demux_Sample_ID", "CaTissueID", "NewClonotypeID", "IgH_gene_ct", "IgLC_gene_ct", "cdr3_aa", "cluster", "Annotation", "Tumor_Type", "Tumor_Type2", "Disease"))
##Read in de-identification matrix
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
sum(plotdf$CaTissueID %in% deid$CaTissueID)/nrow(plotdf)
##1
plotdf$Deidentified_SampleID <- deid$PatientDeIdentified[match(plotdf$CaTissueID, deid$CaTissueID)]
plotdf$Tumor <- ifelse(plotdf$Annotation == "Tumor", paste0(plotdf$Tumor_Type2, "_", plotdf$Deidentified_SampleID), as.character(plotdf$Annotation))
cell_freq <- data.frame(table(plotdf$Tumor))
colnames(cell_freq) <- c("Sample", "Count")
sum(cell_freq$Count)
##48875
cell_freq$Sample <- factor(cell_freq$Sample, levels = cell_freq$Sample[order(cell_freq$Count, decreasing = F)])
write.csv(cell_freq, paste0(output_dir, "TumorCellNumbers_Data.csv"), row.names = F)
png(paste0(output_dir, "TumorCellNumbers.png"), res = 300, units="in", width = 6, height = 5)
ggplot(cell_freq) +
  geom_bar(aes(Count, Sample), stat = "identity", fill = "lightblue", color = "black", width = 0.5) +
  scale_x_log10() +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(size = 14)) +
  ylab("") +
  xlab("# of tumor cells") +
  geom_vline(xintercept = median(cell_freq$Count, na.rm = T), linetype = "dashed", color = "red")
dev.off()
summary(cell_freq$Count, na.rm = T)
##Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##106.0   374.5  1088.5  2036.5  2353.8  8694.0

#Supp. Figure 3C
##Plot normal B cell UMAP
sub2 <- readRDS("results/Subclustering_B/normal_integrated.rds")
annotation <- c('10' = 'TBC', '12' = 'TBC', '18' = 'TBC', '0' = 'NBC', '13' = 'NBC', '1' = 'NCSMBC', '4' = 'NCSMBC', '9' = 'aBC', '6' = 'MZBC', '11' = 'IgE+ MBC', '3' = 'TCF7+ MBC', '5' = 'S100A10+ MBC', '2' = 'TCF4+ MBC', '7' = 'ABC')
renamed <- RenameIdents(sub2, annotation)
renamed_sub <- subset(renamed, idents=unname(annotation))
renamed_sub$Cell_Type <- Idents(renamed_sub)
dim(renamed_sub)
##36601 28991
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
plotdf$Cell_ID <- rownames(plotdf)
centroids <- aggregate(plotdf[, c("UMAP_1","UMAP_2")], list(factor(plotdf$Cell_Type)), median)
colnames(centroids)[1] <- "Cluster"
set.seed(123)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
final_cols <- col_vector[1:length(unique(plotdf$Cell_Type))]

png(paste0(output_dir, "NormalBCells_UMAP.png"), res = 300, units="in",width = 6, height = 5)
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
  coord_cartesian(clip="off") +
  xlim(min(plotdf$UMAP_1, na.rm = T), 6)
dev.off()

#Supp. Figure 3D
##Plot Normal B cell marker heatmap
sub2 <- readRDS("results/Subclustering_B/normal_integrated.rds")
annotation <- c('10' = 'TBC', '12' = 'TBC', '18' = 'TBC', '0' = 'NBC', '13' = 'NBC', '1' = 'NCSMBC', '4' = 'NCSMBC', '9' = 'aBC', '6' = 'MZBC', '11' = 'IgE+ MBC', '3' = 'TCF7+ MBC', '5' = 'S100A10+ MBC', '2' = 'TCF4+ MBC', '7' = 'ABC')
renamed <- RenameIdents(sub2, annotation)
renamed_sub <- subset(renamed, idents=unname(annotation))
renamed_sub$Cell_Type <- Idents(renamed_sub)
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
markers <- c("CD38", "SOX4", "TCL1A", "IGHD", "IL4R", "ABCB1", "FCER2", "PLPP5", "CCR7", "CD27", "TNFRSF13B", "IGHM", "CD1C", "PLD4", "TCF4", "SYK", "RASSF6", "JCHAIN", "IGHG1", "IGHG4", "IGHA1", "IGHA2", "COCH", "TCF7", "TEX9", "PDE4D", "S100A10", "ANXA2", "ANXA4", "LGALS2", "CD99", "ITGB1", "HOMER3", "IGHE", "TBX21", "ITGAX", "FCRL3", "FCRL5", "JUN", "FOS", "FOSB", "DUSP1", "TNFAIP3")
orig_cell_types <- c('TBC', 'NBC', 'NCSMBC', 'MZBC', 'TCF4+ MBC', 'TCF7+ MBC', 'S100A10+ MBC', 'IgE+ MBC', 'ABC', 'aBC')
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

png(paste0(output_dir, "NormalBCells_Annotation_Heatmap.png"), res = 300, units="in", width = 11, height = 5)
ggplot(marker_means_long) +
  geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1, size=12, face="italic", color="black"), axis.text.y=element_text(size=12, color="black")) + 
  xlab("") + 
  ylab("") + 
  labs(fill="Mean Exp")
dev.off()

#Supp. Figure 3E
##Generate & write Normal B cell type count & proportion table
sub2 <- readRDS("results/Subclustering_B/normal_integrated.rds")
annotation <- c('10' = 'TBC', '12' = 'TBC', '18' = 'TBC', '0' = 'NBC', '13' = 'NBC', '1' = 'NCSMBC', '4' = 'NCSMBC', '9' = 'aBC', '6' = 'MZBC', '11' = 'IgE+ MBC', '3' = 'TCF7+ MBC', '5' = 'S100A10+ MBC', '2' = 'TCF4+ MBC', '7' = 'ABC')
renamed <- RenameIdents(sub2, annotation)
renamed_sub <- subset(renamed, idents=unname(annotation))
renamed_sub$Cell_Type <- Idents(renamed_sub)
plotdf <- FetchData(renamed_sub, vars= c("Cell","UMAP_1","UMAP_2", "Demux_Sample_ID", "Cell_Type"))
t_count <- data.frame(table(plotdf$Cell_Type, plotdf$Demux_Sample_ID))
colnames(t_count) <- c("Cell_Type","Sample","Count")
t_count <- spread(t_count,Cell_Type,Count)
rownames(t_count) <- t_count$Sample
t_count$Sample <- NULL
write.csv(t_count, paste0(output_dir, "NormalBCell_count_table.csv"), row.names = T)
t_prop <- t_count/rowSums(t_count, na.rm = T)
write.csv(t_prop, paste0(output_dir, "NormalBCell_prop_table.csv"), row.names = T)

##Compare Normal B cell proportions between AWM patients and HD
threshold <- 50
##Read in normal B cell count table
tmp <- read.csv(paste0(output_dir, "NormalBCell_count_table.csv"), row.names = 1, check.names = F)
sum(tmp, na.rm = T)
##[1] 28991
##Read in metadata
metadata <-read.csv(paste0(input_dir, "wm-cohort-metadata.csv"), check.names = F, na.strings = "")
##Aggregate cell counts by CaTissue ID
tmp$CaTissueID <- metadata$CaTissueID[match(rownames(tmp), metadata$Demux_Sample_ID_Final)]
df <- transform(aggregate(tmp[,!colnames(tmp) %in% c("CaTissueID")], list(tmp$CaTissueID), sum, na.rm=T), row.names=Group.1, Group.1=NULL)
colnames(df) <- colnames(tmp[,!colnames(tmp) %in% c("CaTissueID")])
nrow(df)
##40
##Remove samples with < 50 cells across immune subpopulations
df <- df[rowSums(df, na.rm = T)>= threshold,]
nrow(df)
##27
##Compute proportions
df <- df/rowSums(df, na.rm = T)
##Compare proportions between AWM and NBM using Wilcoxon tests
df$Status <- metadata$Disease[match(rownames(df), metadata$CaTissueID)]
df <- df[df$Status %in% c("AWM", "NBM"),]
nrow(df)
##23
table(df$Status)
##AWM NBM 
##13  10
df$CaTissueID <- rownames(df)
df_long <- gather(df, "Cell_Type", "Proportion", -Status, -CaTissueID)
df_long$Cell_Type <- factor(df_long$Cell_Type)
df_long$Status <- as.character(df_long$Status)
df_long$Status[df_long$Status %in% "NBM"] <- "HD"
df_long$Status <- factor(df_long$Status, levels = c("HD", "AWM"))
df_long$Cell_Type <- factor(df_long$Cell_Type, levels = c('TBC', 'NBC', 'NCSMBC', 'MZBC', 'TCF4+ MBC', 'TCF7+ MBC', 'S100A10+ MBC', 'IgE+ MBC', 'ABC', 'aBC'))
write.csv(df_long, paste0(output_dir, "NormalBCells_Proportions_Data.csv"), row.names = F)
out <- compare_means(Proportion ~ Status, group.by="Cell_Type", data=df_long, p.adjust.method = "BH")
write.csv(out, paste0(output_dir, "NormalBCells_Proportions_AWMvsHD_Stats.csv"), row.names = F)
maxes <- aggregate(df_long$Proportion, list(df_long$Cell_Type), max)
out$ypos <- maxes$x[match(out$Cell_Type, maxes$Group.1)]
png(paste0(output_dir, "NormalBcells_Proportions_AWMvsHD.png"), res = 300, units="in", width = 8, height = 5)
ggplot(df_long) +
  geom_boxplot(aes(Status, Proportion*100, fill = Cell_Type), outlier.size = -1, alpha =0.5, show.legend = F) +
  geom_violin(aes(Status, Proportion*100, fill = Cell_Type), alpha = 0.5, show.legend = F) +
  geom_point(aes(Status, Proportion*100, fill = Cell_Type), shape = 21, size = 3, show.legend = F, position = position_jitter(width = 0.2)) +
  scale_fill_brewer(palette = "Set3", name = "Cell Type") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 14), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "italic")) +
  xlab("") +
  ylab("% of B cells") +
  geom_bracket(aes(xmin = group1, xmax = group2, y.position = ypos*100 + 5, label=paste0("q=", signif(p.adj,2))), data = out[out$p.adj < 0.1,], vjust = 1.5) +
  facet_wrap(.~Cell_Type, scales = "free_y", ncol = 5)
dev.off()

#Supp. Figure 3F
##Plot CXCR4 expression by mutation status
plotdf <- read.csv("results/Figure3/CXCR4_Data.csv", row.names = 1)
sub <- readRDS("results/Subclustering_B/tumor_integrated.rds")
cxcr4 <- data.frame(t(data.frame(GetAssayData(sub["CXCR4",], slot = "data"))))
rownames(cxcr4) <- gsub("^X", "", rownames(cxcr4))
rownames(cxcr4) <- gsub("\\.", "-", rownames(cxcr4))
plotdf$CXCR4_exp <- cxcr4$CXCR4[match(rownames(plotdf), rownames(cxcr4))]
sum(is.na(plotdf$CXCR4_exp))
##0
cxcr4_agg <- aggregate(plotdf$CXCR4_exp, list(plotdf$Tumor), median)
plotdf$Tumor <- factor(plotdf$Tumor, levels = cxcr4_agg$Group.1[order(cxcr4_agg$x, decreasing = T)])
meta <- read.csv("data/wm-cohort-metadata.csv")
sum(plotdf$Demux_Sample_ID %in% meta$Demux_Sample_ID)/nrow(plotdf)
##1
plotdf$BM <- meta$BM_Infiltration[match(plotdf$Demux_Sample_ID, meta$Demux_Sample_ID)]
plotdf$CXCR4_mut <- ifelse(plotdf$CXCR4_mut == "WT" & plotdf$BM < 0.3, "WT & Infiltr. < 30%",
                           ifelse(plotdf$CXCR4_mut == "WT" & plotdf$BM >= 0.3, "WT & Infiltr. >= 30%", as.character(plotdf$CXCR4_mut)))
##Pt5, who has a serial sample should be considered WT & >= 30% based on their baseline sample with a 60% infiltration
deid <- read.csv("data/SampleTable_DeIdentified.csv")
sum(plotdf$CaTissueID %in% deid$CaTissueID)/nrow(plotdf)
##1
plotdf$Deidentified_PatientID <- deid$PatientDeIdentified[match(plotdf$CaTissueID, deid$CaTissueID)]
plotdf$CXCR4_mut[plotdf$Deidentified_PatientID == "Pt5"] <- "WT & Infiltr. >= 30%" 
plotdf$CXCR4_mut <- factor(plotdf$CXCR4_mut, levels = c("Mut", "WT & Infiltr. >= 30%", "WT & Infiltr. < 30%", "Unknown"))
write.csv(plotdf, paste0(output_dir, "CXCR4_expression_by_Mut_Data.csv"), row.names = T)
png(paste0(output_dir, "CXCR4_expression_by_Mut_AtLeast30Infiltration.png"), res = 300, units="in", width = 7, height = 5)
ggplot(plotdf) +
  geom_boxplot(aes(Tumor, CXCR4_exp, fill=CXCR4_mut), show.legend = F) +
  scale_fill_manual(values = c("orange", "steelblue", "lightblue", "gray")) +
  theme(axis.text.x = element_text(angle=65, hjust =1, color = "black", size = 10), 
        axis.text.y = element_text(color = "black", size = 12), 
        strip.background = element_blank(), 
        strip.text = element_text(face = "italic", size = 14), 
        axis.title = element_text(size = 14)) +
  xlab("") +
  ylab("CXCR4 expression") +
  facet_wrap(.~CXCR4_mut, scales = "free_x")
dev.off()  
