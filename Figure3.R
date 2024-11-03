#Import dependencies
library(tidyverse)
library(RColorBrewer)
library(Seurat)
library(cowplot)
library(ggrepel)
library(ggpubr)

#Set input/output directory
input_dir <- "data/"
output_dir <- "results/Figure3/"

#Set seed
set.seed(123)

#Figure 3A
##Plot tumor cell UMAPs
##Read in tumor object
sub <- readRDS("results/Subclustering_B/tumor_integrated.rds")
##Plot annotated tumor UMAP
plotdf <- FetchData(sub, vars= c("Cell", "UMAP_1", "UMAP_2", "Demux_Sample_ID", "CaTissueID", "NewClonotypeID", "IgH_gene_ct", "IgLC_gene_ct", "cdr3_aa", "cluster", "Annotation", "Tumor_Type", "Tumor_Type2", "Disease"))
plotdf$Tumor_Type <- factor(plotdf$Tumor_Type, levels = c("WM", "CLL", "UCL"))
png(paste0(output_dir, "Annotated_UMAP_Tumor_Type.png"), res = 300, units="in", width = 8, height = 5)
ggplot(plotdf) +
  geom_point(aes(UMAP_1, UMAP_2, color = Tumor_Type), stroke = 0.1) + 
  scale_color_manual(values = c("steelblue", "orange", "pink"), name = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.text = element_text(size = 12))
dev.off()

##Plot tumor-specific UMAP
##Read in de-identification matrix
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
sum(plotdf$CaTissueID %in% deid$CaTissueID)/nrow(plotdf)
##1
plotdf$Deidentified_SampleID <- deid$PatientDeIdentified[match(plotdf$CaTissueID, deid$CaTissueID)]
plotdf$Tumor <- ifelse(plotdf$Annotation %in% "Tumor", paste0(plotdf$Tumor_Type2, "_", plotdf$Deidentified_SampleID), as.character(plotdf$Annotation))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
final_cols <- col_vector[c(1, 10, 2, 9, 3, 8, 4, 7, 5, 6, 12:13, 16:18, 25, 27, 29, 31, 33, 40, 44, 46, 14, 15, 28, 30, 32, 36, 39, 42)]
centroids <- aggregate(plotdf[, c("UMAP_1","UMAP_2")], list(factor(plotdf$Tumor)), median)
colnames(centroids)[1] <- "Cluster"
png(paste0(output_dir, "Annotated_UMAP_TumorSpecific.png"), res = 300, units="in", width = 8, height = 5)
ggplot(plotdf) +
  geom_point(aes(UMAP_1, UMAP_2, color = Tumor), alpha = 0.5, stroke = 0.1, size = 0.5, show.legend = F) + 
  scale_color_manual(values = final_cols, name = "Status") +
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank()) +
  geom_text_repel(data = centroids, aes(x=UMAP_1, y=UMAP_2, label= Cluster), size=4, min.segment.length = 0, segment.alpha = 0.7)
dev.off()

#Figure 3B
##Plot diagnostics heatmap
##Read in tumor object
sub <- readRDS("results/Subclustering_B/tumor_integrated.rds")
sub_meta <- sub@meta.data
##Read in de-identification matrix
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
sum(sub_meta$CaTissueID %in% deid$CaTissueID)/nrow(sub_meta)
##1
sub_meta$Deidentified_SampleID <- deid$PatientDeIdentified[match(sub_meta$CaTissueID, deid$CaTissueID)]
sub_meta$Tumor <- ifelse(sub_meta$Annotation == "Tumor", paste0(sub_meta$Tumor_Type2, "_", sub_meta$Deidentified_SampleID), as.character(sub_meta$Annotation))
sub@meta.data <- sub_meta
##Plot marker heatmap
markers <- c("CD5", "FCER2", "SPN", "LEF1", "TCL1A", "MME", "IRF4", "BCL2", "BCL6", "CCND1", "SOX11", "ITGAX", "IL2RA", "CD38", "SDC1", "NCAM1", "IGHD", "IGHM", "CD79B", "FCRLA", "RALGPS2", "CD27", "TNFRSF13B", "ITM2C", "ITM2B", "SYK", "BLNK", "CD1C", "CD9", "DUSP22", "RAC2", "JCHAIN", "CD52", "CD53", "ZNF804A", "ITGB1", "RASSF6", "VPREB3")
order_tumors <- sub_meta[!duplicated(sub_meta$Tumor),] 
order_tumors$Tumor_Type <- factor(order_tumors$Tumor_Type, levels = c("CLL", "UCL", "WM"))
marker_counts <- data.frame(t(data.frame(GetAssayData(sub[markers,], slot = "data"))))
rownames(marker_counts) <- gsub("^X", "", rownames(marker_counts))
rownames(marker_counts) <- gsub("\\.", "-", rownames(marker_counts))
marker_counts <- data.frame(scale(marker_counts))
marker_counts$idents <- sub_meta$Tumor[match(rownames(marker_counts), rownames(sub_meta))]
marker_means <- aggregate(marker_counts[, 1:(ncol(marker_counts)-1)], list(marker_counts$idents), mean)
rownames(marker_means) <- marker_means$Group.1
marker_means$Group.1  <- NULL
marker_means_t <- data.frame(t(marker_means))
colnames(marker_means_t) <- rownames(marker_means)
marker_means_t$Marker <- rownames(marker_means_t)
marker_means_long <- gather(marker_means_t, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels = order_tumors$Tumor[order(order_tumors$Tumor_Type, decreasing = F)])
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)
segment_data_top <- data.frame("x"=c(0.5,5.5,17.5), "x_end"= c(5.5, 17.5, 38.5), "y" = c(5.5,6.5,30.5), "y_end" = c(5.5,6.5,30.5))
segment_data_bottom <- data.frame("x"=c(0.5,5.5,17.5), "x_end"= c(5.5, 17.5, 38.5), "y" = c(0.5,5.5,6.5), "y_end" = c(0.5,5.5,6.5))
segment_data_left <- data.frame("x"=c(0.5, 5.5, 17.5), "x_end"=c(0.5, 5.5, 17.5), "y" = c(0.5, 5.5, 6.5), "y_end" = c(5.5,6.5,30.5))
segment_data_right <- data.frame("x"=c(5.5, 17.5, 38.5), "x_end"=c(5.5, 17.5, 38.5),  "y" = c(0.5, 5.5, 6.5), "y_end" = c(5.5,6.5,30.5))
png(paste0(output_dir, "Diagnostics_Heatmap.png"), res = 300, units="in", width = 11, height = 7)
ggplot(marker_means_long) +
  geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 65, hjust=1, size=12, face="italic", color="black"), 
        axis.text.y = element_text(size=12, color="black")) + 
  xlab("") + 
  ylab("") + 
  labs(fill="Mean Exp") +
  geom_segment(data=segment_data_top, aes(x=x, xend=x_end, y=y, yend=y_end), linetype="dashed", size=0.8) + 
  geom_segment(data=segment_data_bottom, aes(x=x, xend=x_end, y=y, yend=y_end), linetype="dashed", size=0.8) +
  geom_segment(data=segment_data_left, aes(x=x, xend=x_end, y=y, yend=y_end), linetype="dashed", size=0.8) +
  geom_segment(data=segment_data_right, aes(x=x, xend=x_end, y=y, yend=y_end), linetype="dashed", size=0.8)
dev.off()

#Figure 3C
##Plot UMAP colored by CXCR4 mutation status
##Read in tumor object
sub <- readRDS("results/Subclustering_B/tumor_integrated.rds")
plotdf <- FetchData(sub, vars= c("Cell", "UMAP_1", "UMAP_2", "Demux_Sample_ID", "CaTissueID", "NewClonotypeID", "IgH_gene_ct", "IgLC_gene_ct", "cdr3_aa", "cluster", "Annotation", "Tumor_Type", "Tumor_Type2", "Disease"))
metadata <- read.csv(paste0(input_dir, "wm-cohort-metadata.csv"))
sum(plotdf$CaTissueID %in% metadata$CaTissueID)/nrow(plotdf)
##1
plotdf$CXCR4 <- metadata$CXCR4[match(plotdf$CaTissueID, metadata$CaTissueID)]
##Read in de-identification matrix
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
sum(plotdf$CaTissueID %in% deid$CaTissueID)/nrow(plotdf)
##1
plotdf$Deidentified_SampleID <- deid$PatientDeIdentified[match(plotdf$CaTissueID, deid$CaTissueID)]
plotdf$Tumor <- ifelse(plotdf$Annotation %in% "Tumor", paste0(plotdf$Tumor_Type2, "_", plotdf$Deidentified_SampleID), as.character(plotdf$Annotation))
##Only assign CXCR4 status to primary tumors
plotdf$CXCR4[plotdf$Tumor_Type2 %in% c("WM_2", "WM_3")] <- NA
##Keep WM tumors 
plotdf <- plotdf[plotdf$Tumor_Type %in% "WM",]
plotdf$CXCR4_mut <- ifelse(is.na(plotdf$CXCR4), "Unknown", as.character(plotdf$CXCR4))
write.csv(plotdf, paste0(output_dir, "CXCR4_Data.csv"), row.names = T)
plotdf <- read.csv(paste0(output_dir, "CXCR4_Data.csv"))
plotdf$CXCR4_mut <- factor(plotdf$CXCR4_mut, levels = c("Mut", "WT", "Unknown"))
a <- ggplot(plotdf) +
  geom_point(aes(UMAP_1, UMAP_2, color = CXCR4_mut)) +
  scale_color_manual(values = c("orange", "steelblue", "gray"), labels = c("CXCR4-mut", "CXCR4-WT", "Unknown"), name = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        legend.text = element_text(size = 8))
a_leg <- get_legend(a)
b <- ggplot(plotdf) +
  geom_point(aes(UMAP_1, UMAP_2, color = CXCR4_mut), size = 0.01, alpha = 0.5, show.legend = F) +
  scale_color_manual(values = c("orange", "steelblue", "gray"), name = "Status") +
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank())

png(paste0(output_dir, "UMAP_CXCR4_status.png"), res = 300, units="in", width = 7, height = 5)
plot_grid(b, a_leg, align = "h", axis = "tb", rel_widths = c(0.5, 0.1))
dev.off()

#Figure 3D
##Plot CXCR4 expression by mutation status
plotdf <- read.csv(paste0(output_dir, "CXCR4_Data.csv"), row.names = 1)
sub <- readRDS("results/Subclustering_B/tumor_integrated.rds")
cxcr4 <- data.frame(t(data.frame(GetAssayData(sub["CXCR4",], slot = "data"))))
rownames(cxcr4) <- gsub("^X", "", rownames(cxcr4))
rownames(cxcr4) <- gsub("\\.", "-", rownames(cxcr4))
plotdf$CXCR4_exp <- cxcr4$CXCR4[match(rownames(plotdf), rownames(cxcr4))]
sum(is.na(plotdf$CXCR4_exp))
##0
cxcr4_agg <- aggregate(plotdf$CXCR4_exp, list(plotdf$Tumor), median)
plotdf$Tumor <- factor(plotdf$Tumor, levels = cxcr4_agg$Group.1[order(cxcr4_agg$x, decreasing = T)])
plotdf$CXCR4_mut <- factor(plotdf$CXCR4_mut, levels = c("Mut", "WT", "Unknown"))
df <- data.frame(table(plotdf$Tumor))
df$Var1 <- factor(df$Var1, levels = levels(plotdf$Tumor))
df <- df[order(df$Var1),]
df$CXCR4_mut <- plotdf$CXCR4_mut[match(df$Var1, plotdf$Tumor)]
df$CXCR4_mut <- factor(df$CXCR4_mut, levels = c("Mut", "WT", "Unknown"))
df <- df[order(df$CXCR4_mut),]
sums <- data.frame(table(df$CXCR4_mut))
df$x <- c(1:sums$Freq[sums$Var1 == levels(df$CXCR4_mut)[1]], 1:sums$Freq[sums$Var1 == levels(df$CXCR4_mut)[2]], 1:sums$Freq[sums$Var1 == levels(df$CXCR4_mut)[3]])
png(paste0(output_dir, "CXCR4_expression_by_Mut.png"), res = 300, units="in", width = 7, height = 5)
ggplot(plotdf[!is.na(plotdf$CXCR4_mut),]) +
  geom_boxplot(aes(Tumor, CXCR4_exp, fill=CXCR4_mut), show.legend = F) +
  scale_fill_manual(values = c("orange", "steelblue", "gray")) +
  theme(axis.text.x = element_text(angle=65, hjust =1, color = "black", size = 10), 
        axis.text.y = element_text(color = "black", size = 12), 
        strip.background = element_blank(), 
        strip.text = element_text(face = "italic", size = 14), 
        axis.title = element_text(size = 14)) +
  xlab("") +
  ylab("CXCR4 expression") +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(min(plotdf$CXCR4_exp)-0.8, 4.5)) + 
  geom_text(data = df, aes(x = x, y = min(plotdf$CXCR4_exp)-0.15, label = paste0("n=", Freq)), size = 3, angle = 90, hjust = 1) +
  facet_wrap(.~CXCR4_mut, scales = "free_x")
dev.off()  

#Figure 3E
##Compare median CXCR4 expression in paired tumor and normal memory B cells 
##Read in tumor object
sub <- readRDS("results/Subclustering_B/tumor_integrated.rds")
##Read in normal B cell object and subset for patient memory B cells
sub2 <- readRDS("results/Subclustering_B/normal_integrated.rds")
annotation <- c('10' = 'TBC', '12' = 'TBC', '18' = 'TBC', '0' = 'NBC', '13' = 'NBC', '1' = 'NCSMBC', '4' = 'NCSMBC', '9' = 'aBC', '6' = 'MZBC', '11' = 'IgE+ MBC', '3' = 'TCF7+ MBC', '5' = 'S100A10+ MBC', '2' = 'TCF4+ MBC', '7' = 'ABC')
renamed <- RenameIdents(sub2, annotation)
renamed_sub <- subset(renamed, idents = c("MZBC", "NCSMBC", "TCF7+ MBC", "S100A10+ MBC", "TCF4+ MBC", "IgE+ MBC"))
renamed_sub$Cell_Type <- Idents(renamed_sub)
renamed_sub <- subset(renamed_sub, Disease == "AWM")
##Combine two fractions and extract CXCR4 expression info
full <- merge(sub, renamed_sub)
cxcr4 <- data.frame(t(data.frame(GetAssayData(full["CXCR4",], slot = "data"))))
rownames(cxcr4) <- gsub("^X", "", rownames(cxcr4))
rownames(cxcr4) <- gsub("\\.", "-", rownames(cxcr4))
plotdf <- FetchData(full, vars= c("Cell", "Demux_Sample_ID", "CaTissueID", "NewClonotypeID", "IgH_gene_ct", "IgLC_gene_ct", "cdr3_aa", "cluster", "Annotation", "Tumor_Type", "Tumor_Type2", "Disease"))
plotdf$CXCR4_exp <- cxcr4$CXCR4[match(rownames(plotdf), rownames(cxcr4))]
##Annotate with CXCR4 mutation status
metadata <- read.csv(paste0(input_dir, "wm-cohort-metadata.csv"))
sum(plotdf$CaTissueID %in% metadata$CaTissueID)/nrow(plotdf)
##1
plotdf$CXCR4_mut <- metadata$CXCR4[match(plotdf$CaTissueID, metadata$CaTissueID)]
##Remove CLL, UCL 
plotdf <- plotdf[!plotdf$Tumor_Type %in% c("CLL", "UCL"),]
##Only assign CXCR4 status to primary tumors
plotdf$CXCR4_mut[plotdf$Tumor_Type2 %in% c("WM_2", "WM_3")] <- NA
##Read in de-identification matrix
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
sum(plotdf$CaTissueID %in% deid$CaTissueID)/nrow(plotdf)
##1
plotdf$Deidentified_SampleID <- deid$PatientDeIdentified[match(plotdf$CaTissueID, deid$CaTissueID)]
plotdf$Tumor <- paste0(plotdf$Tumor_Type2, "_", plotdf$Deidentified_SampleID)
##Remove secondary tumors and patient with no CXCR4 status annotation
table(plotdf$Tumor[is.na(plotdf$CXCR4_mut)])
##Normal_Pt7  WM_2_Pt13   WM_2_Pt3   WM_2_Pt8  WM_3_Pt13     WM_Pt7 
##2        572        688       1481        106        307
plotdf <- plotdf[!is.na(plotdf$CXCR4_mut),]
##Compute median CXCR4 expression by patient and compartment
cxcr4_agg <- aggregate(plotdf$CXCR4_exp, list(plotdf$Tumor), median)
cxcr4_agg$status <- plotdf$CXCR4_mut[match(cxcr4_agg$Group.1, plotdf$Tumor)]
cxcr4_agg$Deidentified_SampleID <- plotdf$Deidentified_SampleID[match(cxcr4_agg$Group.1, plotdf$Tumor)]
cxcr4_agg$Annotation <- plotdf$Annotation[match(cxcr4_agg$Group.1, plotdf$Tumor)]
cxcr4_agg$Group.1 <- NULL
write.csv(cxcr4_agg, paste0(output_dir, "Median_CXCR4_By_Compartment_Data.csv"), row.names = F)
cxcr4_agg <- read.csv(paste0(output_dir, "Median_CXCR4_By_Compartment_Data.csv"))
paired_df <- spread(cxcr4_agg, "Annotation", "x")
nrow(paired_df)
##20
##Remove patients missing a fraction (n=5)
paired_df <- paired_df[!is.na(paired_df$Normal) & !is.na(paired_df$Tumor),]
nrow(paired_df)
##15
wilcox.test(paired_df$Tumor, paired_df$Normal, paired = T)
##p-value = 0.02801
stats <- data.frame("group1" = "Normal", "group2" = "Tumor", "p" = wilcox.test(paired_df$Tumor, paired_df$Normal, paired = T)$p.value)
cxcr4_agg$status <- factor(cxcr4_agg$status, levels = c("Mut", "WT"))
cxcr4_agg <- cxcr4_agg[cxcr4_agg$Deidentified_SampleID %in% paired_df$Deidentified_SampleID,]
png(paste0(output_dir, "Median_CXCR4_expression_TumorVSNormal.png"), res = 300, units="in", width = 5, height = 5)
ggplot(cxcr4_agg) + 
  geom_boxplot(aes(Annotation, x, fill = Annotation), alpha = 0.5, outlier.size = -1, show.legend = F) + 
  geom_violin(aes(Annotation, x, fill = Annotation), alpha = 0.5, show.legend = F) + 
  geom_point(aes(Annotation, x, color = status), size = 4, alpha = 0.8, position = position_jitter(width = 0.3), show.legend = T) +
  scale_fill_manual(values = c("lemonchiffon1", "tomato3"), name = "", guide = "none") +
  scale_color_manual(values = c("orange", "steelblue", "gray"), name = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.x = element_text(color = "black", size = 14), 
        axis.text.y = element_text(color = "black", size = 12), 
        axis.title = element_text(size = 14)) +
  xlab("") +
  ylab("CXCR4 expression") +
  scale_x_discrete(labels = c("Normal MBCs", "Tumor cells")) +
  geom_bracket(aes(xmin = group1, xmax = group2, label = paste0("p=", signif(p, 2))), y.position = 3.2, vjust = 1.5, label.size = 5, step.increase = 0.1, data = stats)
dev.off()

#Figure 3F
##Plot heatmap of CNV status
##List files
files <- list.files(paste0(input_dir, "CNV/"))
##Read in de-identification matrix
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
dict <- data.frame("sample" = files, "CaTissueID" = gsub("_.*", "", files), "Tumor" = "WM")
dict$Tumor[grep("WM_2", dict$sample)] <- "WM_2"
dict$DeIdentifiedSampleID <- deid$PatientDeIdentified[match(dict$CaTissueID, deid$CaTissueID)]
dict$deidentified_ID <- paste0(dict$Tumor, "_", dict$DeIdentifiedSampleID)
##Prepare seg file
for(i in seq_along(files)){
  seg <- read.delim(paste0(input_dir, "CNV/", files[i], "/segs_consensus_2.tsv"))
  #remove events with length < 100
  seg <- seg[seg$seg_length > 100,]
  #create state variable
  seg$state <- seg$cnv_state
  seg$state <- ifelse(seg$state %in% "bamp", "amp", as.character(seg$state))
  #set ymin, ymax
  seg$ymin = 0
  seg$ymax = 1
  #create sample variable
  seg$sample <- dict[i, "deidentified_ID"]
  if(i==1){
    final_seg <- seg
  } else {
    final_seg <- rbind(final_seg, seg)
  }
}
final_seg$sample <- factor(final_seg$sample, levels = dict$deidentified_ID[order(dict$DeIdentifiedSampleID)])
write.csv(final_seg, paste0(output_dir, "Tri4_FinalSeg.csv"), row.names = F)
##Plot
final_seg <- read.csv(paste0(output_dir, "Tri4_FinalSeg.csv"))
png(paste0(output_dir, "Tri4_Heatmap.png"), res=300, unit="in", width=11, height=4)
ggplot(final_seg, 
       aes(xmin = seg_start, xmax = seg_end, ymin = ymin, ymax = ymax)) +
  geom_rect(aes(fill = state), colour = NA, size=0) +
  scale_fill_manual(values = c("#B2182B", "#2166AC", "orange"), limits = c("amp", "del", "loh"), labels = c("Amp", "Del", "LOH"), na.value = "white") +
  labs(fill = "CNV") +
  facet_grid(rows = vars(sample),
             cols = vars(CHROM),
             scales = "free", 
             space = "free",
             switch = "y") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(0.2, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(size=.3, fill = NA),
        strip.text.y.left = element_text(angle = 0, size = 14, color = "black"),
        strip.text.x.top = element_text(angle = 0, size = 7, color = "black"),
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0.05,0))
dev.off()

#Figure 3G
##Plot UMAP colored by Trisomy 4
##Read in tumor object
sub <- readRDS("results/Subclustering_B/tumor_integrated.rds")
plotdf <- FetchData(sub, vars= c("Cell", "UMAP_1", "UMAP_2", "Demux_Sample_ID", "CaTissueID", "NewClonotypeID", "IgH_gene_ct", "IgLC_gene_ct", "cdr3_aa", "cluster", "Annotation", "Tumor_Type", "Tumor_Type2", "Disease"))
##Read in de-identification matrix
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
sum(plotdf$CaTissueID %in% deid$CaTissueID)/nrow(plotdf)
##1
plotdf$Deidentified_SampleID <- deid$PatientDeIdentified[match(plotdf$CaTissueID, deid$CaTissueID)]
plotdf$Tumor <- ifelse(plotdf$Annotation %in% "Tumor", paste0(plotdf$Tumor_Type2, "_", plotdf$Deidentified_SampleID), as.character(plotdf$Annotation))
##Annotate tumors with Trisomy 4 status
final_seg <- read.csv(paste0(output_dir, "Tri4_FinalSeg.csv"))
plotdf$Status <- ifelse(plotdf$Tumor %in% final_seg$sample[final_seg$seg_cons %in% "4a" & final_seg$state %in% "amp"], "Trisomy 4", "No")
##Keep WM only
plotdf <- plotdf[plotdf$Tumor_Type %in% "WM",]
##Plot
a <- ggplot(plotdf) +
  geom_point(aes(UMAP_1, UMAP_2, color = Status), show.legend = T) + 
  scale_color_manual(values = c("grey", "tomato2"), name = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.text = element_text(size = 8))
a_leg <- get_legend(a)
b <- ggplot(plotdf) +
  geom_point(aes(UMAP_1, UMAP_2, color = Status), alpha = 0.5, size = 0.05, show.legend = F) + 
  scale_color_manual(values = c("grey", "tomato2"), name = "Status") +
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())

png(paste0(output_dir, "UMAP_Tri4_status.png"), res = 300, units="in", width = 7, height = 5)
plot_grid(b, a_leg, align = "h", axis = "tb", rel_widths = c(0.5, 0.1))
dev.off()
