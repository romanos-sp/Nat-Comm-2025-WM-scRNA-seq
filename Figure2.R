#Import dependencies
library(tidyverse)
library(Seurat)
library(ggpubr)
library(ggrepel)

#Set input/output directories
input_dir <- "data/"
output_dir <- "results/Figure2/"

#Set seed
set.seed(123)

#Figure 2B
##Plot UMAP of T & NK subtypes before and after stimulation
##Read integrated object as .rds object
integrated <- readRDS("results/IFN_experiment/integrated.rds")
##Subset T and NK cells
sub_integrated <- subset(integrated, subset = seurat_clusters %in% c("0", "1", "2", "4", "5" , "8", "10", "11", "15"))
##Create annotation vector, rename clusters and remove un-annotated clusters
annotation <- c('4' = 'CD4+ TN', '2' = 'IFN+ CD4+ TN', '1' = 'CD4+ TCM', '0' = 'IFN+ CD4+ TCM', '8' = 'CD8+ TCM/TEM', '5' = 'IFN+ CD8+ TCM/TEM', '10' = 'CD56dim NK', '11' = 'IFN+ CD56dim NK')
renamed <- RenameIdents(sub_integrated, annotation)
renamed_sub <- subset(renamed, idents=unname(annotation))
renamed_sub <- AddModuleScore(renamed_sub, features = list(c("ISG15", "ISG20", "IFI6", "IFI44L", "MX1", "STAT1", "MX2", "OAS1", "OASL", "IFIT1", "IFIT2", "IFIT3", "IFI27", "EIF2AK2")), name = c("IFN_sig"))
renamed_sub$Cell_Type <- Idents(renamed_sub)
saveRDS(renamed_sub, paste0(output_dir, "T_NK_Object.rds"))
##Plot UMAP
renamed_sub <- readRDS(paste0(output_dir, "T_NK_Object.rds"))
plotdf <- FetchData(renamed_sub, vars= c("Cell", "UMAP_1", "UMAP_2", "Demux_Sample_ID", "Status", "Condition", "TissueType", "Condition_Status", "Condition_Status_Tissue", "Cell_Type", "IFN_sig1"))
plotdf$Cell_ID <- rownames(plotdf)
centroids <- aggregate(plotdf[, c("UMAP_1","UMAP_2")], list(factor(plotdf$Cell_Type)), median)
colnames(centroids)[1] <- "Cluster"
plotdf$condition <- ifelse(plotdf$Condition == "M", "Media only", "IFN+")
plotdf$tissue <- ifelse(plotdf$TissueType == "BM", "BMMCs", "PBMCs")
plotdf$Lab <- paste0(plotdf$Status, " ", plotdf$tissue, ", ", plotdf$condition)
plotdf$Lab <- factor(plotdf$Lab, levels = c("AWM BMMCs, Media only", "AWM BMMCs, IFN+", "AWM PBMCs, Media only", "AWM PBMCs, IFN+", "HD BMMCs, Media only", "HD BMMCs, IFN+"))
write.csv(plotdf, paste0(output_dir, "T_NK_Data.csv"), row.names = T)
png(paste0(output_dir, "UMAP_Annotated_Facet.png"), res = 300, units="in",width = 5, height = 5)
ggplot(plotdf[plotdf$UMAP_1 < 3 & plotdf$UMAP_2 > -7.5,],aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=factor(Cell_Type)), show.legend = F) +
  geom_text(data = centroids, aes(x =  UMAP_1, y = UMAP_2, label = Cluster), size = 2.5) +
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), strip.background = element_blank(), strip.text = element_text(size=14, face = "italic")) +
  scale_color_manual(values=c("lemonchiffon1", "yellow", "orange", "tomato2", "lightblue", "steelblue", "lightgreen", "darkgreen")) +
  facet_wrap(.~Lab, ncol = 2)
dev.off()

#Figure 2C
##Plot UMAP of T & NK cells, colored by the activity of an IFN signature
plotdf <- read.csv(paste0(output_dir, "T_NK_Data.csv"), row.names = 1)
plotdf$Lab <- factor(plotdf$Lab, levels = c("AWM BMMCs, Media only", "AWM BMMCs, IFN+", "AWM PBMCs, Media only", "AWM PBMCs, IFN+", "HD BMMCs, Media only", "HD BMMCs, IFN+"))
png(paste0(output_dir, "UMAP_IFNsig_Facet.png"), res = 300, units="in",width = 6, height = 5)
ggplot(plotdf[plotdf$UMAP_1 < 3 & plotdf$UMAP_2 > -7.5,],aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=IFN_sig1), show.legend = T) +
  geom_text(data = centroids, aes(x =  UMAP_1, y = UMAP_2, label = Cluster), size = 2.5) +
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), strip.background = element_blank(), strip.text = element_text(size=14, face = "italic")) +
  scale_color_gradient(low = "lightblue", high = "tomato2", name = "IFN signature") +
  facet_wrap(.~Lab, ncol = 2)
dev.off()

#Figure 2D
##Plot marker heatmap
renamed_sub <- readRDS(paste0(output_dir, "T_NK_Object.rds"))
plotdf <- read.csv(paste0(output_dir, "T_NK_Data.csv"), row.names = 1)
markers <- c("CD3D", "CD4", "LEF1", "TCF7", "SELL", "CCR7", "IL7R", "HNRNPLL", "FAS", "SESN3", "GATA3", "CD8A", "CD8B", "LINC02446", "CXCR3", "GZMA", "GZMK", "LYAR", "CCL5", "GZMH", "GZMB", "PRF1", "FGFBP2", "GNLY", "NKG7","MX1", "ISG15", "IFI6", "IFI44L", "STAT1","IFIT1", "IFIT2", "IFIT3")
orig_cell_types <- c('CD4+ TN', 'IFN+ CD4+ TN', 'CD4+ TCM', 'IFN+ CD4+ TCM', 'CD8+ TN', 'CD8+ TCM/TEM', 'IFN+ CD8+ TCM/TEM', 'CD56dim NK', 'IFN+ CD56dim NK')
marker_counts <- data.frame(t(data.frame(GetAssayData(renamed_sub[markers,], slot = "data"))))
marker_counts <- data.frame(scale(marker_counts))
rownames(marker_counts) <- gsub("\\.", "-", rownames(marker_counts))
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
png(paste0(output_dir, "Annotation_Heatmap.png"), res = 300, units="in",width = 12, height = 5)
ggplot(marker_means_long) +
  geom_tile(aes(Marker, Cluster, fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 65, hjust=1, size=12, face="italic", color="black"), 
        axis.text.y=element_text(size=12, color="black")) + 
  xlab("") + 
  ylab("") + 
  labs(fill="Mean Exp")
dev.off()

#Figure 2E
##Compare IFN expression levels between diagnoses in median-only condition
plotdf <- read.csv(paste0(output_dir, "T_NK_Data.csv"), row.names = 1)
plotdf$Lab <- factor(plotdf$Lab, levels = c("AWM BMMCs, Media only", "AWM BMMCs, IFN+", "AWM PBMCs, Media only", "AWM PBMCs, IFN+", "HD BMMCs, Media only", "HD BMMCs, IFN+"))
mediaonly <- plotdf[plotdf$Condition %in% "M",]
mediaonly$NewCellType <- gsub("IFN\\+ ", "", mediaonly$Cell_Type)
mediaonly$NewCellType <- factor(mediaonly$NewCellType, levels = c("CD4+ TN", "CD4+ TCM", "CD8+ TCM/TEM", "CD56dim NK"))
stats <- compare_means(IFN_sig1 ~ Lab, group.by = c("NewCellType"), data = mediaonly, p.adjust.method = "BH")
stats$mid <-as.numeric(stats$NewCellType)
stats$xmin <- ifelse(stats$group1=="AWM BMMCs, Media only" & stats$group2=="AWM PBMCs, Media only", stats$mid-1/3.8, ifelse(stats$group1=="AWM PBMCs, Media only" & stats$group2=="HD BMMCs, Media only", stats$mid, ifelse(stats$group1=="AWM BMMCs, Media only" & stats$group2=="HD BMMCs, Media only", stats$mid-1/3.8, NA)))
stats$xmax <- ifelse(stats$group1=="AWM BMMCs, Media only" & stats$group2=="AWM PBMCs, Media only", stats$mid, ifelse(stats$group1=="AWM PBMCs, Media only" & stats$group2=="HD BMMCs, Media only", stats$mid + 1/3.8, ifelse(stats$group1=="AWM BMMCs, Media only" & stats$group2=="HD BMMCs, Media only", stats$mid+1/3.8, NA)))
ymaxes <- data.frame(aggregate(mediaonly$IFN_sig1, list(mediaonly$NewCellType), max, na.rm = T))
colnames(ymaxes) <- c("NewCellType", "y.max")
stats$y.max <- ymaxes$y.max[match(stats$NewCellType, ymaxes$NewCellType)]
stats$y.position <- ifelse(stats$group1=="AWM BMMCs, Media only" & stats$group2=="AWM PBMCs, Media only", stats$y.max + 0.2, ifelse(stats$group1=="AWM PBMCs, Media only" & stats$group2=="HD BMMCs, Media only", stats$y.max + 0.7, ifelse(stats$group1=="AWM BMMCs, Media only" & stats$group2=="HD BMMCs, Media only", stats$y.max + 1.2, NA)))
stats$p.format <- ifelse(stats$p.adj < 2e-16, "q<2e-16", paste0("q=", signif(stats$p.adj, 2)))
mediaonly$label <- paste0(mediaonly$NewCellType, ";", mediaonly$Lab)
table(mediaonly$label)
##CD4+ TCM;AWM BMMCs, Media only     CD4+ TCM;AWM PBMCs, Media only      CD4+ TCM;HD BMMCs, Media only 
##587                               3309                                696 
##CD4+ TN;AWM BMMCs, Media only      CD4+ TN;AWM PBMCs, Media only       CD4+ TN;HD BMMCs, Media only 
##608                               2836                                539 
##CD56dim NK;AWM BMMCs, Media only   CD56dim NK;AWM PBMCs, Media only    CD56dim NK;HD BMMCs, Media only 
##114                               1700                                128 
##CD8+ TCM/TEM;AWM BMMCs, Media only CD8+ TCM/TEM;AWM PBMCs, Media only  CD8+ TCM/TEM;HD BMMCs, Media only 
##468                               2085                                676 

df <- data.frame(table(mediaonly$label))
df$NewCellType <- gsub(";.*", "", df$Var1)
df$Lab <- gsub("^.*;", "", df$Var1)
df$NewCellType <- factor(df$NewCellType, levels = levels(mediaonly$NewCellType))
df$Lab <- factor(df$Lab, levels = levels(mediaonly$Lab))
df <- df[order(df$NewCellType, decreasing = F),]
df$x <- ifelse(df$Lab == levels(df$Lab)[1] & df$NewCellType == levels(df$NewCellType)[1], as.numeric(df$NewCellType)[1]*2/3, 
               ifelse(df$Lab == levels(df$Lab)[3] & df$NewCellType == levels(df$NewCellType)[1], as.numeric(df$NewCellType)[1]*3/3, 
                      ifelse(df$Lab == levels(df$Lab)[5] & df$NewCellType == levels(df$NewCellType)[1], as.numeric(df$NewCellType)[1]*4/3, 
                             ifelse(df$Lab == levels(df$Lab)[1] & df$NewCellType == levels(df$NewCellType)[2], 1+ as.numeric(df$NewCellType)[1]*2/3, 
                                    ifelse(df$Lab == levels(df$Lab)[3] & df$NewCellType == levels(df$NewCellType)[2], 1+as.numeric(df$NewCellType)[1]*3/3, 
                                           ifelse(df$Lab == levels(df$Lab)[5] & df$NewCellType == levels(df$NewCellType)[2], 1+as.numeric(df$NewCellType)[1]*4/3, 
                                                  ifelse(df$Lab == levels(df$Lab)[1] & df$NewCellType == levels(df$NewCellType)[3], 2+as.numeric(df$NewCellType)[1]*2/3, 
                                                         ifelse(df$Lab == levels(df$Lab)[3] & df$NewCellType == levels(df$NewCellType)[3], 2+as.numeric(df$NewCellType)[1]*3/3, 
                                                                ifelse(df$Lab == levels(df$Lab)[5] & df$NewCellType == levels(df$NewCellType)[3], 2+as.numeric(df$NewCellType)[1]*4/3, 
                                                                       ifelse(df$Lab == levels(df$Lab)[1] & df$NewCellType == levels(df$NewCellType)[4], 3+as.numeric(df$NewCellType)[1]*2/3, 
                                                                              ifelse(df$Lab == levels(df$Lab)[3] & df$NewCellType == levels(df$NewCellType)[4], 3+as.numeric(df$NewCellType)[1]*3/3, 
                                                                                     ifelse(df$Lab == levels(df$Lab)[5] & df$NewCellType == levels(df$NewCellType)[4], 3+as.numeric(df$NewCellType)[1]*4/3, NA
                                                                                     )
                                                                              )
                                                                       )
                                                                )
                                                         )
                                                  )
                                           )
                                    )
                             )
                      )))



png(paste0(output_dir, "IFN_Sig_MediaOnly.png"), res = 300, units="in", width = 6, height = 5)
ggplot(mediaonly) +
  geom_boxplot(aes(NewCellType, IFN_sig1, fill = Lab), alpha = 0.5, position = position_dodge(width = 0.8)) +
  geom_violin(aes(NewCellType, IFN_sig1, fill = Lab), alpha = 0.5, position = position_dodge(width = 0.8), show.legend = F) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.title = element_text(size = 14), 
        legend.position = "top", 
        plot.title = element_text(hjust = 0.5, size = 16), 
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("tomato2", "orange", "lightblue"), labels = c("AWM BMMCs", "AWM PBMCs", "HD BMMCs"), name = "") +
  ylab("IFN signature") +
  xlab("") +
  geom_bracket(aes(xmin = xmin, xmax = xmax, label = p.format, y.position = y.position), data = stats[stats$p.adj < 0.05,], size = 0.4, label.size = 4) +
  ylim(min(mediaonly$IFN_sig1)-0.2, 4.5) +
  ggtitle("IFN stimulation (-)") +
  geom_text(data = df, aes(x = x, y = min(mediaonly$IFN_sig1)-0.15, label = paste0("n=", Freq)), size = 3)
dev.off()

#Figure 2F
##Compare IFN expression levels between diagnoses after IFN stimulation
plotdf <- read.csv(paste0(output_dir, "T_NK_Data.csv"), row.names = 1)
plotdf$Lab <- factor(plotdf$Lab, levels = c("AWM BMMCs, Media only", "AWM BMMCs, IFN+", "AWM PBMCs, Media only", "AWM PBMCs, IFN+", "HD BMMCs, Media only", "HD BMMCs, IFN+"))
ifnonly <- plotdf[plotdf$Condition %in% "IFN",]
ifnonly$NewCellType <- gsub("IFN\\+ ", "", ifnonly$Cell_Type)
ifnonly$NewCellType <- factor(ifnonly$NewCellType, levels = c("CD4+ TN", "CD4+ TCM", "CD8+ TCM/TEM", "CD56dim NK"))
stats <- compare_means(IFN_sig1 ~ Lab, group.by = c("NewCellType"), data = ifnonly, p.adjust.method = "BH")
stats$mid <-as.numeric(stats$NewCellType)
stats$xmin <- ifelse(stats$group1=="AWM BMMCs, IFN+" & stats$group2=="AWM PBMCs, IFN+", stats$mid-1/3.8, ifelse(stats$group1=="AWM PBMCs, IFN+" & stats$group2=="HD BMMCs, IFN+", stats$mid, ifelse(stats$group1=="AWM BMMCs, IFN+" & stats$group2=="HD BMMCs, IFN+", stats$mid-1/3.8, NA)))
stats$xmax <- ifelse(stats$group1=="AWM BMMCs, IFN+" & stats$group2=="AWM PBMCs, IFN+", stats$mid, ifelse(stats$group1=="AWM PBMCs, IFN+" & stats$group2=="HD BMMCs, IFN+", stats$mid + 1/3.8, ifelse(stats$group1=="AWM BMMCs, IFN+" & stats$group2=="HD BMMCs, IFN+", stats$mid+1/3.8, NA)))
ymaxes <- data.frame(aggregate(ifnonly$IFN_sig1, list(ifnonly$NewCellType), max, na.rm = T))
colnames(ymaxes) <- c("NewCellType", "y.max")
stats$y.max <- ymaxes$y.max[match(stats$NewCellType, ymaxes$NewCellType)]
stats$y.position <- ifelse(stats$group1=="AWM BMMCs, IFN+" & stats$group2=="AWM PBMCs, IFN+", stats$y.max + 0.2, ifelse(stats$group1=="AWM PBMCs, IFN+" & stats$group2=="HD BMMCs, IFN+", stats$y.max + 0.7, ifelse(stats$group1=="AWM BMMCs, IFN+" & stats$group2=="HD BMMCs, IFN+", stats$y.max + 1.2, NA)))
stats$p.format <- ifelse(stats$p.adj < 2e-16, "q<2e-16", paste0("q=", signif(stats$p.adj, 2)))
ifnonly$label <- paste0(ifnonly$NewCellType, ";", ifnonly$Lab)
table(ifnonly$label)
##CD4+ TCM;AWM BMMCs, IFN+     CD4+ TCM;AWM PBMCs, IFN+      CD4+ TCM;HD BMMCs, IFN+      CD4+ TN;AWM BMMCs, IFN+ 
##815                         2990                          814                          980 
##CD4+ TN;AWM PBMCs, IFN+       CD4+ TN;HD BMMCs, IFN+   CD56dim NK;AWM BMMCs, IFN+   CD56dim NK;AWM PBMCs, IFN+ 
##  2290                          671                          155                         1408 
##CD56dim NK;HD BMMCs, IFN+ CD8+ TCM/TEM;AWM BMMCs, IFN+ CD8+ TCM/TEM;AWM PBMCs, IFN+  CD8+ TCM/TEM;HD BMMCs, IFN+ 
##  153                          793                         1622                          802

df <- data.frame(table(ifnonly$label))
df$NewCellType <- gsub(";.*", "", df$Var1)
df$Lab <- gsub("^.*;", "", df$Var1)
df$NewCellType <- factor(df$NewCellType, levels = levels(ifnonly$NewCellType))
df$Lab <- factor(df$Lab, levels = levels(ifnonly$Lab))
df <- df[order(df$NewCellType, decreasing = F),]
df$x <- ifelse(df$Lab == levels(df$Lab)[2] & df$NewCellType == levels(df$NewCellType)[1], as.numeric(df$NewCellType)[1]*2/3, 
               ifelse(df$Lab == levels(df$Lab)[4] & df$NewCellType == levels(df$NewCellType)[1], as.numeric(df$NewCellType)[1]*3/3, 
                      ifelse(df$Lab == levels(df$Lab)[6] & df$NewCellType == levels(df$NewCellType)[1], as.numeric(df$NewCellType)[1]*4/3, 
                             ifelse(df$Lab == levels(df$Lab)[2] & df$NewCellType == levels(df$NewCellType)[2], 1+ as.numeric(df$NewCellType)[1]*2/3, 
                                    ifelse(df$Lab == levels(df$Lab)[4] & df$NewCellType == levels(df$NewCellType)[2], 1+as.numeric(df$NewCellType)[1]*3/3, 
                                           ifelse(df$Lab == levels(df$Lab)[6] & df$NewCellType == levels(df$NewCellType)[2], 1+as.numeric(df$NewCellType)[1]*4/3, 
                                                  ifelse(df$Lab == levels(df$Lab)[2] & df$NewCellType == levels(df$NewCellType)[3], 2+as.numeric(df$NewCellType)[1]*2/3, 
                                                         ifelse(df$Lab == levels(df$Lab)[4] & df$NewCellType == levels(df$NewCellType)[3], 2+as.numeric(df$NewCellType)[1]*3/3, 
                                                                ifelse(df$Lab == levels(df$Lab)[6] & df$NewCellType == levels(df$NewCellType)[3], 2+as.numeric(df$NewCellType)[1]*4/3, 
                                                                       ifelse(df$Lab == levels(df$Lab)[2] & df$NewCellType == levels(df$NewCellType)[4], 3+as.numeric(df$NewCellType)[1]*2/3, 
                                                                              ifelse(df$Lab == levels(df$Lab)[4] & df$NewCellType == levels(df$NewCellType)[4], 3+as.numeric(df$NewCellType)[1]*3/3, 
                                                                                     ifelse(df$Lab == levels(df$Lab)[6] & df$NewCellType == levels(df$NewCellType)[4], 3+as.numeric(df$NewCellType)[1]*4/3, NA
                                                                                     )
                                                                              )
                                                                       )
                                                                )
                                                         )
                                                  )
                                           )
                                    )
                             )
                      )))

png(paste0(output_dir, "IFN_Sig_IFN+.png"), res = 300, units="in", width = 6, height = 5)
ggplot(ifnonly) +
  geom_boxplot(aes(NewCellType, IFN_sig1, fill = Lab), alpha = 0.5, position = position_dodge(width = 0.8)) +
  geom_violin(aes(NewCellType, IFN_sig1, fill = Lab), alpha = 0.5, position = position_dodge(width = 0.8), show.legend = F) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_text(size = 12, color = "black"), axis.text.x = element_text(angle = 65, hjust = 1), axis.title = element_text(size = 14), legend.position = "top", plot.title = element_text(hjust = 0.5, size = 16), legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("tomato2", "orange", "lightblue"), labels = c("AWM BMMCs", "AWM PBMCs", "HD BMMCs"), name = "") +
  ylab("IFN signature") +
  xlab("") +
  geom_bracket(aes(xmin = xmin, xmax = xmax, label = p.format, y.position = y.position), data = stats[stats$p.adj < 0.05,], size = 0.4, label.size = 4) +
  ylim(min(ifnonly$IFN_sig1)-0.2, 4.5) +
  ggtitle("IFN stimulation (+)") +
  geom_text(data = df, aes(x = x, y = min(ifnonly$IFN_sig1)-0.15, label = paste0("n=", Freq)), size = 3)
dev.off()

#Figure 2G
##Visualize DE results in monocytes from HD
res <- read.csv("results/IFN_experiment/IFNvsHD_CD14+ Monocytes_DE_Wilcoxon.csv")
up <- res[res$Log2FC > 0 & res$padj < 0.05,]
topup <- up$Gene[order(up$Log2FC, decreasing = T)][1:40]
down <- res[res$Log2FC < 0 & res$padj < 0.05,]
topdown <- down$Gene[order(down$Log2FC, decreasing = F)][1:40]
write.csv(topup, paste0(output_dir, "IFNstim_up.csv"), row.names = F)
write.csv(topdown, paste0(output_dir, "IFNstim_down.csv"), row.names = F)
res$lab <- ifelse(res$Gene %in% c(topup, topdown), as.character(res$Gene), NA)
res$flag <- factor(ifelse(res$Log2FC > 0 & res$padj >= 0.05, 0, ifelse(res$Log2FC > 0 & res$padj < 0.05, 1, ifelse(res$Log2FC < 0 & res$padj >= 0.05, 2, ifelse(res$Log2FC < 0 & res$padj < 0.05, 3, NA)))), levels = c("0", "1", "2", "3"))
res$lab2 <- ifelse(res$Gene %in% c("MNDA"), "bolditalic('MNDA')", ifelse(res$Gene %in% c("IL1B"), "bolditalic('IL1B')", ifelse(res$Gene %in% c("CXCL8"), "bolditalic('CXCL8')", paste0("italic('", as.character(res$lab), "')"))))
res$lab2[is.na(res$lab)] <- NA
png(paste0(output_dir, "DE_HD_CD14+ Monocytes.png"), res = 300, units="in", width = 6, height = 5)
ggplot(res) +
  geom_point(aes(Log2FC, -log10(padj), color = flag), show.legend = F) +
  geom_text_repel(aes(Log2FC, -log10(padj), label = lab2), parse = T, segment.alpha = 0.2, segment.curvature = 0.5, max.overlaps = Inf, size = 3, force = 5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5, color = "black") +
  scale_color_manual(values = c("pink", "orchid4", "lemonchiffon1", "orange")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA,  color = "black"), 
        axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 16)) +
  ylab(expression("-Log"[10]*"(q-value)")) +
  xlab(expression("Log"[2]*" fold-change")) +
  scale_x_continuous(limits = c(-1.5, 3), breaks = c(-1, 0, 1)) +
  ggtitle("HD BM CD14+ Monocytes: IFN(-) VS IFN (+)")
dev.off()

#Figure 2H
##Read in IFN stimulation signatures and compare activity in patients
topup <- read.csv(paste0(output_dir, "IFNstim_up.csv"))
topdown <- read.csv(paste0(output_dir, "IFNstim_down.csv"))
integrated <- readRDS("results/Figure1/Myeloid_Object.rds")
integrated <- AddModuleScore(integrated, features = list(topup$x), name = c("IFNstim_up"))
integrated <- AddModuleScore(integrated, features = list(topdown$x), name = c("IFNstim_down"))
integrated <- AddModuleScore(integrated, features = list(c("ISG15", "ISG20", "IFI6", "IFI44L", "MX1", "STAT1", "MX2", "OAS1", "OASL", "IFIT1", "IFIT2", "IFIT3", "IFI27", "EIF2AK2")), name = c("IFN_sig"))
plotdf <- FetchData(integrated, vars= c("Cell", "UMAP_1", "UMAP_2", "Demux_Sample_ID", "Cell_Type", "IFN_sig1", "Disease", "IFNstim_up1", "IFNstim_down1"))
plotdf$Disease <- as.character(plotdf$Disease)
plotdf$Disease[plotdf$Disease %in% "NBM"] <- "HD"
plotdf$Disease <- factor(plotdf$Disease, levels = c("HD", "SMM", "AWM", "WM"))
write.csv(plotdf, paste0(output_dir, "Myeloid_IFN_Data.csv"), row.names = T)
png(paste0(output_dir, "IFN_stim_PatientMonocytes.png"), res = 300, units="in",width = 6, height = 5)
ggplot(plotdf[plotdf$Cell_Type %in% c("IL1B+ CD14+ Mono", "IFN+ CD14+ Mono"),]) +
  geom_point(aes(IFNstim_up1, IFNstim_down1, color = Cell_Type), alpha = 1) +
  stat_density2d(aes(IFNstim_up1, IFNstim_down1, alpha=..level../10), geom="polygon", fill = "lemonchiffon1", show.legend = F) +
  scale_color_manual(values = c("steelblue", "tomato2"), name = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 14, face = "italic"), 
        legend.position = "top", 
        legend.text = element_text(size = 14)) +
  xlab(paste0("Genes ", sprintf('\u2191') ," with IFN stimulation")) +
  ylab(paste0("Genes ", sprintf('\u2193') ," with IFN stimulation"))+
  facet_wrap(.~Disease)
dev.off()

#Figure 2I
##Compare levels of the IFN signature in IFN+ CD14+ Monocytes
plotdf <- read.csv(paste0(output_dir, "SourceData_Fig2I.csv"))
plotdf$Disease <- factor(plotdf$Disease, levels = c("HD", "SMM", "AWM", "WM"))
stats <- compare_means(IFN_sig1 ~ Disease, data = plotdf[plotdf$Cell_Type %in% "IFN+ CD14+ Mono",], p.adjust.method = "BH")
stats$p.format <- ifelse(stats$p.adj < 2e-16, "q<2e-16", paste0("q=", signif(stats$p.adj, 2)))
df <- data.frame(table(plotdf$Disease))
df$x <- 1:nrow(df)
png(paste0(output_dir, "IFN_sig_IFN+ CD14+ Mono.png"), res = 300, units="in", width = 5, height = 5)
ggplot(plotdf[plotdf$Cell_Type %in% "IFN+ CD14+ Mono",]) +
  geom_boxplot(aes(Disease, IFN_sig1, fill = Disease), show.legend = F, alpha = 0.5) +
  geom_violin(aes(Disease, IFN_sig1, fill = Disease), show.legend = F, alpha = 0.5) +
  scale_fill_manual(values = c("steelblue", "tomato2", "orange", "orange3")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(size = 14), 
        axis.title = element_text(size = 14)) +
  xlab("") +
  ylab("IFN sig activity in IFN+ CD14+ Mono") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = p.format), data = stats[stats$p.adj < 0.05,], y.position = 2, step.increase = 0.07, size = 0.4, label.size = 5, vjust = 1.5) +
  scale_y_continuous(limits = c(min(plotdf$IFN_sig1)-0.2, 2.7), breaks = c(0, 1, 2)) +
  geom_text(data = df, aes(x = x, y = min(plotdf$IFN_sig1)-0.15, label = paste0("n=", Freq)), size = 5)
dev.off()

#Figure 2J
##Read in Olink data and plot interferon gamma levels
olink <- read.csv(paste0(input_dir, "Olink_Data.csv"))
ifn <- olink[olink$Cytokine %in% "IFNG",]
ifn$Diagnosis <- ifelse(as.character(ifn$Diagnosis) %in% "MGUS", "Non-IgM MGUS", ifelse(as.character(ifn$Diagnosis) %in% "Healthy", "HD", ifelse(as.character(ifn$Diagnosis) %in% "IgM-MGUS", "IgM MGUS", as.character(ifn$Diagnosis))))
ifn$Diagnosis <- factor(ifn$Diagnosis, levels = c("HD", "IgM MGUS", "Non-IgM MGUS", "SMM"))
table(ifn$Diagnosis)
##HD     IgM MGUS Non-IgM MGUS          SMM 
##13            5           15           20
stats <- compare_means(Level ~ Diagnosis, data = ifn, ref.group = "IgM MGUS", p.adjust.method = "BH")
write.csv(stats, paste0(output_dir, "Olink_IFNg_Stats.csv"), row.names = F)
png(paste0(output_dir, "Olink_IFNg.png"), res=300, unit="in", width=6, height=5)
ggplot(ifn) +
  geom_boxplot(aes(Diagnosis, Level, fill = Diagnosis), alpha = 0.5, outlier.size = -1, show.legend = F) +
  geom_violin(aes(Diagnosis, Level, fill = Diagnosis), alpha = 0.5, show.legend = F) +
  geom_point(aes(Diagnosis, Level, fill = Diagnosis), position = position_jitter(width = 0.2), shape = 21, size = 4, show.legend = F) +
  scale_fill_manual(values = c("steelblue", "lemonchiffon1", "lightcoral", "orangered3")) +
  scale_y_log10(limits = c(0.05, 3)) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 14)) +
  ylab("PB plasma Interferon-gamma (pg/mL)") +
  xlab("") +
  geom_bracket(aes(xmin = group1, xmax = group2, label = paste0("p=", signif(p, 2))), y.position = 0.1, data = stats, step.increase = 0.05, label.size = 5, vjust = 1.5)
dev.off()

#Figure 2K
##Import NK cell functionality experiment data
##K562 cancer cells were co-cultured with HD vs AWM PB-derived NK cells and the frequency of apoptotic cells was measured via Flow Cytometry
data <- read.csv(paste0(input_dir, "NK_Functionality_Experiment/NK_Functionality_Data.csv"))
##Compare NK cell functionality before and after stimulation
data$Condition <- ifelse(data$SampleID %in% c("HD1", "HD2", "HD3", "HD4"), "HD",
                         ifelse(data$SampleID %in% c("AWM1", "AWM2", "AWM3", "AWM4", "AWM5"), "AWM", "K562"))
data <- data[, c("SampleID", "Condition", "Baseline", "IFN")]
data_long <- gather(data, "Status", "Rate", -SampleID, -Condition)
data_long$Status <- as.character(data_long$Status)
data_long$Status <- factor(ifelse(data_long$Status %in% "Baseline", "w/o IFN", "w/ IFN"), levels = c("w/o IFN", "w/ IFN"))
data_long$Condition <- factor(data_long$Condition, levels = c("K562", "HD", "AWM"))
stats <- data.frame("group1" = c("w/o IFN", "w/o IFN"), 
                    "group2" = c("w/ IFN", "w/ IFN"), 
                    "Condition" = c("HD", "AWM"),
                    "p" = c(t.test(data$IFN[data$Condition == "HD"], data$Baseline[data$Condition == "HD"], paired = T)$p.value, t.test(data$IFN[data$Condition == "AWM"], data$Baseline[data$Condition == "AWM"], paired = T)$p.value))
agg <- aggregate(data_long$Rate, list(data_long$Condition), max)
stats$y.position <- agg$x[match(stats$Condition, agg$Group.1)]
stats$Condition <- factor(stats$Condition, levels = c("HD", "AWM"))
png(paste0(output_dir, "NKFunctionality_PrevsPost.png"), res = 300, units="in", width = 5, height = 5)
ggplot(data_long[data_long$Condition %in% c("HD", "AWM"),]) + 
  geom_boxplot(aes(Status, Rate, fill = Status), outlier.size = -1, alpha = 0.5, show.legend = F) + 
  geom_point(aes(Status, Rate, fill = Status), shape = 21, size = 4, show.legend = F, alpha = 0.8) +
  geom_line(aes(Status, Rate, group = SampleID), linetype = "dashed") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "italic")) +
  scale_fill_manual(values = c("lemonchiffon1", "darkblue")) +
  ylab("(%) Apoptotic K562") +
  xlab("") +
  geom_bracket(data = stats, aes(xmin = group1, xmax = group2, y.position = y.position + 4, label = paste0("p=", signif(p, 2))), label.size = 4, vjust = 1.8) +
  facet_wrap(.~Condition)
dev.off()

#Figure 2L
##Import NK cell functionality experiment data
##K562 cancer cells were co-cultured with HD vs AWM PB-derived NK cells and the frequency of apoptotic cells was measured via Flow Cytometry
data <- read.csv(paste0(input_dir, "NK_Functionality_Experiment/NK_Functionality_Data.csv"))
##Compare K562 apoptosis rate between patients with AWM and HD at baseline and post-stimulation with IFN-I
data$Condition <- ifelse(data$SampleID %in% c("HD1", "HD2", "HD3", "HD4"), "HD",
                         ifelse(data$SampleID %in% c("AWM1", "AWM2", "AWM3", "AWM4", "AWM5"), "AWM", "K562"))
data <- data[, c("SampleID", "Condition", "Baseline", "IFN")]
data_long <- gather(data, "Status", "Rate", -SampleID, -Condition)
data_long$Status <- as.character(data_long$Status)
data_long$Status <- factor(ifelse(data_long$Status %in% "Baseline", "w/o IFN", "w/ IFN"), levels = c("w/o IFN", "w/ IFN"))
data_long$Condition <- factor(data_long$Condition, levels = c("K562", "HD", "AWM"))
stats <- compare_means(Rate ~ Condition, group.by = "Status", data = data_long[data_long$Condition %in% c("HD", "AWM"),], p.adjust.method = "BH")
agg <- aggregate(data_long$Rate, list(data_long$Status), max)
stats$y.position <- agg$x[match(stats$Status, agg$Group.1)]

png(paste0(output_dir, "NKFunctionality_HDvsAWM.png"), res = 300, units="in", width = 5, height = 5)
ggplot(data_long) + 
  geom_boxplot(aes(Condition, Rate, fill = Condition), outlier.size = -1, alpha = 0.5, show.legend = F) + 
  geom_point(aes(Condition, Rate, fill = Condition), shape = 21, position = position_jitter(), size = 4, show.legend = F, alpha = 0.8) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "italic")) +
  scale_fill_manual(values = c("lightblue", "darkblue", "orange")) +
  ylab("(%) Apoptotic K562") +
  xlab("") +
  geom_bracket(data = stats, aes(xmin = group1, xmax = group2, y.position = y.position + 6, label = paste0("p=", p.format)), label.size = 4, vjust = 1.8) +
  facet_wrap(.~Status)
dev.off()
