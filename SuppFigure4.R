#Import dependencies
library(tidyverse)
library(ggpubr)

#Set input/output directories
input_dir <- "data/"
output_dir <- "results/SuppFigure4/"

#Supp. Figure 4A
##Plot GEX-7 activity by CXCR4 mutation status
h_clust <- read.csv("results/Figure4/Merged_Data.csv", row.names = 1)
nsigs <- 8
h_clust_long <- h_clust[, c("UMAP_1", "UMAP_2", "NewCellID", "CaTissueID", "Annotation", "Tumor_Type2", c(paste0("S", 1:nsigs)))]
h_clust_long <- gather(h_clust_long, "Signature", "Expression", -UMAP_1, -UMAP_2, -NewCellID, -CaTissueID, -Annotation, -Tumor_Type2)
h_clust_long$Signature <- factor(gsub("S", "GEX-", h_clust_long$Signature), levels = paste0("GEX-", 1:nsigs))
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
h_clust_long$Deidentified_SampleID <- deid$PatientDeIdentified[match(h_clust_long$CaTissueID, deid$CaTissueID)]
metadata <- read.csv(paste0(input_dir, "wm-cohort-metadata.csv"))
h_clust_long$CXCR4_mut <- metadata$CXCR4[match(h_clust_long$CaTissueID, metadata$CaTissueID)]
h_clust_long$Tumor <- ifelse(h_clust_long$Annotation == "Tumor", paste0(h_clust_long$Tumor_Type2, "_", h_clust_long$Deidentified_SampleID), as.character(h_clust_long$Annotation))
##Only assign CXCR4 mutant status to primary WM tumors
h_clust_long$CXCR4_mut[h_clust_long$Tumor_Type2 %in% c("WM_2", "WM_3")] <- NA
##Remove tumor with NA CXCR4 status
h_clust_long <- h_clust_long[!is.na(h_clust_long$CXCR4_mut),]
##Keep only GEX-7
gex7 <- h_clust_long[h_clust_long$Signature %in% "GEX-7",]
cxcr4_agg <- aggregate(gex7$Expression, list(gex7$Tumor), median)
h_clust_long$Tumor <- factor(h_clust_long$Tumor, levels = cxcr4_agg$Group.1[order(cxcr4_agg$x, decreasing = T)])
write.csv(h_clust_long, paste0(output_dir, "NMF_Tumor_Poisson_HL1_WL1_n30_GEX7_ByCXCR4Status_Data.csv"), row.names = F)
plotdf <- h_clust_long[h_clust_long$Signature %in% c("GEX-7"),]
df <- data.frame(table(plotdf$Tumor))
df$Var1 <- factor(df$Var1, levels = levels(plotdf$Tumor))
df <- df[order(df$Var1),]
df$CXCR4_mut <- plotdf$CXCR4_mut[match(df$Var1, plotdf$Tumor)]
df$CXCR4_mut <- factor(df$CXCR4_mut, levels = c("Mut", "WT"))
df <- df[order(df$CXCR4_mut),]
sums <- data.frame(table(df$CXCR4_mut))
df$x <- c(1:sums$Freq[sums$Var1 == levels(df$CXCR4_mut)[1]], 1:sums$Freq[sums$Var1 == levels(df$CXCR4_mut)[2]])

png(paste0(output_dir, "NMF_Tumor_Poisson_HL1_WL1_n30_GEX7_ByCXCR4Status.png"), res = 300, units="in", width = 8, height = 5)
ggplot(plotdf) +
  geom_boxplot(aes(Tumor, log10(Expression + 1), fill = CXCR4_mut), show.legend = F) +
  scale_fill_manual(values = c("orange", "steelblue", "lightgrey")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.y = element_text(color = "black", size =12), 
        axis.text.x = element_text(angle = 65, hjust = 1, color = "black", size = 14), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 14, face = "italic"), 
        axis.title = element_text(size = 14)) +
  xlab("") +
  ylab(expression("Log"[10]*"(GEX-7 + 1)")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(min(plotdf$Expression)-0.8, 3)) + 
  geom_text(data = df, aes(x = x, y = min(plotdf$Expression)-0.15, label = paste0("n=", Freq)), size = 3, angle = 90, hjust = 1) +
  facet_wrap(.~CXCR4_mut, scales = "free_x")
dev.off()

#Supp. Figure 4B
##Plot evolution of signature GEX-7 in patient with serial samples (Pt5)
deid <- read.csv(paste0(input_dir, "SampleTable_DeIdentified.csv"))
H_df <- read.csv("results/Figure4/H.csv", row.names = 1)
wm <- readRDS("results/Subclustering_B/WMonly/integrated.rds")
meta <- wm@meta.data
H_df$CaTissueID <- meta$CaTissueID[match(rownames(H_df), rownames(meta))]
sum(H_df$CaTissueID %in% deid$CaTissueID)/nrow(H_df)
##1
H_df$Deidentified_SampleID <- deid$SampleDeIdentified[match(H_df$CaTissueID, deid$CaTissueID)]
serial <- H_df[H_df$Deidentified_SampleID %in% c("S5", "S12"),]
serial$new_max_id <- apply(serial[, c("S2", "S7")], 1, function(x) {colnames(serial[, c("S2", "S7")])[which.max(x)]})
prop_s7 <- data.frame(prop.table(table(serial$new_max_id, serial$Deidentified_SampleID), margin = 2))
colnames(prop_s7) <- c("Signature", "Sample", "Prop")
prop_s7$Signature <- gsub("S", "GEX-", prop_s7$Signature)
prop_s7 <- prop_s7[prop_s7$Signature %in% "GEX-7",]
count_s7 <- data.frame(table(serial$new_max_id, serial$Deidentified_SampleID))
colnames(count_s7) <- c("Signature", "Sample", "Count")
prop_s7$low <- ifelse(prop_s7$Sample == "S5", prop.test(x = count_s7$Count[count_s7$Signature == "S7" & count_s7$Sample == "S5"], n = sum(count_s7$Count[count_s7$Sample == "S5"]))$conf.int[1], prop.test(x = count_s7$Count[count_s7$Signature == "S7" & count_s7$Sample == "S12"], n = sum(count_s7$Count[count_s7$Sample == "S12"]))$conf.int[1])
prop_s7$high <- ifelse(prop_s7$Sample == "S5", prop.test(x = count_s7$Count[count_s7$Signature == "S7" & count_s7$Sample == "S5"], n = sum(count_s7$Count[count_s7$Sample == "S5"]))$conf.int[2], prop.test(x = count_s7$Count[count_s7$Signature == "S7" & count_s7$Sample == "S12"], n = sum(count_s7$Count[count_s7$Sample == "S12"]))$conf.int[2])
write.csv(prop_s7, paste0(output_dir, "Serial_GEX7_Data.csv"), row.names = F)
stats <- data.frame("group1" = c(1), "group2" = c(2), "p" = fisher.test(table(serial$new_max_id, serial$Deidentified_SampleID))$p.value, "ypos" = max(prop_s7$high))

png(paste0(output_dir, "Serial_GEX7.png"), res = 300, units="in", width = 5, height = 5)
ggplot(prop_s7) +
  geom_bar(aes(Sample, Prop*100, fill = Sample), stat = "identity", width = 0.5, show.legend = F) +
  geom_errorbar(aes(x = Sample, ymin = low*100, ymax = high*100), width = 0.1) +
  geom_bracket(data = stats, aes(xmin = group1, xmax = group2, y.position = ypos*100 + 0.02*100, label = paste0("p=", signif(p, 2))), label.size = 5, vjust = 1.8) +
  scale_fill_manual(values = c("orange", "tomato2")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14)) +
  ylab("(%) cells assigned to GEX-7") + 
  xlab("Serial samples") +
  scale_x_discrete(labels = c("T1", "T2"))
dev.off()
