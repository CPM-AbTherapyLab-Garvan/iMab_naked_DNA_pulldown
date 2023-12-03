#set the directory path into where the files are: 
setwd("/Users/cris/Desktop/chip")

# Load the necessary libraries
library(dplyr)
library(readr)
library(ggplot2)
library(org.Hs.eg.db)

# read the peaks values
MCF7_1_MACS <- read.table("MCF7_1_MACS_peaks.narrowPeak", header=FALSE, sep="\t")

#rename columns
colnames(MCF7_1_MACS) <- c("Chr", "Start", "End", "name", "UCSC_score", "strand", "fold_change", "minuslog10pvalue", "minuslog10qvalue", "relative_summit")

#be sure to have an annotation of the peak files using Homer 
# or run in command the following:$ annotatePeaks.pl MCF7_1_MACS_peaks.narrowPeak hg19 > MCF7_1_MACS_peaks_annotated.narrowPeak
#read the annotated file
MCF7_1_MACS_annotated <- read_tsv("MCF7_1_MACS_peaks_annotated.narrowPeak")

#change the peak name to $name and $"Entrez ID" column names into a single word
colnames(MCF7_1_MACS_annotated)[colnames(MCF7_1_MACS_annotated) == 'PeakID (cmd=annotatePeaks.pl MCF7_1_MACS_peaks.narrowPeak hg19)'] <- 'name'
colnames(MCF7_1_MACS_annotated)[colnames(MCF7_1_MACS_annotated) == "Entrez ID"] <- "ENTREZID"

#merge the annotated file and the peak values by name
MCF7_1_merged_annotated <- merge(MCF7_1_MACS, MCF7_1_MACS_annotated, by = "name", all = TRUE)

# Read RNA file the file into R
RNA <- read.delim("RNA_WT/MCF7_WT_counts.txt", header = TRUE, sep = "\t", fill = TRUE, skip = 1)

#rename columns
colnames(RNA) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", "Rep1", "Rep2", "Rep3")

# Get the ENTREZIDs
RNA$ENTREZID <- mapIds(org.Hs.eg.db, 
                       keys=RNA$Geneid, 
                       column="ENTREZID", 
                       keytype="SYMBOL", 
                       multiVals="first")

#Add AVG and Log2 columns to the RNA dataframe
RNA$AVG <- rowMeans(RNA[,c("Rep1", "Rep2", "Rep3")], na.rm = TRUE)

# Adding 1 to avoid log(0) scenarios
RNA$Log2 <- log2(RNA$AVG + 1)  

# Convert log2 to numeric
RNA$log2 <- as.numeric(as.character(RNA$Log2))

# graphical distribution of RNA-seq data (optional)
hist(RNA$log2, main = "Histogram of AVG", xlab = "log2", col = "lightblue", border = "black")

#Q-Q plot
qqnorm(RNA$AVG, main = "Q-Q Plot of AVG")
qqline(RNA$AVG, col = "red")

#Shapiro-Wilk Test on the first 5000 RNA-seq reads
set.seed(123) 

# Setting seed to make your results reproducible
sample_rows <- sample(nrow(RNA), 5000)
subset_df <- RNA[sample_rows, ]

# Run the Shapiro-Wilk test on the AVG column
shapiro.test(subset_df$AVG)
#Data is not normally distributed in this case
#	Shapiro-Wilk normality test
# data:  subset_df$AVG
# W = 0.16404, p-value < 2.2e-16

# Merge both RNA data with the peak file 
final_MCF7_1_peaks_RNA_merged <- merge(RNA, MCF7_1_merged_annotated, by = "ENTREZID", all = TRUE)

# check the unique annotations in the $Annotation column
unique(final_MCF7_1_peaks_RNA_merged$Annotation)

#deleting the de description between parenthesis from the $Annotation column 
final_MCF7_1_peaks_RNA_merged$Annotation <- gsub("\\(.*\\)", "", final_MCF7_1_peaks_RNA_merged$Annotation)

#Also deleting any (dot + number) after the annotationexample "TTS .7" "5' UTR .6"
final_MCF7_1_peaks_RNA_merged$Annotation <- gsub("\\..*", "", final_MCF7_1_peaks_RNA_merged$Annotation)

#Check $Annotation identifiers
unique(final_MCF7_1_peaks_RNA_merged$Annotation)

#omit the NAs in the $Annotation column
final_MCF7_1_peaks_RNA_merged <- final_MCF7_1_peaks_RNA_merged %>%
  filter(!is.na(Annotation))

# Calculate medians
medians <- final_MCF7_1_peaks_RNA_merged %>% 
  group_by(Annotation) %>%
  summarise(median = median(Log2, na.rm = TRUE))

# Add median to original data
final_MCF7_1_peaks_RNA_merged <- final_MCF7_1_peaks_RNA_merged %>%
  left_join(medians, by = "Annotation")

# Convert Annotation to a factor, ordered by median
final_MCF7_1_peaks_RNA_merged$Annotation <- factor(final_MCF7_1_peaks_RNA_merged$Annotation, 
                                                   levels = unique(final_MCF7_1_peaks_RNA_merged$Annotation[order(final_MCF7_1_peaks_RNA_merged$median)]))

# Create plot  geom_boxplot defines the range and whiskers, stat_boxplot adds error bars (horizontal lines) at the ends of the whiskers
Fig3_a <- ggplot(final_MCF7_1_peaks_RNA_merged, aes(x = Annotation, y = Log2, fill = median)) +
  geom_boxplot(lwd=1.5, outlier.shape = NA) +  
  stat_boxplot(geom ='errorbar', lwd=1.5, width=0.2) +  
  scale_fill_gradient(low = "beige", high = "darkseagreen4", limits = c(1.5, 8)) +
  theme_minimal() +
  labs(x = "Annotation", y = "Log2(average mRNA expression + 1)") +
  ggtitle("Boxplot of Log2 for each Annotation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0, 20)

#saveplot
ggsave("Fig3_a.pdf", plot = Fig3_a, width = 7, height = 5)

# Perform pairwise Wilcoxon tests adjusting p-values for multiple testing using Benjamini-Hochberg method
pairwise.wilcox.test(final_MCF7_1_peaks_RNA_merged$Log2, final_MCF7_1_peaks_RNA_merged$Annotation, 
                     p.adjust.method = "BH")  

# Perform Kruskal-Wallis test
kruskal.test(Log2 ~ Annotation, data = final_MCF7_1_peaks_RNA_merged)

# cummulative distribution ECDF plot
Fig3_b <- ggplot(final_MCF7_1_peaks_RNA_merged, aes(x = Log2, color = Annotation)) +
  stat_ecdf(size = 1) +
  theme_minimal() +
  labs(x = "Log2(average mRNA expression)", y = "Cumulative Distribution") +
  ggtitle("ECDF of Log2 for each Annotation") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_brewer(name = "Annotation", palette = "RdYlBu")

#saveplot
ggsave("Fig3_b.pdf", plot = Fig3_b, width = 8, height = 5)
