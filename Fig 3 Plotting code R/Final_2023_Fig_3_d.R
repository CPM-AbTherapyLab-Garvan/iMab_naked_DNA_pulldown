# load packages for differential analysis
library(edgeR)
library(limma)
library(Glimma)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Get published GROseq data
# download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94479
# file GSE94479_GROCountNorm.txt.gz
# $gunzip GSE94479_GROCountNorm.txt.gz
# read file
GROCountNorm <- read.table("GSE94479_GROCountNorm.txt", header = TRUE, stringsAsFactors = FALSE)

# read Sampleinfo:
sampleinfo <- read.delim("SampleInfo.txt", stringsAsFactors = TRUE)
#sampleinfo is a txt file as below
#FileName SampleName CellType Status
#1 GSE94479_GROCountNorm     G0_G11     MCF7  G0_G1
#2 GSE94479_GROCountNorm     G0_G12     MCF7  G0_G1
#3 GSE94479_GROCountNorm         S1     MCF7      S
#4 GSE94479_GROCountNorm         S2     MCF7      S
#5 GSE94479_GROCountNorm      G2_M1     MCF7   G2_M
#6 GSE94479_GROCountNorm      G2_M2     MCF7   G2_M

# Remove non count data
countdata <- GROCountNorm[,-(7:8)]
rownames(countdata) <- NULL

# Store RefseqIDs as rownames
rownames(countdata) <- GROCountNorm[,7]

#Check that sample info and the countdata colnames are the same
table(colnames(countdata)==sampleinfo$SampleName)

#counts into DGEList
y <- DGEList(countdata)

#check that y element contains all the counts and samples
#DGEList group as a factor
group <- paste(sampleinfo$Status)
group <- factor(group)
group

#[1] G0_G1 G0_G1 S     S     G2_M  G2_M 
#Levels: G0_G1 G2_M S
# Add the group information into the DGEList
y$samples$group <- group

#Add the annotation using org.Hs.eg.db
ann <- select(org.Hs.eg.db, keys = rownames(y$counts), columns = c("REFSEQ", "ENTREZID"), keytype = "REFSEQ")

#Check refseq matched y$counts
table(ann$REFSEQ==rownames(y$counts))

#Add the annotation in the genes slot of y
y$genes <- ann

# CPM for low expressed filtering
myCPM <- cpm(countdata)

# threshhold at 0.5 CPM
thresh <- myCPM > 0.5
#check how many grather than 0.5 
table(rowSums(thresh))

#Keeping genes at least 2 meeting the thresh based on all samples
keep <- rowSums(thresh) >= 2
summary(keep)

#Filter the DGEList object based on filtering
y <- y[keep, keep.lib.sizes=FALSE]

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)

# just to distributions
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)

# adding a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs")

# MDS plot
plotMDS(y)

# applying composition bias normalisation to DGEList object just in case
y <- calcNormFactors(y)
y$samples

#Design matrix for differential expression
#Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design

#voom transform the data based on the norm.factors from previous step
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)

# using limma to test for differential expression
# Fit the linear model
fit <- lmFit(v)
names(fit)
#[1] "coefficients"     "stdev.unscaled"   "sigma"            "df.residual"      "cov.coefficients"
#[6] "pivot"            "rank"             "Amean"            "method"           "design"         

# make contrasts
cont.matrix <- makeContrasts(G0_G1VsG2_M=G0_G1 - G2_M, G0_G1VsS=G0_G1 - S,G2_MVsS=G2_M - S, levels=design)

# apply the contrast matrix to the fit object
fit.cont <- contrasts.fit(fit, cont.matrix)

# call the eBayers function to performer Bayers shrinkage on the variances and get s-statistics and p-values
fit.cont <- eBayes(fit.cont)

# using limma to get a default summary, decideTests for the DE genes for the contrasts
summa.fit <- decideTests(fit.cont)

summary(summa.fit)
#G0_G1VsG2_M G0_G1VsS G2_MVsS
#Down          4775     5360    3530
#NotSig        6648     5919    9239
#Up            4689     4833    3343

#Write them as tables:
G0_G1VsG2_M <- topTable(fit.cont, coef=1, n=Inf)
G0_G1VsS <- topTable(fit.cont, coef=2, n=Inf)
G2_MVsS <- topTable(fit.cont, coef=3, n=Inf)

#test was decided to better proceed at log fold change of 0.5 to reduce numbers
fit.treat <- treat(fit.cont,lfc=0.5)
res.treat <- decideTests(fit.treat)
summary(res.treat)
#summary(res.treat)
#G0_G1VsG2_M G0_G1VsS G2_MVsS
#Down          1587     1951     734
#NotSig       13095    12538   14807
#Up            1430     1623     571

# Get the upregulated genes for each contrast
upregulated_genes_G0_G1VsG2_M <- rownames(res.treat)[res.treat[,"G0_G1VsG2_M"] == 1]
upregulated_genes_G0_G1VsS <- rownames(res.treat)[res.treat[,"G0_G1VsS"] == 1]
upregulated_genes_G2_MVsS <- rownames(res.treat)[res.treat[,"G2_MVsS"] == 1]

# Get the downregulated genes for each contrast
downregulated_genes_G0_G1VsG2_M <- rownames(res.treat)[res.treat[,"G0_G1VsG2_M"] == -1]
downregulated_genes_G0_G1VsS <- rownames(res.treat)[res.treat[,"G0_G1VsS"] == -1]
downregulated_genes_G2_MVsS <- rownames(res.treat)[res.treat[,"G2_MVsS"] == -1]

# remapped with the control
# in command 
# cd to the place conatining the bigwig files,
# bigwig compare the ratio, log2 computation will be calculated later after getting the average
# bigwigCompare -p 10 --bigwig1 DIPMCF7_1_sorted_dupmrkd_normalized.bw --bigwig2 INPUT_DP1MCF7_1_sorted_dupmrkd_normalized.bw --skipZeroOverZero --operation ratio --outFileName DIPMCF7_1_minus_INPUT_ratio.bw
# in command also: 
# computeMatrix reference-point -S DIPMCF7_1_minus_INPUT_ratio.bw -R /Users/cris/Desktop/chip/refGene_hg19_TSS.bed --referencePoint center -a 1000 -b 1000 -out DIPMCF7_1_minus_INPUT_ratio_on_TSS_refseq_1000bp.tab.gz -p 10 --skipZeros --missingDataAsZero &
# read the  compressed MATLAB file into a dataframe from the results in Deeptools-ComputeMatrix
DIPMCF7on_TSS_1000 <- read.table("DIP_RNA_integration_plots/DIPMCF7_1_minus_INPUT_ratio_on_TSS_refseq_1000bp.tab.gz"
                                 , header=FALSE, sep="\t", skip=1)
# dimensions
dim(DIPMCF7on_TSS_1000)

#rows values mean rowMeans(df)
DIPMCF7on_TSS_1000_means <- DIPMCF7on_TSS_1000[,1:4]
DIPMCF7on_TSS_1000_means$mean_signal <- rowMeans(DIPMCF7on_TSS_1000[,7:206])

# Alternatively with max signal
#DIPMCF7on_TSS_1000_max_signal <- DIPMCF7on_TSS_1000[,1:4]
#DIPMCF7on_TSS_1000_max_signal$max_signal <- apply(DIPMCF7on_TSS_1000[,7:206], 1, max)
#rename columns
colnames(DIPMCF7on_TSS_1000_means) <- c("Chr", "Start", "End", "Coordinate", "Mean_signal")

# Read in the file for refGene with TSS coordinates
refGene_hg19_TSS <- read.table("refGene_hg19_TSS.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(refGene_hg19_TSS) <- c("Chr", "Start", "End", "refseq_mrna", "GeneSymbol")

#remove the version from each refseq identifier:
refGene_hg19_TSS$refseq_mrna <- gsub("\\..*","",refGene_hg19_TSS$refseq_mrna)

# map identifiers using mapIds
entrez_ids <- mapIds(org.Hs.eg.db, keys=refGene_hg19_TSS$refseq_mrna, column="ENTREZID", keytype="REFSEQ")

# add the ENTREZIDs to your data frame
refGene_hg19_TSS$ENTREZID <- entrez_ids

# create common ID in refGene_hg19_TSS
refGene_hg19_TSS$coor <- paste(refGene_hg19_TSS$Chr, refGene_hg19_TSS$Start, refGene_hg19_TSS$End, sep="_")

# create common coordinate ID in DIPMCF7on_TSS_1000_means
DIPMCF7on_TSS_1000_means$coor <- paste(DIPMCF7on_TSS_1000_means$Chr, DIPMCF7on_TSS_1000_means$Start, DIPMCF7on_TSS_1000_means$End, sep="_")
DIPMCF7on_TSS_1000_means$coor <- sub("chr", "", DIPMCF7on_TSS_1000_means$coor)

# create common coordinate ID in refGene_hg19_TSS
refGene_hg19_TSS$coor <- paste(refGene_hg19_TSS$Chr, refGene_hg19_TSS$Start, refGene_hg19_TSS$End, sep="_")

# merge 'refGene_hg19_TSS' and 'DIPMCF7on_TSS_1000_means' based on 'coor'
merged_TSS_DIPMCF7_1000_means <- merge(refGene_hg19_TSS, DIPMCF7on_TSS_1000_means[c("coor", "Mean_signal")], by="coor")

# Read Bulk RNA file the file into R
RNA <- read.delim("RNA_WT/MCF7_WT_counts.txt", header = TRUE, sep = "\t", fill = TRUE, skip = 1)

#rename columns
colnames(RNA) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", "Rep1", "Rep2", "Rep3")

# Get the ENTREZIDs
RNA$ENTREZID <- mapIds(org.Hs.eg.db, 
                       keys=RNA$Geneid, 
                       column="ENTREZID", 
                       keytype="SYMBOL", 
                       multiVals="first")

# read the peaks values
MCF7_1_MACS <- read.table("MCF7_1_MACS_peaks.narrowPeak", header=FALSE, sep="\t")
#rename columns
colnames(MCF7_1_MACS) <- c("Chr", "Start", "End", "name", "UCSC_score", "strand", "fold_change", "minuslog10pvalue", "minuslog10qvalue", "relative_summit")

#read the annotated file
MCF7_1_MACS_annotated <- read_tsv("MCF7_1_MACS_peaks_annotated.narrowPeak")

#change the peak name to $name and $"Entrez ID" to a single word
colnames(MCF7_1_MACS_annotated)[colnames(MCF7_1_MACS_annotated) == 'PeakID (cmd=annotatePeaks.pl MCF7_1_MACS_peaks.narrowPeak hg19)'] <- 'name'
colnames(MCF7_1_MACS_annotated)[colnames(MCF7_1_MACS_annotated) == "Entrez ID"] <- "ENTREZID"

#merge the annotated file and the peak values by name
MCF7_1_merged_annotated <- merge(MCF7_1_MACS, MCF7_1_MACS_annotated, by = "name", all = TRUE)

# Merge both RNA seq data with the peak file 
RNA_MCF7_1_merged_annotated <- merge(MCF7_1_merged_annotated, RNA, by = "ENTREZID", all = TRUE)

# Rename columns in merged_TSS_DIPMCF7_1000_means
colnames(merged_TSS_DIPMCF7_1000_means) <- c("coor", "Chr", "Start1000pos", "End1000pos", "refseq_mrna", "GeneSymbol", "ENTREZID", "Mean_signal")

# Deduplicate by the mean of Mean_signal for each ENTREZID and ungroup
deduplicated_TSS_DIPMCF7_1000_means <- merged_TSS_DIPMCF7_1000_means %>%
  group_by(ENTREZID) %>% summarize(Mean_signal = mean(Mean_signal, na.rm = TRUE),
    refseq_mrna = dplyr::first(refseq_mrna)) %>%  ungroup()

# Check dimensions
dim(deduplicated_TSS_DIPMCF7_1000_means)
colnames(deduplicated_TSS_DIPMCF7_1000_means)

# final merge with the MCF7 RNAseq data
final_RNA_TSS_DIPMCF7_1000_means <- merge(deduplicated_TSS_DIPMCF7_1000_means, RNA_MCF7_1_merged_annotated, by = "ENTREZID", all.x = TRUE)

# check the unique annotations in the $Annotation column
unique(final_RNA_TSS_DIPMCF7_1000_means$Annotation)

#deleting the de description between parenthesis from the $Annotation column 
final_RNA_TSS_DIPMCF7_1000_means$Annotation <- gsub("\\(.*\\)", "", final_RNA_TSS_DIPMCF7_1000_means$Annotation)

#Also deleting any (dot + number) after the annotationexample "TTS .7" "5' UTR .6"
final_RNA_TSS_DIPMCF7_1000_means$Annotation <- gsub("\\..*", "", final_RNA_TSS_DIPMCF7_1000_means$Annotation)

#Check $Annotation identifiers
unique(final_RNA_TSS_DIPMCF7_1000_means$Annotation)

# combine all genes into a single dataframe with a column to indicate group
upregulated_G0_G1VsG2_M_group <- final_RNA_TSS_DIPMCF7_1000_means[final_RNA_TSS_DIPMCF7_1000_means$refseq_mrna %in% upregulated_genes_G0_G1VsG2_M, c("refseq_mrna", "Mean_signal")]
upregulated_G0_G1VsG2_M_group$Group <- "Up G0_G1VsG2_M"
upregulated_G0_G1VsS_group <- final_RNA_TSS_DIPMCF7_1000_means[final_RNA_TSS_DIPMCF7_1000_means$refseq_mrna %in% upregulated_genes_G0_G1VsS, c("refseq_mrna", "Mean_signal")]
upregulated_G0_G1VsS_group$Group <- "Up G0_G1VsS"
upregulated_G2_MVsS_group <- final_RNA_TSS_DIPMCF7_1000_means[final_RNA_TSS_DIPMCF7_1000_means$refseq_mrna %in% upregulated_genes_G2_MVsS, c("refseq_mrna", "Mean_signal")]
upregulated_G2_MVsS_group$Group <- "Up G2_MVsS"
downregulated_G0_G1VsG2_M_group <- final_RNA_TSS_DIPMCF7_1000_means[final_RNA_TSS_DIPMCF7_1000_means$refseq_mrna %in% downregulated_genes_G0_G1VsG2_M, c("refseq_mrna", "Mean_signal")]
downregulated_G0_G1VsG2_M_group$Group <- "Down G0_G1VsG2_M"
downregulated_G0_G1VsS_group <- final_RNA_TSS_DIPMCF7_1000_means[final_RNA_TSS_DIPMCF7_1000_means$refseq_mrna %in% downregulated_genes_G0_G1VsS, c("refseq_mrna", "Mean_signal")]
downregulated_G0_G1VsS_group$Group <- "Down G0_G1VsS"
downregulated_G2_MVsS_group <- final_RNA_TSS_DIPMCF7_1000_means[final_RNA_TSS_DIPMCF7_1000_means$refseq_mrna %in% downregulated_genes_G2_MVsS, c("refseq_mrna", "Mean_signal")]
downregulated_G2_MVsS_group$Group <- "Down G2_MVsS"

#combine data
genes_grouped <- rbind(upregulated_G0_G1VsG2_M_group, upregulated_G0_G1VsS_group, upregulated_G2_MVsS_group, downregulated_G0_G1VsG2_M_group, downregulated_G0_G1VsS_group, downregulated_G2_MVsS_group)

# get unique refseq_mrna and Comparison rows so duplicated data doesn't get counted several times
genes_grouped_unique <- genes_grouped %>%
  distinct(refseq_mrna, Group, .keep_all = TRUE)

# Create boxplot
#ggplot(genes_grouped_unique, aes(x = Group, y = log2(Mean_signal), fill = Group)) +
#  geom_boxplot() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#  labs(x = "Group", y = "Log2(Mean IP Signal TSS +- 1kbp)", fill = "Group") +
#  theme_minimal()
#SEPARATED
# Add a new column for comparison
genes_grouped_unique$Comparison <- ifelse(genes_grouped_unique$Group %in% c("Up G0_G1VsG2_M", "Down G0_G1VsG2_M"), "G0_G1VsG2_M",
                                   ifelse(genes_grouped_unique$Group %in% c("Up G0_G1VsS", "Down G0_G1VsS"), "G0_G1VsS",
                                          "G2_MVsS"))

# Create boxplot
Fig3_d_box <- ggplot(genes_grouped_unique, aes(x = Group, y = log2(Mean_signal), fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Up G0_G1VsG2_M" = "brown2", 
                               "Down G0_G1VsG2_M" = "deepskyblue3", 
                               "Up G0_G1VsS" = "brown2", 
                               "Down G0_G1VsS" = "deepskyblue3",
                               "Up G2_MVsS" = "brown2",
                               "Down G2_MVsS" = "deepskyblue3"))  +
  labs(x = "Group", y = "log2(Mean Signal)") + 
  theme_minimal() +
  facet_wrap(~Comparison, scales = "free_x", ncol = 3) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  guides(fill=FALSE) +
  coord_cartesian( ylim = c(-6, 2)) 

#plot
Fig3_d_box

#saveplot
ggsave("Fig3_d_box.pdf", plot = Fig3_d_box, width = 4.5, height = 3.5)

# List of comparisons
comparisons <- c("G0_G1VsG2_M", "G0_G1VsS", "G2_MVsS")

# Test if the mean signal is normally distributed
shapiro.test(log2(genes_grouped_unique$Mean_signal))
#	Shapiro-Wilk normality test
#data:  log2(subset_genes_grouped_unique$Mean_signal)
#W = 0.93818, p-value < 2.2e-16
#Mean_signal data is also not normally distributed

# Function to perform a Wilcoxon rank-sum test
perform_wilcoxon_test <- function(comp) {
  # Subset the data for each comparison
  data <- subset(genes_grouped_unique, Comparison == comp)
  # Split the data into two groups based on up- and down-regulation
  group1 <- subset(data, Group == paste0("Up ", comp))$Mean_signal
  group2 <- subset(data, Group == paste0("Down ", comp))$Mean_signal
  # Perform the Wilcoxon test and return the result
  return(wilcox.test(group1, group2))
}

# Apply the function to each comparison
results <- lapply(comparisons, perform_wilcoxon_test)

# Print the results
names(results) <- comparisons
print(results)

#for volcano plots with genes associated to gene transcript iMotifs
# Add a new column to the dataframe to specify whether a gene is upregulated or downregulated
G0_G1VsG2_M$differential <- ifelse(G0_G1VsG2_M$REFSEQ %in% upregulated_genes_G0_G1VsG2_M, "Upregulated",
                                   ifelse(G0_G1VsG2_M$REFSEQ %in% downregulated_genes_G0_G1VsG2_M, "Downregulated", "Others"))
G0_G1VsS$differential <- ifelse(G0_G1VsS$REFSEQ %in% upregulated_genes_G0_G1VsS, "Upregulated",
                                ifelse(G0_G1VsS$REFSEQ %in% downregulated_genes_G0_G1VsS, "Downregulated", "Others"))
G2_MVsS$differential <- ifelse(G2_MVsS$REFSEQ %in% upregulated_genes_G2_MVsS, "Upregulated",
                               ifelse(G2_MVsS$REFSEQ %in% downregulated_genes_G2_MVsS, "Downregulated", "Others"))

# Create a new column for -Log10(P.Value)
G0_G1VsG2_M$minusLog10PValue <- -log10(G0_G1VsG2_M$`P.Value`)
G0_G1VsS$minusLog10PValue <- -log10(G0_G1VsS$`P.Value`)
G2_MVsS$minusLog10PValue <- -log10(G2_MVsS$`P.Value`)

# Create volcano plot
#ggplot(data = G0_G1VsG2_M, aes(x = logFC, y = minusLog10PValue, color = differential)) +  geom_point(alpha = 0.5) +  scale_color_manual(values = c("Upregulated" = "brown2", "Downregulated" = "deepskyblue3", "Others" = "black")) + theme_minimal() +  labs(x = "logFC", y = "-Log10(P-Value)") + theme(legend.position = "right")

# create the other 2 volcano plots as well
#ggplot(data = G0_G1VsS, aes(x = logFC, y = minusLog10PValue, color = differential)) +  geom_point(alpha = 0.5) +  scale_color_manual(values = c("Upregulated" = "brown2", "Downregulated" = "deepskyblue3", "Others" = "black")) +  theme_minimal() +  labs(x = "logFC", y = "-Log10(P-Value)") +  theme(legend.position = "right")
#ggplot(data = G2_MVsS, aes(x = logFC, y = minusLog10PValue, color = differential)) +  geom_point(alpha = 0.5) +  scale_color_manual(values = c("Upregulated" = "brown2", "Downregulated" = "deepskyblue3", "Others" = "black")) +  theme_minimal() +  labs(x = "logFC", y = "-Log10(P-Value)") +  theme(legend.position = "right")

# Create gene_related dataframe to color plots
gene_related <- final_RNA_TSS_DIPMCF7_1000_means %>%
  mutate(Annotation = recode(Annotation,
                             "Intergenic" = "non-related",
                             "intron " = "non-related",
                             "exon " = "gene-related",
                             "3' UTR " = "gene-related",
                             "non-coding " = "non-related",
                             "promoter-TSS " = "gene-related",
                             "5' UTR " = "gene-related",
                             "TTS " = "gene-related")) %>%
  filter(Annotation == "gene-related")

# Get the indices of upregulated and downregulated genes in gene_related dataframe
idx_upregulatedG0_G1VsG2_M <- which(gene_related$`Nearest Refseq` %in% upregulated_genes_G0_G1VsG2_M)
idx_downregulatedG0_G1VsG2_M <- which(gene_related$`Nearest Refseq` %in% downregulated_genes_G0_G1VsG2_M)

# Compute the ratio of upregulated and downregulated genes in gene_related
ratio_upregulatedG0_G1VsG2_M <- length(idx_upregulatedG0_G1VsG2_M) / length(upregulated_genes_G0_G1VsG2_M)
ratio_downregulatedG0_G1VsG2_M <- length(idx_downregulatedG0_G1VsG2_M) / length(downregulated_genes_G0_G1VsG2_M)

# Print the results
cat("Ratio of upregulated genes in gene_related: ", ratio_upregulatedG0_G1VsG2_M, "\n")
cat("Ratio of downregulated genes in gene_related: ", ratio_downregulatedG0_G1VsG2_M)

# Get the indices of upregulated and downregulated genes in gene_related dataframe
idx_upregulatedG0_G1VsS <- which(gene_related$`Nearest Refseq` %in% upregulated_genes_G0_G1VsS)
idx_downregulatedG0_G1VsS <- which(gene_related$`Nearest Refseq` %in% downregulated_genes_G0_G1VsS)

# Compute the ratio of upregulated and downregulated genes in gene_related
ratio_upregulatedG0_G1VsS <- length(idx_upregulatedG0_G1VsS) / length(upregulated_genes_G0_G1VsS)
ratio_downregulatedG0_G1VsS <- length(idx_downregulatedG0_G1VsS) / length(downregulated_genes_G0_G1VsS)

# Print the result
cat("Ratio of upregulated genes in gene_related: ", ratio_upregulatedG0_G1VsS, "\n")
cat("Ratio of downregulated genes in gene_related: ", ratio_downregulatedG0_G1VsS)

# Get the indices of upregulated and downregulated genes in gene_related dataframe
idx_upregulatedG2_MVsS <- which(gene_related$`Nearest Refseq` %in% upregulated_genes_G2_MVsS)
idx_downregulatedG2_MVsS <- which(gene_related$`Nearest Refseq` %in% downregulated_genes_G2_MVsS)

# Compute the ratio of upregulated and downregulated genes in gene_related
ratio_upregulatedG2_MVsS <- length(idx_upregulatedG2_MVsS) / length(upregulated_genes_G2_MVsS)
ratio_downregulatedG2_MVsS <- length(idx_downregulatedG2_MVsS) / length(downregulated_genes_G2_MVsS)

# Print the results
cat("Ratio of upregulated genes in gene_related: ", ratio_upregulatedG2_MVsS, "\n")
cat("Ratio of downregulated genes in gene_related: ", ratio_downregulatedG2_MVsS)

#coloring orange for &Related
# Create a new column 'gene_related_type' in your dataframe
G0_G1VsG2_M$gene_related_type <- G0_G1VsG2_M$differential

# Find indices of the gene_related dataframe in G0_G1VsG2_M
idx_upregulated_related_G0_G1VsG2_M <- which(G0_G1VsG2_M$REFSEQ %in% gene_related$`Nearest Refseq` & G0_G1VsG2_M$differential == "Upregulated")
idx_downregulated_related_G0_G1VsG2_M <- which(G0_G1VsG2_M$REFSEQ %in% gene_related$`Nearest Refseq` & G0_G1VsG2_M$differential == "Downregulated")

# Change 'gene_related_type' of these indices
G0_G1VsG2_M$gene_related_type[idx_upregulated_related_G0_G1VsG2_M] <- "Upregulated & gene related"
G0_G1VsG2_M$gene_related_type[idx_downregulated_related_G0_G1VsG2_M] <- "Downregulated & gene related"

# Create a volcano plot
Fig3_d_G0_G1VsG2_M <- ggplot(data = G0_G1VsG2_M, aes(x = logFC, y = minusLog10PValue, color = gene_related_type)) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = c("Upregulated" = "salmon", "Downregulated" = "skyblue", "Upregulated & gene related" = "purple3", "Downregulated & gene related" = "purple3", "Others" = "black")) +
  theme_minimal() +
  labs(x = "logFC", y = "-Log10(P-Value)") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(-6, 6), ylim = c(0, 30)) 

#plot
Fig3_d_G0_G1VsG2_M

#saveplot
ggsave("Fig3_d_G0_G1VsG2_M.pdf", plot = Fig3_d_G0_G1VsG2_M, width = 3.5, height = 3.5)

#coloring orange for &Related
# Create a new column 'gene_related_type' in your dataframe
G0_G1VsS$gene_related_type <- G0_G1VsS$differential

# Find indices of the gene_related dataframe in G0_G1VsS
idx_upregulated_related_G0_G1VsS <- which(G0_G1VsS$REFSEQ %in% gene_related$`Nearest Refseq` & G0_G1VsS$differential == "Upregulated")
idx_downregulated_related_G0_G1VsS <- which(G0_G1VsS$REFSEQ %in% gene_related$`Nearest Refseq` & G0_G1VsS$differential == "Downregulated")

# Change 'gene_related_type' of these indices
G0_G1VsS$gene_related_type[idx_upregulated_related_G0_G1VsS] <- "Upregulated & gene related"
G0_G1VsS$gene_related_type[idx_downregulated_related_G0_G1VsS] <- "Downregulated & gene related"

# Create a volcano plot
Fig3_d_G0_G1VsS <- ggplot(data = G0_G1VsS, aes(x = logFC, y = minusLog10PValue, color = gene_related_type)) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = c("Upregulated" = "salmon", "Downregulated" = "skyblue", "Upregulated & gene related" = "purple3", "Downregulated & gene related" = "purple3", "Others" = "black")) +
  theme_minimal() +
  labs(x = "logFC", y = "-Log10(P-Value)") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(-6, 6), ylim = c(0, 30)) 

#plot
Fig3_d_G0_G1VsS

#saveplot
ggsave("Fig3_d_G0_G1VsS.pdf", plot = Fig3_d_G0_G1VsS, width = 3.5, height = 3.5)

#coloring orange for &Related
# Create a new column 'gene_related_type' in your dataframe
G2_MVsS$gene_related_type <- G2_MVsS$differential

# Find indices of the gene_related dataframe in G2_MVsS
idx_upregulated_related_G2_MVsS <- which(G2_MVsS$REFSEQ %in% gene_related$`Nearest Refseq` & G2_MVsS$differential == "Upregulated")
idx_downregulated_related_G2_MVsS <- which(G2_MVsS$REFSEQ %in% gene_related$`Nearest Refseq` & G2_MVsS$differential == "Downregulated")

# Change 'gene_related_type' of these indices
G2_MVsS$gene_related_type[idx_upregulated_related_G2_MVsS] <- "Upregulated & gene related"
G2_MVsS$gene_related_type[idx_downregulated_related_G2_MVsS] <- "Downregulated & gene related"

# Create a volcano plot
Fig3_d_G2_MVsS <- ggplot(data = G2_MVsS, aes(x = logFC, y = minusLog10PValue, color = gene_related_type)) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = c("Upregulated" = "salmon", "Downregulated" = "skyblue", "Upregulated & gene related" = "purple3", "Downregulated & gene related" = "purple3", "Others" = "black")) +
  theme_minimal() +
  labs(x = "logFC", y = "-Log10(P-Value)") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(-6, 6), ylim = c(0, 30)) 

#plot
Fig3_d_G2_MVsS

#saveplot
ggsave("Fig3_d_G2_MVsS.pdf", plot = Fig3_d_G2_MVsS, width = 3.5, height = 3.5)