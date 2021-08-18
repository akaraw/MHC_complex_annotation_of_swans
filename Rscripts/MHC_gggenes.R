#install.packages("gggenes")
library(ggplot2)
library(gggenes)

setwd("D:/05.OneDrive/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/")
df <- read.table("data/MHC_loci_in_MS_BS.txt", sep = "\t", header = T)
head(df)
df$exon_no <- paste0("exon_", df$exon_no)
df$exon_no <- paste0(df$exon_no, " of ", df$MHC_class)
df$exon_no
df$y_axis <- paste0(df$seqname, " ", df$MHC_class)
df

dim(df)

dfI <- df[df$MHC_class == "MHC I",]
dfII <- df[df$MHC_class == "MHC II",]

dfI
bsMHCI_min_loc <- 565331 - 10000
bsMHCI_max_loc <- 576889 + 10000
bsMHCI_max_loc
dfbs <- df2[df2$seqname == "HiC_scaffold_31",]
dfbs$end
dfbs$start
dfbsMHC_ClassI <- dfbs[dfbs$start > bsMHCI_min_loc & dfbs$end < bsMHCI_max_loc,]
dfbsMHC_ClassI

msMHCI_min_loc <- 242591 - 10000
msMHCI_max_loc <- 251768 + 10000
dfms <- df2[df2$seqnames == "chr33",]
dfmsMHC_classI_loc1 <- dfbs[dfms$start > msMHCI_min_loc & dfms$end < msMHCI_max_loc,]
dfmsMHC_classI_loc1

dfII
bsMHCII_min_loc1 <- 695817 - 10000
bsMHCII_max_loc1 <- 703901 + 10000
msMHCII_max_loc1 <- 131501 + 10000
msMHCII_min_loc1 <- 130558 - 10000
msMHCII_min_loc2 <- 137914 - 10000
msMHCII_max_loc2 <- 138857 + 10000 

(dfbs[dfbs$start >= bsMHCI_min_loc & dfbs$end <= bsMHCII_max_loc1,])


dfbs
dfms

dfbs$type <- ifelse(startsWith(dfbs$product,"Killer cell lectin-like"), "killer cell lectin", "NA")
dfbs$type <- ifelse(startsWith(dfbs$product, "C-type lectin"), "C-type lectin", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "Zinc finger protein"), "Zinc finger protein", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "Butyrophilin subfamily"), "Butyrophilin subfamily", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "HLA class II histocompatibility antigen, DM beta"), "DMB", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "Mamu class II histocompatibility antigen, DR alpha"), "DMA/DRA", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "Bromodomain-containing"), "BRD", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "Class II histocompatibility antigen, M alpha"), "M-alpha", dfbs$type)
dfbs$type <- ifelse(endsWith(dfbs$product, "TRIM7") | endsWith(dfbs$product, "RFP") | endsWith(dfbs$product, "TRIM39"), "TRIM Family", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "Mamu class II histocompatibility antigen, DR alpha"), "DMA/DRA", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "Antigen peptide transporter 1"), "TAP1", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "Antigen peptide transporter 2"), "TAP2", dfbs$type)
colnames(dfbs)
dfII
#write.table(dfbs, file = "data/dfbs_MHC_anno.tsv", quote = F, row.names = F, sep = "\t")
dfbs <- read.table(file = "data/dfbs_MHC_anno.tsv", header = T, sep = "\t")
dim(dfbs)
dim(dfbs <- dfbs[!dfbs$product == "hypothetical protein",])
dfbs$type
dfbs$type <- ifelse(is.na(dfbs$type), "Other genes", dfbs$type)



dfms
dfms$type <- ifelse(startsWith(dfms$product,"Killer cell lectin-like"), "killer cell lectin", "NA")
dfms$type <- ifelse(startsWith(dfms$product, "C-type lectin"), "C-type lectin", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "Zinc finger protein"), "Zinc finger protein", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "Butyrophilin subfamily"), "Butyrophilin subfamily", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "HLA class II histocompatibility antigen, DM beta"), "DMB", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "Mamu class II histocompatibility antigen, DR alpha"), "DMA/DRA", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "Bromodomain-containing"), "BRD", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "Class II histocompatibility antigen, M alpha"), "M-alpha", dfms$type)
dfms$type <- ifelse(endsWith(dfms$product, "TRIM7") | endsWith(dfms$product, "RFP") | endsWith(dfms$product, "TRIM39") | endsWith(dfms$product, "TRIM41"), 
                    "TRIM Family", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "Mamu class II histocompatibility antigen, DR alpha"), "DMA/DRA", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "Antigen peptide transporter 1"), "TAP1", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "Antigen peptide transporter 2"), "TAP2", dfms$type)
colnames(dfms)
dfms$type
dfms

#write.table(dfms, file = "data/dfms_MHC_anno.tsv", quote = F, row.names = F, sep = "\t")
dfms <- read.table(file = "data/dfms_MHC_anno.tsv", header = T, sep = "\t")
dim(dfms)
dim(dfms <- dfms[!dfms$product == "hypothetical protein",])
dfms$type
dfms$type <- ifelse(is.na(dfms$type), "Other genes", dfms$type)
dfms$product

svg(filename = "figures/MHC_class_I_swans.svg", width = 6, height = 4.5)
ggplot(dfI, aes(xmin = start, xmax = end, y = y_axis, 
                fill = exon_no)) +
  geom_gene_arrow() +
  facet_wrap(~ seqname, scales = "free", nrow = 4) +
  theme_genes() + geom_gene_arrow()
#scale_fill_brewer(palette = "Pastel2")
dev.off()
svg(filename = "figures/MHC_class_II_swans.svg", width = 6, height = 4.5)
ggplot(dfII, aes(xmin = start, xmax = end, y = y_axis, 
                fill = exon_no)) +
  geom_gene_arrow() +
  facet_wrap(~ seqname, scales = "free", nrow = 4) +
  theme_genes() + geom_gene_arrow()
#scale_fill_brewer(palette = "Pastel2")
dev.off()

svg(filename = "figures/MHC_class_bsmsswans.svg", width = 9.5, height = 7.5)
dfmerged_MHC <- rbind(dfbs, dfms)
ggplot(dfmerged_MHC, aes(xmin = start, xmax = end, y = seqnames, forward = (strand == "+"),
                fill = type)) +
  geom_gene_arrow() +
  facet_wrap(~ seqnames, scales = "free", nrow = 4) +
  theme_genes() + geom_gene_arrow()
#scale_fill_brewer(palette = "Pastel2")
dev.off()

