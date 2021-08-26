library(circlize)
library(dplyr)
library(grid)
library(ComplexHeatmap)
library(ggplot2)
library(gggenes)
require(rtracklayer)

#setwd('C:/Users/uqakaraw/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/')
setwd('D:/05.OneDrive/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/')
kar <- read.table('data/chr.kar', header = F, sep = ' ')
head(kar)
kar <- select(kar, 'V4', 'V5', 'V6', 'V7')
head(kar, 10)
kar <- data.frame(lapply(kar, function(x){
  gsub(" ", "\t", x)
}))
head(kar)
tail(kar)

kar[,1] = gsub('HiC_scaffold_', 'bs_chr', kar[,1])
tail(kar)
colnames(kar) <- c('chr', 'start', 'end', 'col')
head(kar)

kar_bs = kar[kar$chr == "bs_chr31",]
kar_bs
kar_ms = kar[kar$chr == 'chr33',]
kar_ms

kar_bs$end <- as.numeric(kar_bs$end)
kar_bs =kar_bs[order(kar_bs$end, decreasing = T),]
kar_bs

kar <- rbind(kar_bs, kar_ms)
kar
head(kar,20)
tail(kar)

getwd()



######################################################################################
gtfms <- rtracklayer::readGFF('data/ms.mhc.gff3')
head(gtfms)
gtfms <- as.data.frame(gtfms)
gtfms <- gtfms[gtfms$seqid == "chr33",]


gff <- gtfms[gtfms$type == "mRNA",]
head(gff)
gff
df <- gff
df$gene <- gsub("MS_gene-", "",df$ref.gene)
df$gene<- gsub("BS_gene-", "", df$gene)
df$gene
dim(df)
head(df)
df1 <- df
df <- dplyr::select(df, 'seqid', 'start', 'end', 'strand', 'ID', 'gene')
head(df)
rownames(df) <- NULL
df$len <- as.numeric(df$end - df$start)
dim(df)
dfms <- df %>%
  group_by(gene) %>%
  dplyr::slice(which.max(len)) %>%
  as.data.frame()
  
dfms
###########################################################################################
gtfbs <- rtracklayer::readGFF('data/bs.mhc.gff3')
head(gtfbs)
gtfbs <- as.data.frame(gtfbs)
gtfbs <- gtfbs[gtfbs$seqid == "HiC_scaffold_31",]
head(gtfbs)
gff <- gtfbs[gtfbs$type == "mRNA",]
head(gff)
df <- gff
df$gene <- gsub("MS_gene-", "",df$ref.gene)
df$gene<- gsub("BS_gene-", "", df$gene)
df$gene
dim(df)
head(df)
df1 <- df
df <- dplyr::select(df, 'seqid', 'start', 'end', 'strand', 'ID', 'gene')
head(df)
rownames(df) <- NULL
df$len <- as.numeric(df$end - df$start)
dim(df)
dfbs <- df %>%
  group_by(gene) %>%
  dplyr::slice(which.max(len)) %>%
  as.data.frame()

dfmerged = rbind(dfbs,dfms)
dim(dfmerged)
df2 <- dfmerged
head(df2)
tail(df2)

ori_bs = as.data.frame(readGFF("data/bs_ori.gff"))
head(ori_bs)
ori_ms = as.data.frame(readGFF("data/ms_chr33.gff"))
head(dfmerged)

ori_ms <- select(ori_ms, gene, product)
head(ori_ms)
ori_bs = select(ori_bs, gene, product)
head(ori_bs)
ori <- rbind(ori_ms, ori_bs)
ori = ori[!is.na(ori$product),]

ori <- ori %>% distinct(gene, .keep_all = T)
ori_merged <- merge(dfmerged, ori, by = "gene", all = F, all.y =F)
head(ori_merged)

############################ Data prep over ##############################################################
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
dfII

dfbs <- ori_merged[ori_merged$seqid == "HiC_scaffold_31",]
dfbs
dfbs$type <- ifelse(startsWith(dfbs$product,"Killer cell lectin-like"), "killer cell lectin", "NA")
dfbs$type <- ifelse(startsWith(dfbs$product, "c-type lectin"), "C-type lectin", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "zinc finger protein"), "Zinc finger protein", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "butyrophilin subfamily"), "Butyrophilin subfamily", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "HLA class II histocompatibility antigen, DM beta"), "DMB", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "HLA class II histocompatibility antigen, DR alpha"), "DMA/DRA", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "bromodomain containing"), "BRD", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "class II histocompatibility antigen, M alpha"), "MA", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "class II histocompatibility antigen, M beta"), "MB", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "E3 ubiquitin-protein ligase TRIM"), "Trim Family", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "class II histocompatibility antigen, B-L beta chain"), "BLB", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "ring finger"), "Ring Finger", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "receptor for activated C kinase 1"), "RACK1", dfbs$type)
dfbs$type <- ifelse(endsWith(dfbs$product, "TRIM7-like") | endsWith(dfbs$product, "RFP") | endsWith(dfbs$product, "TRIM39-like") | startsWith(dfbs$product, "tripartite motif-containing"), "Trim Family", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "class II histocompatibility antigen, DR alpha"), "DMA/DRA", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "transporter 1"), "TAP1", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "transporter 2"), "TAP2", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "tenascin XB"), "TNXB", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "tubulin beta class I"), "TUBB", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "class I histocompatibility antigen, F10 alpha chain-like"), "F10A", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "centromere protein A"), "CENPA", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "flotillin 1"), "FLOT1", dfbs$type)
dfbs$type <- ifelse(dfbs$type == "NA", "Other genes", dfbs$type)
colnames(dfbs)
dfbs$type
dfbs
dfII

#write.table(dfbs, file = "data/dfbs_MHC_anno.tsv", quote = F, row.names = F, sep = "\t")
dfbs <- read.table(file = "data/dfbs_MHC_anno.tsv", header = T, sep = "\t")
dim(dfbs)
dfbs$type <- ifelse(startsWith(dfbs$product, "C-type lectin"), "C-type lectin", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "complement C4-"), "C4", dfbs$type)
dfbs[order(dfbs$start),]

dfms <- ori_merged[ori_merged$seqid == "chr33",]
dfms
dfms$type <- ifelse(startsWith(dfms$product,"Killer cell lectin-like"), "killer cell lectin", "NA")
dfms$type <- ifelse(startsWith(dfms$product, "C-type lectin"), "C-type lectin", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "zinc finger protein"), "Zinc finger protein", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "butyrophilin subfamily"), "Butyrophilin subfamily", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "HLA class II histocompatibility antigen, DM beta"), "DMB", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "HLA class II histocompatibility antigen, DR alpha"), "DMA/DRA", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "bromodomain containing"), "BRD", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "class II histocompatibility antigen, M alpha"), "MA", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "class II histocompatibility antigen, M beta"), "MB", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "E3 ubiquitin-protein ligase TRIM"), "Trim Family", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "class II histocompatibility antigen, B-L beta chain"), "BLB", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "ring finger"), "Ring Finger", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "receptor for activated C kinase 1"), "RACK1", dfms$type)
dfms$type <- ifelse(endsWith(dfms$product, "TRIM7-like") | endsWith(dfms$product, "RFP") | endsWith(dfms$product, "TRIM39-like") | startsWith(dfms$product, "tripartite motif-containing"), "Trim Family", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "class II histocompatibility antigen, DR alpha"), "DMA/DRA", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "transporter 1"), "TAP1", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "transporter 2"), "TAP2", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "tenascin XB"), "TNXB", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "tubulin beta class I"), "TUBB", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "class I histocompatibility antigen, F10 alpha chain-like"), "F10A", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "centromere protein A"), "CENPA", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "flotillin 1"), "FLOT1", dfms$type)
dfms$type <- ifelse(dfms$type == "NA", "Other genes", dfms$type)

colnames(dfms)
dfms$type

#write.table(dfms, file = "data/dfms_MHC_anno.tsv", quote = F, row.names = F, sep = "\t")
dfms <- read.table(file = "data/dfms_MHC_anno.tsv", header = T, sep = "\t")
dfms$type <- ifelse(startsWith(dfms$product, "complement C4-like"), "C4", dfms$type)
dim(dfms)
dfms[order(dfms$start),]

################## Seperate exon of MHC class molecules visualization ####################
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
###################################### End exon visualization ####################################

###################################### Starting complete gene visua. #############################
svg(filename = "figures/MHC_class_bsmsswans.svg", width = 9.5, height = 7.5)
dfmerged_MHC <- rbind(dfbs, dfms)
head(dfmerged_MHC)
ggplot(dfmerged_MHC, aes(xmin = start, xmax = end, y = seqid, forward = (strand == "+"),
                         fill = type)) +
  geom_gene_arrow() +
  facet_wrap(~ seqid, scales = "free", nrow = 4) +
  theme_genes() + geom_gene_arrow()
#scale_fill_brewer(palette = "Pastel2")
dev.off()



df <- dfmerged_MHC
rownames(df) <- NULL
df
df <- select(df, seqid, start, end, type, product, ID)
write.table(df, file = "data/MHC_annotation.tsv", quote = F, row.names = F, sep = "\t")
graphics.off()
head(df)
pdf('data/circlize_ms_bs_mhc.pdf', width = 6, height = 5)
circos.clear()
chromosome.index = c(paste0("chr", c(33)), 
                     rev(paste0("bs_chr", c(31))))
kar
chromosome.index
(kar[1,2] <- 0)
kar[1,3] <- 806914 


##### Start ####
circos.clear()
circos.par(gap.after = c(2,10), cell.padding=c(0,0,0,0), start.degree=85)

circos.initializeWithIdeogram(kar, plotType = NULL, chromosome.index = chromosome.index)
col_text <- 'grey40'
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr=gsub(".*chr", "", CELL_META$sector.index)
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim),mean(ylim),chr,cex=0.6,col=col_text,
              facing="bending.inside",niceFacing=TRUE)
},bg.col="grey90",bg.border=F,track.height=0.06)

# genomes x axis
brk <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,20)*10^6
brk
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h="top",major.at=brk,labels=round(brk/10^6,1),labels.cex=0.4,
              col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
},bg.border=F)

df[,1] <- gsub('HiC_scaffold_', 'bs_chr', df[,1])


df$end - df$start
head(df)
circos.genomicTrack(df, ylim = c(0,0.3), track.height = 0.01, panel.fun = function(region, value, ...){
  circos.genomicRect(region, value, col = as.numeric(factor(value[[1]])), border = 'grey90', ybottom = 0.2, ytop = -0.4)
}, bg.border = F)

#Labels
circos.genomicLabels(df,labels.column=4,cex=0.6,col=as.numeric(factor(df$type)),line_lwd=0.5,
                     line_col=as.numeric(factor(df$type)),
                     side="inside",connection_height=0.08,labels_height=0.025, font = 8, niceFacing = T, padding = 0.01)


#Color legenes for the type of genes
sort(unique(factor(df$type)))
sort(unique(as.numeric(factor(df$type))))
gpar(col = as.numeric(unique(factor(df$type))))
lgd_genes = ComplexHeatmap::Legend(at=sort(unique(factor(df$type))), type = "points", 
                                   legend_gp = gpar(col=sort(as.numeric(unique(factor(df$type))))),
                                   title_position = "lefttop", title = "Genes")


lgd_genes
lgd <- packLegend(lgd_genes)
lgd
draw(lgd, x = unit(6, "mm"), y = unit(15, "mm"), just = c("left", "bottom"))

circos.clear()
dev.off()
############################### The end ###################################################################

graphics.off()
rm(list = ls())

##########################################MHC expression - black swan #####################################
flnc.bed <- genomicDensity(flnc.bed, window.size = 5e5)





