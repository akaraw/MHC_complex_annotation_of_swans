########################################################## Loading librariries ##################################################
library(circlize)
library(dplyr)
library(grid)
library(ComplexHeatmap)
library(ggplot2)
library(gggenes)
require(rtracklayer)
library(GenomicRanges)
library(IRanges)
require(viridis)

############################################################## Set up envioronment #############################################
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
head(kar,20)
tail(kar)

########################################################## Import annotations - Mute swan ############################################
gtfms <- rtracklayer::readGFF('data/ms.mhc.gff3')
head(gtfms)
gtfms <- as.data.frame(gtfms)
gtfms <- gtfms[gtfms$seqid == "chr33",]


gff <- gtfms[gtfms$type == "mRNA",]
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
dfms <- df %>%
  group_by(gene) %>%
  dplyr::slice(which.max(len)) %>%
  as.data.frame()

########################################################## Import annotations - Black swan ##########################################
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
################################################### Data prep over ##############################################################

################################################### Data type selection #########################################################
df <- read.table("data/MHC_loci_in_MS_BS.txt", sep = "\t", header = T)
head(df)
df$exon_no <- paste0("exon_", df$exon_no)
df$exon_no <- paste0(df$exon_no, " of ", df$MHC_class)
df$y_axis <- paste0(df$seqname, " ", df$MHC_class)
dim(df)

dfI <- df[df$MHC_class == "MHC I",]
dfII <- df[df$MHC_class == "MHC II",]

dfbs <- ori_merged[ori_merged$seqid == "HiC_scaffold_31",]
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
dfbs$type <- ifelse(startsWith(dfbs$product, "C-type lectin"), "C-type lectin", dfbs$type)
dfbs$type <- ifelse(startsWith(dfbs$product, "complement C4-"), "C4", dfbs$type)
#write.table(dfbs, file = "data/dfbs_MHC_anno.tsv", quote = F, row.names = F, sep = "\t")
dfbs <- read.table(file = "data/dfbs_MHC_anno.tsv", header = T, sep = "\t")
dim(dfbs)

dfms <- ori_merged[ori_merged$seqid == "chr33",]
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
dfms$type <- ifelse(startsWith(dfms$product, "complement C4-like"), "C4", dfms$type)
dfms$type <- ifelse(startsWith(dfms$product, "complement C4-like"), "C4", dfms$type)
#write.table(dfms, file = "data/dfms_MHC_anno.tsv", quote = F, row.names = F, sep = "\t") #After writing manually enter the MHC loci
dfms <- read.table(file = "data/dfms_MHC_anno.tsv", header = T, sep = "\t")
dim(dfms)

###################################################### Seperate exon of MHC class molecules visualization ####################
#svg(filename = "figures/MHC_class_I_swans.svg", width = 6, height = 4.5)
ggplot(dfI, aes(xmin = start, xmax = end, y = y_axis, 
                fill = exon_no)) +
  geom_gene_arrow() +
  facet_wrap(~ seqname, scales = "free", nrow = 4) +
  theme_genes() + geom_gene_arrow()
#scale_fill_brewer(palette = "Pastel2")
#dev.off()
#svg(filename = "figures/MHC_class_II_swans.svg", width = 6, height = 4.5)
ggplot(dfII, aes(xmin = start, xmax = end, y = y_axis, 
                 fill = exon_no)) +
  geom_gene_arrow() +
  facet_wrap(~ seqname, scales = "free", nrow = 4) +
  theme_genes() + geom_gene_arrow()
#dev.off()
################################################################# End exon visualization ####################################

################################################################# Starting complete gene visualization ######################
#svg(filename = "figures/MHC_class_bsmsswans.svg", width = 9.5, height = 7.5)
dfmerged_MHC <- rbind(dfbs, dfms)
head(dfmerged_MHC)
ggplot(dfmerged_MHC, aes(xmin = start, xmax = end, y = seqid, forward = (strand == "+"),
                         fill = type)) +
  geom_gene_arrow() +
  facet_wrap(~ seqid, scales = "free", nrow = 4) +
  theme_genes() + geom_gene_arrow()
#scale_fill_brewer(palette = "Pastel2")
#dev.off()

df <- dfmerged_MHC
rownames(df) <- NULL
df <- select(df, seqid, start, end, type, product, ID)
write.table(df, file = "data/MHC_annotation.tsv", quote = F, row.names = F, sep = "\t")
#graphics.off()
head(df)
#Circos - 01
#pdf('data/circlize_ms_bs_mhc.pdf', width = 6, height = 5)
circos.clear()
chromosome.index = c(paste0("chr", c(33)), 
                     rev(paste0("bs_chr", c(31))))
chromosome.index
(kar[1,2] <- 0)
kar[1,3] <- 806914 
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
#dev.off()
########################################## The end - first pipeline###################################################################

########################################## MHC expression - swans ####################################################################
dfbk <- df
df2 <- dfbk

########################################## Start MHC expression analysis  Black swan #################################################
flnc.bs.bed <- read.table("data/bs2flnc.bed", header = F)
head(flnc.bs.bed)
flnc.bs.bed$V1 <- gsub('HiC_scaffold_', 'bs_chr', flnc.bs.bed[,1])
flnc.bs.bed <- flnc.bs.bed[flnc.bs.bed$V1 == "bs_chr31",]
flnc.bs.bed <- genomicDensity(flnc.bs.bed, window.size = 500, count_by = "percent")
write.table(flnc.bs.bed, "data/flnc.bs.bed", sep = "\t", quote = F, row.names = F)
head(flnc.bs.bed)
tail(flnc.bs.bed)

tagments <- data.frame(tagments=paste(df2$seqid, df2$start, df2$end, sep = "-"), start=df2$start, end=df2$end, chr=df2$seqid)
head(tagments)
correspondance <- cbind(dplyr::select(df2, seqid, start,end), tagments)
correspondance <-correspondance[,-7]
tail(tagments)
head(correspondance)

isnps <- with(flnc.bs.bed, IRanges(start, end, names=chr, value=value))
head(tagments)
tagmentsbs <- tagments[tagments$chr == "bs_chr31",]
head(tagmentsbs)
igenes <- with(tagmentsbs, IRanges(start, end, names=tagments, chr=chr))
olaps <- findOverlaps(isnps, igenes)
bs <- cbind(flnc.bs.bed[queryHits(olaps),], tagmentsbs[subjectHits(olaps),])
tail(bs)
colnames(bs) <- c("chr","start1","end1","value","tagments","start","end","chr")
bs <- dplyr::select(bs, tagments, start1, end1, value)
head(bs)
chr_bg_color = rand_color(2, transparency = 0.8)
names(chr_bg_color) = c("bs_chr31", "chr33")

head(correspondance)
kar
kar_bs
head(df2) 

#check the MHC locations
#bs_chr31 541642 562672(1083 - ) bs_chr31 622721 625612 (- 1252) - MHC I
#MHC class I 82     bs_chr31 695164 697066  bs_chr31 707341 712959
rownames(bs) <- NULL
bs[bs$tagments == "bs_chr31-541642-562672",]
bs[bs$tagments == "bs_chr31-622721-625612",]
bs[bs$tagments == "bs_chr31-695164-697066",]
bs[bs$tagments == "bs_chr31-707341-712959",]

head(bs)
bs_MHC_I <- bs[c(1469:1702),]
bs_MHC_II <- bs[c(1861:1904),]
head(bs_MHC_I)
head(bs_MHC_II)

chro <- kar_bs
tagment <- tagmentsbs
tagment <- tagment[order(tagment$start),]
rownames(tagment) <- NULL
tagment <- tagment[c(47:52, 62:64),]
tagment
mytrack <- rbind(bs_MHC_I, bs_MHC_II)
chromosome.index <- "bs_chr31"
head(mytrack)

dfbs <- read.table(file = "data/dfbs_MHC_anno.tsv", header = T, sep = "\t")
dfbs <- dfbs[order(dfbs$start),]
rownames(dfbs) <- NULL
dfbs$seqid <- gsub("HiC_scaffold_", "bs_chr", dfbs$seqid)
head(dfbs)
dfbs$chr <- paste(dfbs$seqid, dfbs$start, dfbs$end, sep = "-")
dfbs$chr == "bs_chr31-541642-562672" #Testing
tagment <- merge(tagment, dfbs, by.x = "tagments", by.y = "chr", all.x = T, all.y = F)
head(tagment)
col = viridis_pal(option = "D")(9)
col

f1 = function() {
  circos.par(gap.after = 2, start.degree = 90)
  circos.initializeWithIdeogram(chro, plotType = NULL, chromosome.index = chromosome.index)
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr=gsub(".*chr", "", CELL_META$sector.index)
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim),mean(ylim),chr,cex=0.6,col=col_text,
                facing="bending.inside",niceFacing=TRUE)
  },bg.col="grey90",bg.border=F,track.height=0.06)
}

f2 = function() {
  circos.par(cell.padding = c(0, 0, 0, 0), gap.after = c(rep(1, nrow(tagment)-1), 10)) #track.margin=c(0,0))
  circos.genomicInitialize(tagment, plotType = NULL)
  circos.genomicTrack(mytrack, ylim = c(-0.04, 1.1),
                        panel.fun = function(region, value, ...) {
                        for(h in seq(0, 1, by = 0.25)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = col)
                        }
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 3, col = col) 
                        circos.genomicLines(region, value, lwd = 0.1, lty = 1, area = T, 
                                           col = add_transparency("grey8", 0.4))
                        #circos.genomicPoints(region, value, 
                        #                     col = ifelse(value[[1]] > 0.51, "#E41A1C", "#377EB8"), 
                        #                     pch = 16, cex = 0.5)
                      }, bg.col = add_transparency(col, 0.5)
                      , track.margin = c(0.02, 0))
  
  # genomes x axis
  brk <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,20)*10^6
  brk
  col_text
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.axis(h="top",major.at=brk,labels=round(brk/10^6,1),labels.cex=0.4,
                col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
  },bg.border=F)
  circos.yaxis(side = "left", at = seq(0, 100, by = 0.3), 
               sector.index = get.all.sector.index()[1], labels.cex = 0.4)
  
  circos.track(ylim = c(0, 1), track.height = mm_h(2), 
               bg.col = add_transparency(col, 0.5))
  #Labels
  circos.genomicLabels(tagment,labels.column=13,cex=0.6,col=col,line_lwd=0.5,
                       line_col=col,
                      side="inside",connection_height=0.08,labels_height=0.025, font = 8, niceFacing = T, padding = 0.01)
  

}

circos.clear()
correspondance <- correspondance[correspondance$seqid == "bs_chr31",]
circos.nested(f1, f2, correspondance, 
              connection_col = "#add8e6",
              connection_border = "grey40" , adjust_start_degree = 90)

#graphics.off()
pdf("figures/bs_MHC_expresion.pdf", width = 6, height = 5)
circos.nested(f1, f2, correspondance, 
              connection_col = "#add8e6",
              connection_border = "grey40" , adjust_start_degree = 90)
dev.off()

####################################### MHC expression - Mute swan ###########################################################################
dfms <- dfbk[dfbk$seqid == "chr33",]
dfms <- dfms[order(dfms$start),]
head(dfms)
flnc.ms.bed <- read.table("data/ms_iso_to_genome.bed", header = F)
head(flnc.ms.bed)
flnc.ms.bed <- flnc.ms.bed[flnc.ms.bed == "chr33",]
head(flnc.ms.bed)
flnc.ms.bed <- genomicDensity(flnc.ms.bed, window.size = 500, count_by = "percent")
head(flnc.ms.bed)
isnps <- with(flnc.ms.bed, IRanges(start, end, names=chr, value=value))
head(tagments)
tagmentsms <- tagments[tagments$chr == "chr33",]
head(tagmentsms)
igenes <- with(tagmentsms, IRanges(start, end, names=tagments, chr=chr))
olaps <- findOverlaps(isnps, igenes)
ms <- cbind(flnc.ms.bed[queryHits(olaps),], tagmentsms[subjectHits(olaps),])
head(ms)
colnames(ms) <- c("chr","start1","end1","value","tagments","start","end","chr")
ms <- dplyr::select(ms, tagments, start1, end1, value)
head(ms)

tagments <- data.frame(tagments=paste(df2$seqid, df2$start, df2$end, sep = "-"), start=df2$start, end=df2$end, chr=df2$seqid)
head(tagments)
correspondance <- cbind(dplyr::select(df2, seqid, start,end), tagments)
head(correspondance)
correspondance <-correspondance[,-7]
tail(ms)

chr_bg_color = rand_color(2, transparency = 0.8)
names(chr_bg_color) = c("bs_chr31", "chr33")
kar
kar_ms <- "chr33"

#check the MHC locations - mute swans
#chr33 230455 236534(1083 - ) bs_chr31 622721 625612 (- 1252) - MHC I
#MHC class I 82     bs_chr31 695164 697066  bs_chr31 707341 712959
rownames(ms) <- NULL
ms[ms$tagments == "chr33-230455-236534",]
ms[ms$tagments == "chr33-250813-251768",]
ms[ms$tagments == "chr33-129930-131805",]
ms[ms$tagments == "chr33-142275-147824",]

ms_MHC_I <- ms[c(385:451),]
ms_MHC_II <- ms[c(237:280),]

head(ms_MHC_I)
head(ms_MHC_II)

chro <- kar[2,]
tagment <- tagmentsms
tagment <- tagment[order(tagment$start),]
rownames(tagment) <- NULL
head(tagment)
tagment <- tagment[c(21:25, 11:13),]
head(tagment)
mytrack <- rbind(ms_MHC_I, ms_MHC_II)
chromosome.index <- "chr33"
head(mytrack)

rownames(dfms) <- NULL
head(dfms)
dfms$chr <- paste(dfms$seqid, dfms$start, dfms$end, sep = "-")
dfms$chr == "chr33-250813-251768"
tagment <- merge(tagment, dfms, by.x = "tagments", by.y = "chr", all.x = T, all.y = F)
head(tagment)

col = viridis_pal(option = "D")(9)
chro
chromosome.index
mytrack <- mytrack[order(mytrack$start1),]
mytrack

f1 = function() {
  circos.par(gap.after = 2, start.degree = 90)
  circos.initializeWithIdeogram(chro, plotType = NULL, chromosome.index = chromosome.index)
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr=gsub("chr", "", CELL_META$sector.index)
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim),mean(ylim),chr,cex=0.6,col=col_text,
                facing="bending.inside",niceFacing=TRUE)
  },bg.col="grey90",bg.border=F,track.height=0.06)
}

f2 = function() {
  circos.par(cell.padding = c(0, 0, 0, 0), gap.after = c(rep(1, nrow(tagment)-1), 10)) #track.margin=c(0,0))
  circos.genomicInitialize(tagment, plotType = NULL)
  circos.genomicTrack(mytrack, ylim = c(-0.04, 1.1),
                      panel.fun = function(region, value, ...) {
                        for(h in seq(0, 1, by = 0.25)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "grey40")
                        }
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 3, col = "grey40") 
                        circos.genomicLines(region, value, lwd = 0.1, lty = 1, area = T, 
                                            col = add_transparency("grey8", 0.4))
                    }, bg.col = add_transparency(col, 0.5)
                      , track.margin = c(0.02, 0))
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.axis(h="top",major.at=brk,labels=round(brk/10^6,1),labels.cex=0.4,
                col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
  },bg.border=F)
  circos.yaxis(side = "left", at = seq(0, 100, by = 0.3), 
               sector.index = get.all.sector.index()[1], labels.cex = 0.4)
  
  circos.track(ylim = c(0, 1), track.height = mm_h(2), 
               bg.col = add_transparency(col, 0.5))
  #Labels
  circos.genomicLabels(tagment,labels.column=8,cex=0.6,col=col,line_lwd=0.5,
                       line_col=col,
                       side="inside",connection_height=0.08,labels_height=0.025, font = 8, niceFacing = T, padding = 0.01)
}

correspondance <- correspondance[correspondance$seqid == "chr33",]
circos.clear()
circos.nested(f1, f2, correspondance, 
              connection_col = "#add8e6",
              connection_border = "grey40" , adjust_start_degree = 90)

pdf("figures/ms_MHC_expresion.pdf", width = 6, height = 5)
circos.nested(f1, f2, correspondance, 
              connection_col = "#add8e6",
              connection_border = "grey40" , adjust_start_degree = 90)
dev.off()

graphics.off()
rm(list=all())
