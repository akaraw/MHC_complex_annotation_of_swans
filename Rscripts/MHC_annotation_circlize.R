library(circlize)
library(dplyr)

#setwd('C:/Users/uqakaraw/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/')
setwd('D:/05.OneDrive/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/')
source('black_swan/MHC_class_annotation.R')
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

#flnc.bed <- read.table('mute_swan/coverage.flnctogenome.bed', sep = '\t', header = F)
#head(flnc.bed)
#flnc.bed[,1] = gsub('HiC_scaffold_', 'bs_chr', flnc.bed[,1])
#head(flnc.bed)
#tail(flnc.bed)
#colnames(flnc.bed) <- c('chr', 'start', 'end', 'value')
#flnc.bed<- flnc.bed[flnc.bed$chr%in%kar$chr,]
#flnc.bed<- flnc.bed[,1:4]
#flnc.bed <- genomicDensity(flnc.bed, window.size = 1000)
#head(flnc.bed)
#dim(flnc.bed)
#flnc.bed

##### A track for the genome GC ####
getwd()
write.table(df, file = "MHC_annotation.tsv", quote = F, row.names = F)
graphics.off()
pdf('circlize_ms_bs_mhc.pdf', width = 6, height = 5)
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

#FLNC reads mapping
#circos.track(factors=flnc.bed$chr,x=flnc.bed$end,y=flnc.bed$value,panel.fun=function(x,y) {
#  circos.lines(x,y,col="#F07575",lwd=0.4, type = 'l', area = T, baseline = "bottom", border = "#F07575")
#  circos.yaxis(side = 'left', at = brk, sector.index = "chr33", lwd = 0.05, labels = T, labels.cex = 0.4)
#},ylim=range(flnc.bed$value),track.height=0.1,bg.border=T)

#Genes
circos.genomicTrack(df, ylim = c(0,0.3), track.height = 0.01, panel.fun = function(region, value, ...){
  circos.genomicRect(region, value, col = value[[1]], border = 'grey90', ybottom = 0.2, ytop = -0.4)
}, bg.border = F)
col_text
head(df)
tail(df)
ifelse(df$col == "#bd6959", "#bd6959", "#5d5cb8")
#Labels
circos.genomicLabels(df,labels.column=6,cex=0.3,col=ifelse(df$col == "#bd6959", "#bd6959", "#5d5cb8"),line_lwd=0.5,
                     line_col=ifelse(df$col == "#bd6959", "#bd6959", "#5d5cb8"),
                     side="inside",connection_height=0.08,labels_height=0.03, font = 8, niceFacing = T, padding = 0.01)


circos.clear()




