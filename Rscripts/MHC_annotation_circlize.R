library(circlize)
library(dplyr)
library(grid)
library(ComplexHeatmap)
#setwd('C:/Users/uqakaraw/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/')
setwd('D:/05.OneDrive/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/')
#source('black_swan/MHC_class_annotation.R')
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

df <- dfmerged_MHC
rownames(df) <- NULL
write.table(df, file = "data/MHC_annotation.tsv", quote = F, row.names = F)
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


head(df)
df[97,]
df$end - df$start
df$end[97] <- 243906
df$start[97] <- 242858


circos.genomicTrack(df, ylim = c(0,0.3), track.height = 0.01, panel.fun = function(region, value, ...){
  circos.genomicRect(region, value, col = as.numeric(factor(value[[7]])), border = 'grey90', ybottom = 0.2, ytop = -0.4)
}, bg.border = F)
col_text
head(df)
tail(df)
ifelse(df$col == "#bd6959", "#bd6959", "#5d5cb8")
#Labels
circos.genomicLabels(df,labels.column=5,cex=0.6,col=as.numeric(factor(df$type)),line_lwd=0.5,
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







