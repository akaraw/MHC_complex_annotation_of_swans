require(rtracklayer)
require(circlize)
require(dplyr)
############################################################################################
gtfbs <- rtracklayer::readGFFAsGRanges('D:/05.OneDrive/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/data/bs_liftover.gff')
head(gtfbs)
gtfbs <- as.data.frame(gtfbs)
gtfbs <- gtfbs[gtfbs$seqnames == "HiC_scaffold_31",]

head(gtfbs)
tail(gtfbs)
gff <- gtfbs[gtfbs$type == "mRNA",]
head(gff)

df <- gff
dim(df)
df$col <- ifelse(gff$partial == "true" & !is.na(gff$partial), "blue", "blue")
df$col
df$col <- ifelse(df$valid_ORF == "False", "#5d5cb8", "#5d5cb8")
head(df)
df$Name <- ifelse(is.na(df$Name), df$ID, df$Name)
head(df)
df$Name
df1 <- df
df <- dplyr::select(df, 'seqnames', 'start', 'end', 'col', 'gene')
head(df)
rownames(df) <- NULL
df <- df[!duplicated(df$gene),]
dim(df)

dim(df <- distinct(df, gene, .keep_all = T))
dfbs <- df
dfcom <- read.table('D:/05.OneDrive/OneDrive - The University of Queensland/Black swan genome/CIRCLIZE/data/MHC_annotation_Swans.txt', sep = '\t', header = T)
dfcom$ID.MS <- trimws(dfcom$ID.MS)
dfcom
dfcom <- as.data.frame(apply(dfcom,2,function(x)gsub('\\s+', '',x)))
dfbs$gene
(dfbs[dfcom$ID.BS%in%dfbs$gene,])

dfmerged <- merge(dfbs, dfcom, by.x = "gene", by.y = "ID.BS", all.y = T)
head(dfmerged, 100)
dfmerged$gene
dfmerged$gene = ifelse(dfmerged$gene == "Not_found", dfmerged$Common_name, dfmerged$gene)


######################################################################################
gtfms <- rtracklayer::readGFFAsGRanges('C:/assemblyub/ms.gff')
head(gtfms)
gtfms <- as.data.frame(gtfms)
gtfms <- gtfms[gtfms$seqnames == "chr33",]

head(gtfms)
gff <- gtfms[gtfms$type == "mRNA",]
head(gff)
gff
df <- gff
df$col <- ifelse(gff$partial == "true" & !is.na(gff$partial), "red", "red")
df$col <- ifelse(df$transcript_id == "False", "#bd6959", "#bd6959")

head(df)
df$gene
df1 <- df
df <- dplyr::select(df, 'seqnames', 'start', 'end', 'col', 'gene', 'product')
head(df)
rownames(df) <- NULL

df <- df[!duplicated(df$gene),]
dim(df)
library(dplyr)
dim(df <- distinct(df, gene, .keep_all = T))
dfms <- df
dfms
dim(dfms)
dim(dfbs)

###########################################################################################
gtfbs <- rtracklayer::readGFFAsGRanges('C:/Users/ackar/OneDrive/Desktop/ms_exclude_partial.gff')
head(gtfbs)
gtfbs <- as.data.frame(gtfbs)
gtfbs <- gtfbs[gtfbs$seqnames == "HiC_scaffold_31",]

head(gtfbs)
gff <- gtfbs[gtfbs$type == "mRNA",]
head(gff)

df <- gff
df$col <- ifelse(gff$partial == "true" & !is.na(gff$partial), "blue", "blue")
df$col <- ifelse(df$valid_ORF == "False", "#bd6959", "#5d5cb8")

head(df)
df$gene
df1 <- df
df <- dplyr::select(df, 'seqnames', 'start', 'end', 'col', 'gene', 'product')
head(df)
rownames(df) <- NULL

df <- df[!duplicated(df$gene),]
dim(df)

dim(df <- distinct(df, gene, .keep_all = T))
dfbs <- df

dfbs$gene
dim(dfmerged[is.na(dfmerged$seqnames),])
dim(t <- dfbs[dfbs$gene %in%dfmerged$ID.MS[is.na(dfmerged$seqnames)],])
head(t)
head(dfmerged)
d <- merge(dfmerged, t, by.x = "ID.MS", by.y = "gene")
head(d)
r1 <- select(merge(dfmerged, t, by.x = "ID.MS", by.y = "gene"), seqnames.y, start.y, end.y, col.y, gene, product, Common_name)
dim(r1)
tail(dfmerged)
r2 <- select(dfmerged[!is.na(dfmerged$seqnames),], seqnames, start, end, col, gene, Common_name)
r1 <- select(r1, -"product")
head(r2)
head(r1)
colnames(r1) <- c("seqnames", "start", "end", "col", "gene", "Common_name")
dfbs = rbind(r1,r2)
dim(dfbs)
dfbs



dfcom
dfmerged
dfmerged$ID.MS 
dfms1 <- dfms[dfms$gene%in%dfcom$ID.MS,]
dim(dfms1)
head(dfms1)

(merge(dfms1, dfmerged, by.x = "gene", by.y = "ID.MS"))

dfms2 <- select(merge(dfms1, dfmerged, by.x = "gene", by.y = "ID.MS"), seqnames.x, start.x, end.x, col.x, gene, Common_name)
colnames(dfms2) <- c("seqnames", "start", "end", "col", "gene", "Common_name")
head(dfms2)
head(dfbs)
df <- rbind(dfbs, dfms2)
dim(df)
df
#######################################################################################################





