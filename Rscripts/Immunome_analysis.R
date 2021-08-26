library(PANTHER.db)
library(ggplot2)
library(topGO)
#reading files

setwd("D:/05.OneDrive/OneDrive - The University of Queensland/Orthofinder_final_reslts/")

list <- list.files("PANTHER/", pattern = ".list")
list
all(exists(list))
#list <- list[-c(16,17)]
df <- read.table("PANTHER/data/panther_sf_acc.list", header = F)
head(df)
dim(df)
length(list)
list[1]

for (i in 1:length(list)) {
  x <- read.table(paste0("PANTHER/", list[i]), header = F)
  df <- merge(df, x, all.x = T, by = "V1", )
}

head(df)
dim(df)
names <- gsub(".list", "", list)
names
colnames(df) <- c("Orthogroup", names)
colnames(df)
tail(df)

df[is.na(df)] <- 0
class(df$Anas_platyrhynchos)
head(df)
df <- df[rowSums(df[2:18]) != 0,]
dim(df)
null <- data.frame(Desc=rep("(null)", 23378))
head(null)
df <- cbind(null, df)
mdf <- df

write.table(mdf, file = "PANTHER/data/panther_sf_to_cafe.txt", quote = F, row.names = F, sep = "\t")
head(mdf)
dim(mdf)
####PANTHER input file for CAFE5 completed########

###Gene families of Black swans that have lesser members than ducks and mute swans
rownames(mdf) <- mdf$Orthogroup
(rownames(mdf)[ mdf$Cygnus_atratus < mdf$Cygnus_olor & mdf$Cygnus_atratus < mdf$Anas_platyrhynchos])
###Common Gene families that have members in both black swan and chicken lesser than both mute swan and duck
(chi_bs_common <- rownames(mdf)[mdf$Cygnus_atratus < mdf$Cygnus_olor & mdf$Cygnus_atratus < mdf$Anas_platyrhynchos & mdf$Gallus_gallus < mdf$Anas_platyrhynchos & mdf$Gallus_gallus < mdf$Cygnus_olor])
l = length(chi_bs_common)
l
ids <- data.frame(IDs=paste0("myid", c(1:l)))
ids
chi_bs_common <- cbind(ids, chi_bs_common)
write.table(chi_bs_common, file = "PANTHER/data/chi_bscommon_sf.txt", quote = F, row.names = F, sep = "\t")
mdf[mdf$Orthogroup == "PTHR24100:SF116",]
mdf[mdf$Orthogroup == "PTHR12080:SF55",]

##########Black swan fast evolving gene families ##################################################
cygA <- read.table("PANTHER/BS_data/Cygg_at_significantfam.list")
dim(cygA)

ids <- data.frame(ids=paste0("my_id", c(1:294)))
ids
cygA <- cbind(ids, cygA)
colnames(cygA) <- c("V1", "V2")
colnames(cygA)
cygA
write.table(cygA, file = "PANTHER/BS_data/Cygnus_atratus_fast_evolving_upload.txt", quote = F, col.names = F, row.names = F, sep = "\t")
grp <- as.data.frame(mdf$Orthogroup)
grp
write.table(grp, file = "PANTHER/BS_data/Background_list.txt", quote = F, col.names = F, row.names = F, sep = "\t")

##########Mute swan fast evolving gene families ##################################################
cygO <- read.table("PANTHER/BS_data/cygol_only_significantfam.list")
dim(cygO)

ids <- data.frame(ids=paste0("my_id", c(1:143)))
ids
cygO <- cbind(ids, cygO)
colnames(cygO) <- c("V1", "V2")
colnames(cygO)
cygO
write.table(cygO, file = "PANTHER/BS_data/cygnus_olor_fast_evol_upload.txt", quote = F, col.names = F, row.names = F, sep = "\t")

##########Duck fast evolving gene families ##################################################
duck <- read.table("PANTHER/BS_data/anaplat_only_significantfam.list")
dim(duck)
duck
ids <- data.frame(ids=paste0("my_id", c(1:141)))
ids
duck <- cbind(ids, duck)
colnames(duck) <- c("V1", "V2")
colnames(duck)
duck
write.table(duck, file = "PANTHER/data/duck_fast_evol_upload.txt", quote = F, col.names = F, row.names = F, sep = "\t")

##########Chick fast evolving gene families ##################################################
chick <- read.table("PANTHER/data/chick_fast_evolving_fam.list")
dim(chick)
chick
ids <- data.frame(ids=paste0("my_id", c(1:159)))
ids
chick <- cbind(ids, chick)
colnames(chick) <- c("V1", "V2")
colnames(chick)
chick
write.table(chick, file = "PANTHER/data/chick_fast_evol_upload.txt", quote = F, col.names = F, row.names = F, sep = "\t")

######Immunome analysis############################################################################
immunome <- read.table("PANTHER/data/immunome.list", sep = "\t", header = F)
head(immunome)
dim(immunome)
dim(mdf)
immunome_spp <- mdf[mdf$Orthogroup%in%immunome$V1,]
dim(immunome)
dim(immunome_spp)
immunome_spp_sp <- dplyr::select(immunome_spp, -c("Orthogroup", "Desc"))
colSums(immunome_spp_sp)

bs_immune <- immunome_spp[immunome_spp$Cygnus_atratus < immunome_spp$Cygnus_olor & immunome_spp$Cygnus_atratus < immunome_spp$Anas_platyrhynchos,]
bs_immune


#######################GO slim category mapping####################################################
keytypes(PANTHER.db)
columns(PANTHER.db)

#################################################TOPGO analysis#################################################
fam_ids <- cygA$V2
fam_ids
cols <- "GOSLIM_ID"
de.fams <- mapIds(PANTHER.db, keys=fam_ids, column=cols, keytype="FAMILY_ID", multiVals="list")
de.fams <- de.fams[!is.na(names(de.fams))]
lengths(de.fams)

#Creating the gene vector
assayed.fam <- read.table("PANTHER/BS_data/Base_count.tab", header = T, sep = "\t")
head(assayed.fam)
colnames(assayed.fam)
assayed.fam <- dplyr::select(assayed.fam, "FamilyID", "Cygnus_atratus.0.", "Cygnus_olor.1.", "X.2.")
head(assayed.fam)
assayed.fam <- assayed.fam[assayed.fam$X.2. != 0,]
dim(assayed.fam)


fam_ids
source("PANTHER/rscripts/mytopgo.R")
mytopgo(fam_ids = fam_ids, assayed.fam = assayed.fam, name = "all_bs")
mytopgo(fam_ids = cygO$V2, assayed.fam = assayed.fam, name = "all_ms")
mytopgo(fam_ids = duck$V2, assayed.fam = assayed.fam, name = "all_duck")
mytopgo(fam_ids = chick$V2, assayed.fam = assayed.fam, name = "all_chick")

#Expand Black swan
assayed.fam <- read.table("PANTHER/BS_data/Base_count.tab", header = T, sep = "\t")
head(assayed.fam)
colnames(assayed.fam)
assayed.fam <- dplyr::select(assayed.fam, "FamilyID", "Cygnus_atratus.0.", "Cygnus_olor.1.", "X.2.")
head(assayed.fam)
assayed.fam <- assayed.fam[assayed.fam$X.2. != 0,]
fam_ids <- cygA$V2
expand <- assayed.fam[assayed.fam$Cygnus_atratus.0. > assayed.fam$X.2.,]
dim(expand)
expand <- expand[expand$FamilyID%in%fam_ids,]
mytopgo(fam_ids = expand$FamilyID, assayed.fam = assayed.fam, name = "bs_expand")

#Contract Black swan
contract <- assayed.fam[assayed.fam$Cygnus_atratus.0. < assayed.fam$X.2.,]
dim(contract)
contract <- contract[contract$FamilyID%in%fam_ids,]
mytopgo(fam_ids = contract$FamilyID, assayed.fam = assayed.fam, name = "bs_contract")

#Expand Mute swan
fam_ids <- cygO$V2
expand <- assayed.fam[assayed.fam$Cygnus_olor.1. > assayed.fam$X.2.,]
dim(expand)
expand <- expand[expand$FamilyID%in%fam_ids,]
mytopgo(fam_ids = expand$FamilyID, assayed.fam = assayed.fam, name = "ms_expand")

#Contract Mute swan
contract <- assayed.fam[assayed.fam$Cygnus_olor.1. < assayed.fam$X.2.,]
dim(contract)
contract <- contract[contract$FamilyID%in%fam_ids,]
mytopgo(fam_ids = contract$FamilyID, assayed.fam = assayed.fam, name = "ms_contract")



#####################topGO with Galloanserinae####################################################
#Creating the gene vector
assayed.fam <- read.table("PANTHER/BS_data/Base_count.tab", header = T, sep = "\t")
head(assayed.fam)
colnames(assayed.fam)
assayed.fam <- dplyr::select(assayed.fam, "FamilyID", "Anas_platyrhynchos.3.", "X.10.")
head(assayed.fam)
assayed.fam <- assayed.fam[assayed.fam$X.10. != 0,]
dim(assayed.fam)

#Expand = duck
fam_ids <- duck$V2
expand <- assayed.fam[assayed.fam$Anas_platyrhynchos.3. > assayed.fam$X.10.,]
dim(expand)
expand <- expand[expand$FamilyID%in%fam_ids,]
mytopgo(fam_ids = expand$FamilyID, assayed.fam = assayed.fam,name = "duck_expand")

#Contract - duck
contract <- assayed.fam[assayed.fam$Anas_platyrhynchos.3. < assayed.fam$X.10.,]
dim(contract)
contract <- contract[contract$FamilyID%in%fam_ids,]
mytopgo(fam_ids = contract$FamilyID, assayed.fam = assayed.fam, name = "duck_contract")

#Creating the gene vector
assayed.fam <- read.table("PANTHER/BS_data/Base_count.tab", header = T, sep = "\t")
head(assayed.fam)
colnames(assayed.fam)
assayed.fam <- dplyr::select(assayed.fam, "FamilyID", "Gallus_gallus.5.", "X.11.")
head(assayed.fam)
assayed.fam <- assayed.fam[assayed.fam$X.11. != 0,]
dim(assayed.fam)

#Expand = chicken
fam_ids <- chick$V2
expand <- assayed.fam[assayed.fam$Gallus_gallus.5. > assayed.fam$X.11.,]
dim(expand)
expand <- expand[expand$FamilyID%in%fam_ids,]
mytopgo(fam_ids = expand$FamilyID, assayed.fam = assayed.fam, name = "chick_expand")

#Contract - chicken
contract <- assayed.fam[assayed.fam$Gallus_gallus.5. < assayed.fam$X.11.,]
dim(contract)
contract <- contract[contract$FamilyID%in%fam_ids,]
mytopgo(fam_ids = contract$FamilyID, assayed.fam = assayed.fam, name = "chick_contract")
