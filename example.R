## Simulate some pedigrees and example chromsomes
## source("~/pedigree/code/libs/pedigree_lib.R")
source("~/pedigree/code/pedigree_lib.R")

## par(mfrow=c(1,2))
## ped <- read.table("~/pedigree/code/peds/uncle_neice.txt", as.is=T, header=T)
## founders<-ped$IND[ped$MAT==0 & ped$PAT==0]
## cols <- brewer.pal(length(founders), "Set1")
## names(cols) <- paste0("FND", founders)
## plot.ped(ped, cols=cols)
## map <- simualate.genome(ped, len)
## plot.chromosome(ped, map, len)

ped <- read.table("~/pedigree/code/peds/uncle_neice.txt", as.is=T, header=T)
chrs <- read.table("~/pedigree/code/chromosome_lengths.txt", as.is=T)

maps <- simulate.genome(ped, chrs[,2])

founders<-ped$IND[ped$MAT==0 & ped$PAT==0]
cols <- brewer.pal(length(founders), "Set1")
names(cols) <- paste0("FND", founders)
par(mfrow=c(1,2))
plot.ped(ped, cols=cols)
plot.genome(ped, maps, chr.lens, cols=cols)
