## Simulate some pedigrees and example chromsomes
source("~/pedigree/code/pedigree_lib.R")
library("RColorBrewer")

## what <- "uncle_niece"
## what <- "half_sibs"
## what <- "grandfather_granddaughter"
## what <- "double_parallel_cousins"
what <- "double_cross_cousins"
## what <- "grandmother_grandson"
## what <- "aunt_nephew"
## what <- "first_cousins"

ped <- read.table(paste0("~/pedigree/code/peds/",what,".txt"), as.is=T, header=T)
chrs <- read.table("~/pedigree/code/chromosome_lengths.txt", as.is=T)

maps <- simulate.genome(ped, chrs[,2])

founders<-ped$IND[ped$MAT==0 & ped$PAT==0]
cols <- brewer.pal(length(founders), "Set1")
names(cols) <- paste0("FND", founders)
pdf(paste0("example_", what, ".pdf", 6, 12))
par(mfrow=c(1,2))
plot.ped(ped, cols=cols)
plot.genome(ped, maps, chrs[,2], cols=cols)
ibd.chunks.genome(maps, chrs[,2])
dev.off()
