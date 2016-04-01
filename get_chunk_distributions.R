## Simulate the autozygous chunk distirbutions under different pedigrees
source("~/pedigree/code/pedigree_lib.R")
library(ggplot2)
library(RColorBrewer)

set.seed(12345)

pedigrees <- c(  "half_sibs", "uncle_niece", "double_first_cousins", "grandfather_granddaughter")
nsims <- 1000

chrs <- read.table("~/pedigree/code/chromosome_lengths.txt", as.is=T)

results <- data.frame(ped=NULL, value=NULL, stringsAsFactors=FALSE)
for(what in pedigrees){
    cat(paste0("\r", what))
    ped <- read.table(paste0("~/pedigree/code/peds/",what,".txt"), as.is=T, header=T)
    check.ped(ped)
    for(i in 1:nsims){
        cat(paste0("\r", what, "\t", i))
        maps <- simulate.genome(ped, chrs[,2])
        chunks <- ibd.chunks.genome(maps, chrs[,2])
        results <- rbind(results, data.frame(ped=rep(what, NROW(chunks)), value=chunks$LEN, stringsAsFactors=FALSE))
    }
}
cat("\n")

cols <- brewer.pal(length(pedigrees), "Set1")
names(cols) <- pedigrees

pdf("Length_distirbution.pdf", height=6, width=9)
p <- ggplot(results, aes(value, fill = ped, col=ped)) + geom_density(alpha=0.5)+scale_color_manual(values=cols)+scale_fill_manual(values=cols)+xlab("Autozygous chunk length (Mb)")
print(p)
dev.off()

average.total.length <- aggregate(results$value, by=list(results$ped), FUN=sum)
average.total.length$x <- average.total.length$x/sum(chrs[,2])/nsims
average.chunk.count <- table(results$ped)/nsims

results.10mb <- results[results$value>10,]
average.10mb.total.length <- aggregate(results.10mb$value, by=list(results.10mb$ped), FUN=sum)
average.10mb.total.length$x <- average.10mb.total.length$x/sum(chrs[,2])/nsims
average.10mb.chunk.count <- table(results.10mb$ped)/nsims

cutoff <- 10

data<-scan("~/pedigree/Altai_auto_chunks.txt")
data <- data[data>cutoff]

res.10<-results[results$value>cutoff,]
for(what in pedigrees){
    if(all(which(pedigrees==what)==1)){
        qqplot(wt, data, col=cols[what], type="l", bty="n", xlab="Model", ylab="Observed", xlim=c(cutoff,100))
    }else{
       pt <- qqplot( wt, data,plot.it=FALSE)
       lines(pt, col=cols[what])
   }
    wt <- res.10[res.10$ped==what,"value"]
    ks <- ks.test(data, wt)
    cat(paste0(what, " ", ks$p.value, "\n"))
}
abline(0,1, lwd=2, lty=2)
legend("bottomright", pedigrees, col=cols[pedigrees], lwd=2, bty="n")
