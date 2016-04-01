## Pedigree manipulation/simulation functions. 

library(kinship2)
library(RColorBrewer)

############################################################################
##
## Check ped
##
## Check that the pedigree is well-specified for this analysis.
## In particular, each founder must have exactly one child.
## 
############################################################################

check.ped <- function(ped){
    ## kinship2 checks
    pedAll <- pedigree(id=ped$IND, dadid=ped$PAT, momid=ped$MAT, sex=ped$SEX)
    
    founders<-ped$IND[ped$MAT==0 & ped$PAT==0]
    children.m <- ped[ped$MAT %in% founders,"IND"]
    children.p <- ped[ped$PAT %in% founders,"IND"]
    n.m <- table(children.m)
    n.p <- table(children.p)
    if(!all(n.m==1)){
        stop(paste0("Founders ", names(n.m)[n.m!=1], "do not have exactly one child"))
    }
    if(!all(n.p==1)){
        stop(paste0("Founders ", names(n.p)[n.p!=1], "do not have exactly one child"))
    }
}

############################################################################
##
## Simulate the whole (autosomal) genome
## 
############################################################################

simulate.genome <- function(ped, chr.lens){
    maps <- vector("list", length(chr.lens))
    for(i in 1:length(chr.lens)){
        maps[[i]] <- simulate.chromosome(ped, chr.lens[i])
    }
    return(maps)
}

############################################################################
##
## Simulate the two chromosomes of individual 1
##
## ped: pedigree (data frame with IND, MAT, PAT columns)
## len: chromosome length in mb
##
############################################################################

simulate.chromosome <- function(ped, len){

    ## parents of individual 1
    parents <- ped[ped$IND==1,c("MAT", "PAT")]
    ## Generate recombination
    recs <- generate.recombinations(ped, len, 1)
    ## Get the maternal and paternal chromosome maps
    mat.map <- decode.vectors(ped, recs, parents$MAT, 1)[,1:2]
    pat.map <- decode.vectors(ped, recs, parents$PAT, 1)[,1:2]
    colnames(mat.map) <- c("POS", "MAT")
    colnames(pat.map) <- c("POS", "PAT")
    mat.map$PAT <- NA
    pat.map$MAT <- NA
    
    map <- rbind(mat.map, pat.map)
    map <- map[order(map$POS),]
    for(i in 2:NROW(map)){
        for(what in c("MAT", "PAT")){
            if(is.na(map[i,what])){
                map[i,what] <- map[i-1,what]
            }
        }
    }
    map <- map[2:NROW(map),]            #First row is duplicated 0 and superfluous
    map <- map[!duplicated(map$POS),]
    return(map)
}

############################################################################
##
## Decode vectors
## from a list of recombinations, recursively get a vector of the chromosome
## painted by founder.
## 
############################################################################

decode.vectors <- function(ped, recs, ind, trans, vec=data.frame(POS=0, FND=ind, TRANS=trans)){
    nonfounders <- ped$IND[ped$MAT!=0 & ped$PAT!=0]
    parents <- ped[ped$IND==ind,c("MAT", "PAT")]

    current.paints <- unique(vec$FND)

    chr.segments <- unique(vec[vec$FND %in% nonfounders,c("FND","TRANS")])

    ## Until all the painted segments are founders
    while(NROW(chr.segments)){    
        for(i in 1:NROW(chr.segments)){
            ind <- chr.segments[i,"FND"]
            trans <- chr.segments[i,"TRANS"]
        
            ## Set the current chunk to one of its parents
            parents <- ped[ped$IND==ind,c("MAT", "PAT")]
            current <- ifelse(runif(1)<0.5, parents$MAT, parents$PAT)
            vec[vec$FND==ind & vec$TRANS==trans,c(2,3)] <- c(current, ind)

            ## Every time we hit a recombination in the current gamete (ind->trans), switch parent. 
            rec.pos <- recs[recs$REC==ind & recs$TRANS==trans,"pos"]
            for(rp in rec.pos){
                other <- parents[parents!=current]
                ## Switch all the chromosomes to the right
                vec[vec$FND==current & vec$TRANS==ind & vec$pos>rp,2] <- other

                ## Add an additional break, if the recombination happens in an ancestral segment
                segment.i <- max(which(vec$POS<rp)) #This is the ancestral segment at the recomb
                if(all(vec[segment.i,c(2,3)]==c(current,ind))){
                    vec <- rbind(vec, c(rp, other, ind))
                }
                current <- other
            }
            vec <- vec[order(vec$POS),]
            
            ## Merge adjacent regions with the same transmission
            ## are you the same as the row above?
            if(NROW(vec)>1){
                same <- c(FALSE, apply(vec[2:NROW(vec),c(2,3)]== vec[1:(NROW(vec)-1),c(2,3)], 1, all))
                vec <- vec[!same,]
            }
        }
        chr.segments <- unique(vec[vec$FND %in% nonfounders,c(2,3)])
    }

    return(vec)
}

############################################################################
##
## recursively generate recombinations in the history of individual ind
##
## Output:
## REC: parent
## TRANS: child
## pos: position of recombination
## 
############################################################################

generate.recombinations <- function(ped, len, ind){
    nonfounders <- ped$IND[ped$MAT!=0 & ped$PAT!=0]

    parents <- ped[ped$IND==ind,c("MAT", "PAT")]
    mum <- parents$MAT
    dad <- parents$PAT
    ## maternal
    if(ind %in% nonfounders){
        ## maternal
        mrate <- 0.5+0.0104*len
        ## mrate <- 0.013*len
        ## mrate <- 1
        nmrec <- rpois(1, mrate)
        pos <- round(runif(nmrec)*len,6)
        mrecs <- data.frame(REC=rep(mum,nmrec), TRANS=rep(ind,nmrec), pos=pos)
        ## recurse
        amrecs <- generate.recombinations(ped, len, mum)

        ## paternal
        prate <- 0.5+0.0032*len
        ## prate <- 0.013*len
        ## prate <- 1
        pmrec <- rpois(1, prate)
        pos <- round(runif(pmrec)*len,6)
        precs <- data.frame(REC=rep(dad,pmrec), TRANS=rep(ind,pmrec), pos=pos)
        ## recurse
        aprecs <- generate.recombinations(ped, len, dad)
    }else{
        mrecs <- amrecs <- data.frame(REC=NULL, TRANS=NULL, pos=NULL)
        precs <- aprecs <- data.frame(REC=NULL, TRANS=NULL, pos=NULL)
    }

    all.recs <- rbind(mrecs, precs, amrecs, aprecs)
    return(all.recs)                      
}

############################################################################
##
## plot the simulated chromsomes
##
############################################################################

plot.chromosome <- function(ped, map, len, cols=NA){
    if(all(is.na(cols))){
        founders<-ped$IND[ped$MAT==0 & ped$PAT==0]
        cols <- brewer.pal(length(founders), "Set1")
        names(cols) <- paste0("FND", founders)
    }

    map$POS <- map$POS/len
    map$MAT <- paste0("FND", map$MAT)
    map$PAT <- paste0("FND", map$PAT)
    map <- rbind(map, c(1,NA, NA))
    plot(c(0,2.75), c(0,1), col="white", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
    for(i in 1:(NROW(map)-1)){
        rect(0,map[i,"POS"],1,map[i+1,"POS"], col=cols[map[i,"MAT"]], border=NA)
        rect(1.5,map[i,"POS"],2.5,map[i+1,"POS"], col=cols[map[i,"PAT"]], border=NA)
    }
    rect(0,0,1,1, border="black")
    rect(1.5,0,2.5,1, border="black")
    
}

############################################################################
##
## plot the simulated genome
##
############################################################################

plot.genome <- function(ped, maps, chr.lens, cols=NA){
    if(all(is.na(cols))){
        founders<-ped$IND[ped$MAT==0 & ped$PAT==0]
        cols <- brewer.pal(length(founders), "Set1")
        names(cols) <- paste0("FND", founders)
    }

    nchr <- length(maps)
    plot(c(0,nchr), c(0,max(chr.lens)), col="white", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
    for(i in 1:nchr){
        map <- maps[[i]]
        map <- rbind(map, c(chr.lens[i],NA, NA))
        map$MAT <- paste0("FND", map$MAT)
        map$PAT <- paste0("FND", map$PAT)
        x1 <- i-1
        x2 <- i-0.7
        wd <- 0.3
        for(j in 1:(NROW(map)-1)){
            rect(x1,map[j,"POS"],x1+wd,map[j+1,"POS"], col=cols[map[j,"MAT"]], border=NA)
            rect(x2,map[j,"POS"],x2+wd,map[j+1,"POS"], col=cols[map[j,"PAT"]], border=NA)
        }
        rect(x1,0,x1+wd,chr.lens[i], border="black")
        rect(x2,0,x2+wd,chr.lens[i], border="black")
    }
}

############################################################################
##
## plot the pedigree
##
############################################################################

plot.ped <- function(ped,  cols=NA, color.founders=FALSE){
    if(all(is.na(cols)) & color.founders){
        founders<-ped$IND[ped$MAT==0 & ped$PAT==0]
        cols <- brewer.pal(length(founders), "Set1")
        names(cols) <- paste0("FND", founders)
    }

    pedAll <- pedigree(id=ped$IND, dadid=ped$PAT, momid=ped$MAT, sex=ped$SEX)
    cc <- cols[paste0("FND", pedAll$id)]
    plot(pedAll, col=ifelse(is.na(cc), "black", cc), affected=ifelse(is.na(cc), 0, 1))
    
}

############################################################################
##
## get autozygous (i.e. IBD) chunks from map
##
############################################################################

ibd.chunks <- function(map, len){
    if(NROW(map)==1){
        map$END <- len
    }else{
        map$END <- c(map$POS[2:NROW(map)], len)
    }
    map$LEN <- map$END-map$POS
    ibd.map <- map[map$PAT==map$MAT,]
}

############################################################################
##
## get autozygous (i.e. IBD) chunks from list of maps
##
############################################################################

ibd.chunks.genome <- function(maps, chr.lens){
    nchr=length(maps)
    chunks <- ibd.chunks(maps[[1]], chr.lens[1])
    if(nchr>1){
        for(i in 2:nchr){
            chunks <- rbind(chunks, ibd.chunks(maps[[i]], chr.lens[i]))
        }
    }
    return(chunks)
}
