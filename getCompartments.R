library(GRanges)
library(pbapply)
library(mixOmics)

mat.bin <- binMatrix(mat=mat, gr.mat=gr.mat, chr="chr14", chr.max=XX)
cor.bin <- corrMatrix(mat.bin)
pc.bin  <- getABSignal(cor.bin)
mybarplot(pc.bin$pc)


# mat: p x n matrices where p represents loci and n cells
# gr.mat: granges object corresponding to the genomic locations of the p loci
# chr: which chromosome should be analyzed? 
# chr.min: what should be the starting pos of the chr?
# chr.max: what should be the endind pos of the chr?
# res: what should be the binning resolution?
# verbose: should the function be verbose?
# FUN: what function should be used to summarize information within a bin?
binMatrix <- function(mat, gr.mat, chr="chr14", chr.min = 0, chr.max,
    res = 100*1000, verbose = TRUE, FUN=sum){

    if (any(is.na(mat))){
        stop("Matrix must not contain NAs. ")
    }
    if (nrow(mat)!=length(gr.mat)){
        stop("Provided granges must have length equal to the matrix number of rows")
    }
    if (!verbose){
        pboptions(type="none")   
    }
    start <- seq(chr.min,chr.max,res)
    end <- c(start[-1],chr.max) -1L
    gr.bin <- GRanges(seqnames = chr,
                       ranges = IRanges(start = start, end = end))
    ids <- findOverlaps(gr.mat, gr.bin, select="first")
    n <- length(gr.bin)
    if (verbose) cat(paste0(n, " bins are created. \n"))

    mat.bin <- apply(mat, 2, function(x){
        temp <- rep(0, n)
        a <- tapply(x, INDEX=ids, FUN=FUN)
        temp[as.numeric(names(a))] <- a
        temp
    })
    colnames(mat.bin) <- colnames(mat)
    wh <- rowSums(mat.bin)!=0
    mat.bin <- mat.bin[wh,]
    gr.bin  <- gr.bin[wh]

    return(list(gr=gr.bin, mat=mat.bin))
}

# To be applied on an object returned by binMatrix
corrMatrix <- function(object){
    mat.cor <- cor(t(object$mat))
    gr.cor  <- object$gr
    return(list(gr=gr.cor, mat=mat.cor))
}


# To be applied on an object returned by corrMatrix
getABSignal <- function(object, k=2, iter=2){
    getFirstPC <- function (matrix, ncomp = 1){  
        matrix <- t(scale(t(matrix), center = TRUE, scale = FALSE))
        if (ncomp > 1) {
            a <- nipals(matrix, ncomp = ncomp)$p
        }
        else {
            a <- nipals(matrix, ncomp = ncomp)$p[, 1]
            sums <- colSums(matrix, na.rm=TRUE)
            if (cor(sums,a) <0){
                a <- -a
            }
        }
        a <- a*sqrt(length(a)) # Chr lenght normalization.
        return(a)
    }

    pc <- getFirstPC(object$mat)
    pc <- meanSmoother(pc, k=k, iter=iter)
    return(list(pc=pc, gr=object$gr))
}

# To visualize the AB signal
mybarplot <- function(x, main="",ylim=c(-1, 1), unitarize=FALSE, reverse=FALSE,
 top.col = "deeppink4", bot.col = "grey50"){
    if (unitarize){
      x <- unitarize(x)
    }
    x <- as.numeric(x)
    if (reverse){
        x <- -x
    }
    
    n <- length(x)
    col <- rep(top.col, n)
    col[x<0] <- bot.col
    barplot(x, ylim=ylim, 
        bty="n", xlab="", ylab="",border=col, col=col, main=main, yaxt='n')
}

