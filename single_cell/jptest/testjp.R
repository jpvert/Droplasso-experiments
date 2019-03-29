# Test droplasso on single-cell RNA-seq from conquer
# Compare glmnet (lasso, ridge, elastic net) to droplasso (dropout, droplasso)
# JP Vert
# March 23, 2019

#### Initialization ###

# Load packages
library(MultiAssayExperiment)
library(parallel)
library(glmnet)
library(droplasso)
library(ROCR)
library(xtable)

# Random number initialization
set.seed(1234)

# Number of cores to use
numCores <- detectCores()-1

# Data to extract from conquer
conquerdata = list(
  list(
    name = "EMTAB2805",
    classid = "cell_cycle_stage",
    class = c("G1",
              "G2M"),
    displayname = "EMTAB2805, G1 vs G2M"
  ),
  list(
    name = "EMTAB2805",
    classid = "cell_cycle_stage",
    class = c("G1",
              "S"),
    displayname = "EMTAB2805, G1 vs S"
  ),
  list(
    name = "EMTAB2805",
    classid = "cell_cycle_stage",
    class = c("S",
              "G2M"),
    displayname = "EMTAB2805, S vs G2M"
  ),
  list(
    name = "GSE45719",
    classid = "source_name_ch1",
    class = c("16-cell stage blastomere",
              "Mid blastocyst cell (92-94h post-fertilization)"),
    displayname = "GSE45719, 16-cell vs Mid blastocyst"
  ),
  list(
    name = "GSE45719",
    classid = "source_name_ch1",
    class = c("16-cell stage blastomere",
              "8-cell stage blastomere"),
    displayname = "GSE45719, 16-cell vs 8-cell"
  ),
  list(
    name = "GSE48968-GPL13112",
    classid = "source_name_ch1",
    class = c("BMDC (1h LPS Stimulation)",
              "BMDC (4h LPS Stimulation)"),
    displayname = "GSE48968, 1h vs 4h"
  ),
  list(
    name = "GSE48968-GPL13112",
    classid = "source_name_ch1",
    class = c("BMDC (4h LPS Stimulation)",
              "BMDC (6h LPS Stimulation)"),
    displayname = "GSE48968, 4h vs 6h"
  ),
  list(
    name = "GSE74596",
    classid = "source_name_ch1",
    class = c("Single_cell_RNA-seq_NKT0",
              "Single_cell_RNA-seq_NKT17"),
    displayname = "GSE74596, NKT0 vs NKT17"
  ),
  list(
    name = "GSE74596",
    classid = "source_name_ch1",
    class = c("Single_cell_RNA-seq_NKT0",
              "Single_cell_RNA-seq_NKT1"),
    displayname = "GSE74596, NKT0 vs NKT1"
  ),
  list(
    name = "GSE74596",
    classid = "source_name_ch1",
    class = c("Single_cell_RNA-seq_NKT1",
              "Single_cell_RNA-seq_NKT2"),
    displayname = "GSE74596, NKT1 vs NKT2"
  )
)

# Fraction of cells that need to have nonzero counts to keep a gene
keepgenethreshold = 0.1

# Threshold to count non-zero weights in final models
epsilon <- 1e-8

# Number of repeats and folds for cross-validation
nrepeats <- 2
nfolds <- 5

# Number of lambda values to test (regularization parameter for glmnet and droplasso)
nlambda <- 20


# Methods to test
classifmethods <- list(
  list(displayname="LASSO",
       fcall="glmnet",
       fargs=list(alpha=1, intercept = F,standardize = F, nlambda=nlambda)),
  list(displayname="Elastic Net",
       fcall="glmnet",
       fargs=list(alpha=0.5, intercept = F,standardize = F, nlambda=nlambda)),
  list(displayname="Ridge",
       fcall="glmnet",
       fargs=list(alpha=0, intercept = F,standardize = F, nlambda=nlambda)),
  list(displayname="Dropout 0.5",
       fcall="droplasso",
       fargs=list(keep_prob=0.5, lambda=0)),
  list(displayname="Dropout 0.1",
       fcall="droplasso",
       fargs=list(keep_prob=0.1, lambda=0)),
  list(displayname="Droplasso 0.5",
       fcall="droplasso",
       fargs=list(keep_prob=0.5, nlambda=nlambda)),
  list(displayname="Droplasso 0.1",
       fcall="droplasso",
       fargs=list(keep_prob=0.1, nlambda=nlambda))
)


### Main loop to run experiments ###
result = list()
idata = 1

for (dataset in conquerdata) {
  print(
    paste(
      "Now working on dataset ",
      dataset$name,
      " with classes ",
      dataset$class[1],
      " and ",
      dataset$class[2],
      sep = ""
    )
  )
  
  # Read data from conquer
  print("Loading data")
  d = readRDS(gzcon(url(
    paste(
      "http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/",
      dataset$name,
      ".rds",
      sep = ""
    )
  )))
  
  # Extract count matrix
  cts = assays(experiments(d)[["gene"]])[["count_lstpm"]]
  print(paste(ncol(cts), " cells, ", nrow(cts), " genes.", sep = ""))
  
  # Extract positive and negative cells
  pdata <- colData(d)
  pos = which(pdata[[dataset$classid]] == dataset$class[1])
  neg = which(pdata[[dataset$classid]] == dataset$class[2])
  npos = length(pos)
  nneg = length(neg)
  n = npos + nneg
  x <- t(cts)[c(pos, neg),]
  # Labels
  y = c(rep(0, npos), rep(1, nneg))
  # Remove genes with positive counts in enough cells
  nposcell <- apply(x, 2, function(i) {
    sum(i > 0)
  })
  x = x[, nposcell >= n * keepgenethreshold]
  # Take log(count+1)
  x = log(x + 1)
  print(paste(
    npos,
    " positives, ",
    nneg,
    " negatives, ",
    ncol(x),
    " genes kept.",
    sep = ""
  ))
  
  # Make folds
  folds <- list()
  for (i in seq(nrepeats)) {
    folds <- c(folds, split(sample(seq(n)), rep(1:nfolds, length = n)))
  }
  
  # Main loop over methods
  resultloc <-
    parallel::mclapply(seq(length(folds)), function(ifold) {
      cat('.')
      itrain <- folds[[ifold]]
      ivalid <- folds[[ifold + 1 - nfolds*(ifold%%nfolds==0)]]
      itest <- seq(n)[-c(itrain,ivalid)]
      
      auc <- list()
      nfeat <- list()
      
      sharedargs <- list(x=x[itrain,], y=y[itrain], family="binomial")
      for (met in classifmethods) {
        # Train model on train set
        m <- do.call(met$fcall,c(sharedargs,met$fargs))
        # Parameter selection on the validation set
        yval <- predict(m, x[ivalid,])
        aucval <- c()
        for (ipar in seq(ncol(yval))) {
          pred <- prediction(yval[,ipar], y[ivalid])
          aucval[ipar] <- performance(pred, "auc")@y.values[[1]]
        }
        bestparam <- m$lambda[which.max(aucval)]
        # Performance on the test set
        ytest <- predict(m, x[itest,], s=bestparam)
        pred <- prediction(ytest, y[itest])
        auc[[met$displayname]] <- performance(pred, "auc")@y.values[[1]]
        nfeat[[met$displayname]] <- sum(abs(predict(m,s=bestparam,type="coefficients"))>epsilon)
        }
      
      return(list(auc=auc, nfeat=nfeat))
    }, mc.cores = numCores)
  
  # Compute mean AUC for each method
  auc <- sapply(resultloc, function(v) unlist(v[["auc"]]))
  print("Mean AUC")
  print(apply(auc, 1, mean))
  print("Standard deviation")
  print(apply(auc, 1, sd))
  # Compute P-values
  nmethods <- nrow(auc)
  pval <- matrix(NA,nrow=nmethods,ncol=nmethods,dimnames = list(rownames(auc),rownames(auc)))
  for (i in seq(nmethods)) {
    for (j in seq(nmethods)) {
      if (sd(auc[i,])+sd(auc[j,])>0) {
        pval[i,j] = t.test(auc[i,],auc[j,],alternative="greater")[["p.value"]]
      }
    }
  }
  # Number of selected features
  nfeat <- sapply(resultloc, function(v) unlist(v[["nfeat"]]))
  print("Mean Nfeat")
  print(apply(nfeat, 1, mean))
  print("Standard deviation")
  print(apply(nfeat, 1, sd))
  
  result[[idata]] = list(dataset=dataset, npos=npos, nneg=nneg, ngenes=ncol(x), auc=auc, auc.mean=apply(auc, 1, mean), auc.sd=apply(auc, 1, sd), pval=pval, nfeat=nfeat, nfeat.mean=apply(nfeat, 1, mean), nfeat.sd=apply(nfeat, 1, sd))
  idata = idata+1
  save(result, file="result.Rdata")
}


### Print and plot results ###

datasetnames <- sapply(conquerdata, function(v) v$displayname)
# Print mean AUC table
meanAUC <- t(sapply(result, function(v) v[["auc.mean"]]))
rownames(meanAUC) <- datasetnames
print(xtable(meanAUC), file="meanAUC.tex")
# Print mean nfeat table
meanfeat <- t(sapply(result, function(v) v[["nfeat.mean"]]))
rownames(meanfeat) <- datasetnames
print(xtable(meanAUC), file="meanNfeat.tex")
# Print pvalue table
smallpval <- Reduce('+', lapply(result, function(v) {pp <- v[["pval"]]; pp[is.na(pp)]=0.5; return(pp<0.05)}))
print(xtable(smallpval, caption=paste("Number of experiments, out of a total of ",length(result),", where each method (in row) is significantly better than each other method (in column).",sep="")),file="pval.tex")
# Plot how many times each method is significantly better than another
pdf("pvalwins.pdf",width=6,height=4)
par(mar = c(4, 7, 0, 0) + 0.2)
barplot(apply(smallpval, 1, sum),horiz=T,las=2,xlab="Number of wins")
dev.off()
# Table to summarize the datasets
datasettable <- t(sapply(result, function(v) c(v$dataset$displayname,v$npos,v$nneg,v$ngenes)))
colnames(datasettable) <- c("Data set","N pos","N neg","N genes")
print(xtable(datasettable), file="dataset.tex")
