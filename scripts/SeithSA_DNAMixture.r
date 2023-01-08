# Code in progress to optimize DNA Mixture likelihood function
# S.A. Sethi, suresh.sethi@cornell.edu
# v. April 16, 2020

library(foreach)
library(doParallel)
library(doRNG) # needed for 'set seed' across parallelization capabilities
library(RcppAlgos) # faster combinations?
library(forensim)
library(gtools)

# This code is work under progress to optimize DNA mixture likelihood functions.  Currently implementing CPP modules, 
# and vectorizing to speed up. Also, 'parsing' tasks to dump RAM and avoid memory overflow.  Works, but needs
# further vetting and cleaning, e.g. to include error traps.

# Development notes:
# Chokepoint 1 is combinations
# chokepoint 2 is sum over large number of combinations
# see:https://stackoverflow.com/questions/30239278/r-combinations-looking-for-faster-and-more-efficient-waypackage-code-parallel
# Check Rcpp for sum portion
# https://teuder.github.io/rcpp4everyone_en/030_basic_usage.html#writing-your-rcpp-code


###### Original version from Sethi et al.: Mixture likelihood function
MixtureLikelihood <- function(x,p.v){
  require(gtools)
  counter <- 0
  if (length(p.v)>0){
    sum.p <- sum(p.v)^(2*x)
    if(length(p.v)>1){ # only conduct if > 1 allele is observed
      for(i in 1:(length(p.v)-1)){
        counter <- counter+1
        temp.combo <- combinations(length(p.v),length(p.v)-i)
        temp.sum <- apply(temp.combo,MARGIN=1,function(z){sum(p.v[z])^(2*x)})
        sum.p <- sum.p+((-1)^counter)*sum(temp.sum)
      }# end i loop
    } # end if
  } else sum.p <- NA # if no alleles observed, assign NA
  return(sum.p)
} # end function


###### Optimized mixture function under progress (works, but includes a cheezy means of parsing tasks and dumping RAM to avoid
# memory overflow)
# This function version does a few things, first combinations are optimized and parallelized via RcppAlgos.
# Second, the sum over rows of the combinations is vectorized which is like a 6X speedup. 
MixtureLikelihood_Optimized <- function(x,p.v){
  # require(gtools)
  require(RcppAlgos)
  counter <- 0
  if (length(p.v)>0){
    sum.p <- sum(p.v)^(2*x)
    if(length(p.v)>1){ # only conduct if > 1 allele is observed
      for(i in 1:(length(p.v)-1)){
        counter <- counter+1
        # Combinatorics
        # original
        # temp.combo <- combinations(length(p.v),length(p.v)-i) # this gets very large when many alleles and i  
        # new, Rcpp based combinations function from RcppAlgos
        temp.combo <- comboGeneral(v=length(p.v), m = length(p.v)-i, repetition = FALSE, Parallel=TRUE,nThreads=10) #; i; dim(temp.combo);i<-i+1
        # Row sums across allele combinations  
        # original
        # temp.sum <- apply(temp.combo,MARGIN=1,function(z){sum(p.v[z])^(2*x)})  # this slows way down too with large combos
        # sum.p <- sum.p+((-1)^counter)*sum(temp.sum)
        # new, vectorized version
        temp.combo.v <- c(temp.combo)
        p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.combo)[1],nc=dim(temp.combo)[2],byrow=F)
        row.sums.pow2x <- rowSums(p.v.m)^(2*x) # reference as folows
        sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x) 
      }# end i loop
    } # end if
  } else sum.p <- NA # if no alleles observed, assign NA
  return(sum.p)
} # end function

###### Optimized mixture function with memory parsing to avoid memory overflow, under progress 
# (works, but includes a cheezy means of parsing tasks and dumping RAM to avoid memory overflow)
# This function version does a few things, first combinations are optimized and parallelized via RcppAlgos.
# Second, the sum over rows of the combinations is vectorized which is like a 6X speedup.
# Finally, to deal with memory overflow, I included a loop--which is slow--to parse sums over the huge
# combinations matrix.  This code calculates the row sums calculations over chunks of 999,999 rows at a time
# which should avoid memory problems. !!!BUT NOTE, there could be an error here if it so happens that a tempcombo is exactly
# divisible by 999999.  Further error trapping to avoid that would be important, but I don't have time right now.

MixtureLikelihood_Optimized_Parsed <- function(x,p.v){
  # require(gtools)
  require(RcppAlgos)
  counter <- 0
  if (length(p.v)>0){
    sum.p <- sum(p.v)^(2*x)
    if(length(p.v)>1){ # only conduct if > 1 allele is observed
      for(i in 1:(length(p.v)-1)){
        counter <- counter+1
        # Combinatorics
        # original
        # temp.combo <- combinations(length(p.v),length(p.v)-i) # this gets very large when many alleles and i  
        # new, Rcpp based combinations function from RcppAlgos
        temp.combo <- comboGeneral(v=length(p.v), m = length(p.v)-i, repetition = FALSE, Parallel=TRUE,nThreads=10) #; i; dim(temp.combo);i<-i+1
        # Row sums across allele combinations  
        # original
        # temp.sum <- apply(temp.combo,MARGIN=1,function(z){sum(p.v[z])^(2*x)})  # this slows way down too with large combos
        # sum.p <- sum.p+((-1)^counter)*sum(temp.sum)
        # new, vectorized version that parses row sums to save on memory
        if(nrow(temp.combo)>15e6){
          ix <- 999999		
          row.sums.pow2x <- 1:(floor(nrow(temp.combo)/ix)+1)
          for(m in 1:length(row.sums.pow2x)){
            temp.m <- temp.combo[ ((m-1)*ix+1) : min(m*ix,nrow(temp.combo)), ]
            temp.combo.v <- c(temp.m)
            p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.m)[1],nc=dim(temp.m)[2],byrow=F)
            row.sums.pow2x[m] <- sum(rowSums(p.v.m)^(2*x)) # reference as folows
          } # end m loop
          sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x) } else 
          {
            # non-parsed
            temp.combo.v <- c(temp.combo)
            p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.combo)[1],nc=dim(temp.combo)[2],byrow=F)
            row.sums.pow2x <- rowSums(p.v.m)^(2*x) # reference as folows
            sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x) 
          } # end parsing option
      }# end i loop
    } # end if
  } else sum.p <- NA # if no alleles observed, assign NA
  return(sum.p)
} # end function



#### Multi-locus functions across 1:x number of putative contributors

# Using Haned et al. function from the forensim{} package
dataL_multi <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- dataL(x=j,p=unlist(loci[i]),theta=0) # function
    }  # end j loop
  } # end i loop
  Likelihood_df$Product <- apply(Likelihood_df, 1, function(x) prod(x, na.rm = TRUE)) # product across all loci, exclude NAs
  Likelihood_df$Contributors <- c(1:y) # putative numbers of contributors
  return(Likelihood_df[Likelihood_df$Product==max(Likelihood_df$Product),]) # find row with max likelihood across all loci
} # end function

# Original version
MixtureLikelihood_multi <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- MixtureLikelihood(x=j,p.v=unlist(loci[i])) # function
    }  # end j loop
  } # end i loop
  Likelihood_df$Product <- apply(Likelihood_df, 1, function(x) prod(x, na.rm = TRUE)) # product across all loci, exclude NAs
  Likelihood_df$Contributors <- c(1:y) # putative numbers of contributors
  return(Likelihood_df[Likelihood_df$Product==max(Likelihood_df$Product),]) # find row with max likelihood across all loci
} # end function

# Optimized version
MixtureLikelihood_Optimized_multi <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- MixtureLikelihood_Optimized(x=j,p.v=unlist(loci[i])) # function
    }  # end j loop
  } # end i loop
  Likelihood_df$Product <- apply(Likelihood_df, 1, function(x) prod(x, na.rm = TRUE)) # product across all loci, exclude NAs
  Likelihood_df$Contributors <- c(1:y) # putative numbers of contributors
  return(Likelihood_df[Likelihood_df$Product==max(Likelihood_df$Product),]) # find row with max likelihood across all loci
} # end function

# Optimized and task parsed version to avoid memory overflow
MixtureLikelihood_Optimized_Parsed_multi <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- MixtureLikelihood_Optimized_Parsed(x=j,p.v=unlist(loci[i])) # function
    }  # end j loop
  } # end i loop
  Likelihood_df$Product <- apply(Likelihood_df, 1, function(x) prod(x, na.rm = TRUE)) # product across all loci, exclude NAs
  Likelihood_df$Contributors <- c(1:y) # putative numbers of contributors
  return(Likelihood_df[Likelihood_df$Product==max(Likelihood_df$Product),]) # find row with max likelihood across all loci
} # end function

