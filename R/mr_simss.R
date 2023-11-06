#'Main function for MR-SimSS
#'
#'\code{mr_simss} is a function which is based on simulated sample splitting in
#' order to alleviate Winner's Curse bias in MR causal effect estimates. It also
#' takes into account sample overlap between the exposure and outcome GWASs. It
#' uses GWAS summary statistics together with knowledge of the number of
#' overlapping samples and the correlation between the exposure and the outcome.
#' It works in combination with existing MR methods, such as IVW and MR-RAPS.
#'
#'@param data A data frame to be inputted by the user containing summary
#' statistics from the exposure and outcome GWASs. It must have at least five
#' columns with column names \code{SNP}, \code{beta.exposure},
#' \code{beta.outcome}, \code{se.exposure} and \code{se.outcome}. Each row must correspond to a unique SNP, identified by
#' \code{SNP}.
#'@param est.lambda A logical value.
#'@param lambda.val A numerical value.
#'@param n.exposure A numerical value to be specified by the user which is equal
#' to the number of individuals that were in the exposure GWAS.
#'@param n.outcome A numerical value to be specified by the user which is equal
#' to the number of individuals that were in the outcome GWAS.
#'@param n.overlap A numerical value to be specified by the user which is equal
#' to the number of individuals that were in both the exposure and outcome GWAS.
#'@param cor.xy A numerical value to be specified by the user which is
#' equal to the observed correlation between the exposure and the outcome. This
#' value must be between -1 and 1.
#'@param n.iter A numerical value which specifies the number of iterations of
#' the method, i.e. the number of times sample splits are randomly simulated.
#' The default setting is \code{n.iter=1000}.
#'@param splits A numerical value that must be equal to 2 or 3, indicating
#' whether splits of 2 or 3 should be simulated. It is recommended that in the
#' case of no overlap between the two GWASs that splits of 2 should be used
#' while in the presence of overlap, especially full overlap, splits of 3 should
#' be used. The default setting is \code{splits=2}.
#'@param pi A numerical value . The default setting is \code{pi=0.5}.
#'@param pi2 A numerical value . The default setting is \code{pi2=0.5}.
#'@param threshold A numerical value. The default setting is
#'\code{threshold=5e-8}.
#'@param mr_method A string. The default setting is \code{mr_method="mr_ivw"}.
#'@param parallel A logical value . The default setting is \code{parallel=TRUE}.
#'@param n.cores A numerical value . The default setting is \code{n.cores=NULL}.
#'
#' @return A list
#' @importFrom foreach %dopar%
#' @export
#'

mr_simss <- function(data,est.lambda=FALSE,lambda.val=0,n.exposure=1,n.outcome=1,n.overlap=1,cor.xy=0,
                     n.iter=1000,splits=2,pi=0.5,pi2=0.5,threshold=5e-8,mr_method="mr_ivw",
                     parallel=TRUE,n.cores=NULL){

  ## ensuring correct use of function
  stopifnot(all(c("SNP", "beta.exposure","beta.outcome","se.exposure","se.outcome") %in% names(data)))
  stopifnot(splits ==  2 | splits == 3)
  stopifnot(cor.xy >= -1 && cor.xy <= 1)
  stopifnot(pi > 0 && pi < 1 && pi2 > 0 && pi2 < 1)
  stopifnot(mr_method %in% TwoSampleMR::mr_method_list()$obj)
  stopifnot(n.overlap <= min(n.outcome,n.exposure))
  stopifnot(threshold >= 0 && threshold <= 1)
  stopifnot(est.lambda == FALSE && n.exposure > 1 && n.outcome > 1)

  if(parallel == TRUE){
    if(is.null(n.cores)){n_cores <- parallel::detectCores()-1}else{n_cores <- n.cores}
    my.cluster <- parallel::makeCluster(n_cores,type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)

    if(splits==2){
      results <- foreach::foreach(i = 1:n.iter,.packages=c('tidyverse','mr.simss'),.combine = 'rbind') %dopar% {mr.simss::split2(data,est.lambda,lambda.val,n.exposure,n.outcome,n.overlap,cor.xy,pi,mr_method,threshold)}
    }else{
      results <- foreach::foreach(i = 1:n.iter,.packages=c('tidyverse','mr.simss'),.combine = 'rbind') %dopar% {mr.simss::split3(data,est.lambda,lambda.val,n.exposure,n.outcome,n.overlap,cor.xy,pi,pi2,mr_method,threshold)}
    }
    parallel::stopCluster(cl = my.cluster)

  }else{
    results <- c()
    if(splits==2){
      for (i in 1:n.iter){
        wc_remove <- mr.simss::split2(data,est.lambda,lambda.val,n.exposure,n.outcome,n.overlap,cor.xy,pi,mr_method,threshold)
        if(is.null(wc_remove) == FALSE){results <- rbind(results,wc_remove)}
      }
    }else{
      for (i in 1:n.iter){
        wc_remove <- mr.simss::split3(data,est.lambda,lambda.val,n.exposure,n.outcome,n.overlap,cor.xy,pi,pi2,mr_method,threshold)
        if(is.null(wc_remove) == FALSE){results <- rbind(results,wc_remove)}
      }
    }
  }

  if(length(results) == 0){return(NULL)}else{
    est_se <- sqrt((sum(results$se^2)-sum((results$b-mean(results$b))^2))/(nrow(results)))
    summary <- data.frame(method=c(mr_method), nsnp = c(mean(results$nsnp)), niter= nrow(results), b = c(mean(results$b)), se=c(est_se), pval=c(2*stats::pnorm(mean(results$b)/est_se, lower.tail=FALSE)))
    return(list("summary" = summary, "results" = results))
  }
}

