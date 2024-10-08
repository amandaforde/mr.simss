#'MR-SimSS main function
#'
#'\code{mr_simss} is the main function for the method, \strong{MR-SimSS}, which
#'is a method based on simulated sample splitting in order to alleviate
#'\emph{Winner's Curse} bias in MR causal effect estimates. It also takes into
#'account sample overlap between the exposure and outcome GWASs. It uses GWAS
#'summary statistics and works in combination with existing MR methods, such as
#'IVW and MR-RAPS.
#'
#'@param data A data frame to be inputted by the user containing summary
#'  statistics from the exposure and outcome GWASs. It must have at least five
#'  columns with column names \code{SNP}, \code{beta.exposure},
#'  \code{beta.outcome}, \code{se.exposure} and \code{se.outcome}. Each row must
#'  correspond to a unique SNP, identified by \code{SNP}.
#'@param subset A logical which permits the user to perform this method with
#'  either the original complete set of SNPs or a subset of SNPs in order to
#'  reduce computational time. The default setting is \code{subset=FALSE}.
#'@param sub.cut A numerical value required if \code{subset=TRUE}, which ensures
#'  that for a single iteration of our method, the number of instruments
#'  selected if the full set of SNPs is used and the number of instruments if
#'  merely the subset is used will be equal with probability at least
#'  \code{1-sub.cut}.
#'@param est.lambda A logical which allows the user to specify if they want to use
#'  the function, \code{est_lambda}, to obtain an estimate for \emph{lambda}, a
#'  term used to describe the correlation between the SNP-outcome and
#'  SNP-exposure effect sizes, with the data stored in \code{data}. This correlation is affected by the number of
#'  overlapping samples between the two GWASs and the correlation between the
#'  exposure and the outcome. Thus, it is recommended to use \code{est_lambda}
#'  if the fraction of overlap and the correlation between exposure and outcome
#'  are unknown. The default setting is \code{est.lambda=TRUE}.
#'@param lambda A numerical value which should be specified by the user if they have used the function \code{est_lambda} with different data than that inputted into this function in order to to obtain an estimate for \emph{lambda}.
#' The default setting is \code{lambda=NULL}.
#'@param n.exposure A numerical value to be specified by the user which is equal
#'  to the number of individuals that were in the exposure GWAS. It should be
#'  specified by the user if \code{est.lambda=FALSE}. The default setting is
#'  \code{n.exposure=1}.
#'@param n.outcome A numerical value to be specified by the user which is equal
#'  to the number of individuals that were in the outcome GWAS. It should be
#'  specified by the user if \code{est.lambda=FALSE}. The default setting is
#'  \code{n.outcome=1}.
#'@param n.overlap A numerical value to be specified by the user which is equal
#'  to the number of individuals that were in both the exposure and outcome
#'  GWAS. It should be specified by the user if \code{est.lambda=FALSE}. The
#'  default setting is \code{n.overlap=1}. The function requires that this value
#'  is less than or equal to the minimum of \code{n.exposure} and
#'  \code{n.outcome}.
#'@param cor.xy A numerical value to be specified by the user which is equal to
#'  the observed correlation between the exposure and the outcome. This value
#'  must be between -1 and 1. It should be specified by the user if
#'  \code{est.lambda=FALSE}. The default setting is \code{cor.xy=0}. If this
#'  value is unknown, the user is encouraged to use the function
#'  \code{est_lambda}.
#'@param n.iter A numerical value which specifies the number of iterations of
#'  the method, i.e. the number of times sample splits are randomly simulated.
#'  The default setting is \code{n.iter=1000}.
#'@param splits A numerical value that must be equal to 2 or 3, indicating
#'  whether splits of 2 or 3 should be simulated. It is recommended that in the
#'  case of no overlap between the two GWASs that splits of 2 should be used
#'  while in the presence of overlap, especially full overlap, splits of 3
#'  should be used. The default setting is \code{splits=2}.
#'@param pi A numerical value which determines the fraction of the first split
#'  in both the 2 and 3 split approaches. This is the fraction that will be used
#'  for SNP selection. The default setting is \code{pi=0.5}. This value must be
#'  between 0 and 1.
#'@param pi2 A numerical value which determines the fraction of the second split
#'  in the 3 split approach. The default setting is \code{pi2=0.5}. This value
#'  must be between 0 and 1.
#'@param threshold A numerical value which specifies the threshold used to
#'  select instrument SNPs for MR at each iteration. The default setting is
#'  \code{threshold=5e-8}. This value must be between 0 and 1.
#'@param mr_method A string which specifies the MR method that MR-SimSS works in
#'  combination with. It is possible to use any method outputted in the list
#'  \code{TwoSampleMR::mr_method_list()$obj}. However, it is currently advised
#'  that the user chooses \code{"mr_ivw"} or \code{"mr_raps"}. The default
#'  setting is \code{mr_method="mr_ivw"}.
#'@param parallel A logical value which allows the user to specify if they wish
#'  to use this function in parallel or in series. The default setting is
#'  \code{parallel=TRUE}. It is advisable to use this default, especially when
#'  \code{n.iter} is large.
#'@param n.cores A numerical value which determines how many cores will be used
#'  if \code{parallel=TRUE}. This value should be supplied by the user if they
#'  wish to use less cores than the output of \code{parallel::detectCores()-1}.
#'  The default setting is \code{n.cores=NULL}.
#'@param lambda.thresh A value which is used when estimating \emph{lambda} to
#'  obtain a subset of SNPs which have absolute \emph{z}-statistics for both exposure and outcome GWASs less than
#'  this value. The method then assumes that both of the true SNP-outcome and
#'  SNP-exposure effect sizes of each SNP in this subset are approximately 0.
#'  The default setting is \code{lambda.thresh=0.5}.
#'
#'@return A list containing two elements, \code{summary} and \code{results}.
#'  \code{summary} is a data frame with one row which outputs \code{b}, the
#'  estimated causal effect of exposure on outcome obtained using the
#'  \strong{MR-SimSS} method, as well as \code{se}, the associated standard
#'  error of this estimate and \code{pval}, corresponding \emph{p}-value. It
#'  also contains the MR method used, the average number of instrument SNPs used
#'  in each iteration and the number of iterations used. \code{results} is a
#'  data frame which contains the output from each iteration. It is in a similar
#'  style as the output from using the function \code{mr} from the
#'  \code{TwoSampleMR} R package.
#'@seealso
#'\url{https://amandaforde.github.io/mr.simss/articles/perform-MR-SimSS.html}
#'for illustration of the use of \code{mr_simss} with a toy data set and further
#'information regarding this MR method.
#'@importFrom foreach %dopar%
#'@export
#'

mr_simss <- function(data,subset=FALSE,sub.cut=0.05,est.lambda=TRUE,lambda = NULL,n.exposure=1,n.outcome=1,n.overlap=1,cor.xy=0,
                     n.iter=1000,splits=2,pi=0.5,pi2=0.5,threshold=5e-8,mr_method="mr_ivw",
                     parallel=TRUE,n.cores=NULL,lambda.thresh=0.5){

  ## ensuring correct use of function
  stopifnot(all(c("SNP", "beta.exposure","beta.outcome","se.exposure","se.outcome") %in% names(data)))
  stopifnot(splits ==  2 | splits == 3)
  stopifnot(cor.xy >= -1 && cor.xy <= 1)
  stopifnot(pi > 0 && pi < 1 && pi2 > 0 && pi2 < 1)
  stopifnot(mr_method %in% TwoSampleMR::mr_method_list()$obj)
  stopifnot(n.overlap <= min(n.outcome,n.exposure))
  stopifnot(threshold >= 0 && threshold <= 1)

  ## work out lambda here - neater ways to write this?
  if(is.null(lambda) == TRUE){
    if(est.lambda==TRUE){
      lambda.val <- mr.simss::est_lambda(data,z.threshold=lambda.thresh)
    }else{
      lambda.val <- (n.overlap*cor.xy)/(sqrt(n.exposure*n.outcome))
    }
  }else{
    lambda.val <- lambda
  }


  if(subset == TRUE){
    thresh1 <- stats::qnorm((threshold)/2, lower.tail=FALSE)
    data$mu.exposure1 <- data$beta.exposure/sqrt(((1)/(pi))*(data$se.exposure)^2)
    data$p.exposure1 <- stats::pnorm((-thresh1 + data$mu.exposure1)/(sqrt(1-pi))) + stats::pnorm((-thresh1 - data$mu.exposure1)/(sqrt(1-pi)))
    data1 <- data[order(data$p.exposure1,decreasing=FALSE),]
    data1$cumul.sum <- cumsum(data1$p.exposure1)
    data <- data1[-which(data1$cumul.sum < sub.cut),]
  }

  if(parallel == TRUE){
    if(is.null(n.cores)){n_cores <- parallel::detectCores()-1}else{n_cores <- n.cores}
    my.cluster <- parallel::makeCluster(n_cores,type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)

    if(splits==2){
      results <- foreach::foreach(i = 1:n.iter,.packages=c('tidyverse','mr.simss'),.combine = 'rbind') %dopar% {mr.simss::split2(data,lambda.val=lambda.val,pi,mr_method,threshold)}
    }else{
      results <- foreach::foreach(i = 1:n.iter,.packages=c('tidyverse','mr.simss'),.combine = 'rbind') %dopar% {mr.simss::split3(data,lambda.val=lambda.val,pi,pi2,mr_method,threshold)}
    }
    parallel::stopCluster(cl = my.cluster)

  }else{
    results <- c()
    if(splits==2){
      for (i in 1:n.iter){
        wc_remove <- mr.simss::split2(data,lambda.val=lambda.val,pi,mr_method,threshold)
        if(is.null(wc_remove) == FALSE){results <- rbind(results,wc_remove)}
      }
    }else{
      for (i in 1:n.iter){
        wc_remove <- mr.simss::split3(data,lambda.val=lambda.val,pi,pi2,mr_method,threshold)
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

