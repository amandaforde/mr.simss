#'Estimating correlation between SNP-exposure and SNP-outcome effect sizes
#'
#'\code{est_lambda} is a function which allows users to obtain an
#'\emph{unbiased} estimate for \emph{lambda}, a term used to describe the
#'correlation between the SNP-outcome and SNP-exposure effect sizes, using a conditional log-likelihood approach. This
#'correlation is affected by the number of overlapping samples between the two
#'GWASs and the correlation between the exposure and the outcome. Thus, when
#'using the function \code{mr_simss}, if the fraction of overlap and the
#'correlation between exposure and outcome are unknown, it is recommended to
#'employ \code{est_lambda} and use the value returned from \code{est_lambda} in
#'\code{mr_simss}.
#'\strong{Note:} For greater accuracy in the estimation of
#'\emph{lambda}, it is advisable to use summary  statistics of the entire set of
#'\emph{unpruned} SNPs from the exposure and outcome GWASs.
#'
#'@param data A data frame to be inputted by the user containing summary
#'  statistics from the exposure and outcome GWASs. It must have at least five
#'  columns with column names \code{SNP}, \code{beta.exposure},
#'  \code{beta.outcome}, \code{se.exposure}, and \code{se.outcome}. Each row
#'  must correspond to a unique SNP, identified by \code{SNP}.
#'@param z.threshold A value which is used to obtain a subset of SNPs which have
#'  absolute \emph{z}-statistics for both exposure and outcome GWASs less than
#'  this value. The method then assumes that both of the true SNP-outcome and
#'  SNP-exposure effect sizes of each SNP in this subset are approximately 0.
#'  The default setting is \code{z.threshold=0.5}.
#'
#'@return A value which is an estimate of \emph{lambda}, the correlation between the SNP-outcome and SNP-exposure effect sizes, using a conditional log-likelihood approach. Note that this estimate is unbiased but potentially has a high degree of variance.
#'@seealso
#'\url{https://amandaforde.github.io/mr.simss/articles/perform-MR-SimSS.html}
#'for illustration of the use of \code{est_lambda} with a toy data set and
#'\url{https://amandaforde.github.io/mr.simss/articles/derive-MR-SimSS.html} for
#'the theoretical derivation of this method based on a conditional
#'log-likelihood approach for estimating \emph{lambda}.
#'@export
#'
#'
est_lambda <- function(data,z.threshold=0.5){
  data <- data[abs(data$beta.exposure/data$se.exposure) < z.threshold & abs(data$beta.outcome/data$se.outcome) < z.threshold,]
  n <- nrow(data)
  B <- sum((data$beta.exposure/data$se.exposure)*(data$beta.outcome/data$se.outcome))
  A <- sum((data$beta.exposure/data$se.exposure)^2)
  C <- sum((data$beta.outcome/data$se.outcome)^2)

  fun_obj <- function(lambda){
    - ((1)/(2*(1-lambda^2)))*(A-2*lambda*B+C) - n*log((pracma::integral2(function(x,y)
      exp((-(1)/(2*(1-lambda^2)))*(x^2-2*lambda*x*y+y^2)),
      xmin = -z.threshold,xmax = z.threshold, ymin = -z.threshold, ymax = z.threshold)$Q))
  }

  est.lambda <- stats::optimize(fun_obj, interval=c(-1,1),maximum=TRUE)$maximum
  return(est.lambda)
}


