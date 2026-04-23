#' Simulating 2 splits with MR-SimSS
#'
#' \code{split2} is a function which is used by the main function, \code{mr_simss} in order to perform the 2 split approach of the method,  \strong{MR-SimSS}.
#'
#' @param data A data frame to be inputted by the user containing summary
#'  statistics from the exposure and outcome GWASs. It must have at least five
#'  columns with column names \code{SNP}, \code{beta.exposure},
#'  \code{beta.outcome}, \code{se.exposure} and \code{se.outcome}. Each row must
#'  correspond to a unique SNP, identified by \code{SNP}.
#' @param lambda.val A numerical value which is computed within the main function, \code{mr_simss}.
#'  It is an estimate of \emph{lambda}, a term used to describe the correlation between the SNP-outcome and
#'  SNP-exposure effect sizes. The default setting is \code{lambda.val=0}.
#' @param pi A numerical value which determines the fraction of the first split. This is the fraction that will be used
#'  for SNP selection. The default setting is \code{pi=0.5}.
#' @param mr_method A string which specifies the MR method that MR-SimSS works in
#'  combination with. It is possible to use any method outputted in the list
#'  \code{TwoSampleMR::mr_method_list()$obj}. However, it is currently advised
#'  that the user chooses \code{"mr_ivw"} or \code{"mr_raps"}. The default
#'  setting is \code{mr_method="mr_ivw"}.
#' @param threshold A numerical value which specifies the threshold used to
#'  select instrument SNPs for MR at each iteration. The default setting is
#'  \code{threshold=5e-8}.
#'
#' @return A data frame which contains
#'  the output from one iteration. It is in a similar style as the output from
#'  using the function \code{mr} from the \code{TwoSampleMR} R package. The MR method used, the number of instrument SNPs, the causal effect estimate, it associated standard error and \emph{p}-value are all outputted.
#' @importFrom magrittr %>%
#' @export
#'


split2 <- function(data,lambda.val=0,pi=0.5,mr_method="mr_ivw", threshold=5e-8){

  data$beta.exposure.1 <- stats::rnorm(n=nrow(data),mean=data$beta.exposure,sd=sqrt((1-pi)/pi)*data$se.exposure)

  se.exposure.1 <-  sqrt(((1)/(pi))*((data$se.exposure)^2))
  pval.exposure.1 <- 2*(stats::pnorm(abs(data$beta.exposure.1/se.exposure.1), lower.tail=FALSE))

  data <- data %>% dplyr::filter(pval.exposure.1 < threshold)

  if(nrow(data) < 3){return(NULL)}else{
    ## Use conditional distribution to obtain variant-outcome estimates for selected variants
    ## reduces comp time considerably compared with use of bivariate normals for all variants
    mean.outcome <- data$beta.outcome + ((lambda.val*data$se.outcome)/(data$se.exposure))*(data$beta.exposure.1 - data$beta.exposure)
    sd.outcome <- data$se.outcome*sqrt(1-(((1-pi)/(pi))*((lambda.val)^2)))
    data$beta.outcome.1 <- stats::rnorm(n=nrow(data),mean=mean.outcome, sd=sd.outcome)

    beta.exposure.2 <- (data$beta.exposure - pi*data$beta.exposure.1)/(1-pi)
    beta.outcome.2 <- (data$beta.outcome - pi*data$beta.outcome.1)/(1-pi)
    se.exposure.2 <- sqrt(((1)/(1-pi))*((data$se.exposure)^2))
    se.outcome.2 <- sqrt(((1)/(1-pi))*((data$se.outcome)^2))

    data <- tibble::tibble(
      SNP = data$SNP,
      id.exposure="X",
      id.outcome="Y",
      exposure="X",
      outcome="Y",
      beta.exposure = beta.exposure.2,
      beta.outcome = beta.outcome.2,
      se.exposure = se.exposure.2,
      se.outcome = se.outcome.2,
      mr_keep=TRUE
    )

    if(mr_method=="mr_raps"){
      results <- tryCatch(
        mr.raps::mr.raps.mle(data$beta.exposure,data$beta.outcome,data$se.exposure,data$se.outcome),
        error = function(e) NULL
      )
      if(is.null(results)) return(NULL)
      results <- data.frame(method="mr_raps", nsnp=nrow(data), b=results$beta.hat, se=results$beta.se, pval=results$beta.p.value)
      return(results)
    }else{
      results <- TwoSampleMR::mr(data,method_list=mr_method)
      return(results[,5:9])
    }
  }
}

