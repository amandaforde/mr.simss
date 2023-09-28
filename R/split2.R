#' Simulating 2 splits
#'
#' @param data A data frame to be inputted by the user containing summary
#' statistics from the exposure and outcome GWASs. It must have at least five
#' columns with column names \code{SNP}, \code{beta.exposure},
#' \code{beta.outcome}, \code{se.exposure}, \code{se.outcome} and
#' \code{eaf.exposure}. Each row must correspond to a unique SNP, identified by
#' \code{SNP}.
#' @param n.exposure A numerical value to be specified by the user which is equal
#' to the number of individuals that were in the exposure GWAS.
#' @param n.outcome A numerical value to be specified by the user which is equal
#' to the number of individuals that were in the outcome GWAS.
#' @param n.overlap A numerical value to be specified by the user which is equal
#' to the number of individuals that were in both the exposure and outcome GWAS.
#' @param correlation A numerical value to be specified by the user which is
#' equal to the observed correlation between the exposure and the outcome. This
#' value must be between -1 and 1.
#' @param pi A numerical value. The default setting is \code{pi=0.5}.
#' @param mr_method A string. The default setting is \code{mr_method="mr_ivw"}.
#' @param threshold A numerical value. The default setting is
#'\code{threshold=5e-8}.
#'
#' @return A data frame
#' @importFrom magrittr %>%
#' @export
#'


split2 <- function(data,n.exposure,n.outcome,n.overlap,correlation,pi=0.5,mr_method="mr_ivw", threshold=5e-8){
  data$maf <- data$eaf.exposure

  # create covariance matrix for the conditional distribution of each SNP
  cond_var_gx <- ((1-pi)/(pi))*(1/(n.exposure*2*data$maf*(1-data$maf)))
  cond_var_gy <- ((1-pi)/(pi))*(1/(n.outcome*2*data$maf*(1-data$maf)))
  cond_cov_gx_gy <- ((1-pi)/(pi))*(((n.overlap*correlation)/(n.exposure*n.outcome))*(1/(2*data$maf*(1-data$maf))))

  cond_cov_array <- array(dim=c(2, 2, nrow(data)))
  cond_cov_array[1,1,] <- cond_var_gx
  cond_cov_array[2,1,] <- cond_cov_gx_gy
  cond_cov_array[1,2,] <- cond_cov_array[2,1,]
  cond_cov_array[2,2,] <- cond_var_gy

  summary_stats_sub <- apply(cond_cov_array, 3, function(x) {MASS::mvrnorm(n=1, mu=c(0,0), Sigma=x)})

  summary_stats_sub1 <- t(summary_stats_sub + rbind(data$beta.exposure, data$beta.outcome))
  colnames(summary_stats_sub1) <- c("beta.exposure.1", "beta.outcome.1")
  data <- cbind(data, summary_stats_sub1)

  se.exposure.1 <-  sqrt(((1)/(pi))*(1/(n.exposure*2*data$maf*(1-data$maf))))
  pval.exposure.1 <- 2*(stats::pnorm(abs(data$beta.exposure.1/se.exposure.1), lower.tail=FALSE))

  data <- data %>% dplyr::filter(pval.exposure.1 < threshold)
  if(nrow(data) < 3){return(NULL)}else{
    beta.exposure.2 <- (data$beta.exposure - pi*data$beta.exposure.1)/(1-pi)
    beta.outcome.2 <- (data$beta.outcome - pi*data$beta.outcome.1)/(1-pi)
    se.exposure.2 <- sqrt(((1)/(1-pi))*(1/(n.exposure*2*data$maf*(1-data$maf))))
    se.outcome.2 <- sqrt(((1)/(1-pi))*(1/(n.outcome*2*data$maf*(1-data$maf))))

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
      fval.exposure = (beta.exposure.2/se.exposure.2)^2,
      fval.outcome = (beta.outcome.2/se.outcome.2)^2,
      pval.exposure = 2*(stats::pnorm(abs(beta.exposure.2/se.exposure.2), lower.tail=FALSE)),
      pval.outcome = 2*(stats::pnorm(abs(beta.outcome.2/se.outcome.2), lower.tail=FALSE)),
      eaf.exposure = data$maf,
      eaf.outcome = data$maf,
      mr_keep=TRUE
    )

    if(mr_method=="mr_raps"){
      results <- mr.raps::mr.raps(data$beta.exposure,data$beta.outcome,data$se.exposure,data$se.outcome)
      results <- data.frame(method="mr_raps", nsnp=nrow(data), b=results$beta.hat, se=results$beta.se, pval=results$beta.p.value)
      return(results)
    }else{
      results <- TwoSampleMR::mr(data,method_list=mr_method)
      return(results[,5:9])
    }
  }
}

