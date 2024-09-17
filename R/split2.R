#' Simulating 2 splits with MR-SimSS
#'
#' \code{split2} is a function which is used by the main function, \code{mr_simss} in order to perform the 2 split approach of the method,  \strong{MR-SimSS}.
#'
#'@param data A data frame to be inputted by the user containing summary
#'  statistics from the exposure and outcome GWASs. It must have at least five
#'  columns with column names \code{SNP}, \code{beta.exposure},
#'  \code{beta.outcome}, \code{se.exposure} and \code{se.outcome}. Each row must
#'  correspond to a unique SNP, identified by \code{SNP}.
#'@param lambda.val A numerical value which is computed within the main function, \code{mr_simss}.
#' It is an estimate of \emph{lambda}, a term used to describe the correlation between the SNP-outcome and
#'  SNP-exposure effect sizes. The default setting is \code{lambda.val=0}.
#'@param pi A numerical value which determines the fraction of the first split. This is the fraction that will be used
#'  for SNP selection. The default setting is \code{pi=0.5}.
#'@param mr_method A string which specifies the MR method that MR-SimSS works in
#'  combination with. It is possible to use any method outputted in the list
#'  \code{TwoSampleMR::mr_method_list()$obj}. However, it is currently advised
#'  that the user chooses \code{"mr_ivw"} or \code{"mr_raps"}. The default
#'  setting is \code{mr_method="mr_ivw"}.
#'@param threshold A numerical value which specifies the threshold used to
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

  # create covariance matrix for the conditional distribution of each SNP
  cond_var_gx <- ((1-pi)/(pi))*(data$se.exposure)^2
  cond_var_gy <- ((1-pi)/(pi))*(data$se.outcome)^2
  cond_cov_gx_gy <- ((1-pi)/(pi))*(data$se.exposure)*(data$se.outcome)*(lambda.val)

  cond_cov_array <- array(dim=c(2, 2, nrow(data)))
  cond_cov_array[1,1,] <- cond_var_gx
  cond_cov_array[2,1,] <- cond_cov_gx_gy
  cond_cov_array[1,2,] <- cond_cov_array[2,1,]
  cond_cov_array[2,2,] <- cond_var_gy

  #summary_stats_sub <- apply(cond_cov_array, 3, function(x) {MASS::mvrnorm(n=1, mu=c(0,0), Sigma=x)})
  # mvrnorm replaced by the following:
  tr_vec <- apply(cond_cov_array,3,function(x){x[1,1]+x[2,2]})
  s_vec <- apply(cond_cov_array,3,function(x){sqrt(sum(x[1,1]*x[2,2]-x[1,2]*x[2,1]))})
  t_vec <- sqrt(tr_vec+2*s_vec)
  s_vec <- rep(s_vec,times=rep(4,nrow(data)))
  I_vec <- rep(c(1,0,0,1),times=rep(nrow(data)))
  t_vec <- rep(t_vec,times=rep(4,nrow(data)))
  sqrt_array <- array((as.vector(cond_cov_array)+s_vec*I_vec)/t_vec,dim=c(2,2,nrow(data)))

  Z_array <- array(stats::rnorm(2*nrow(data)),dim=c(1,2,nrow(data)))
  Z_array <- abind::abind(Z_array,Z_array,along=1)
  sqrt_array_normal <- sqrt_array*Z_array
  # rearrange array so to use matrix multiplication
  sqrt_array_normal <- aperm(a=sqrt_array_normal,perm=c(3,1,2))
  dim1 <- sqrt_array_normal[,1,] %*% matrix(c(1,1),nrow=2)
  dim2 <- sqrt_array_normal[,2,] %*% matrix(c(1,1),nrow=2)
  summary_stats_sub <- cbind(dim1,dim2)

  summary_stats_sub1 <- (summary_stats_sub + cbind(data$beta.exposure, data$beta.outcome))

  #summary_stats_sub1 <- t(summary_stats_sub + rbind(data$beta.exposure, data$beta.outcome))
  colnames(summary_stats_sub1) <- c("beta.exposure.1", "beta.outcome.1")
  data <- cbind(data, summary_stats_sub1)

  se.exposure.1 <-  sqrt(((1)/(pi))*((data$se.exposure)^2))
  pval.exposure.1 <- 2*(stats::pnorm(abs(data$beta.exposure.1/se.exposure.1), lower.tail=FALSE))

  data <- data %>% dplyr::filter(pval.exposure.1 < threshold)
  if(nrow(data) < 3){return(NULL)}else{
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
      results <- mr.raps::mr.raps(data$beta.exposure,data$beta.outcome,data$se.exposure,data$se.outcome)
      results <- data.frame(method="mr_raps", nsnp=nrow(data), b=results$beta.hat, se=results$beta.se, pval=results$beta.p.value)
      return(results)
    }else{
      results <- TwoSampleMR::mr(data,method_list=mr_method)
      return(results[,5:9])
    }
  }
}

