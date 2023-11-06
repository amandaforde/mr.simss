#' Estimating lambda
#'
#' @param data A data frame to be inputted by the user containing summary
#' statistics from the exposure and outcome GWASs. It must have at least five
#' columns with column names \code{SNP}, \code{beta.exposure},
#' \code{beta.outcome}, \code{se.exposure}, and \code{se.outcome}. Each row
#' must correspond to a unique SNP, identified by \code{SNP}.
#' @param z.threshold A value
#'
#' @return A value
#' @export
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
      xmin = -0.5,xmax = 0.5, ymin = -0.5, ymax = 0.5)$Q))
  }

  est.lambda <- stats::optimize(fun_obj, interval=c(-1,1),maximum=TRUE)$maximum
  return(est.lambda)
}
