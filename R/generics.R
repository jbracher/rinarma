#' Summarize an INARMA(1, 1) model fit
#'
#' @param model the model fit (an object of class `"inarma"`).
#' @export
summary.inarma <- function(model){
  cat("Summarizing an INARMA(1, 1) model with family =", model$family, "\n \n")
  cat("Estimated model parameters:\n")
  print(model$coefficients)
  cat("\n")
  if(model$fitting_method == "maximum_likelihood"){
    cat("Estimated standard errors:\n")
    print(model$se)
    cat("\n")
    cat("Number of observations used for fitting:", model$nobs, "\n")
    cat("AIC:", model$AIC, "\n")
    convergence_code <- model$optim$convergence
    convergence_text <- ifelse(convergence_code == 0, "successful", "failed")
    cat("\n")
    cat("Convergence of optimizer", convergence_text, paste0("(optim$convergence = ", convergence_code, ")."))
  }
  if(model$fitting_method == "moments"){
    cat("No standard errors or AIC are available for fits via method of moments.\n")
    cat("\n")
    cat("Number of observations used for fitting:", model$nobs, "\n")
  }
}

#' Print an INARMA(1, 1) model fit
#'
#' @param model the model fit (an object of class `"inarma"`).
#' @export
print.inarma <- function(model){
  summary(model) #re-use summary function.
}

#' Get AIC of an INARMA(1, 1) model
#' @param model the model (an object of class `"inarma"`)
#' @export
AIC.inarma <- function(model){
  if(model$fitting_method == "moments"){
    stop("No AIC available for model fitted via method of moments.")
  }else{
    model$AIC
  }
}

#' Get fitted values of an INARMA(1, 1) model
#' @param model the model (an object of class `"inarma"`)
#' @export
fitted.inarma <- function(model){
  if(model$fitting_method == "moments"){
    stop("No fitted values available for model fitted via method of moments.")
  }else{
    model$fitted_values
  }
}