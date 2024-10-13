#' @title Measuring Agreement between a New Device and Multiple Imperfect Reference Standards
#' @param reference a references 'formula'. The expression 'rater1 + rater2 + rater3' indicates that three raters are being used, with the column names being 'rater1', 'rater2', and 'rater3'.
#' @param device the column name for device.
#' @param index the agreement index. Can be "CCC" or "TDI"
#' @param tranformation the transformation function g(x)
#' @param data a dataframe containing raters and CAD
#' @param p available when using TDI; indicating the fraction of items with differences numerically exceeding the TDI
#' @return a list of estimates, standard error, and confidence interval
#' @export
agreement <- function(references = "rater1 + rater2 + rater3", device = "CAD", index = c("CCC", "TDI"), transformation = function(x) 1 / x, data = data,p = 0.05, boot = NULL){
  raters_col <- trimws(unlist(strsplit(references, "\\+")))
  raters <- data[, c(raters_col)]
  CAD <- data[,c(device)]
  index <- match.arg(index)
  result <- switch(index,
                   "CCC" = agreement_CCC(references = raters, device = CAD, transformation = tranformation),
                   "TDI" = agreement_TDI(references = raters, device = CAD, transformation = tranformation, p = p))
  return(result)
}

#' @title Calculating Concordance Correlation Coefficient (CCC)
#' @export
#' @noRd
CCC <- function(y1,y2){
  n <- length(y1)
  m1<-mean(y1)
  m2<-mean(y2)
  S1<-var(y1)
  S2<-var(y2)
  S12<-cov(y1, y2)
  rho_C<-2*S12/(S1+S2+(m1-m2)^2)
  return(rho_C)
}

#' @title Calculate CCC between a New Device and Multiple Imperfect Reference Standards
#' @export
#' @noRd
agreement_CCC <- function(references = references, device = device,
                          transformation = function(x) 1 / (1-x)){
  # number of references
  n_raters <- as.character(ncol(references))
  results <- switch(n_raters,
                    "3" = agreement_CCC3(references, device, transformation),
                    "4" = agreement_CCC4(references, device, transformation))
  return(results)
}

#' @title Get weights 3 observations CCC
#' @export
#' @noRd
get_weight_3obs_CCC <- function(references, transformation = function(x) 1 / (1-x)){
  n <- ncol(references)
  u <- rep(0,n)
  for(i in 1:n) {
    u[i] <- with(as.data.frame(references), CCC( references[,i] , 0.5 * references[,-i][,1] + 0.5* references[,-i][,2]))
  }
  lambda1 <- transformation(u[1]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]))
  lambda2 <- transformation(u[2]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]))
  lambda3 <- transformation(u[3]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]))
  lambda <- c(lambda1,lambda2,lambda3)
  return(lambda)
}

#' @title Get CCC 3 observations
#' @export
#' @noRd
agreement_CCC3 <- function(references = references, device = device, transformation = function(x) 1 / (1-x)){
  references <- as.data.frame(references)
  weights <- get_weight_3obs_CCC(references, transformation)
  ccc_proposed <- CCC(device, as.matrix(references) %*% weights)
  ccc_naive <- CCC(device, as.matrix(references) %*% c(1/3, 1/3, 1/3))
  output <- c(ccc_proposed, ccc_naive, weights)
  return(output)
}

#' @title Get weights 4 observations CCC
#' @export
#' @noRd
get_weight_4obs_CCC <- function(references, transformation = function(x) 1 / (1-x)){
  n <- ncol(references)
  u <- rep(0,n)
  for(i in 1:n) {
    u[i] <- with(as.data.frame(references), CCC( references[,i] , as.matrix(references[,-i]) %*% get_weight_3obs_CCC(references[,-i], transformation) ))
  }
  lambda1 <- transformation(u[1]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]), transformation(u[4]))
  lambda2 <- transformation(u[2]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]), transformation(u[4]))
  lambda3 <- transformation(u[3]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]), transformation(u[4]))
  lambda4 <- transformation(u[4]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]), transformation(u[4]))
  lambda <- c(lambda1,lambda2,lambda3,lambda4)
  return(lambda)
}

#' @title Get CCC 4 observations
#' @export
#' @noRd
agreement_CCC4 <- function(references = references, device = device, transformation = function(x) 1 / (1-x)){
  references <- as.data.frame(references)
  weights <- get_weight_4obs_CCC(references, transformation)
  ccc_proposed <- CCC(device, as.matrix(references) %*% weights)
  ccc_naive <- CCC(device, as.matrix(references) %*% c(1/4, 1/4, 1/4, 1/4))
  output <- c(ccc_proposed, ccc_naive, weights)
  return(output)
}

#' @title Calculating Total Deviation Index (TDI)
#' @export
#' @noRd
TDI <- function( y1, y2, p = 0.05)
{
  if( length(y1) != length(y2) )
    stop( "Lengths of y1 and y2 must be the same!" )

  # Analytical approximation
  mu  <- mean(y1-y2)
  sigma <- sd(y1-y2)
  # TDI.appr <- qnorm(1-p/2)*sqrt(mu^2+sigma^2)

  # Limits of agreement
  LoA <- mu + c(-1,1)*qnorm(1-p/2)*sigma
  names( LoA ) <- c("lower","upper")

  # Numerical calculation
  FF <- function( x ) pnorm( (x-mu)/sigma ) - pnorm( (-x-mu)/sigma ) - (1-p)
  int <- 2 * max(abs(LoA)) * c(-1,1)
  TDI.num <- uniroot( FF, int )$root
  return(TDI.num)
}

#' @title Calculate TDI between a New Device and Multiple Imperfect Reference Standards
#' @export
#' @noRd
agreement_TDI <- function(references = references, device = device,
                          transformation = function(x) 1 / (x), p = 0.05){
  # number of references
  n_raters <- as.character(ncol(references))
  results <- switch(n_raters,
                    "3" = agreement_TDI3(references, device, transformation, p = p),
                    "4" = agreement_TDI4(references, device, transformation, p = p))
  return(results)
}

#' @title Get weights 3 observations TDI
#' @export
#' @noRd
get_weight_3obs_TDI <- function(references, transformation = function(x) 1 / (x), p = 0.05){
  n <- ncol(references)
  u <- rep(0,n)
  for(i in 1:n) {
    u[i] <- with(as.data.frame(references), TDI(references[,i] , 0.5 * references[,-i][,1] + 0.5* references[,-i][,2], p = p))
  }
  lambda1 <- transformation(u[1]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]))
  lambda2 <- transformation(u[2]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]))
  lambda3 <- transformation(u[3]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]))
  lambda <- c(lambda1,lambda2,lambda3)
  return(lambda)
}

#' @title Get TDI 3 observations
#' @export
#' @noRd
agreement_TDI3 <- function(references = references, device = device, transformation = function(x) 1 / (x), p = 0.05){
  references <- as.data.frame(references)
  weights <- get_weight_3obs_TDI(references, transformation, p = p)
  tdi_proposed <- TDI(device, as.matrix(references) %*% weights, p = p)
  tdi_naive <- TDI(device, as.matrix(references) %*% c(1/3, 1/3, 1/3), p = p)
  output <- c(tdi_proposed, tdi_naive, weights)
  return(output)
}


#' @title Get weights 4 observations TDI
#' @export
#' @noRd
get_weight_4obs_TDI <- function(references, transformation = function(x) 1 / (x), p = 0.05){
  n <- ncol(references)
  u <- rep(0,n)
  for(i in 1:n) {
    u[i] <- with(as.data.frame(references), TDI( references[,i] , as.matrix(references[,-i]) %*% get_weight_3obs_TDI(references[,-i], transformation), p = p))
  }
  lambda1 <- transformation(u[1]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]), transformation(u[4]))
  lambda2 <- transformation(u[2]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]), transformation(u[4]))
  lambda3 <- transformation(u[3]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]), transformation(u[4]))
  lambda4 <- transformation(u[4]) / sum(transformation(u[1]), transformation(u[2]), transformation(u[3]), transformation(u[4]))
  lambda <- c(lambda1,lambda2,lambda3,lambda4)
  return(lambda)
}

#' @title Get TDI 4 observations
#' @export
#' @noRd
agreement_TDI4 <- function(references = references, device = device, transformation = function(x) 1 / (x), p = 0.05){
  references <- as.data.frame(references)
  weights <- get_weight_4obs_TDI(references, transformation, p = p)
  tdi_proposed <- TDI(device, as.matrix(references) %*% weights)
  tdi_naive <- TDI(device, as.matrix(references) %*% c(1/4, 1/4, 1/4, 1/4))
  output <- c(tdi_proposed, tdi_naive, weights)
  return(output)
}
