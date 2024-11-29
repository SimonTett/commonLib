# code to call extRemes::fevd and generate summary statistics
fevd_sum <- function( x,data, threshold = NULL, threshold.fun = ~1, location.fun = ~1,
  scale.fun = ~1, shape.fun = ~1, use.phi = FALSE,
  type = c("GEV", "GP", "PP", "Gumbel", "Exponential"),
  method = c("MLE", "GMLE", "Bayesian", "Lmoments"), initial = NULL,
  span, units = NULL, time.units = "days", period.basis = "year",
  #na.action = na.fail,  # commented out because not used as stuffs up the covariance matrix
  optim.args = NULL, priorFun = NULL,
  priorParams = NULL, proposalFun = NULL, proposalParams = NULL,
  iter = 9999, weights = 1, blocks = NULL, verbose = FALSE,
  test_time = FALSE) {

  if (verbose) {

    print("args are ")
    print(paste("x:",x))
    print(paste("data.colnames:",colnames(data)))
    print(paste("data shape:",dim(data)))
    print(paste("location.fun:",location.fun))
    print(paste("scale.fun:",scale.fun))
    print(paste("initial:",initial))
    print(paste("use.phi:",use.phi))

  }
  if (test_time) { # do not run the function, just return some dummy data.
    sum = list(par.names = c("location", "scale", "shape"), par = c(1, 2, 3),
               se.theta=c(0.1,0.2,0.3),
               cov.theta=matrix(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),3,3),nllh=1.0,AIC=2.0,BIC=3.0)
    return(sum)
  }

  fit <- extRemes::fevd(x=x,data,
    threshold = threshold, threshold.fun = threshold.fun,
    location.fun = location.fun, scale.fun = scale.fun, shape.fun = shape.fun,
    use.phi = use.phi, type = type, method = method, initial = initial,
    span = span, units = units, time.units = time.units, period.basis = period.basis,
    #na.action = na.action, # having this means the covariance matrix is not computed...
    optim.args = optim.args, priorFun = priorFun,
    priorParams = priorParams, proposalFun = proposalFun, proposalParams = proposalParams,
    iter = iter, weights = weights,
    blocks = blocks,
    verbose = verbose)
  sum <- summary(fit,silent=TRUE)
  sum$par.names <- names(fit$results$par) # add the parameter estimates to the summary with names



  return(sum)
}