#' Models currently supported by baggr
#' see \code{baggr()}
baggr_models <- c("rubin", "quantiles", "mutau", "individual")

#' Predict method for baggr objects
#' @param x model fitted via baggr()
#' @param newdata new data to use with the fitted model
#' @param allow_new_levels whether to allow new groups or levels to be passed to the model.
#' @param ... other arguments to pass to predict
#' @export
#' @details Predict responses based on the fitted model. Can be used
#' with the data used to fit the model to perform posterior
#' predictive checks, or with new data to generate new predictions.
#' These will include not only the parameter uncertainty but also
#' the variance that comes as a function of the likelihood. Because
#' residual error is incorporated, it can be used for visual model
#' criticism. For linear predictions (e.g. that only incorporate
#' )
predict.baggr <- function(x, newdata = NULL,
                          allow_new_levels = T, nsamples, ...) {
  do.call(paste0("predict_", x$model),
          args = list(
            x = x,
            newdata = newdata,
            allow_new_levels = allow_new_levels,
            nsamples = nsamples,
            ...)
          )
}

#' Make model matrix for the rubin data
#' @importFrom reshape2 dcast
#' @inheritParams predict.baggr
rubin_data <- function(x, newdata = NULL, allow_new_levels = T) {
  if(x$model != "rubin") {
    stop("Model must be type Rubin.")
  }

  group_label <- attr(x$inputs, "group_label")
  group_num <- 1:attr(x$inputs, "n_groups")
  group_col <- names(which(sapply(x$data,
                            function(x, labels) any(labels %in% x),
                            labels = group_label)))

  dat <- x$data

  other_cols <- setdiff(colnames(dat), group_col)

  if(is.null(newdata)) {
    dat <- x$data
  } else {
    dat <- newdata
  }

  dat[,group_col] <-
        factor(dat[,group_col],
               levels = as.character(group_label))


  if(allow_new_levels != T) {
    if(any(is.na(dat[,group_col]))) stop("Data contains new levels. If this behavior is desired, set allow_new_levels to TRUE.")
  }

  predmat <- matrix(nrow = nrow(dat), ncol = x$n_groups)

  for(i in 1:ncol(predmat)) {
    lvl <- as.integer(dat[i,group_col])

    predmat[lvl,i] <- 1
  }

  predmat[which(is.na(predmat))] <- 0
  cbind(1, predmat)

}

#' Predict function for the rubin model
#' @importFrom rstan extract
#' @param x model to predict from
#' @param newdata new data to predict not allowed for rubin model (yet)
#' @param allow_new_levels allow the predictive of new, unobserved groups
#' @param nsamples number of samples to predict
predict_rubin <- function(x,
                          newdata = NULL,
                          allow_new_levels = T,
                          nsamples,
                          ...) {

  if(!is.null(newdata)){
    stop("New data not currently allowed for the rubin model. ")
  }

  pred_data <- rubin_data(x, newdata)

  sigmas <- sapply(x$data$se, rep, times = nsamples)

  params <- rstan::extract(x$fit, c("tau","sigma_tau","tau_k","eta"))

  tau_k <- params$tau_k[1:nsamples, ]
  eta <- params$eta[1:nsamples, ]

  tau <- params$tau[1:nsamples, ]
  sigma_tau <- params$sigma_tau[1:nsamples, ]

  means <- cbind(tau, eta)

  # measurement effect variance + treatment effect variance
  # I think this was wrong before
  full_sigmas <- sqrt(sigmas^2 + sigma_tau^2)

  pred_means <- pred_data %*% t(means)

  epsilon <- rnorm(length(sigmas),
                   mean = 0,
                   full_sigmas)

  pp_dist <- pred_means + epsilon

  t(pp_dist)

}

#' Stop with informative error
stop_not_implemented <- function() {
  stop("Method not implemented.")
}

#' Get Y from various models
#' @param x model to get y
#' @details grabs the relevant Y variable
#' for use with posterior or prior predictive checks.
get_y <- function(x) {
  switch(x$model,
         rubin = "tau",
         mutau = "tau",
         quantiles = stop_not_implemented(),
         individual = stop_not_implemented()
         )
}

#' Posterior predictive checks for baggr model
#' @import bayesplot
#' @param x model fit by baggr()
#' @param type type of pp_check. For a list see \link{bayesplot::available_ppc}
#' @param nsamples number of samples to draw for the comparison
#' @details Posterior predictive checks are a form of model criticism
#' that compares the likelihood based simulations from the posterior
#' predictive distribution to observed data. Given a model fit by
#' baggr, you can compare kernel density estimates, histograms, and much
#' more. For full details, please refer to the bayesplot package.
#' @export
pp_check.baggr <- function(x, type = "dens_overlay", nsamples = 40) {
  pp_fun <- getFromNamespace(paste0("ppc_",type),ns = "bayesplot")
  model_type <- x$model
  y <- x$data[,get_y(x)]
  yrep <- predict(x, nsamples = nsamples)
  pp_fun(y, yrep)
}

