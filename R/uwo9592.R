globalVariables(c("x", "y"))
#' Create CR Plots for GAMLSS objects
#'
#' @param obj A `gamlss` object
#' @param term Name of the variable for the plot
#' @param plot Logical indicating whether the plot should be returned (`TRUE`) or just the data (`FALSE`)
#' @param points Should the points be plotted too (`TRUE`) or just the lines (`FALSE`)
#' @param what Which moment of the distribution should be considered
#' @param ... Other arguments, currently not implemented
#' @import ggplot2
#' @import gamlss
#' @importFrom stats formula get_all_vars median predict quantile residuals update
#' @importFrom utils capture.output
#' @export
cr.gamlss <- function(obj,
                      term,
                      plot=TRUE,
                      points = TRUE,
                      what = c("mu", "sigma", "nu", "tau"),
                      ...){
  wht = match.arg(what)
  trms <- predict(obj, type="terms", wht)
  cpr <- residuals(obj, type="partial", what=wht, term=term)
  tmp <- data.frame(x = trms[,term],
                    y = cpr[,term])
  attr(tmp, "term") <- term
  if(plot){
    g <- ggplot(tmp, aes(x=x, y=y)) +
      geom_smooth(aes(linetype="Linear"), method="lm", formula = y ~ x, se=FALSE, color="black", show.legend=FALSE) +
      geom_smooth(aes(linetype="LOESS"), method="loess", formula = y ~ x, se=FALSE, color="black", show.legend=FALSE) +
      scale_linetype_manual(values=c(3,1)) +
      theme_classic() +
      labs(x=term, y="Component + Residual", linetype="")
    if(points){
      g <- ggplot(tmp, aes(x=x, y=y)) +
        geom_point(shape=1, col="gray50") +
        geom_smooth(aes(linetype="Linear"), method="lm", formula = y ~ x, se=FALSE, color="black", show.legend=FALSE) +
        geom_smooth(aes(linetype="LOESS"), method="loess", formula = y ~ x, se=FALSE, color="black", show.legend=FALSE) +
        scale_linetype_manual(values=c(3,1)) +
        theme_classic() +
        labs(x=term, y="Component + Residual", linetype="")
    }
    g
  }else{
    return(tmp)
  }
}


#' Bootstraps (percentile) confidence intervals for effects from a GAMLSS model.
#'
#' @param obj A fitted `gamlss` object.
#' @param data The original data used to fit `obj`.
#' @param R Number of bootstrap iterations
#' @param what Which moment of the distribution to predict - only one allowed at a time.
#' @param terms Named list of variable names and values used to generate predictions. All
#' possible combinations of variables provided will be expanded.
#' @param progress Logical indicating whether a progress bar should be printed.
#' @param ... Currently not implemented.
#'
#' @details Uses nonparametric bootstrapping to calculate confidence intervals for effects
#' plots.  As the function refits the model for every bootstrap iteration, computation times
#' could be long. One of the main benefits of this approach is it also essentially takes into
#' account variability in the smoothing parameters of the penalized splines.  Use `term.plot` from
#' the `gamlss` package as a quick alternative.
#' @import dplyr
#' @importFrom progress progress_bar
#' @export
boot_eff.gamlss <- function(obj,
                            data,
                            R = 1000,
                            what = c("mu", "sigma", "nu", "tau"),
                            terms,
                            progress=TRUE,
                            ...){
  is_fac_char <- function(x){
    inherits(x, "factor") | inherits(x, "character")
  }
  modal <- function(x){
    if(is.factor(x)){
      levs <- levels(x)
      tab <- table(x)
      mlev <- names(tab)[which.max(tab)]
      factor(mlev, levels=levs)
    }else{
      tab <- table(x)
      mlev <- names(tab)[which.max(tab)]
    }
  }
  w <- match.arg(what)
  newdat <- origdat <- get_all_vars(formula(obj), data)
  newdat <- newdat %>%
    summarise(across(where(is_fac_char), modal),
              across(where(is.numeric), median))
  newdat <- as.list(newdat)
  check <- all(names(terms) %in% names(newdat))
  if(!check)stop("Not all terms in model\n")
  for(i in 1:length(terms)){
    newdat[[names(terms)[i]]] <- terms[[i]]
  }
  newdat <- do.call(expand.grid, newdat)
  preds <- NULL
  if(progress)pb <- progress_bar$new(total = R)
  for(i in 1:R){
    samp <- origdat[sample(1:nrow(origdat), nrow(origdat), replace=TRUE), ]
    capture.output(res <- update(obj, data=samp))
    capture.output(preds <- cbind(preds, predict(res, newdata=newdat, what=w, data = origdat)))
    if(progress)pb$tick()
  }
  out <- as.data.frame(newdat)[,names(terms), drop=FALSE]
  capture.output(out$estimate <- predict(obj, newdata=newdat, data=origdat, what=w))
  out$lwr <- apply(preds, 1, quantile, probs=.025)
  out$upr <- apply(preds, 1, quantile, probs=.975)
  return(out)
}


#' Bootstraps (percentile) confidence intervals for first differences from a GAMLSS model.
#'
#' @param obj A fitted `gamlss` object.
#' @param data The original data used to fit `obj`.
#' @param R Number of bootstrap iterations
#' @param what Which moment of the distribution to predict - only one allowed at a time.
#' @param terms Named list where the name corresponds to a variable in the model.  The
#' values of the list should be either a vector of two numbers or a function.  See details below.
#' @param progress Logical indicating whether a progress bar should be printed.
#' @param ... Currently not implemented.
#'
#' @details If `terms` is two values, then those two values will be used for all observations
#' and the average first difference will be calculated.  If the value is a character string,
#' that should refer to a function whose only required argument is `x` and that `x` should be the
#' variable named in the list.  In that case, the variable - the evaluated function and the variable
#' plust the evaluated function will be used.  For example if `terms = list(x = "sd")`, then
#' predictions will be made for `x-sd(x)` and `x+sd(x)`.
#' @export
boot_fd.gamlss <- function(obj,
                           data,
                           R = 1000,
                           what = c("mu", "sigma", "nu", "tau"),
                           terms,
                           progress=TRUE,
                           ...){
  w <- match.arg(what)
  newdat0 <- newdat1 <- origdat <- get_all_vars(formula(obj), data)
  check <- all(names(terms) %in% names(newdat0))
  if(!check)stop("Not all terms in model\n")
  if(length(terms[[1]]) == 1){
    newdat0[[names(terms)[1]]] <- newdat0[[names(terms)[1]]] - do.call(terms[[1]], list(x=newdat0[[names(terms)[1]]]))
    newdat1[[names(terms)[1]]] <- newdat1[[names(terms)[1]]] + do.call(terms[[1]], list(x=newdat1[[names(terms)[1]]]))
    vals <- c(paste0(names(terms)[1], " - ", terms[[1]]),
              paste0(names(terms)[1], " + ", terms[[1]]))
  }
  if(length(terms[[1]]) == 2){
    newdat0[[names(terms)[1]]] <- terms[[1]][1]
    newdat1[[names(terms)[1]]] <- terms[[1]][2]
    vals <- paste(names(terms)[1], terms[[1]], sep=" = ")
  }
  if(length(terms[[1]]) > 2){
    stop("terms must be a named list that has the name of a function or two values.\n")
  }
  preds0 <- preds1 <- NULL
  if(progress)pb <- progress_bar$new(total = R)
  for(i in 1:R){
    samp <- origdat[sample(1:nrow(origdat), nrow(origdat), replace=TRUE), ]
    capture.output(res <- update(obj, data=samp))
    capture.output(preds0 <- cbind(preds0, predict(res, newdata=newdat0, what=w, data = origdat)))
    capture.output(preds1 <- cbind(preds1, predict(res, newdata=newdat1, what=w, data = origdat)))
    if(progress)pb$tick()
  }
  capture.output(p0_hat <- predict(obj, newdata=newdat0, data=origdat, what=w))
  capture.output(p1_hat <- predict(obj, newdata=newdat1, data=origdat, what=w))
  p0bar <- colMeans(preds0)
  p1bar <- colMeans(preds1)
  diffbar <- colMeans(preds1-preds0)
  out <- data.frame(vals=c(vals, 'diff'),
                    estimate = c(mean(p0_hat), mean(p1_hat), mean(p1_hat-p0_hat)),
                    lwr = c(quantile(p0bar, .025), quantile(p1bar, .025), quantile(diffbar, .025)),
                    upr = c(quantile(p0bar, .975), quantile(p1bar, .975), quantile(diffbar, .975)),
                    mean = c(mean(p0bar), mean(p1bar), mean(diffbar)))
  return(out)
}


#' Brief summary of gamlss objects
#' 
#' @param obj A `gamlss` object
#' @param ... currently not implemented
#' @export
brief <- function(obj, ...){
  capture.output(s <- summary(obj))
  class(s) <- "brief"
  return(s)
}

#' Print method for brief summaries
#' 
#' @param x an object of class `brief`
#' @param ... currently not implemented
#' @method print brief
#' @importFrom stats printCoefmat
#' @export
print.brief <- function(x, ...){
  printCoefmat(x$result$coef.table, ...)
}

#' Tidy method for `comparisons` from `marginaleffects` package
#' 
#' @param x An object of class `comparisons`
#' @param ... currently not implemented
#' @method tidy comparisons
#' @importFrom broom tidy
#' @export
tidy.comparisons <- function(x, ...){
  x %>% select("term", "estimate", "std.error",
                   "p.value", "conf.low", "conf.high")
}

#' Plot layer behind existing layers in ggplot
#' 
#' @param plot A `ggplot`
#' @param layer A layer to be added to the plot
#' @importFrom ggplot2 is.ggplot
#' @export
`-.gg` <- function(plot, layer) {
  if (missing(layer)) {
    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
  }
  if (!is.ggplot(plot)) {
    stop('Need a plot on the left side')
  }
  plot$layers = c(layer, plot$layers)
  plot
}

#' Calculates Values from a Plot
#' 
#' @param plot A `ggplot` 
#' @param scale The fraction of the y-axis to identify.  For example, if `scale = .3`, then the function will return a value 30 percent of the way up the y-axis. 
#' @importFrom ggplot2 ggplot_build
#' @export
calc_ymax <- function(plot, scale=.3) {
  gb = ggplot_build(plot)
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  ymin + (ymax-ymin)*scale
}

#' Calculates Information Criterion Delta Values
#' 
#' @param ... Models whose IC values will be calculated
#' @importFrom AICcmodavg AICc
#' @importFrom stats AIC BIC
#' @export
IC_delta <- function(...){
  mods <- list(...)
  a <- AIC(...)
  a$AICc <- sapply(mods, AICc)
  a$BIC <- BIC(...)[,2]
  out <- a %>% 
    mutate(across(c(AIC, AICc, BIC), ~.x-min(.x)))
  names(out)[2:4] <- paste0("D_", names(out)[2:4])
  out
}


#' Creates a CR Plot using ggplot2
#' 
#' @param obj An object of class `glm`
#' @param term The model term to be plotted
#' @param plot Logical indicating whether a plot should be returned (`TRUE`) or just the data (`FALSE`)
#' @param res_type The type of residuals to be used. 
#' @param ... Currently not implemented
#' @importFrom statmod qresid
#' @importFrom stats rstandard
#' @export
gg_crplot <- function(obj, term, plot=TRUE, 
                      res_type = c("std_deviance", "std_pearson", "deviance", "pearson", "response", "working", "quantile"), 
                      ...){
  trms <- predict(obj, type="terms")
  restyp <- match.arg(res_type)
  if(!term %in% colnames(trms)){
    stop("term must be one of: ", paste(paste('"', colnames(trms), '"', sep=""), collapse=", "))
  }
  if(restyp == "std_deviance")e <- rstandard(obj, type="deviance")
  if(restyp == "std_pearson")e <- rstandard(obj, type="pearson")
  if(restyp == "deviance")e <- residuals(obj, type="deviance")
  if(restyp == "pearson")e <- residuals(obj, type="pearson")
  if(restyp == "response")e <- residuals(obj, type="response")
  if(restyp == "working")e <- residuals(obj, type="working")
  if(restyp == "quantile")e <- qresid(obj)
  
  cpr <- apply(trms, 2, \(x)x+e)
  tmp <- data.frame(x = trms[,term], 
                    y = cpr[,term])
  attr(tmp, "term") <- term
  attr(tmp, "res_type") <- restyp
  if(plot){
    ggplot(tmp, aes(x=x, y=y)) + 
      geom_point(shape=1, col="gray50") + 
      geom_smooth(aes(linetype="Linear"), method="lm", formula = y ~ x, se=FALSE, color="black") + 
      geom_smooth(aes(linetype="LOESS"), method="loess", formula = y ~ x, se=FALSE, color="black") + 
      scale_linetype_manual(values=c(3,1)) + 
      theme_classic() + 
      labs(x=term, y="Component + Residual", linetype="")
    
  }else{
    return(tmp)
  }
}


#' Make Hypothetical Prediction Data for mlogit models
#' 
#' @param obj An object of class `mlogit`
#' @param data The data frame used to fit `obj`
#' @param change A named list with a single element whose values are the values over which the named variable will be varied
#' @importFrom stats aggregate model.frame model.response
#' @export
makeFakeData <- function(obj, data, change){
  rl <- attr(data, "reshapeLong")
  idvar <- rl$idvar
  altvar <- rl$timevar
  varying <- names(rl$varying)
  tmpdata <- data[which(data[[idvar]] == 1), ]
  y <- model.response(model.frame(obj))
  vars <- all.vars(formula(obj))[-1]
  vars <- vars[-which(vars %in% names(change))]
  if(!is.null(varying)){
    vars <- vars[-which(vars %in% varying)]
  }
  if (any(!(vars %in% names(data)))) {
    vars <- vars[-which(!vars %in% names(data))]
  }
  var.classes <- sapply(vars, function(x) 
    case_when(inherits(data[[x]], "factor") ~ "factor", 
              inherits(data[[x]], "numeric") ~ "numeric", 
              TRUE ~ NA_character_))
  meds <- lapply(vars, function(x) NA)
  names(meds) <- vars
  levs <- list()
  for(i in 1:length(var.classes)){
    if(var.classes[i] == "factor"){
      levs[[vars[i]]] <- levels(data[[vars[i]]])
    }
  }
  if (length(levs) > 0) {
    for (i in 1:length(levs)) {
      tmp.tab <- table(data[[names(levs)[i]]])
      tmpdata[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)], 
                                          levels = levs[[i]])
      vars <- vars[-which(vars == names(levs)[i])]
    }
  }
  for (i in 1:length(vars)) {
    tmpdata[[vars[i]]] <- median(data[[vars[i]]], na.rm = T)
  }
  alt <- data[[altvar]]
  if(!is.null(varying)){
    for(i in 1:length(varying)){
      ag <- aggregate(data[[varying[i]]], list(alt), mean)
      tmpdata[[varying[i]]] <- ag[match(tmpdata[[altvar]], ag[,1]), 2]
    }
  }
  tmpdata2 <- data[which(data[[idvar]] %in% 1:length(change[[1]])), ]
  ch <- tmpdata2[[idvar]]
  for(i in 1:length(change[[1]])){
    tmpdata2[which(tmpdata2[[idvar]] == i), ] <- tmpdata
  }
  tmpdata2[[idvar]] <- ch
  for(i in 1:length(change[[1]])){
    tmpdata2[which(tmpdata2[[idvar]] == i), names(change)[1]] <- change[[1]][i]
  }
  attr(tmpdata2, "index") <- attr(tmpdata2, "index")[1:nrow(tmpdata2), ]
  tmpdata2
}


