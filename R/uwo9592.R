globalVariables(c("x", "y", "mu_fit", "mu_se", "obs", "vbl", "wn", "xbar"))
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


#' Brief summary of GAMLSS objects
#' 
#' @param object A `gamlss` object
#' @param ... currently not implemented
#' @importFrom car brief 
#' @importFrom purrr quietly
#' @method brief gamlss
#' @export
brief.gamlss <- function(object, ...){
  qsum <- quietly(summary)
  s <- qsum(object)
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
  x <- x$result
  ints <- which(rownames(x) == "(Intercept)")
  if(length(ints) == 1){
    l <- list(mu = x)
  }
  if(length(ints) == 2){
    l <- list(mu = x[1:(ints[2]-1), , drop=FALSE],
         sigma = x[ints[2]:nrow(x), , drop=FALSE])
  }
  if(length(ints) == 3){
    l <-list(mu = x[1:(ints[2]-1), , drop=FALSE],
         sigma = x[ints[2]:(ints[3]-1), , drop=FALSE], 
         nu = x[ints[3]:nrow(x), , drop=FALSE])
  }
  if(length(ints) == 4){
    l <-  list(mu = x[1:(ints[2]-1), , drop=FALSE],
         sigma = x[ints[2]:(ints[3]-1), , drop=FALSE], 
         nu = x[ints[3]:(ints[4]-1), , drop=FALSE], 
         tau = x[ints[4]:nrow(x), , drop=FALSE])
  }
  out <- lapply(l, printCoefmat)
  out
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

#' Generate Effect Plot Data
#' 
#' Generate linear predictor and standard error values from GAMLSS models
#' 
#' @param obj A `gamlss` object
#' @param data The data frame used to fit `obj`.
#' @param what The moment of the data to be predicted
#' @param ... Currently not implemented
#' @importFrom stats get_all_vars setNames
#' @importFrom tidyr pivot_longer
#' @export
tp_data <- function(obj, 
                    data, 
                    what = c("mu", "sigma", "nu", "tau"), 
                    ...){
  wht <- match.arg(what)
  dat <- get_all_vars(obj, data)[,-1]
  c_plot <- dat %>% 
    mutate(obs=row_number()) %>% 
    pivot_longer(-obs, names_to = "vbl", values_to = "x") 
  pred <- lpred(obj, wht, type="terms", se.fit=TRUE)
  fit <- pred$fit %>% 
    as.data.frame() %>% 
    mutate(obs = row_number()) %>% 
    setNames(c(names(dat), "obs")) %>% 
    pivot_longer(-obs, names_to = "vbl", values_to = "mu_fit") 
  se <- pred$se.fit %>% 
    as.data.frame() %>% 
    mutate(obs = row_number()) %>% 
    setNames(c(names(dat), "obs")) %>% 
    pivot_longer(-obs, names_to = "vbl", values_to = "mu_se") 
  c_plot$mu_fit <- mu_fit$mu_fit
  c_plot$mu_se <- mu_se$mu_se
  c_plot %>% 
    group_by(vbl) %>% 
    arrange(x) %>% 
    ungroup() %>% 
    group_by(vbl, x) %>% 
    slice_head(n=1)
  
}


#' Within Group Transformation Function
#' 
#' Performs within transformation on `x` within the groups of `id`.  If `x` is numeric, 
#' then `means` should either be `NULL` or a named vector of means with the names of the
#' elements have corresponding elements in `id`.  If `x` is a factor, then `means` should
#' be an `id` by levels of `x` matrix giving the proportion of observations in each category
#' of `x` for each different value of `id`.  If `means` is `NULL`, it will be calculated 
#' internally and then passed to the function. 
#' @param x A variable whose within transformation you want to calculate. This should be a vector of values, not a variable name. 
#' @param id The grouping variable within which the transformation be performed. 
#' @param means An optional vector of matrix of means depending on the class of `x`. 
#' @param ... Currently not implemented. 
#' @export
win <- function(x, id, means=NULL, ...){
  if(is.null(means)){
    tmp <- wint(x, id)
    means <- attr(tmp, "means")
    win(x, id, means=means)
  }else{
    wint(x, id, means=means)
  }
}

#' Performs within transformation 
#' 
#' @param x A variable whose within transformation you want to calculate. This should be a vector of values, not a variable name. 
#' @param id The grouping variable within which the transformation be performed. 
#' @param means An optional vector of matrix of means depending on the class of `x`. 
#' @param ... Currently not implemented. 
wint <- function(x, id, means=NULL, ...){
  UseMethod("wint")
}


#' Numeric method for `wint`
#' @param x A variable whose within transformation you want to calculate. This should be a vector of values, not a variable name. 
#' @param id The grouping variable within which the transformation be performed. 
#' @param means An optional vector of matrix of means depending on the class of `x`. 
#' @param ... Currently not implemented. 
#' @method wint numeric
wint.numeric <- function(x, id, means=NULL, ...){
  if(is.factor(id))id <- droplevels(id)
  tmp <- data.frame(x=x, id=id)
  if(!is.null(means)){
    mns <- data.frame(xbar = means, id = names(means))
    if(is.numeric(id)){
      mns$id <- as.numeric(mns$id)
    }
    mvec <- mns$xbar
    names(mvec) <- mns$id
    tmp <- left_join(tmp, mns, by = join_by(id))
    tmp <- tmp %>% mutate(wn = x - xbar)
  }else{
    tmp <- tmp %>% group_by(id) %>% mutate(xbar = mean(x)) %>% ungroup() %>% mutate(wn = x-xbar)
    mns <- tmp %>% group_by(id) %>% slice_head(n=1)
    mvec <- mns$xbar
    names(mvec) <- mns$id
  }
  out <- tmp %>% select(wn) %>% as.matrix()
  attr(out, "means") <- mvec
  attr(out, "class") <- "win"
  out
  #structure(out, class=c("win", "matrix"))
}
#' Factor method for `wint`
#' @param x A variable whose within transformation you want to calculate. This should be a vector of values, not a variable name. 
#' @param id The grouping variable within which the transformation be performed. 
#' @param means An optional vector of matrix of means depending on the class of `x`. 
#' @param ... Currently not implemented. 
#' @importFrom stats lm model.matrix
#' @method wint factor
wint.factor <- function(x, id, means=NULL, ...){
  if(is.factor(id))id <- droplevels(id)
  X <- model.matrix(~x-1)
  if(!is.null(means)){
    tmp <- data.frame(id=id)
    Xb <- left_join(tmp, as_tibble(means, rownames= "id"), by=join_by(id))
    Xb <- Xb %>% select(-id) %>% as.matrix()
    Xw <- X - Xb[,colnames(X), drop=FALSE]
    Xw <- Xw[,-1, drop=FALSE]
    mns <- means
  }else{
    aux <- lm(X ~ id)
    Xw <- aux$residuals[,-1, drop=FALSE]
    Xb <- aux$fitted[, , drop=FALSE]
    mns <- by(X, list(id), colMeans)
    mns <- do.call(rbind, mns)
  }
  out <- Xw
  attr(out, "means") <- mns
  attr(out, "class") <- "win"
  out
  #structure(out, class=c("win", "matrix"))
}

#' Predict method for win
#' 
#' @param object An object of class `win`
#' @param newdata An optional data frame with which to generate the predictions
#' @param ... Other arguments to be passed down
#' @method predict win
predict.win <- function(object, newdata, ...){
  if(missing(newdata))
    object
  else win(newdata, means = attr(object, "means"))
}

#' Makepredictcall method for win
#' @param var A variable
#' @param call A term in the formula, as a call. 
#' @method makepredictcall win
makepredictcall.win <- function(var, call){
  if (as.character(call)[1L] == "win" || (is.call(call) && 
                                          identical(eval(call[[1L]]), win))) 
    call$means <- attr(var, "means")
  call
}

#' Between Group Transformation Function
#' 
#' Performs between transformation on `x` within the groups of `id`.  If `x` is numeric, 
#' then `means` should either be `NULL` or a named vector of means with the names of the
#' elements have corresponding elements in `id`.  If `x` is a factor, then `means` should
#' be an `id` by levels of `x` matrix giving the proportion of observations in each category
#' of `x` for each different value of `id`.  If `means` is `NULL`, it will be calculated 
#' internally and then passed to the function. 
#' @param x A variable whose within transformation you want to calculate. This should be a vector of values, not a variable name. 
#' @param id The grouping variable within which the transformation be performed. 
#' @param means An optional vector of matrix of means depending on the class of `x`. 
#' @param ... Currently not implemented. 
#' @export
bwn <- function(x, id, means=NULL, ...){
  if(is.null(means)){
    tmp <- bwnt(x, id)
    means <- attr(tmp, "means")
    bwn(x, id, means=means)
  }else{
    bwnt(x, id, means=means)
  }
  
}

#' Performs Between Transformation
#' 
#' @param x A variable whose between transformation you want to calculate. This should be a vector of values, not a variable name. 
#' @param id The grouping variable between which the transformation be performed. 
#' @param means An optional vector of matrix of means depending on the class of `x`. 
#' @param ... Currently not implemented. 
#' @export
bwnt <- function(x, id, means=NULL, ...){
  UseMethod("bwnt")
}

#' Numeric method for `bwnt`
#' @param x A variable whose between transformation you want to calculate. This should be a vector of values, not a variable name. 
#' @param id The grouping variable between which the transformation be performed. 
#' @param means An optional vector of matrix of means depending on the class of `x`. 
#' @param ... Currently not implemented. 
#' @method wint numeric
bwnt.numeric <- function(x, id, means=NULL, ...){
  if(is.factor(id))id <- droplevels(id)
  tmp <- data.frame(x=x, id=id)
  if(!is.null(means)){
    mns <- data.frame(bwn = means, id = names(means))
    if(is.numeric(id)){
      mns$id <- as.numeric(mns$id)
    }
    mvec <- mns$bwn
    names(mvec) <- mns$id
    tmp <- left_join(tmp, mns, by = join_by(id))
  }else{
    tmp <- tmp %>% group_by(id) %>% mutate(bwn = mean(x)) %>% ungroup() 
    mns <- tmp %>% group_by(id) %>% slice_head(n=1)
    mvec <- mns$bwn
    names(mvec) <- mns$id
  }
  out <- tmp %>% select(bwn) %>% as.matrix()
  attr(out, "means") <- mvec
  attr(out, "class") <- "bwn"
  out
  #structure(out, class=c("bwn", "matrix"))
}

#' Factor method for `bwnt`
#' @param x A variable whose between transformation you want to calculate. This should be a vector of values, not a variable name. 
#' @param id The grouping variable between which the transformation be performed. 
#' @param means An optional vector of matrix of means depending on the class of `x`. 
#' @param ... Currently not implemented. 
#' @method wint numeric
bwnt.factor <- function(x, id, means=NULL, ...){
  if(is.factor(id))id <- droplevels(id)
  X <- model.matrix(~x-1)
  if(!is.null(means)){
    tmp <- data.frame(id=id)
    Xb <- left_join(tmp, as_tibble(means, rownames= "id"), by=join_by(id))
    Xb <- Xb %>% select(-id) %>% as.matrix()
    Xb <- Xb[,colnames(X), drop=FALSE]
    mns <- means
  }else{
    aux <- lm(X ~ id)
    Xb <- aux$fitted[, , drop=FALSE]
    mns <- by(Xb, list(id), colMeans)
    mns <- do.call(rbind, mns)
  }
  out <- Xb
  attr(out, "means") <- mns
  attr(out, "class") <- "bwn"
  out
  #structure(out, class=c("bwn", "matrix"))
}

#' Predict method for bwn
#' 
#' @param object An object of class `bwn`
#' @param newdata An optional data frame with which to generate the predictions
#' @param ... Other arguments to be passed down
#' @method predict bwn
predict.bwn <- function(object, newdata, ...){
  if(missing(newdata))
    object
  else bwn(newdata, means = attr(object, "means"))
}

#' Makepredictcall method for bwn
#' @param var A variable
#' @param call A term in the formula, as a call. 
#' @method makepredictcall bwn
makepredictcall.bwn <- function(var, call){
  if (as.character(call)[1L] == "bwn" || (is.call(call) && 
                                          identical(eval(call[[1L]]), bwn))) 
    call$means <- attr(var, "means")
  call
}



