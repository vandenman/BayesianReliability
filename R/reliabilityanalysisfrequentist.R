reliabilityFrequentist <- function(jaspResults, dataset, options) {

    # sink("~/Downloads/log_freq.txt")
    # on.exit(sink(NULL))

  dataset <- .frequentistReliabilityReadData(dataset, options)
  
  .frequentistReliabilityCheckErrors(dataset, options)
  
  model <- .frequentistReliabilityMainResults(jaspResults, dataset, options)
  
  .frequentistReliabilityScaleTable(         jaspResults, model, options)
  .frequentistReliabilityItemTable(          jaspResults, model, options)
  .freqentistReliabilitySingleFactorFitTable(jaspResults, model, options)
  return()
  
}

# read data, check errors----
.frequentistReliabilityDerivedOptions <- function(options) {
  
  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimatorsF  = unlist(options[c("mcDonaldScalef","alphaScalef", "guttman2Scalef", "guttman6Scalef", 
                                            "glbScalef", "averageInterItemCor", "meanScale", "sdScale")]),
    itemDroppedSelectedF = unlist(options[c("mcDonaldItemf", "alphaItemf", "guttman2Itemf", "guttman6Itemf",
                                            "glbItemf", "itemRestCor", "meanItem", "sdItem")]),
    namesEstimators     = list(
      tables = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6", 
                 "Greatest Lower Bound", "Average interitem correlation", "mean", "sd"),
      tables_item = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6", 
                      "Greatest Lower Bound", "Item-rest correlation", "mean", "sd"),
      coefficients = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6", 
                       "Greatest Lower Bound", "Item-rest correlation"),
      plots = list(expression("McDonald's"~omega), expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]), 
                   expression("Guttman's"~lambda[6]), "Greatest Lower Bound")
    )
  )
  # order to show in JASP
  derivedOptions[["order"]] <- c(5, 1, 2, 3, 4, 6, 7, 8)
  
  return(derivedOptions)
}

.frequentistReliabilityReadData <- function(dataset, options) {
  
  variables <- unlist(options[["variables"]])
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.numeric = variables, columns.as.factor = NULL, exclude.na.listwise = NULL)
  }
  return(dataset)
}

.frequentistReliabilityCheckErrors <- function(dataset, options) {
  
  .hasErrors(dataset = dataset, perform = "run",
             type = c("infinity", "variance", "observations"),
             observations.amount = " < 3",
             exitAnalysisIfErrors = TRUE)
  
}

# estimate reliability ----
# maybe in the future it would be easier to have one function for every estimator...
.frequentistReliabilityMainResults <- function(jaspResults, dataset, options) {
  samples <- options[["noSamplesf"]]
  p <- ncol(dataset)
  
  model <- jaspResults[["modelObj"]]$object
  relyFit <- model[["relyFit"]]
  if (is.null(model)) {
    
    model <- list()
    variables <- options[["variables"]]
    if (length(variables) > 2L) {
      
      dataset <- as.matrix(dataset) # fails for string factors!
      if (length(options[["reverseScaledItems"]]) > 0L) {
        # nvar <- length(variables)
        # key <- rep(1, nvar)
        # key[match(.v(unlist(options[["reverseScaledItems"]])), nvar)] <- -1
        # dataset <- dataset %*% diag(key, nvar, nvar) # this seems like it does not work
        cols <- match(unlist(options[["reverseScaledItems"]]), .unv(colnames(dataset)))
        total <- min(dataset, na.rm = T) + max(dataset, na.rm = T)
        dataset[ ,cols] = total - dataset[ ,cols]
      }

      # observations for alpha interval need to be speccified: 
      model[["obs"]] <- nrow(dataset)
      
      model[["footnote"]] <- .frequentistReliabilityCheckLoadings(dataset, variables)
      
      if (options[["missingValuesf"]] == "excludeCasesPairwise") {
        missing <- "pairwise"
        use.cases <- "pairwise.complete.obs"
        model[["footnote"]] <- paste0(model[["footnote"]], ". Using pairwise complete cases.")
      } else {
        pos <- which(is.na(dataset), arr.ind = T)[, 1]
        dataset <- dataset[-pos, ] 
        use.cases <- "complete.obs"
        model[["footnote"]] <- paste0(model[["footnote"]], ". Using ", nrow(dataset), 
                                      " complete cases.")
      }
      
      if (options[["alphaInterval"]] == "alphaAnalytic") {
        alphaAna <- TRUE
        alphaSteps <- 0
      } else {
        alphaAna <- FALSE
        alphaSteps <- samples
        if (options[["alphaMethod"]] == "alphaStand") {
          alphaSteps <- alphaSteps + samples
        }
      }
      
      if (options[["omegaEst"]] == "pfa") {
        omegaSteps <- samples
        omegaAna <- TRUE
      } else {
        if (options[["omegaInterval"]] == "omegaAnalytic") {
          omegaAna <- TRUE
          omegaSteps <- 0
        } else {
          omegaAna <- FALSE
          omegaSteps <- 0
          # omegaSteps <- samples  # not working because we cant tick in the bootstrapLavaan function
        }
      }

      if (options[["bootType"]] == "bootNonpara") {
        para <- FALSE
      } else {
        para <- TRUE
      }
      
      startProgressbar(samples * 5 # cov_mat bootstrapping and coefficients (also avg_cor) without alpha and omega
                       + alphaSteps
                       + omegaSteps) # dont need ifitem steps since that is very fast

      if (options[["alphaMethod"]] == "alphaStand") {
        model[["dat_cov"]] <- Bayesrel:::make_symmetric(cov2cor(cov(dataset, use = use.cases)))
        relyFit <- try(Bayesrel::strel(data = dataset, estimates=c("lambda2", "lambda6", "glb", "omega"), 
                                       Bayes = FALSE, n.boot = options[["noSamplesf"]],
                                       item.dropped = TRUE, omega.freq.method = options[["omegaEst"]], 
                                       omega.int.analytic = omegaAna,
                                       para.boot = para,
                                       missing = missing, callback = progressbarTick))
        
        relyFit$freq$est$freq_alpha <- Bayesrel:::applyalpha(model[["dat_cov"]])
        p <- ncol(dataset)
        Ctmp <- array(0, c(p, p - 1, p - 1))
        for (i in 1:p){
          Ctmp[i, , ] <- model[["dat_cov"]][-i, -i]
        }
        relyFit$freq$ifitem$alpha <- apply(Ctmp, 1, Bayesrel:::applyalpha)
        
        if (!alphaAna) { # when standardized alpha, but bootstrapped alpha interval:
          cors <- array(0, c(options[["noSamplesf"]], p, p))
          for (i in 1:options[["noSamplesf"]]) {
            cors[i, , ] <- .cov2cor.callback(relyFit$freq$covsamp[i, , ], progressbarTick)
          }
          relyFit$freq$boot$alpha <- apply(cors, 1, Bayesrel:::applyalpha)
          if (omegaAna) {
            relyFit[["freq"]][["boot"]] <- relyFit[["freq"]][["boot"]][c(4, 1, 2, 3)]
          } else {
            relyFit[["freq"]][["boot"]] <- relyFit[["freq"]][["boot"]][c(5, 1, 2, 3, 4)]
          }
          
        }
        
        relyFit[["freq"]][["est"]] <- relyFit[["freq"]][["est"]][c(5, 1, 2, 3, 4)]
        relyFit[["freq"]][["ifitem"]] <- relyFit[["freq"]][["ifitem"]][c(5, 1, 2, 3, 4)]
        
      } else {
        model[["dat_cov"]] <- Bayesrel:::make_symmetric(cov(dataset, use = use.cases))
        relyFit <- try(Bayesrel::strel(data = dataset, estimates=c("alpha", "lambda2", "lambda6", "glb", "omega"), 
                                       Bayes = FALSE, n.boot = options[["noSamplesf"]],
                                       item.dropped = TRUE, omega.freq.method = options[["omegaEst"]], 
                                       alpha.int.analytic = alphaAna, 
                                       omega.int.analytic = omegaAna,
                                       missing = missing, callback = progressbarTick))
      }


      
      if (!is.null(relyFit[["freq"]][["omega.error"]])) {
        model[["footnote"]] <- paste0(model[["footnote"]], " omega estimation method switched to PFA because the CFA
                                      did not find a solution.")
      }
      
      
      # first the scale statistics
      cordat <- cor(dataset, use = use.cases)
      relyFit$freq$est$avg_cor <- mean(cordat[lower.tri(cordat)])
      relyFit$freq$est$mean <- mean(rowMeans(dataset, na.rm = T))
      relyFit$freq$est$sd <- sd(colMeans(dataset, na.rm = T))

      corsamp <- apply(relyFit$freq$covsamp, c(1), .cov2cor.callback, progressbarTick)
      relyFit$freq$boot$avg_cor <- apply(corsamp, 2, function(x) mean(x[x!=1]))
      relyFit$freq$boot$mean <- c(NA_real_, NA_real_)
      relyFit$freq$boot$sd <- c(NA_real_, NA_real_)

      
      # now the item statistics
      relyFit$freq$ifitem$ircor <- NULL
      
      for (i in 1:ncol(dataset)) {
        relyFit$freq$ifitem$ircor[i] <- cor(dataset[, i], rowMeans(dataset[, -i], na.rm = T), use = use.cases)
      }
      relyFit$freq$ifitem$mean <- colMeans(dataset, na.rm = T)
      relyFit$freq$ifitem$sd <- apply(dataset, 2, sd, na.rm = T)  

      ops <- .BayesianReliabilityDerivedOptions(options)
      order <- ops[["order"]]
      # relyFit[["freq"]][["boot"]] <- relyFit[["freq"]][["boot"]][order]
      relyFit[["freq"]][["est"]] <- relyFit[["freq"]][["est"]][order]
      relyFit[["freq"]][["ifitem"]] <- relyFit[["freq"]][["ifitem"]][order]
      
      # Consider stripping some of the contents of relyFit to reduce memory load
      if (inherits(relyFit, "try-error")) {
        
        model[["error"]] <- paste("The analysis crashed with the following error message:\n", relyFit)
        
      } else {
        
        model[["dataset"]] <- dataset
        
        model[["relyFit"]] <- relyFit
        
        stateObj <- createJaspState(model)
        stateObj$dependOn(options = c("variables", "reverseScaledItems", "noSamplesf", "missingValuesf", "omegaEst", 
                                      "alphaMethod", "alphaInterval", "omegaInterval", "bootType"))
        jaspResults[["modelObj"]] <- stateObj

      }
    }
  } else {
    print("model from state")
  }
  
  if (is.null(model[["error"]])) {
    cfiState <- jaspResults[["cfiObj"]]$object
    if (is.null(cfiState) && !is.null(relyFit)) {
      
      scaleCfi <- .frequentistReliabilityCalcCfi(relyFit[["freq"]][["boot"]],             
                                              options[["confidenceIntervalValue"]])

      # alpha int is analytical, not from the boot sample, so:
      if (options[["alphaInterval"]] == "alphaAnalytic") {

        alphaCfi <- Bayesrel:::ciAlpha(1 - options[["confidenceIntervalValue"]], model[["obs"]], model[["dat_cov"]])
        names(alphaCfi) <- c("lower", "upper")
        scaleCfi$alpha <- alphaCfi
      }

      
      # omega cfa analytic interval:
      if (is.null(relyFit[["freq"]][["omega.pfa"]]) & (options[["omegaInterval"]] == "omegaAnalytic")) { 
        fit <- relyFit[["freq"]][["fit.object"]]
        params <- lavaan::parameterestimates(fit, level = options[["confidenceIntervalValue"]])
        om_low <- params$ci.lower[params$lhs=="omega"]
        om_up <- params$ci.upper[params$lhs=="omega"]
        omegaCfi <- c(om_low, om_up)
        names(omegaCfi) <- c("lower", "upper")
        scaleCfi$omega <- omegaCfi
        if (options[["alphaInterval"]] == "alphaAnalytic") {
          scaleCfi <- scaleCfi[c(8, 7, 1, 2, 3, 4, 5, 6)] # check this when more estimators come in
        } else {
          scaleCfi <- scaleCfi[c(8, 1, 2, 3, 4, 5, 6, 7)] # check this when more estimators come in
        }
      } else {
        if (options[["alphaInterval"]] == "alphaAnalytic") {
          scaleCfi <- scaleCfi[c(4, 8, 1, 2, 3, 5, 6, 7)] # check this when more estimators come in
        } else {
          scaleCfi <- scaleCfi[c(5, 1, 2, 3, 4, 6, 7, 8)] # check this when more estimators come in
        }
      }

      cfiState <- list(scaleCfi = scaleCfi)
      jaspCfiState <- createJaspState(cfiState)
      jaspCfiState$dependOn(options = "confidenceIntervalValue", optionsFromObject = jaspResults[["modelObj"]])
      jaspResults[["cfiObj"]] <- jaspCfiState
      
    }
    model[["cfi"]] <- cfiState
  }
  
  model[["derivedOptions"]] <- .frequentistReliabilityDerivedOptions(options)
  model[["itemsDropped"]] <- .unv(colnames(dataset))
  progressbarTick()
  
  return(model)
}

.frequentistReliabilityCalcCfi <- function(boot, cfiValue) {
  
  cfi <- vector("list", length(boot))
  names(cfi) <- names(boot)
  
  for (nm in names(boot)) {
    if (any(is.na(boot[[nm]]))) {
      cfi[[nm]] <- c(NA_real_, NA_real_)
    } else {
      cfi[[nm]] <- quantile(boot[[nm]], prob = c((1-cfiValue)/2, 1-(1-cfiValue)/2))
    }
    names(cfi[[nm]]) <- c("lower", "upper")
  }
  return(cfi)
}

.frequentistReliabilityCheckLoadings <- function(dataset, variables) {
  
  prin <- psych::principal(dataset)
  idx <- prin[["loadings"]] < 0
  sidx <- sum(idx)
  hasSchar <- if (sidx == 1L) "" else "s"
  footnote <- sprintf("The following item%s correlated negatively with the scale: %s",
                      hasSchar, paste0(variables[idx], collapse = ", "))
  return(footnote)
}


# tables ----


.frequentistReliabilityScaleTable <- function(jaspResults, model, options) {
  
  if (!is.null(jaspResults[["scaleTableF"]])) {
    print("jaspResults[['scaleTableF']] from state")
    return()
  }
  
  scaleTableF <- createJaspTable("Frequentist Scale Reliability Statistics")
  scaleTableF$dependOn(options = c("variables", "mcDonaldScalef", "alphaScalef", "guttman2Scalef", "guttman6Scalef",
                                   "glbScalef", "reverseScaledItems", "confidenceIntervalValue", "noSamplesf", 
                                   "averageInterItemCor", "meanScale", "sdScale", "missingValuesf", "omegaEst", 
                                   "alphaMethod", "alphaInterval", "omegaInterval", "bootType"))
  
  overTitle <- sprintf("%s%% Confidence interval",
                       format(100*options[["confidenceIntervalValue"]], digits = 3, drop0trailing = TRUE))
  scaleTableF$addColumnInfo(name = "statistic", title = "Statistic",        type = "string")
  scaleTableF$addColumnInfo(name = "pointEst",  title = "Point Estimate",   type = "number")
  scaleTableF$addColumnInfo(name = "lower",     title = "Lower",            type = "number", overtitle = overTitle)
  scaleTableF$addColumnInfo(name = "upper",     title = "Upper",            type = "number", overtitle = overTitle)
  
  relyFit <- model[["relyFit"]]
  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  order    <- derivedOptions[["order"]]
  selected <- derivedOptions[["selectedEstimatorsF"]]
  idxSelected <- which(selected)
  
  if (!is.null(relyFit)) {
    allData <- cbind.data.frame(
      statistic = opts,
      pointEst = unlist(relyFit$freq$est, use.names = FALSE),
      do.call(rbind, model[["cfi"]][["scaleCfi"]]))[idxSelected, ]
    
    scaleTableF$setData(allData)
    
    if (!is.null(model[["footnote"]]))
      scaleTableF$addFootnote(model[["footnote"]])
    
  } else if (sum(selected) > 0L) {
    
    scaleTableF[["statistic"]] <- opts[idxSelected]
    
    nvar <- length(options[["variables"]])
    if (nvar > 0L && nvar < 3L)
      scaleTableF$addFootnote("Please enter at least 3 variables to do an analysis.")
  }
  if (!is.null(model[["error"]]))
    scaleTableF$setError(model[["error"]])
  
  if (!is.null(model[["footnote"]]))
    scaleTableF$addFootnote(model[["footnote"]])
  
  jaspResults[["scaleTableF"]] <- scaleTableF
  
  return()
}



.frequentistReliabilityItemTable <- function(jaspResults, model, options) {
  
  if (!is.null(jaspResults[["itemTableF"]]) || !any(model[["derivedOptions"]][["itemDroppedSelectedF"]])) {
    print("jaspResults[['hasItemTableF']] from state or not wanted")
    return()
  }
  
  derivedOptions <- model[["derivedOptions"]]
  itemDroppedSelectedF <- derivedOptions[["itemDroppedSelectedF"]]
  # order <- derivedOptions[["order_item"]]
  estimators <- derivedOptions[["namesEstimators"]][["tables_item"]]
  overTitle <- "If item dropped"
  
  itemTableF <- createJaspTable("Frequentist Individual Item Reliability Statistics")
  itemTableF$dependOn(options = c("variables",
                                  "mcDonaldScalef", "alphaScalef", "guttman2Scalef", "guttman6Scalef", "glbScalef", 
                                  "averageInterItemCor", "meanScale", "sdScale",
                                  "mcDonaldItemf",  "alphaItemf",  "guttman2Itemf", "guttman6Itemf", "glbItemf",
                                  "reverseScaledItems", "meanItem", "sdItem", "itemRestCor", "missingValuesf", 
                                  "omegaEst", "alphaMethod"))
  itemTableF$addColumnInfo(name = "variable", title = "Item", type = "string")
  
  idxSelectedF <- which(itemDroppedSelectedF)
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]
  for (i in idxSelectedF) {
    if (estimators[i] %in% coefficients) {
      itemTableF$addColumnInfo(name = paste0("pointEst", i), title = estimators[i], type = "number", 
                               overtitle = overTitle)
    } else {
      itemTableF$addColumnInfo(name = paste0("pointEst", i), title = estimators[i], type = "number")
    }
  }
  
  relyFit <- model[["relyFit"]]
  if (!is.null(relyFit)) {
    tb <- data.frame(variable = model[["itemsDropped"]])
    for (i in idxSelectedF) {
      newtb <- cbind(pointEst = relyFit$freq$ifitem[[i]])
      colnames(newtb) <- paste0(colnames(newtb), i)
      tb <- cbind(tb, newtb)
      
    }
    itemTableF$setData(tb)
    
  } else if (length(model[["itemsDropped"]]) > 0) {
    itemTableF[["variables"]] <- model[["itemsDropped"]]
  }
  
  jaspResults[["itemTableF"]] <- itemTableF
  return()
}

# once the package is updated check this again and apply:
.freqentistReliabilitySingleFactorFitTable <- function(jaspResults, model, options) {

  if (!options[["fitMeasures"]] || !options[["mcDonaldScalef"]])
    return()
  if (!is.null(jaspResults[["fitTable"]]) || !options[["fitMeasures"]]) {
    print("jaspResults[['fitTable']] from state or not wanted")
    return()
  }

  fitTable <- createJaspTable(sprintf("Fit Measures of Single Factor Model Fit"))
  fitTable$dependOn(options = c("variables", "mcDonaldScalef", "reverseScaledItems", "fitMeasures", "missingValuesf", 
                                "omegaEst"))
  fitTable$addColumnInfo(name = "measure", title = "Fit Measure",   type = "string")
  fitTable$addColumnInfo(name = "value",     title = "Value", type = "number")

  relyFit <- model[["relyFit"]]
  derivedOptions <- model[["derivedOptions"]]
  opts     <- names(relyFit$freq$omega_fit)

  if (!is.null(relyFit)) {
    if (is.null(opts)) {
      allData <- data.frame(
        measure = NA_real_,
        value = NA_real_
      )
      if (!is.null(relyFit[["freq"]][["omega.error"]])) {
          fitTable$addFootnote("Fit measures cannot be displayed because the omega estimation method switched 
                                to PFA as the CFA did not find a solution.")
      }
    } else {
      allData <- data.frame(
        measure = opts,
        value = as.vector(unlist(relyFit$freq$omega_fit, use.names = FALSE))
      )
    }

    fitTable$setData(allData)
  }
  if (!is.null(model[["error"]]))
    fitTable$setError(model[["error"]])


  jaspResults[["fitTable"]] <- fitTable
}



# get bootstrapped sample for omega with cfa
.applyomega_cfa_cov <- function(cov, n){
  data <- MASS::mvrnorm(n, numeric(ncol(cov)), cov)
  out <- Bayesrel:::omegaFreqData(data, pairwise = F)
  om <- out$omega
  return(om)
}
