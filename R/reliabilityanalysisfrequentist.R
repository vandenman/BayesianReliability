reliabilityFrequentist <- function(jaspResults, dataset, options) {
  sink("~/Downloads/log_freq.txt")
  on.exit(sink(NULL))
  
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
    selectedEstimatorsF  = unlist(options[c("alphaScalef", "guttman2Scalef", "guttman6Scalef", "glbScalef", 
                                            "mcDonaldScalef", "averageInterItemCor", "meanScale", "sdScale")]),
    itemDroppedSelectedF = unlist(options[c("mcDonaldItemf", "alphaItemf", "guttman2Itemf", "guttman6Itemf",
                                            "glbItemf", "itemRestCor", "meanItem", "sdItem")]),
    namesEstimators     = list(
      tables = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6", "Greatest Lower Bound", 
                 "McDonald's \u03C9", "Average interitem correlation", "mean", "sd"),
      tables_item = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6", "Greatest Lower Bound", 
                 "McDonald's \u03C9", "Item-rest correlation", "mean", "sd"),
      coefficients = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6", "Greatest Lower Bound", 
                       "McDonald's \u03C9"),
      plots = list(expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]), 
                   expression("Guttman's"~lambda[6]), "Greatest Lower Bound", expression("McDonald's"~omega))
    )
  )
  # order to show in JASP
  derivedOptions[["order"]] <- match(c("McDonald's \u03C9", "Cronbach's \u03B1", 
                                       "Guttman's \u03BB2", "Guttman's \u03BB6", "Greatest Lower Bound",
                                       "Average interitem correlation", "mean", "sd"), 
                                     derivedOptions[["namesEstimators"]][["tables"]])
  derivedOptions[["order_item"]] <- match(c("McDonald's \u03C9", "Cronbach's \u03B1", 
                                       "Guttman's \u03BB2", "Guttman's \u03BB6", "Greatest Lower Bound",
                                       "Item-rest correlation", "mean", "sd"), 
                                     derivedOptions[["namesEstimators"]][["tables_item"]])
  
  
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
.frequentistReliabilityMainResults <- function(jaspResults, dataset, options) {
  
  model <- jaspResults[["modelObj"]]$object
  relyFit <- model[["relyFit"]]
  if (is.null(model)) {
    
    model <- list()
    variables <- options[["variables"]]
    if (length(variables) > 2L) {
      
      dataset <- as.matrix(dataset) # fails for string factors!
      if (length(options[["reverseScaledItems"]]) > 0L) {
        nvar <- length(variables)
        key <- rep(1, nvar)
        key[match(.v(unlist(options[["reverseScaledItems"]])), nvar)] <- -1
        dataset <- dataset %*% diag(key, nvar, nvar)
      }
      
      if (options[["missingValuesF"]] == "excludeCasesPairwise") {missing <- "pairwise"}
      else if (options[["missingValuesF"]] == "excludeCasesListwise") {missing <- "listwise"}
      
      model[["footnote"]] <- .frequentistReliabilityCheckLoadings(dataset, variables)
      relyFit <- try(Bayesrel::strel(x = dataset, estimates=c("alpha", "lambda2", "lambda6", "glb", "omega"), 
                                     Bayes = FALSE, n.boot = options[["noSamplesf"]],
                                     item.dropped = TRUE, omega.freq.method = "cfa", 
                                     alpha.int.analytic = TRUE, 
                                     missing = missing))
      
      if (any(is.na(dataset))) {
        if (!is.null(relyFit[["miss_pairwise"]])) {
          model[["footnote"]] <- paste0(model[["footnote"]], ". Using pairwise complete cases.")
          use.cases <- "pairwise.complete.obs"
        } 
        if (!is.null(relyFit[["complete"]])) {
          model[["footnote"]] <- paste0(model[["footnote"]], ". Using ", relyFit[["complete"]], 
                                        " listwise complete cases.")
          use.cases <- "complete.obs"
        } 
      } else {
        use.cases <- "everything"
        }
      
      # first the scale statistics
      cordat <- cor(dataset, use = use.cases)
      relyFit$freq$est$avg_cor <- mean(cordat[lower.tri(cordat)])
      relyFit$freq$est$mean <- mean(rowMeans(dataset, na.rm = T))
      relyFit$freq$est$sd <- sd(colMeans(dataset, na.rm = T))

      
      corsamp <- apply(relyFit$freq$covsamp, c(1), cov2cor)
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
      
      
      # Consider stripping some of the contents of relyFit to reduce memory load
      if (inherits(relyFit, "try-error")) {
        
        model[["error"]] <- paste("The analysis crashed with the following error message:\n", relyFit)
        
      } else {
        
        model[["dataset"]] <- dataset
        
        model[["relyFit"]] <- relyFit
        
        stateObj <- createJaspState(model)
        stateObj$dependOn(options = c("variables", "reverseScaledItems", "noSamplesf", "missingValues"))
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
      
      # omega doesnt work with bootstrapping:
      tmp <- Bayesrel::strel(dataset, estimates = c("alpha", "omega"), Bayes = F, item.dropped = T, 
                             alpha.int.analytic = TRUE, omega.freq.method = "cfa", 
                             interval = options[["confidenceIntervalValue"]]) 
      # this will work with the github package, so far the interval for omega doesnt change.

      alphaCfi <- as.vector(unlist(tmp$freq$conf))[c(1, 3)]
      names(alphaCfi) <- c("lower", "upper")
      scaleCfi$alpha <- alphaCfi
      omegaCfi <- as.vector(unlist(tmp$freq$conf))[c(2, 4)]
      names(omegaCfi) <- c("lower", "upper")
      scaleCfi$omega <- omegaCfi
      scaleCfi <- scaleCfi[c(7, 1, 2, 3, 8, 4, 5, 6)] # check this when more estimators come in

      cfiState <- list(scaleCfi = scaleCfi)
      jaspCfiState <- createJaspState(cfiState)
      jaspCfiState$dependOn(options = "confidenceIntervalValue", optionsFromObject = jaspResults[["modelObj"]])
      jaspResults[["cfiObj"]] <- jaspCfiState
      
    }
    model[["cfi"]] <- cfiState
  }
  
  model[["derivedOptions"]] <- .frequentistReliabilityDerivedOptions(options)
  model[["itemsDropped"]] <- .unv(colnames(dataset))
  
  return(model)
}

.frequentistReliabilityCalcCfi <- function(boot, cfiValue) {
  
  cfi <- vector("list", length(boot))
  names(cfi) <- names(boot)
  
  for (nm in names(boot)) {
    if (any(is.na(boot[[nm]]))) {
      cfi[[nm]] <- c(NA_real_, NA_real_)
    } else {
      cfi[[nm]] <- quantile(boot[[nm]], prob = c(0+(1-cfiValue)/2, 1-(1-cfiValue)/2))
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
                                   "averageInterItemCor", "meanScale", "sdScale"))
  
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
  

  if (!is.null(relyFit)) {
    allData <- cbind.data.frame(
      statistic = opts,
      pointEst = unlist(relyFit$freq$est, use.names = FALSE),
      do.call(rbind, model[["cfi"]][["scaleCfi"]])
    )[order, ][selected[order], ] # TODO: <- simplify this
    
    scaleTableF$setData(allData)
    
    if (!is.null(model[["footnotes"]]))
      scaleTableF$addFootnote(model[["footnotes"]])
    
  } else if (sum(selected) > 0L) {
    
    scaleTableF[["statistic"]] <- opts[order][selected[order]]
    
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
  order <- derivedOptions[["order_item"]]
  estimators <- derivedOptions[["namesEstimators"]][["tables_item"]][order]
  overTitle <- "If item dropped"
  
  itemTableF <- createJaspTable("Frequentist Individual Item Reliability Statistics")
  itemTableF$dependOn(options = c("variables",
                                  "mcDonaldScalef", "alphaScalef", "guttman2Scalef", "guttman6ScaleF", "glbScalef", 
                                  "averageInterItemCor", "meanScale", "sdScale",
                                  "mcDonaldItemf",  "alphaItemf",  "guttman2Itemf", "guttman6ItemF", "glbItemf",
                                  "reverseScaledItems", "meanItem", "sdItem", "itemRestCor"))
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
      idx <- order[i]
      newtb <- cbind(pointEst = relyFit$freq$ifitem[[idx]])
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
  fitTable$dependOn(options = c("variables", "mcDonaldScalef", "reverseScaledItems", "fitMeasures"))
  fitTable$addColumnInfo(name = "measure", title = "Fit Measure",   type = "string")
  fitTable$addColumnInfo(name = "value",     title = "Value", type = "number")

  relyFit <- model[["relyFit"]]
  derivedOptions <- model[["derivedOptions"]]
  opts     <- names(relyFit$freq$omega_fit)

  if (!is.null(relyFit)) {
    allData <- data.frame(
      measure = opts,
      value = as.vector(unlist(relyFit$freq$omega_fit, use.names = FALSE))
    )

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
