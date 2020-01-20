freqSTReliability <- function(jaspResults, dataset, options) {
  
  dataset <- .frequentistReliabilityReadData(dataset, options)
  
  .frequentistReliabilityCheckErrors(dataset, options)
  
  model <- .frequentistReliabilityMainResults(jaspResults, dataset, options)
  
  .frequentistReliabilityScaleTable(         jaspResults, model, options)
  .frequentistReliabilityItemTable(          jaspResults, model, options)
  # .freqentistReliabilitySingleFactorFitTable(jaspResults, model, options)
  return()
  
}

# read data, check errors----
.frequentistReliabilityDerivedOptions <- function(options) {
  
  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimatorsF  = unlist(options[c("alphaScalef", "guttman2Scalef", "glbScalef", "mcDonaldScalef", 
                                            "averageInterItemCor", "meanScale", "sdScale")]),
    itemDroppedSelectedF = unlist(options[c("mcDonaldItemf", "alphaItemf", "guttman2Itemf", "glbItemf",
                                            "itemRestCor", "meanItem", "sdItem")]),
    namesEstimators     = list(
      tables = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Greatest Lower Bound", 
                 "McDonald's \u03C9", "Average interitem correlation", "mean", "sd"),
      tables_item = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Greatest Lower Bound", 
                 "McDonald's \u03C9", "Item-rest correlation", "mean", "sd"),
      coefficients = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Greatest Lower Bound", 
                       "McDonald's \u03C9"),
      plots = list(expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]), 
                  "Greatest Lower Bound", expression("McDonald's"~omega))
    )
  )
  # order to show in JASP
  derivedOptions[["order"]] <- match(c("McDonald's \u03C9", "Cronbach's \u03B1", 
                                       "Guttman's \u03BB2", "Greatest Lower Bound",
                                       "Average interitem correlation", "mean", "sd"), 
                                     derivedOptions[["namesEstimators"]][["tables"]])
  derivedOptions[["order_item"]] <- match(c("McDonald's \u03C9", "Cronbach's \u03B1", 
                                       "Guttman's \u03BB2", "Greatest Lower Bound",
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
      
      model[["footnote"]] <- .frequentistReliabilityCheckLoadings(dataset, variables)
      relyFit <- try(Bayesrel::strel(x = dataset, estimates=c("alpha", "lambda2", "glb", "omega"), 
                                     n.iter = 250, boot.n = options[["noSamplesf"]],
                                     item.dropped = TRUE, omega.freq.method = "pa", ))
      
      # add the old funtionality here, beware of missings!!!
      # first the scale statistics
      cordat <- cor(dataset)
      relyFit$freq$est$avg_cor <- mean(cordat[lower.tri(cordat)])
      relyFit$freq$est$mean <- mean(dataset)
      relyFit$freq$est$sd <- sd(apply(dataset, 2, mean))
      
      corsamp <- apply(relyFit$freq$covsamp, c(1), cov2cor)
      relyFit$freq$boot$avg_cor <- apply(corsamp, 2, function(x) mean(x[x!=1]))
      relyFit$freq$boot$mean <- 0
      relyFit$freq$boot$sd <- 0
      
      
      # now the item statistics
      relyFit$freq$ifitem$ircor <- NULL
      
      for (i in 1:ncol(dataset)) {
        idx <- seq(1, ncol(dataset))
        idx <- idx[idx!=i]
        relyFit$freq$ifitem$ircor[i] <- cor(dataset[, i], apply(dataset[, idx], 1, mean))
      }
      relyFit$freq$ifitem$mean <- apply(dataset, 2, mean)
      relyFit$freq$ifitem$sd <- apply(dataset, 2, sd)  
      
      
      # Consider stripping some of the contents of relyFit to reduce memory load
      if (inherits(relyFit, "try-error")) {
        
        model[["error"]] <- paste("The analysis crashed with the following error message:\n", relyFit)
        
      } else {
        
        model[["dataset"]] <- dataset
        
        model[["relyFit"]] <- relyFit
        
        stateObj <- createJaspState(model)
        stateObj$dependOn(options = c("variables", "reverseScaledItems", "noSamplesf"))
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
                                              options[["ConfidenceIntervalValue"]])
      
      cfiState <- list(scaleCfi = scaleCfi)
      jaspCfiState <- createJaspState(cfiState)
      jaspCfiState$dependOn(options = "ConfidenceIntervalValue", optionsFromObject = jaspResults[["modelObj"]])
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

    cfi[[nm]] <- quantile(boot[[nm]], prob = c(0+(1-cfiValue)/2, 1-(1-cfiValue)/2))
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
  scaleTableF$dependOn(options = c("variables", "mcDonaldScalef", "alphaScalef", "guttman2Scalef",
                                   "glbScalef", "reverseScaledItems", "ConfidenceIntervalValue", "noSamplesf", 
                                   "averageInterItemCor", "meanScale", "sdScale"))
  
  overTitle <- sprintf("%s%% Confidence interval",
                       format(100*options[["ConfidenceIntervalValue"]], digits = 3, drop0trailing = TRUE))
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
                                  "mcDonaldScalef", "alphaScalef", "guttman2Scalef", "glbScalef", 
                                  "averageInterItemCor", "meanScale", "sdScale",
                                  "mcDonaldItemf",  "alphaItemf",  "guttman2Itemf", "glbItemf",
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

# # once the package is updated check this again and apply: 
# .freqentistReliabilitySingleFactorFitTable <- function(jaspResults, model, options) {
#   
#   if (!options[["fitMeasures"]] || !options[["mcDonaldScalef"]])
#     return()
#   if (!is.null(jaspResults[["fitTable"]]) || !options[["fitMeasures"]]) {
#     print("jaspResults[['fitTable']] from state or not wanted")
#     return()
#   }
#   
#   fitTable <- createJaspTable(sprintf("Fit Measures of Single Factor Model Fit"))
#   fitTable$dependOn(options = c("variables", "mcDonaldScalef", "reverseScaledItems", "fitMeasures"))
#   fitTable$addColumnInfo(name = "measure", title = "Fit Measure",   type = "string")
#   fitTable$addColumnInfo(name = "value",     title = "Value", type = "number")
#   
#   relyFit <- model[["relyFit"]]
#   derivedOptions <- model[["derivedOptions"]]
#   opts     <- names(relyFit$freq$fit$omega)
# 
#   if (!is.null(relyFit)) {
#     allData <- data.frame(
#       measure = opts,
#       value = unlist(relyFit$freq$fit$omega, use.names = FALSE)
#     )
#     
#     fitTable$setData(allData)
#   }
#   if (!is.null(model[["error"]]))
#     fitTable$setError(model[["error"]])
#   
#   jaspResults[["fitTable"]] <- fitTable
# }
# 
