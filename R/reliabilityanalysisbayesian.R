reliabilityBayesian <- function(jaspResults, dataset, options) {
# 
#   sink("~/Downloads/log_Bay.txt")
#   on.exit(sink(NULL))
  

	dataset <- .BayesianReliabilityReadData(dataset, options)

	.BayesianReliabilityCheckErrors(dataset, options)

	model <- .BayesianReliabilityMainResults(jaspResults, dataset, options)

	.BayesianReliabilityScaleTable(         jaspResults, model, options)
	.BayesianReliabilityItemTable(          jaspResults, model, options)
	.BayesianReliabilityProbTable(          jaspResults, model, options)
	.BayesianReliabilityPosteriorPlot(      jaspResults, model, options)
	.BayesianReliabilityPosteriorPredictive(jaspResults, model, options)
	.BayesianReliabilityIfItemPlot(         jaspResults, model, options)
	.BayesianReliabilityConvergenceTable(   jaspResults, model, options)
	.BayesianReliabilityTracePlot(          jaspResults, model, options)
	  

	return()

}

# read data, check errors----
.BayesianReliabilityDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("alphaScale", "guttman2Scale", "guttman6Scale", "glbScale", "mcDonaldScale",
                                           "averageInterItemCor", "meanScale", "sdScale")]),
    selectedEstimatorsPlots  = unlist(options[c("alphaScale", "guttman2Scale", "guttman6Scale", 
                                                "glbScale", "mcDonaldScale")]),
    itemDroppedSelected = unlist(options[c("mcDonaldItem", "alphaItem", "guttman2Item", "guttman6Item","glbItem",
                                           "itemRestCor", "meanItem", "sdItem")]),
    itemDroppedSelectedItem = unlist(options[c("mcDonaldItem", "alphaItem", "guttman2Item", "guttman6Item", 
                                               "glbItem")]),

    namesEstimators     = list(
      tables = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6", "Greatest Lower Bound", 
                 "McDonald's \u03C9", "Average interitem correlation", "mean", "sd"),
      tables_item = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6", "Greatest Lower Bound", 
                 "McDonald's \u03C9", "Item-rest correlation", "mean", "sd"),
      coefficients = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6", "Greatest Lower Bound", 
                       "McDonald's \u03C9", "Item-rest correlation"),
      plots = list(expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]), 
                   expression("Guttman's"~lambda[6]), "Greatest Lower Bound", expression("McDonald's"~omega))
    )

  )
  # order to show in JASP
  derivedOptions[["order"]] <- match(c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", 
                                       "Guttman's \u03BB6",  "Greatest Lower Bound", "Average interitem correlation", 
                                       "mean", "sd"), 
                                     derivedOptions[["namesEstimators"]][["tables"]])
  derivedOptions[["order_item"]] <- match(c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", 
                                            "Guttman's \u03BB6", "Greatest Lower Bound", "Item-rest correlation", 
                                            "mean", "sd"), 
                                     derivedOptions[["namesEstimators"]][["tables_item"]])


  return(derivedOptions)
}

.BayesianReliabilityReadData <- function(dataset, options) {

  variables <- unlist(options[["variables"]])
	if (is.null(dataset)) {
		dataset <- .readDataSetToEnd(columns.as.numeric = variables, columns.as.factor = NULL, exclude.na.listwise = NULL)
	}
  return(dataset)
}

.BayesianReliabilityCheckErrors <- function(dataset, options) {

  .hasErrors(dataset = dataset, perform = "run",
             type = c("infinity", "variance", "observations"),
             observations.amount = " < 3",
             exitAnalysisIfErrors = TRUE)
}

# estimate reliability ----
.BayesianReliabilityMainResults <- function(jaspResults, dataset, options) {
  
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
        # dataset <- dataset %*% diag(key, nvar, nvar) # seems to be not working
        cols <- match(unlist(options[["reverseScaledItems"]]), .unv(colnames(dataset)))
        total <- min(dataset, na.rm = T) + max(dataset, na.rm = T)
        dataset[ ,cols] = total - dataset[ ,cols]
      }
      
      missing <- "none" 
      options[["missings"]] <- "everything"
      if (any(is.na(dataset))) {
        if (options[["missingValues"]] == "excludeCasesPairwise") {
          missing <- "pairwise"
          options[["missings"]] <- "pairwise.complete.obs"
        } else if (options[["missingValues"]] == "excludeCasesListwise") {
          missing <- "listwise"
          options[["missings"]] <- "complete.obs"
          }
      }
      model[["footnote"]] <- .BayesianReliabilityCheckLoadings(dataset, variables)
      relyFit <- try(Bayesrel::strel(data = dataset, estimates=c("alpha", "lambda2", "lambda6", "glb", "omega"), 
                                     n.iter = options[["noSamples"]], n.burnin = options[["noBurnin"]], 
                                     # n.chains = options[["noChains"]], thin = options[["noThin"]],
                                     freq = F, item.dropped = TRUE, missing = missing))
      
      if (any(is.na(dataset))) {
        if (!is.null(relyFit[["miss_pairwise"]])) {
          model[["footnote"]] <- paste0(model[["footnote"]], ". Using pairwise complete cases.")
        }
        if (!is.null(relyFit[["complete"]])) {
          model[["footnote"]] <- paste0(model[["footnote"]], ". Using ", relyFit[["complete"]], 
                                        " listwise complete cases.")
        }
      } 
      
      # add the scale info
      corsamp <- apply(relyFit$Bayes$covsamp, 1, cov2cor)
      relyFit$Bayes$samp$avg_cor <- coda::mcmc(apply(corsamp, 2, function(x) mean(x[x!=1])))
      relyFit$Bayes$est$avg_cor <- median(relyFit$Bayes$samp$avg_cor)

      relyFit$Bayes$samp$mean <- c(NA_real_, NA_real_)
      relyFit$Bayes$est$mean <- mean(rowMeans(dataset, na.rm = T))
      relyFit$Bayes$samp$sd <- c(NA_real_, NA_real_)
      relyFit$Bayes$est$sd <- sd(colMeans(dataset, na.rm = T))


      # now the item statistics
      relyFit$Bayes$ifitem$samp$ircor <- .reliabilityItemRestCor(dataset, options[["noSamples"]], options[["noBurnin"]], 
                                                                 missing)
      relyFit$Bayes$ifitem$est$ircor <- apply(relyFit$Bayes$ifitem$samp$ircor, 2, median)

      relyFit$Bayes$ifitem$est$mean <- colMeans(dataset, na.rm = T)
      relyFit$Bayes$ifitem$est$sd <- apply(dataset, 2, sd, na.rm = T)
      relyFit$Bayes$ifitem$samp$mean <- (matrix(NA_real_, ncol(dataset), 2))
      relyFit$Bayes$ifitem$samp$sd <- (matrix(NA_real_, ncol(dataset), 2))

      
      # Consider stripping some of the contents of relyFit to reduce memory load
      if (inherits(relyFit, "try-error")) {

        model[["error"]] <- paste("The analysis crashed with the following error message:\n", relyFit)

      } else {
        
        model[["dataset"]] <- dataset

        model[["relyFit"]] <- relyFit
        
        model[["options"]] <- options

        stateObj <- createJaspState(model)
        stateObj$dependOn(options = c("variables", "reverseScaledItems", "noSamples", "noBurnin", "noChains", "noThin",
                                      "missingValues"))
        jaspResults[["modelObj"]] <- stateObj
      }
    }
  } else {
    print("model from state")
  }

	if (is.null(model[["error"]])) {
	  criState <- jaspResults[["criObj"]]$object
	  if (is.null(criState) && !is.null(relyFit)) {

	    scaleCri <- .BayesianReliabilityCalcCri(relyFit[["Bayes"]][["samp"]],             
	                                            options[["credibleIntervalValue"]])
	    itemCri  <- .BayesianReliabilityCalcCri(relyFit[["Bayes"]][["ifitem"]][["samp"]], 
	                                            options[["credibleIntervalValue"]])
      
	    
	    criState <- list(scaleCri = scaleCri, itemCri  = itemCri)
	    jaspCriState <- createJaspState(criState)
	    jaspCriState$dependOn(options = "credibleIntervalValue", optionsFromObject = jaspResults[["modelObj"]])
	    jaspResults[["criObj"]] <- jaspCriState

	  }
	  model[["cri"]] <- criState
	}

  model[["derivedOptions"]] <- .BayesianReliabilityDerivedOptions(options)
  model[["itemsDropped"]] <- .unv(colnames(dataset))

	return(model)
}

.BayesianReliabilityCalcCri <- function(samps, criValue) {

  cri <- vector("list", length(samps))
  names(cri) <- names(samps)
  
  for (nm in names(samps)) {
    if (any(is.na(samps[[nm]]))) {
      cri[[nm]] <- NA_real_
    } else {
      cri[[nm]] <- coda::HPDinterval(coda::mcmc(samps[[nm]]), prob = criValue)
    }
  }
  return(cri)
}


.BayesianReliabilityCheckLoadings <- function(dataset, variables) {

	prin <- psych::principal(dataset)
	idx <- prin[["loadings"]] < 0
	sidx <- sum(idx)
	hasSchar <- if (sidx == 1L) "" else "s"
	footnote <- sprintf("The following item%s correlated negatively with the scale: %s",
	                    hasSchar, paste0(variables[idx], collapse = ", "))
	return(footnote)
}


# ----------------------------- tables ------------------------------------
.BayesianReliabilityScaleTable <- function(jaspResults, model, options) {

  if (!is.null(jaspResults[["scaleTable"]])) {
    print("jaspResults[['scaleTable']] from state")
    return()
  }

  scaleTable <- createJaspTable("Bayesian Scale Reliability Statistics")
  scaleTable$dependOn(options = c("variables", "mcDonaldScale", "alphaScale", "guttman2Scale",  "guttman6Scale",
                                  "glbScale", "reverseScaledItems", "credibleIntervalValue", "noSamples", "noBurnin",
                                  "noChains", "noThin",
                                  "averageInterItemCor", "meanScale", "sdScale", "missingValues"))

  overTitle <- sprintf("%s%% Credible interval",
                       format(100*options[["credibleIntervalValue"]], digits = 3, drop0trailing = TRUE))
	scaleTable$addColumnInfo(name = "statistic", title = "Statistic",        type = "string")
	scaleTable$addColumnInfo(name = "postMean",  title = "Posterior Median", type = "number")
	scaleTable$addColumnInfo(name = "lower",     title = "Lower",            type = "number", overtitle = overTitle)
	scaleTable$addColumnInfo(name = "upper",     title = "Upper",            type = "number", overtitle = overTitle)

	relyFit <- model[["relyFit"]]
	derivedOptions <- model[["derivedOptions"]]
	opts     <- derivedOptions[["namesEstimators"]][["tables"]]
	order    <- derivedOptions[["order"]]
	selected <- derivedOptions[["selectedEstimators"]]

	if (!is.null(relyFit)) {
	  allData <- cbind.data.frame(
	    statistic = opts,
	    postMean = unlist(relyFit$Bayes$est, use.names = FALSE),
	    do.call(rbind, model[["cri"]][["scaleCri"]])
	  )[order, ][selected[order], ] # TODO: <- simplify this

		scaleTable$setData(allData)

		if (!is.null(model[["footnotes"]]))
		  scaleTable$addFootnote(model[["footnotes"]])

	} else if (sum(selected) > 0L) {

    scaleTable[["statistic"]] <- opts[order][selected[order]]

    nvar <- length(options[["variables"]])
    if (nvar > 0L && nvar < 3L)
      scaleTable$addFootnote("Please enter at least 3 variables to do an analysis.")
	}
	if (!is.null(model[["error"]]))
	  scaleTable$setError(model[["error"]])

	if (!is.null(model[["footnote"]]))
	  scaleTable$addFootnote(model[["footnote"]])

	jaspResults[["scaleTable"]] <- scaleTable

	return()
}


.BayesianReliabilityItemTable <- function(jaspResults, model, options) {

	if (!is.null(jaspResults[["itemTable"]]) || !any(model[["derivedOptions"]][["itemDroppedSelected"]])) {
	  print("jaspResults[['hasItemTable']] from state or not wanted")
		return()
	}

  derivedOptions <- model[["derivedOptions"]]
  itemDroppedSelected <- derivedOptions[["itemDroppedSelected"]]
  order <- derivedOptions[["order_item"]]
  overTitles <- format(derivedOptions[["namesEstimators"]][["tables_item"]][order], digits = 3, drop0trailing = T)
  overTitles <- paste0(overTitles, " (if item dropped)")
  
  cred <- format(100*options[["credibleIntervalValue"]], digits = 3, drop0trailing = TRUE)
  itemTable <- createJaspTable("Bayesian Individual Item Reliability Statistics")
  itemTable$dependOn(options = c("variables",
                                 "mcDonaldScale", "alphaScale", "guttman2Scale",  "guttman6Scale", "glbScale", 
                                 "averageInterItemCor", "meanScale", "sdScale",
                                 "mcDonaldItem",  "alphaItem",  "guttman2Item",  "guttman6Item", "glbItem",
                                 "reverseScaledItems", "credibleIntervalValue", 
                                 "itemRestCor", "meanItem", "sdItem", "missingValues"))
  itemTable$addColumnInfo(name = "variable", title = "Item", type = "string")

  idxSelected <- which(itemDroppedSelected)
  estimators <- derivedOptions[["namesEstimators"]][["tables_item"]][order]
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]
  
  for (i in idxSelected) {
    if (estimators[i] %in% coefficients) {
      if (estimators[i] == "Item-rest correlation") {
        itemTable$addColumnInfo(name = paste0("postMean", i), title = "Posterior Median", type = "number", 
                                overtitle = "Item-rest correlation")
        itemTable$addColumnInfo(name = paste0("lower", i), title = paste0("Lower ", cred, "%"), type = "number", 
                                overtitle = "Item-rest correlation")
        itemTable$addColumnInfo(name = paste0("upper", i), title = paste0("Upper ", cred, "%"), type = "number", 
                                overtitle = "Item-rest correlation")
      } else {
        itemTable$addColumnInfo(name = paste0("postMean", i), title = "Posterior Median", type = "number", 
                                overtitle = overTitles[i])
        itemTable$addColumnInfo(name = paste0("lower", i), title = paste0("Lower ", cred, "%"), type = "number", 
                                overtitle = overTitles[i])
        itemTable$addColumnInfo(name = paste0("upper", i), title = paste0("Upper ", cred, "%"), type = "number", 
                                overtitle = overTitles[i])
      }
    } else {
      itemTable$addColumnInfo(name = paste0("postMean", i), title = estimators[i], type = "number")
    }

  }

  relyFit <- model[["relyFit"]]
  if (!is.null(relyFit)) {
    cris <- model[["cri"]][["itemCri"]]
    tb <- data.frame(variable = model[["itemsDropped"]])
    for (i in idxSelected) {
      idx <- order[i]
      if (idx %in% c(1:6)) { # check this when more estimators are included !!!!!!!!!!!!!!!!!!!!!
        newtb <- cbind(postMean = relyFit$Bayes$ifitem$est[[idx]], cris[[idx]])
      } else {
        newtb <- cbind(postMean = relyFit$Bayes$ifitem$est[[idx]])
      }
      colnames(newtb) <- paste0(colnames(newtb), i)
      tb <- cbind(tb, newtb)
    }
		itemTable$setData(tb)

  } else if (length(model[["itemsDropped"]]) > 0) {
    itemTable[["variables"]] <- model[["itemsDropped"]]
  }

  jaspResults[["itemTable"]] <- itemTable
	return()
}


.BayesianReliabilityProbTable <- function(jaspResults, model, options) {

  if (!is.null(jaspResults[["probTable"]]) || !options[["probTable"]]) {
    print("jaspResults[['probTable']] from state or not wanted")
		return()
  }

  probTable <- createJaspTable(
    sprintf("Probability that Reliability Statistic is Larger than %.2f and Smaller than %.2f", 
            options[["probTableValueLow"]], options[["probTableValueHigh"]]))
  probTable$dependOn(options = c("variables", "mcDonaldScale", "alphaScale", "guttman2Scale", "guttman6Scale",
                                 "glbScale", "reverseScaledItems", "probTableValueLow", "probTable",
                                 "probTableValueHigh", "missingValues"))
  overTitle <- format("Probability",
                      digits = 3, drop0trailing = T)
  probTable$addColumnInfo(name = "statistic", title = "Statistic",   type = "string")
	probTable$addColumnInfo(name = "prior",     title = "Prior", type = "number", overtitle = overTitle )
	probTable$addColumnInfo(name = "posterior", title = "Posterior", type = "number", overtitle = overTitle )
	

  relyFit <- model[["relyFit"]]
	derivedOptions <- model[["derivedOptions"]]
	opts     <- derivedOptions[["namesEstimators"]][["tables"]]
	order    <- derivedOptions[["order"]]
	selected <- derivedOptions[["selectedEstimators"]]
	idx      <- which(selected)
  
	n.item <- dim(relyFit$Bayes$covsamp)[2]
	prior <- priors[[as.character(n.item)]] 
  end <- length(prior[[1]][["x"]])
  poslow <- end - sum(prior[[1]][["x"]] > options[["probTableValueLow"]]) 
  poshigh <- end - sum(prior[[1]][["x"]] > options[["probTableValueHigh"]]) 
  # since the priors are only available in density form, the prior probability for the estimator being larger than
  # a cutoff is given by caculating the relative probability of the density from the cutoff to 1.
  # maybe check this with Don though
  
	if (!is.null(relyFit)) {
    probsPost <- numeric(sum(selected))
    probsPrior <- numeric(sum(selected))
    for (i in seq_along(idx)) {
      probsPost[i] <- mean(relyFit[["Bayes"]][["samp"]][[idx[i]]] > options[["probTableValueLow"]]) -
                      mean(relyFit[["Bayes"]][["samp"]][[idx[i]]] > options[["probTableValueHigh"]])
      probsPrior[i] <- sum(prior[[idx[i]]][["y"]][poslow:end]) / sum(prior[[idx[i]]][["y"]]) - 
                       sum(prior[[idx[i]]][["y"]][poshigh:end]) / sum(prior[[idx[i]]][["y"]])
    }
    df <- data.frame(statistic = opts[idx], prior = probsPrior, posterior = probsPost)
    probTable$setData(df)
  } else if (sum(selected) > 0) {
    probTable[["statistic"]] <- opts[idx]
  }

  jaspResults[["probTable"]] <- probTable
  return()
}


.BayesianReliabilityConvergenceTable <- function(jaspResults, model, options) {
  
  if (!is.null(jaspResults[["convTable"]]) || !options[["convTable"]]) {
    print("jaspResults[['convTable']] from state or not wanted")
    return()
  }
  convTable <- createJaspTable(
    sprintf("Convergence Diagnostics"))
  convTable$dependOn(options = c("variables", "mcDonaldScale", "alphaScale", "guttman2Scale", "guttman6Scale",
                                 "glbScale", "reverseScaledItems", "rHat", "missingValues"))
  convTable$addColumnInfo(name = "statistic", title = "Statistic",   type = "string")
  convTable$addColumnInfo(name = "rhat",     title = "R-hat", type = "number")
  
  relyFit <- model[["relyFit"]]
  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  order    <- derivedOptions[["order"]]
  selected <- derivedOptions[["selectedEstimators"]]
  idx      <- which(selected)
  
  if (!is.null(relyFit)) {
    rhat <- numeric(sum(selected))
    for (i in seq_along(idx)) {
      tmp <- lapply(as.data.frame(t(relyFit[["Bayes"]][["samp"]][[idx[i]]])), coda::mcmc)
      rhat[i] <- coda::gelman.diag(coda::as.mcmc.list(tmp))[["psrf"]][, 1]
    }
    df <- data.frame(statistic = opts[idx], rhat = rhat)
    convTable$setData(df)
  } else if (sum(selected) > 0) {
    convTable[["statistic"]] <- opts[idx]
  }
  jaspResults[["convTable"]] <- convTable
  return()
}





# -------------------------------------------- plots ---------------------------------
.BayesianReliabilityPosteriorPlot <- function(jaspResults, model, options) {

	if (!options[["plotPosterior"]])
		return()

  plotContainer <- jaspResults[["plotContainer"]]
  if (is.null(plotContainer)) {
    print("Plotcontainer remade")
    plotContainer <- createJaspContainer("Posterior Plots")
    plotContainer$dependOn(options = c("variables", "reverseScaledItems", "plotPosterior", "shadePlots",
                                       "probTable", "probTableValueLow", "probTableValueHigh", "fixXRange", 
                                       "dispPrior", "noSamples", "noBurnin", "noChains", "noThin",
                                       "credibleIntervalValue","alphaScale", "guttman2Scale", "guttman6Scale", 
                                       "glbScale", "mcDonaldScale", "missingValues"))
    jaspResults[["plotContainer"]] <- plotContainer
  } else {
    print("Plotcontainer from state")
  }

	derivedOptions <- model[["derivedOptions"]]
	order     <- derivedOptions[["order"]]
	indices   <- which(derivedOptions[["selectedEstimatorsPlots"]][order])
	nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]][order]
	nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]][order]

	relyFit  <- model[["relyFit"]]
	# probably better to do this once directly after computation!
	scaleCri <- model[["cri"]][["scaleCri"]][order]
	relyFit[["Bayes"]][["samp"]] <- relyFit[["Bayes"]][["samp"]][order]
	n.item <- dim(relyFit$Bayes$covsamp)[2]
	prior <- priors[[as.character(n.item)]][c(5, 1, 2, 3, 4)] ##### change this when more estimators are included!!!


	if (options[["shadePlots"]] && options[["probTable"]]) {
	  shadePlots <- c(options[["probTableValueLow"]], options[["probTableValueHigh"]])
	} else {
	  shadePlots <- NULL
	}


	if (!is.null(relyFit)) {
	  for (i in indices) {
	    if (is.null(plotContainer[[nmsObjs[i]]])) {

	      p <- .BayesianReliabilityMakeSinglePosteriorPlot(relyFit, scaleCri, i, nmsLabs[[i]], options[["fixXRange"]],
	                                                       shadePlots, options[["dispPrior"]], prior)
	      plotObj <- createJaspPlot(plot = p, title = nmsObjs[i])
	      plotObj$dependOn(options = names(indices[i]))
	      plotObj$position <- i
	      plotContainer[[nmsObjs[i]]] <- plotObj

	    } else {
	      print(sprintf("plotContainer[[%s]] from state", nmsObjs[i]))
	    }
	  }
	} else if (length(indices) > 0) {
	  for (i in indices) {
	    plotObj <- createJaspPlot(title = nmsObjs[i])
	    plotObj$dependOn(options = names(indices[i]))
	    plotObj$position <- i
	    plotContainer[[nmsObjs[i]]] <- plotObj
	  }
	} else {
	  plotContainer[["Posterior Plots"]] <- createJaspPlot()
	}

	return()
}

.BayesianReliabilityMakeSinglePosteriorPlot <- function(relyFit, scaleCri, i, nms, fixXRange, 
                                                        shade = NULL, priorTrue, priorSample) {

  # TODO: consider precomputing all densities (maybe with kernsmooth?) and reducing memory that way
  pr <- priorSample[[i]]
  if (fixXRange) {
    d <- stats::density(relyFit$Bayes$samp[[i]], from = 0, to = 1, n = 2^10)
  } else {
    d <- stats::density(relyFit$Bayes$samp[[i]], n = 2^10)
  }
	datDens <- data.frame(x = d$x, y = d$y)
	datPrior <- data.frame(x = pr$x, y = pr$y)
	


	xBreaks <- JASPgraphs::getPrettyAxisBreaks(datDens$x)
	# max height posterior is at 90% of plot area; remainder is for credible interval
	ymax <- max(d$y) / .9
	yBreaks <- JASPgraphs::getPrettyAxisBreaks(c(0, ymax))
	ymax <- max(yBreaks)
	datCri <- data.frame(xmin = scaleCri[[i]][1L], xmax = scaleCri[[i]][2L], y = .925 * ymax)
	height <- (ymax - .925 * ymax) / 2
	if (fixXRange) {
	  datTxt <- data.frame(x = c(datCri$xmin*.9, datCri$xmax*1.08),
	                       y = 0.985 * ymax,
	                       label = sapply(c(datCri$xmin, datCri$xmax), format, digits = 3, scientific = -1),
	                       stringsAsFactors = FALSE)
	} else {
	  datTxt <- data.frame(x = c(datCri$xmin, datCri$xmax),
	                       y = 0.985 * ymax,
	                       label = sapply(c(datCri$xmin, datCri$xmax), format, digits = 3, scientific = -1),
	                       stringsAsFactors = FALSE)
	}


	if (datCri$xmin[1L] < 0) {
	  datCri$xmin[1L] <- 0
	  datTxt$x[1L] <- 0
	  datTxt$label[1L] <- "< 0"
	}

	# if the bounds are less than 0.05 away from 0 or 1, expand the axis by 0.1 so the credible interval text does not
	# get chopped off.
	xExpand <- .1 * ((c(0, 1) - datTxt$x) <= 0.05)

	g <- ggplot2::ggplot(data = datDens, mapping = ggplot2::aes(x = x, y = y)) +
		ggplot2::geom_line(size = .85) +
		ggplot2::geom_errorbarh(data = datCri, mapping = ggplot2::aes(xmin = xmin, xmax = xmax, y = y),
		                        height = height, inherit.aes = FALSE) +
		# ggrepel::geom_text_repel(data = datTxt, mapping = ggplot2::aes(x = x, y = y, label = label), 
	  #                          inherit.aes = FALSE, segment.alpha = 0) +
	ggplot2::geom_text(data = datTxt, mapping = ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE) +
	ggplot2::scale_y_continuous(name = "Density", breaks = yBreaks, limits = range(yBreaks)) +
	ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, expand = xExpand)
	
	if (!is.null(shade)) {
	  datFilter <- datDens[datDens[["x"]] >= shade[1] & datDens[["x"]] <= shade[2], ]
	  g <- g + ggplot2::geom_ribbon(data = datFilter, mapping = ggplot2::aes(ymin = 0, ymax = y), 
	                                fill = "grey", alpha = 0.95) +
	           ggplot2::geom_line(size = .85)
	}
	
	if (priorTrue) {
	  # g <- g + ggplot2::geom_hline(yintercept = 1)
	  g <- g + ggplot2::geom_line(data = datPrior, mapping = ggplot2::aes(x = x, y = y),
	                              linetype = "dashed", size = .85) +
	           ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, expand = xExpand,
	                                       limits= c(min(xBreaks), max(xBreaks)))
	           # ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, expand = xExpand, 
	           #                             limits= c(min(d$x), max(d$x)))
	  

	  
	}

	# if (!is.null(cutoffs)) {
	#   cut1_peak <- d$y[findInterval(cutoffs[1], d$x)]
	#   cut2_peak <- d$y[findInterval(cutoffs[2], d$x)]
	#   if (length(cut1_peak) != 0 & findInterval(cutoffs[1], d$x) != 2^10) {
	#     g <- g +
	#       ggplot2::geom_segment(ggplot2::aes(x = cutoffs[1], y = 0, xend = cutoffs[1], yend = cut1_peak), 
	#                             color = "grey60", linetype = 1, alpha = .5, size = .3) +
	#       ggplot2::geom_line(size = .85)
	#     
	#   } 
	#   if (length(cut2_peak) != 0 & findInterval(cutoffs[2], d$x) != 2^10) {
	#     g <- g +
	#       ggplot2::geom_segment(ggplot2::aes(x = cutoffs[2], y = 0, xend = cutoffs[2], yend = cut2_peak), 
	#                           color = "grey60", linetype = 1, alpha = .5, size = .3) +
	#       ggplot2::geom_line(size = .85)
	#   }
	# }
	


	return(JASPgraphs::themeJasp(g))

}


.BayesianReliabilityPosteriorPredictive <- function(jaspResults, model, options) {

  if (!options[["dispPPC"]] || !options[["mcDonaldScale"]])
    return()
  
  relyFit <- model[["relyFit"]]
  dataset <- model[["dataset"]]
  
  if (is.null(relyFit)) {
    g <- NULL
  } else {
    ll <- relyFit[["Bayes"]][["loadings"]]
    rr <- relyFit[["Bayes"]][["resid_var"]]
    print(options[["missings"]])
    cimpl <- ll %*% t(ll) + diag(rr)
    cobs <- cov(dataset, use = model[["options"]][["missings"]])
    k <- ncol(cobs)
    eframe <- data.frame(number = seq(1, k), eigen_value = eigen(cobs)$values)
    ee_impl <- matrix(0, 1e3, k)
    for (i in 1:1e3) {
      dtmp <- MASS::mvrnorm(nrow(dataset), rep(0, k), cimpl)
      ee_impl[i, ] <- eigen(cov(dtmp))$values
    }
    eframe$eigen_sim_low <- apply(ee_impl, 2, quantile, prob = .025)
    eframe$eigen_sim_up<- apply(ee_impl, 2, quantile, prob = .975)
    leg_pos <- (max(eframe$eigen_value) + min(eframe$eigen_value)) * .75
    yBreaks <- JASPgraphs::getPrettyAxisBreaks(c(0, max(eframe$eigen_sim_up)))
    
    
    g <- ggplot2::ggplot(eframe, mapping = ggplot2::aes(x = number, y = eigen_value)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = eigen_sim_low, ymax = eigen_sim_up), fill = "grey80") +
      ggplot2::geom_line(ggplot2::aes(x = number, y = eigen_sim_low), colour = "grey", linetype = 2) +
      ggplot2::geom_line(ggplot2::aes(x = number, y = eigen_sim_up), colour = "grey", linetype = 2) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      # ggplot2::geom_segment(ggplot2::aes(x = k/2.5, y = leg_pos, xend = k/2.5*1.1, yend = leg_pos), size = .4) +
      # ggplot2::geom_segment(ggplot2::aes(x = k/2.5, y = leg_pos*.9, xend = k/2.5*1.1, yend = leg_pos*.9), size = .4,
      #              colour = "gray", linetype = 2) +
      # ggplot2::annotate(geom = "text", x = k/2.5, y = leg_pos, 
      #                   label = "data set covariance matrix", hjust = -.2) +
      # ggplot2::annotate(geom = "text", x = k/2.5, y = leg_pos*.9, 
      #                   label = "model implied covariance matrix", hjust = -.164) +
      ggplot2::scale_y_continuous(name = "Eigenvalue", breaks = yBreaks, limits = range(yBreaks)) +
      ggplot2::xlab("Factors")
    
    g <- JASPgraphs::themeJasp(g)
  }
  plot <- createJaspPlot(plot = g, title = "Posterior Predictive Check Omega", width = 400)
  plot$dependOn(options = c("variables", "reverseScaledItems", "noSamples", "noBurnin", "noChains", "noThin",
                            "credibleIntervalValue", "dispPPC", "mcDonaldScale", "missingValues"))
  jaspResults[["OmegaPosteriorPredictive"]] <- plot
}


.BayesianReliabilityIfItemPlot <- function(jaspResults, model, options) {
  
  if (!options[["plotItem"]])
    return()
  
  plotContainerItem <- jaspResults[["plotContainerItem"]]
  if (is.null(plotContainerItem)) {
    print("Plotcontainer remade")
    plotContainerItem <- createJaspContainer("If Item Dropped Posterior Plots")
    plotContainerItem$dependOn(options = c("variables", "plotItem", "noSamples", "noBurnin", "noChains", "noThin",
                                           "credibleIntervalValue", 
                                          "orderItemKL", "orderItemKS", "reverseScaledItems", "missingValues"))
    jaspResults[["plotContainerItem"]] <- plotContainerItem
  } else {
    print("Plotcontainer from state")
  }
  
  derivedOptions <- model[["derivedOptions"]]
  order     <- derivedOptions[["order"]]
  order_item  <- derivedOptions[["order_item"]]
  indices   <- which(derivedOptions[["itemDroppedSelectedItem"]])
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]][order_item]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables_item"]][order_item]
  
  relyFit  <- model[["relyFit"]]
  # probably better to do this once directly after computation!
  relyFit[["Bayes"]][["samp"]] <- relyFit[["Bayes"]][["samp"]][order]
  relyFit[["Bayes"]][["ifitem"]][["samp"]] <- relyFit[["Bayes"]][["ifitem"]][["samp"]][order_item]
  
  if (options[["orderItemKL"]]) {ordering <- "kl"}
  else if (options[["orderItemKS"]]) {ordering <- "ks"}
  else {ordering <-  NULL}

  if (!is.null(relyFit)) {
    for (i in indices) {
      if (is.null(plotContainerItem[[nmsObjs[i]]])) {
        
        p <- .BayesianReliabilityMakeIfItemPlot(relyFit, i, nmsLabs[[i]], options[["credibleIntervalValue"]], 
                                                ordering = ordering)
        plotObjItem <- createJaspPlot(plot = p, title = nmsObjs[i])
        plotObjItem$dependOn(options = names(indices[i]))
        plotObjItem$position <- i
        plotContainerItem[[nmsObjs[i]]] <- plotObjItem
        
      } else {
        print(sprintf("plotContainerItem[[%s]] from state", nmsObjs[i]))
      }
    }
  } else if (length(indices) > 0) {
    for (i in indices) {
      plotObjItem <- createJaspPlot(title = nmsObjs[i])
      plotObjItem$dependOn(options = names(indices[i]))
      plotObjItem$position <- i
      plotContainerItem[[nmsObjs[i]]] <- plotObjItem
    }
  } else {
    plotContainerItem[["If Item Dropped Posterior Plots"]] <- createJaspPlot()
  }
  
  return()
}

.BayesianReliabilityMakeIfItemPlot <- function(relyFit, i, nms, int, ordering) {
  
  n_row <- length(unlist(relyFit$Bayes$ifitem$est[1]))
  lower <- (1-int)/2
  upper <- int + (1-int)/2

  dat <- data.frame(as.matrix(unlist(relyFit$Bayes$samp[[i]])), row.names =  NULL)
  names(dat) <- "value"
  dat$colos <- "1"
  dat$var <- "original"
  
  dat_del <- t(as.matrix(as.data.frame(relyFit$Bayes$ifitem$samp[[i]])))
  
  names <- NULL
  for(j in 1:(n_row)){
    names[j] <- paste0("x", j)
  }
  
  for (j in 1:n_row){
    tmp <- as.data.frame(dat_del[j, ])
    colnames(tmp) <- "value"
    tmp$var <- names[j]
    tmp$colos <- "2"
    dat <- rbind(dat, tmp)
  }
  
  dat$var <- factor(dat$var, levels = unique(dat$var))
  
  if (!is.null(ordering)) {
    if (ordering == "kl") {
      est <- as.data.frame(relyFit$Bayes$ifitem$est[[i]])
      est[n_row + 1, ] <- 1
      colnames(est) <- "value"
      est$name <- c(names, "original")
      samps <- relyFit$Bayes$ifitem$samp[[i]]
      og_samp <- relyFit$Bayes$samp[[i]]
      dists <- apply(samps, 2, .KLD.statistic, y = og_samp) # kl divergence
      dists[length(dists)+1] <- 0
      est <- est[order(dists), ]
      dat$var <- factor(dat$var, levels = c(est$name))
    
    } else if (ordering == "ks") {
      est <- as.data.frame(relyFit$Bayes$ifitem$est[[i]])
      est[n_row + 1, ] <- 1
      colnames(est) <- "value"
      est$name <- c(names, "original")
      samps <- relyFit$Bayes$ifitem$samp[[i]]
      og_samp <- relyFit$Bayes$samp[[i]]
      dists <- apply(samps, 2, .ks.test.statistic, y = og_samp) # ks distance
      dists[length(dists)+1] <- 0
      est <- est[order(dists), ]
      dat$var <- factor(dat$var, levels = c(est$name))
    }
  }
  
  
  g <- ggplot2::ggplot(dat, ggplot2::aes(x = value, y = var, fill = colos)) +
    ggridges::stat_density_ridges(quantile_lines = T, quantiles = c(lower, 0.5, upper),
                                  alpha = .85, show.legend = F, scale =1) +
    # ggplot2::theme_linedraw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                   strip.text = ggplot2::element_text(colour = "black")) +
    ggplot2::xlab(nms) +
    ggplot2::ylab("Item Dropped") +
    ggplot2::scale_fill_grey()
    # ggplot2::scale_y_discrete(expand = ggplot2::expand_scale(add = c(0.25, 1.5))) 
    # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust = 4, size = 20),
    #                axis.title = ggplot2::element_text(size = 16),
    #                axis.text = ggplot2::element_text(size = 12))
  
  
  return(JASPgraphs::themeJasp(g))
  
}

  

.BayesianReliabilityTracePlot <- function(jaspResults, model, options) {
  
  if (!options[["tracePlot"]])
    return()
  
  plotContainerTP <- jaspResults[["plotContainerTP"]]
  if (is.null(plotContainerTP)) {
    print("PlotcontainerTP remade")
    plotContainerTP <- createJaspContainer("Convergence Traceplot")
    plotContainerTP$dependOn(options = c("variables", "tracePlot", "noSamples", "noBurnin", "noChains", "noThin",
                                         "missingValues", "reverseScaledItems"))
    jaspResults[["plotContainerTP"]] <- plotContainerTP
    
  } else {
    print("PlotcontainerTP from state")
  }
  
  derivedOptions <- model[["derivedOptions"]]
  order     <- derivedOptions[["order"]]
  indices   <- which(derivedOptions[["selectedEstimatorsPlots"]][order])
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]][order]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]][order]
  
  relyFit  <- model[["relyFit"]]
  # probably better to do this once directly after computation!
  relyFit[["Bayes"]][["samp"]] <- relyFit[["Bayes"]][["samp"]][order]
  
  if (!is.null(relyFit)) {
    for (i in indices) {
      if (is.null(plotContainerTP[[nmsObjs[i]]])) {
        
        p <- .BayesianReliabilityMakeTracePlot(relyFit, i, nmsLabs[[i]])
        plotObjTP <- createJaspPlot(plot = p, title = nmsObjs[i], width = 400)
        plotObjTP$dependOn(options = names(indices[i]))
        plotObjTP$position <- i
        plotContainerTP[[nmsObjs[i]]] <- plotObjTP
        
      } else {
        print(sprintf("plotContainerTP[[%s]] from state", nmsObjs[i]))
      }
    }
  } else if (length(indices) > 0) {
    for (i in indices) {
      plotObjTP <- createJaspPlot(title = nmsObjs[i])
      plotObjTP$dependOn(options = names(indices[i]))
      plotObjTP$position <- i
      plotContainerTP[[nmsObjs[i]]] <- plotObjTP
    }
  } else {
    plotContainerItem[["If Item Dropped Posterior Plots"]] <- createJaspPlot()
  }
  
  return()
}


.BayesianReliabilityMakeTracePlot <- function(relyFit, i, nms) {
  
  dt <- data.frame(relyFit$Bayes$samp[[i]])
  names(dt)[1] <- "Value"
  dt$Iterations <- seq(1, nrow(dt))
  
  g <- ggplot2::ggplot(dt, ggplot2::aes(x = Iterations, y = Value)) +
    ggplot2::geom_line() +
    ggplot2::ylab(nms)
  
  return(JASPgraphs::themeJasp(g))
  
}




# ----- some other functions -----------------
.reliabilityItemRestCor <- function(dataset, n.iter, n.burnin, missing) {
  if (missing == "listwise") {
    pos <- which(is.na(dataset), arr.ind = T)[, 1]
    dataset <- dataset[-pos, ]  
  } else {missing <- "pairwise"}
  
  help_dat <- array(0, c(ncol(dataset), nrow(dataset), 2))
  
  for (i in 1:ncol(dataset)) {
    help_dat[i, , ] <- cbind(dataset[, i], rowMeans(dataset[, -i], na.rm = T))
  }
  
  ircor_samp <- apply(help_dat, c(1), .WishartCorTransform, n.iter = n.iter, n.burnin = n.burnin, missing)
  return(ircor_samp)
}



.WishartCorTransform <- function(x, n.iter, n.burnin, missing) {
  pairwise <- FALSE
  if (missing == "pairwise") {pairwise <- TRUE}
  tmp_cov <- Bayesrel:::covSamp(x, n.iter, n.burnin, pairwise)
  tmp_cor <- apply(tmp_cov, 1, cov2cor)
  out <- apply(tmp_cor, 2, function(x) mean(x[x!=1]))
  return(out)
}


# calculate the kolomogorov smirnov distances between some samples and the original sample
.ks.test.statistic <- function(x, y) {
  t <- stats::ks.test(x, y)
  t$statistic
}

# calculate the kublack leibler distance between two samples
.KLD.statistic <- function(x, y) {
  t <- LaplacesDemon::KLD(x, y)
  t$sum.KLD.py.px
}

