STReliability <- function(jaspResults, dataset, options) {

  #.libPaths("/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
  # print(.libPaths())
#  print('hoi')
#  Bayesrel::strel(matrix(rnorm(30), 10, 3))
#  print('hoi2')
  
    # if (!require("Bayesrel"))
  #   install.packages("Bayesrel")

	dataset <- .BayesianReliabilityReadData(dataset, options)

	.BayesianReliabilityCheckErrors(dataset, options)

	model <- .BayesianReliabilityMainResults(jaspResults, dataset, options)

	.BayesianReliabilityScaleTable(         jaspResults, model, options)
	.BayesianReliabilityScaleTableF(         jaspResults, model, options)
	.BayesianReliabilityItemTable(          jaspResults, model, options)
	.BayesianReliabilityItemTableF(          jaspResults, model, options)
	.BayesianReliabilityProbTable(          jaspResults, model, options)
	.BayesianReliabilityPosteriorPlot(      jaspResults, model, options)
	.BayesianReliabilityPosteriorPredictive(jaspResults, model, options)

	return()

}

# read data, check errors----
.BayesianReliabilityDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("alphaScale", "guttman2Scale", "glbScale", "mcDonaldScale")]),
    itemDroppedSelected = unlist(options[c("mcDonaldItem", "alphaItem", "guttman2Item", "glbItem")]),
    selectedEstimatorsF  = unlist(options[c("alphaScalef", "guttman2Scalef", "glbScalef", "mcDonaldScalef")]),
    itemDroppedSelectedF = unlist(options[c("mcDonaldItemf", "alphaItemf", "guttman2Itemf", "glbItemf")]),
    namesEstimators     = list(
      tables = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Greatest Lower Bound", "McDonald's \u03C9"),
      plots = list(expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]), "Greatest Lower Bound",
                   expression("McDonald's"~omega))
    )
  )
  # order to show in JASP
  derivedOptions[["order"]] <- match(c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2",
                                       "Greatest Lower Bound"), derivedOptions[["namesEstimators"]][["tables"]])


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
        nvar <- length(variables)
        key <- rep(1, nvar)
        key[match(.v(unlist(options[["reverseScaledItems"]])), nvar)] <- -1
        dataset <- dataset %*% diag(key, nvar, nvar)
      }

      model[["footnote"]] <- .BayesianReliabilityCheckLoadings(dataset, variables)
      relyFit <- try(Bayesrel::strel(x = dataset, estimates=c("alpha", "lambda2", "glb", "omega"), 
                                     n.iter = options[["noSamples"]], boot.n = options[["noSamplesf"]],
                                     item.dropped = TRUE, omega.freq.method = "pa"))
      # Consider stripping some of the contents of relyFit to reduce memory load
      if (inherits(relyFit, "try-error")) {

        model[["error"]] <- paste("The analysis crashed with the following error message:\n", relyFit)

      } else {

        # # should be done inside of Bayesrel?
        for (nm in names(relyFit[["bay"]][["ifitem"]][["samp"]])) {
          relyFit[["bay"]][["ifitem"]][["samp"]][[nm]] <- coda::mcmc(t(relyFit[["bay"]][["ifitem"]][["samp"]][[nm]]))
        }
        
        model[["dataset"]] <- dataset

        model[["relyFit"]] <- relyFit

        stateObj <- createJaspState(model)
        stateObj$dependOn(options = c("variables", "reverseScaledItems", "noSamples", "noSamplesf"))
        jaspResults[["modelObj"]] <- stateObj
      }
    }
  } else {
    print("model from state")
  }

	if (is.null(model[["error"]])) {
	  criState <- jaspResults[["criObj"]]$object
	  cfiState <- jaspResults[["cfiObj"]]$object
	  if (is.null(criState) && is.null(cfiState) && !is.null(relyFit)) {

	    scaleCri <- .BayesianReliabilityCalcCri(relyFit[["bay"]][["samp"]],             
	                                            options[["CredibleIntervalValue"]])
	    itemCri  <- .BayesianReliabilityCalcCri(relyFit[["bay"]][["ifitem"]][["samp"]], 
	                                            options[["CredibleIntervalValue"]])
	    
	    scaleCfi <- .BayesianReliabilityCalcCfi(relyFit[["freq"]][["boot"]],             
	                                            options[["ConfidenceIntervalValue"]])
	    
	    criState <- list(scaleCri = scaleCri, itemCri  = itemCri)
	    cfiState <- list(scaleCfi = scaleCfi)
	    jaspCriState <- createJaspState(criState)
	    jaspCfiState <- createJaspState(cfiState)
	    jaspCriState$dependOn(options = "CredibleIntervalValue", optionsFromObject = jaspResults[["modelObj"]])
	    jaspCfiState$dependOn(options = "ConfidenceIntervalValue", optionsFromObject = jaspResults[["modelObj"]])
	    jaspResults[["criObj"]] <- jaspCriState
	    jaspResults[["cfiObj"]] <- jaspCfiState
	    
	  }
	  model[["cri"]] <- criState
	  model[["cfi"]] <- cfiState
	}

  model[["derivedOptions"]] <- .BayesianReliabilityDerivedOptions(options)
  model[["itemsDropped"]] <- .unv(colnames(dataset))

	return(model)
}

.BayesianReliabilityCalcCri <- function(samps, criValue) {
  
  cri <- vector("list", length(samps))
  names(cri) <- names(samps)
  
  for (nm in names(samps)) {
    cri[[nm]] <- coda::HPDinterval(samps[[nm]], prob = criValue)
  }
  return(cri)
}

.BayesianReliabilityCalcCfi <- function(boot, cfiValue) {
  
  cfi <- vector("list", length(boot))
  names(cfi) <- names(boot)
  
  for (nm in names(boot)) {
    cfi[[nm]] <- quantile(boot[[nm]], prob = c(0+(1-cfiValue)/2, 1-(1-cfiValue)/2))
    names(cfi[[nm]]) <- c("lower", "upper")
  }
  return(cfi)
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


# tables ----
.BayesianReliabilityScaleTable <- function(jaspResults, model, options) {

  if (!is.null(jaspResults[["scaleTable"]])) {
    print("jaspResults[['scaleTable']] from state")
    return()
  }

  scaleTable <- createJaspTable("Bayesian Scale Reliability Statistics")
  scaleTable$dependOn(options = c("variables", "mcDonaldScale", "alphaScale", "guttman2Scale", "glbScale", "meanScale",
                                  "sdScale", "reverseScaledItems", "CredibleIntervalValue", "noSamples"))

  overTitle <- sprintf("%s%% Credible interval",
                       format(100*options[["CredibleIntervalValue"]], digits = 3, drop0trailing = TRUE))
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
	    postMean = unlist(relyFit$bay$est, use.names = FALSE),
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

.BayesianReliabilityScaleTableF <- function(jaspResults, model, options) {
  
  if (!is.null(jaspResults[["scaleTableF"]])) {
    print("jaspResults[['scaleTable']] from state")
    return()
  }
  
  scaleTableF <- createJaspTable("Frequentist Scale Reliability Statistics")
  scaleTableF$dependOn(options = c("variables", "mcDonaldScalef", "alphaScalef", "guttman2Scalef", "glbScalef", "meanScale",
                                  "sdScale", "reverseScaledItems", "ConfidenceIntervalValue", "noSamplesf"))
  
  overTitle <- sprintf("%s%% Confidence interval",
                       format(100*options[["ConfidenceIntervalValue"]], digits = 3, drop0trailing = TRUE))
  scaleTableF$addColumnInfo(name = "statistic", title = "Statistic",        type = "string")
  scaleTableF$addColumnInfo(name = "pointEst",  title = "Point Estimate", type = "number")
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

.BayesianReliabilityItemTable <- function(jaspResults, model, options) {

	if (!is.null(jaspResults[["itemTable"]]) || !any(model[["derivedOptions"]][["itemDroppedSelected"]])) {
	  print("jaspResults[['hasItemTable']] from state or not wanted")
		return()
	}

  derivedOptions <- model[["derivedOptions"]]
  itemDroppedSelected <- derivedOptions[["itemDroppedSelected"]]
  order <- derivedOptions[["order"]]
  overTitles <- derivedOptions[["namesEstimators"]][["tables"]][order]

  itemTable <- createJaspTable("Bayesian If Item Dropped Scale Reliability Statistics")
  itemTable$dependOn(options = c("variables",
                                 "mcDonaldScale", "alphaScale", "guttman2Scale", "glbScale", "meanScale", "sdScale",
                                 "mcDonaldItem",  "alphaItem",  "guttman2Item",  "glbItem",  "meanItem",  "sdItem",
                                 "reverseScaledItems"))
  itemTable$addColumnInfo(name = "variable", title = "if item dropped", type = "string")

  idxSelected <- which(itemDroppedSelected)

    for (i in idxSelected) {
    itemTable$addColumnInfo(name = paste0("postMean", i), title = "Posterior Median", type = "number", 
                            overtitle = overTitles[i])
    itemTable$addColumnInfo(name = paste0("lower", i),    title = "Lower",      type = "number", 
                            overtitle = overTitles[i])
    itemTable$addColumnInfo(name = paste0("upper", i),    title = "Upper",      type = "number", 
                            overtitle = overTitles[i])
  }

  relyFit <- model[["relyFit"]]
  if (!is.null(relyFit)) {
    cris <- model[["cri"]][["itemCri"]]
    tb <- data.frame(variable = model[["itemsDropped"]])
    for (i in idxSelected) {
      idx <- order[i]
      newtb <- cbind(postMean = relyFit$bay$ifitem$est[[idx]], cris[[idx]])
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

.BayesianReliabilityItemTableF <- function(jaspResults, model, options) {
  
  if (!is.null(jaspResults[["itemTableF"]]) || !any(model[["derivedOptions"]][["itemDroppedSelectedF"]])) {
    print("jaspResults[['hasItemTable']] from state or not wanted")
    return()
  }
  
  derivedOptions <- model[["derivedOptions"]]
  itemDroppedSelectedF <- derivedOptions[["itemDroppedSelectedF"]]
  order <- derivedOptions[["order"]]
  overTitles <- derivedOptions[["namesEstimators"]][["tables"]][order]
  
  itemTableF <- createJaspTable("Frequentist If Item Dropped Scale Reliability Statistics")
  itemTableF$dependOn(options = c("variables",
                                 "mcDonaldScalef", "alphaScalef", "guttman2Scalef", "glbScalef", "meanScalef", "sdScale",
                                 "mcDonaldItemf",  "alphaItemf",  "guttman2Itemf",  "glbItemf",  "meanItem",  "sdItem",
                                 "reverseScaledItems"))
  itemTableF$addColumnInfo(name = "variable", title = "if item dropped", type = "string")
  
  idxSelectedF <- which(itemDroppedSelectedF)
  
  for (i in idxSelectedF) {
    itemTableF$addColumnInfo(name = paste0("pointEst", i), title = "Point Estimate", type = "number", 
                            overtitle = overTitles[i])

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


.BayesianReliabilityProbTable <- function(jaspResults, model, options) {

  if (!is.null(jaspResults[["probTable"]]) || !options[["probTable"]]) {
    print("jaspResults[['probTable']] from state or not wanted")
		return()
  }

  probTable <- createJaspTable(sprintf("Probability that Reliability Statistic is Larger than %.2f", 
                                       options[["probTableValue"]]))
  probTable$dependOn(options = c("variables", "mcDonaldScale", "alphaScale", "guttman2Scale", "glbScale", "meanScale",
                                  "sdScale", "reverseScaledItems", "probTableValue", "probTable"))
  probTable$addColumnInfo(name = "statistic", title = "Statistic",   type = "string")
	probTable$addColumnInfo(name = "probs",     title = "Probability", type = "number")

  relyFit <- model[["relyFit"]]
	derivedOptions <- model[["derivedOptions"]]
	opts     <- derivedOptions[["namesEstimators"]][["tables"]]
	order    <- derivedOptions[["order"]]
	selected <- derivedOptions[["selectedEstimators"]]
	idx      <- which(selected)

  if (!is.null(relyFit)) {
    probs <- numeric(sum(selected))

    for (i in seq_along(idx)) {
      probs[i] <- mean(relyFit[["bay"]][["samp"]][[idx[i]]] > options[["probTableValue"]])
    }
    df <- data.frame(statistic = opts[idx], probs = probs)
    probTable$setData(df)
  } else if (sum(selected) > 0) {
    probTable[["statistic"]] <- opts[idx]
  }

  jaspResults[["probTable"]] <- probTable
  return()
}

# plots ----
.BayesianReliabilityPosteriorPlot <- function(jaspResults, model, options) {

	if (!options[["plotPosterior"]])
		return()

  plotContainer <- jaspResults[["plotContainer"]]
  if (is.null(plotContainer)) {
    print("Plotcontainer remade")
    plotContainer <- createJaspContainer("Posteriors Plots")
    plotContainer$dependOn(options = c("variables", "reverseScaledItems", "plotPosterior", "shadePlots",
                                       "probTable", "probTableValue", "cutoff", "fixXRange", 
                                       "cutoffValue1", "cutoffValue2", "dispPrior", "noSamples"))
    jaspResults[["plotContainer"]] <- plotContainer
  } else {
    print("Plotcontainer from state")
  }

	derivedOptions <- model[["derivedOptions"]]
	order     <- derivedOptions[["order"]]
	indices   <- which(derivedOptions[["selectedEstimators"]][order])
	nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]][order]
	nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]][order]

	relyFit  <- model[["relyFit"]]
	# probably better to do this once directly after computation!
	scaleCri <- model[["cri"]][["scaleCri"]][order]
	relyFit[["bay"]][["samp"]] <- relyFit[["bay"]][["samp"]][order]
	n.item <- dim(relyFit$bay$covsamp)[2]
	prior <- priors[[as.character(n.item)]][c(4, 1, 2, 3)]


	if (options[["shadePlots"]] && options[["probTable"]]) {
	  plotContainer$dependOn("probTableValue")
	  shadePlots <- options[["probTableValue"]]
	} else {
	  shadePlots <- NULL
	}

	if (options[["cutoff"]]) {
    plotContainer$dependOn("cutoffValue1")
    plotContainer$dependOn("cutoffValue2")
    cutoffs <- c(options[["cutoffValue1"]], options[["cutoffValue2"]])
	} else {
	  cutoffs <- NULL
	}

	if (!is.null(relyFit)) {
	  for (i in indices) {
	    if (is.null(plotContainer[[nmsObjs[i]]])) {

	      p <- .BayesianReliabilityMakeSinglePosteriorPlot(relyFit, scaleCri, i, nmsLabs[[i]], options[["fixXRange"]],
	                                                       shadePlots, cutoffs, options[["dispPrior"]], prior)
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
	  plotContainer[["Posterior Distribution"]] <- createJaspPlot()
	}

	return()
}

.BayesianReliabilityMakeSinglePosteriorPlot <- function(relyFit, scaleCri, i, nms, fixXRange, 
                                                        shade = NULL, cutoffs = NULL, priorTrue, priorSample) {

  # TODO: consider precomputing all densities (maybe with kernsmooth?) and reducing memory that way
  pr <- priorSample[[i]]
  if (fixXRange) {
    d <- stats::density(relyFit$bay$samp[[i]], from = 0, to = 1, n = 2^10)
  } else {
    d <- stats::density(relyFit$bay$samp[[i]], n = 2^10)
  }
	datDens <- data.frame(x = d$x, y = d$y)
	datPrior <- data.frame(x = pr$x, y = pr$y)
	


	xBreaks <- JASPgraphs::getPrettyAxisBreaks(datDens$x)
	print(xBreaks)
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
		ggplot2::scale_x_continuous(name = "Reliability", breaks = xBreaks, expand = xExpand)
	
	if (!is.null(shade)) {
	  datFilter <- datDens[datDens[["x"]] >= shade, ]
	  g <- g + ggplot2::geom_ribbon(data = datFilter, mapping = ggplot2::aes(ymin = 0,ymax = y), 
	                                fill = "grey", alpha = 0.95) +
	           ggplot2::geom_line(size = .85)
	}
	
	if (priorTrue) {
	  # g <- g + ggplot2::geom_hline(yintercept = 1)
	  g <- g + ggplot2::geom_line(data = datPrior, mapping = ggplot2::aes(x = x, y = y),
	                              linetype = 3, size = .85) +
	           ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, expand = xExpand,
	                                       limits= c(min(xBreaks), max(xBreaks)))
	           # ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, expand = xExpand, 
	           #                             limits= c(min(d$x), max(d$x)))
	  

	  
	}
	
	if (!is.null(cutoffs)) {
	  g <- g +
	    ggplot2::geom_segment(ggplot2::aes(x = cutoffs[1], y = 0, xend = cutoffs[1], yend = max(d$y)), 
	                          color = "grey60", linetype = 1, alpha = .5, size = .3) +
	    ggplot2::geom_segment(ggplot2::aes(x = cutoffs[2], y = 0, xend = cutoffs[2], yend = max(d$y)), 
	                          color = "grey60", linetype = 1, alpha = .5, size = .3) +
	    ggplot2::geom_line(size = .85)
	}
	


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
    ll <- relyFit[["bay"]][["loadings"]]
    rr <- relyFit[["bay"]][["resid.var"]]
    
    cimpl <- ll %*% t(ll) + diag(rr)
    cobs <- cov(dataset)
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
      ggplot2::geom_ribbon(ggplot2::aes(ymin = eigen_sim_low, ymax = eigen_sim_up), fill = "gray80") +
      ggplot2::geom_line(ggplot2::aes(x = number, y = eigen_sim_low), colour = "gray", linetype = 2) +
      ggplot2::geom_line(ggplot2::aes(x = number, y = eigen_sim_up), colour = "gray", linetype = 2) +
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
  plot <- createJaspPlot(plot = g, title = "Posterior Predictive Check Omega")
  plot$dependOn(options = c("dispPPC", "mcDonaldScale", "reverseScaledItems"))
  jaspResults[["OmegaPosteriorPredictive"]] <- plot
}
