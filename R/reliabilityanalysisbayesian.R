BayesianReliability <- function(jaspResults, dataset, options) {

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

	.BayesianReliabilityScaleTable(   jaspResults, model, options)
	.BayesianReliabilityItemTable(    jaspResults, model, options)
	.BayesianReliabilityProbTable(    jaspResults, model, options)
	.BayesianReliabilityPosteriorPlot(jaspResults, model, options)

	return()

}

# read data, check errors----
.BayesianReliabilityDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("alphaScale", "guttman2Scale", "glbScale", "mcDonaldScale")]),
    itemDroppedSelected = unlist(options[c("mcDonaldItem", "alphaItem", "guttman2Item", "glbItem")]),
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
      relyFit <- try(Bayesrel::strel(x = dataset, estimates=c("alpha", "lambda2", "glb", "omega"), freq = FALSE, n.iter = options[["noSamples"]],
                                     item.dropped = TRUE, prior.samp = TRUE))
      # Consider stripping some of the contents of relyFit to reduce memory load
      if (inherits(relyFit, "try-error")) {

        model[["error"]] <- paste("The analysis crashed with the following error message:\n", relyFit)

      } else {

        # # should be done inside of Bayesrel?
        # for (nm in names(relyFit[["bay"]][["ifitem"]][["samp"]])) {
        #   relyFit[["bay"]][["ifitem"]][["samp"]][[nm]] <- coda::mcmc(t(relyFit[["bay"]][["ifitem"]][["samp"]][[nm]]))
        # }
        
        # temporary workarounds until glb is fixed
        # relyFit[["bay"]][["est"]] <- c(relyFit[["bay"]][["est"]], list(bayes.glb = 0))[c(1, 2, 4, 3)]
        # 
        # relyFit[["bay"]][["samp"]] <- c(relyFit[["bay"]][["samp"]], 
        #                                 list(bayes.glb = coda::mcmc(numeric(options[["noSamples"]]-50))))[c(1, 2, 4, 3)]
        # relyFit[["bay"]][["ifitem"]][["samp"]] <- c(relyFit[["bay"]][["ifitem"]][["samp"]], 
        #                                             list(glb = coda::mcmc(matrix(0, 
        #                                                                          options[["noSamples"]]-50, 
        #                                                                          dim(relyFit[["bay"]][["covsamp"]])[2])
        #                                                                   )))[c(1, 2, 4, 3)]
        # relyFit[["priors"]] <- c(relyFit[["priors"]], list(glb = numeric(2000)))[c(1, 2, 4, 3)]
        
        ####

        model[["relyFit"]] <- relyFit

        stateObj <- createJaspState(model)
        stateObj$dependOn(options = c("variables", "reverseScaledItems"))
        jaspResults[["modelObj"]] <- stateObj
      }
    }
  } else {
    print("model from state")
  }

	if (is.null(model[["error"]])) {
	  criState <- jaspResults[["criObj"]]$object
	  if (is.null(criState) && !is.null(relyFit)) {

	    scaleCri <- .BayesianReliabilityCalcCri(relyFit[["bay"]][["samp"]],             options[["CredibleIntervalValue"]])
	    itemCri  <- .BayesianReliabilityCalcCri(relyFit[["bay"]][["ifitem"]][["samp"]], options[["CredibleIntervalValue"]])
	    
	    criState <- list(scaleCri = scaleCri, itemCri  = itemCri)
	    jaspCriState <- createJaspState(criState)
	    jaspCriState$dependOn(options = "CredibleIntervalValue", optionsFromObject = jaspResults[["modelObj"]])
	    jaspResults[["criObj"]] <- jaspCriState

	  }
	  model[["cri"]] <- criState
	}

  model[["derivedOptions"]] <- .BayesianReliabilityDerivedOptions(options)
  model[["itemsDropped"]] <- .unv(colnames(dataset))

	return(model)
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

.BayesianReliabilityCalcCri <- function(samps, criValue) {

  cri <- vector("list", length(samps))
  names(cri) <- names(samps)

  for (nm in names(samps)) {
    cri[[nm]] <- coda::HPDinterval(samps[[nm]], prob = criValue)
  }
  return(cri)
}

# tables ----
.BayesianReliabilityScaleTable <- function(jaspResults, model, options) {

  if (!is.null(jaspResults[["scaleTable"]])) {
    print("jaspResults[['scaleTable']] from state")
    return()
  }

  scaleTable <- createJaspTable("Scale Reliability Statistics")
  scaleTable$dependOn(options = c("variables", "mcDonaldScale", "alphaScale", "guttman2Scale", "glbScale", "meanScale",
                                  "sdScale", "reverseScaledItems", "CredibleIntervalValue"))

  overTitle <- sprintf("%s%% Credible interval",
                       format(100*options[["CredibleIntervalValue"]], digits = 3, drop0trailing = TRUE))
	scaleTable$addColumnInfo(name = "statistic", title = "Statistic",      type = "string")
	scaleTable$addColumnInfo(name = "postMean",  title = "PosteriorMedian", type = "number")
	scaleTable$addColumnInfo(name = "lower",     title = "Lower",          type = "number", overtitle = overTitle)
	scaleTable$addColumnInfo(name = "upper",     title = "Upper",          type = "number", overtitle = overTitle)

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

.BayesianReliabilityItemTable <- function(jaspResults, model, options) {

	if (!is.null(jaspResults[["itemTable"]]) || !any(model[["derivedOptions"]][["itemDroppedSelected"]])) {
	  print("jaspResults[['hasItemTable']] from state or not wanted")
		return()
	}

  derivedOptions <- model[["derivedOptions"]]
  itemDroppedSelected <- derivedOptions[["itemDroppedSelected"]]
  order <- derivedOptions[["order"]]
  overTitles <- derivedOptions[["namesEstimators"]][["tables"]][order]

  itemTable <- createJaspTable("Scale Reliability Statistics")
  itemTable$dependOn(options = c("variables",
                                 "mcDonaldScale", "alphaScale", "guttman2Scale", "glbScale", "meanScale", "sdScale",
                                 "mcDonaldItem",  "alphaItem",  "guttman2Item",  "glbItem",  "meanItem",  "sdItem",
                                 "reverseScaledItems"))
  itemTable$addColumnInfo(name = "variable", title = "if item dropped", type = "string")

  idxSelected <- which(itemDroppedSelected)

    for (i in idxSelected) {
    itemTable$addColumnInfo(name = paste0("postMean", i), title = "Posterior Median", type = "number", overtitle = overTitles[i])
    itemTable$addColumnInfo(name = paste0("lower", i),    title = "Lower Cri",      type = "number", overtitle = overTitles[i])
    itemTable$addColumnInfo(name = paste0("upper", i),    title = "Upper Cri",      type = "number", overtitle = overTitles[i])
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

.BayesianReliabilityProbTable <- function(jaspResults, model, options) {

  if (!is.null(jaspResults[["probTable"]]) || !options[["probTable"]]) {
    print("jaspResults[['probTable']] from state or not wanted")
		return()
  }

  probTable <- createJaspTable(sprintf("Probability that Reliability Statistic is Larger than %.2f", options[["probTableValue"]]))
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
                                       "cutoffValue1", "cutoffValue2", "dispPrior"))
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
	                                                       shadePlots, cutoffs, options[["dispPrior"]])
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

.BayesianReliabilityMakeSinglePosteriorPlot <- function(relyFit, scaleCri, i, nms, fixXRange, shade = NULL, cutoffs = NULL, prior) {

  # TODO: consider precomputing all densities (maybe with kernsmooth?) and reducing memory that way
  if (fixXRange) {
    d <- stats::density(relyFit$bay$samp[[i]], from = 0, to = 1, n = 2^11, adjust = 1.5)
    pr <- stats::density(relyFit$priors[[i]], from = 0, to = 1, n = 2^11)
  } else {
    d <- stats::density(relyFit$bay$samp[[i]], n = 2^11, adjust = 1.5)
    pr <- stats::density(relyFit$priors[[i]], n = 2^11)
  }
	datDens <- data.frame(x = d$x, y = d$y)
	datPrior <- data.frame(x = pr$x, y = pr$y)

	xBreaks <- JASPgraphs::getPrettyAxisBreaks(datDens$x)
	# max height posterior is at 90% of plot area; remained is for credible interval
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
		# ggrepel::geom_text_repel(data = datTxt, mapping = ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE, segment.alpha = 0) +
	  ggplot2::geom_text(data = datTxt, mapping = ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE) +
		ggplot2::scale_y_continuous(name = "Density", breaks = yBreaks, limits = range(yBreaks)) +
		ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, expand = xExpand)
	
	if (!is.null(shade)) {
	  datFilter <- datDens[datDens[["x"]] >= shade, ]
	  g <- g + ggplot2::geom_ribbon(data = datFilter, mapping = ggplot2::aes(ymin = 0,ymax = y), fill = "grey",
	                                alpha = 0.95) +
	           ggplot2::geom_line(size = .85)
	}
	
	if (prior) {
	  g <- g + ggplot2::geom_hline(yintercept = 1)
	  # g <- g + ggplot2::ggplot(data = datPrior, mapping = ggplot2::aes(x = x, y = y))
	}
	
	if (!is.null(cutoffs)) {
	  g <- g +
	    ggplot2::geom_segment(ggplot2::aes(x = cutoffs[1], y = 0, xend = cutoffs[1], yend = ymax*.875), color = "grey60", linetype = 2) +
	    ggplot2::geom_segment(ggplot2::aes(x = cutoffs[2], y = 0, xend = cutoffs[2], yend = ymax*.875), color = "grey60", linetype = 2) 
	}
	


	return(JASPgraphs::themeJasp(g))

}

.BayesianReliabilityPosteriorPredictive <- function(jaspResults, model, options) {

	if (!options[["plotPosterior"]] || !is.null(jaspResults[["OmegaPosteriorPredictive"]]))
		return()

  # browser()

  relyfit <- model[["relyFit"]]
  if (is.null(relyfit)) {
    g <- NULL
  } else {
    # df <- data.frame(
    #   x = ,
    #   y =
    # ),
    g <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x, y = y)) +
      geom_line() +
      geom_point() +
      geom_ribbon()
  }
  plot <- createJaspPlot(plot = g, title = "Posterior Predictive Check Omega")
  plot$dependOn(options = c("plotPosterior", "variables", "reverseScaledItems", "CredibleIntervalValue"))
  jaspResults[["OmegaPosteriorPredictive"]] <- plot
}
