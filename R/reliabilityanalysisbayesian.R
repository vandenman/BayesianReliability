BaySTReliability <- function(jaspResults, dataset, options) {

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
	.BayesianReliabilityItemTable(          jaspResults, model, options)
	.BayesianReliabilityProbTable(          jaspResults, model, options)
	.BayesianReliabilityPosteriorPlot(      jaspResults, model, options)
	.BayesianReliabilityPosteriorPredictive(jaspResults, model, options)

	return()

}

# read data, check errors----
.BayesianReliabilityDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("alphaScale", "guttman2Scale", "glbScale", "mcDonaldScale",
                                           "averageInterItemCor", "meanScale", "sdScale")]),
    itemDroppedSelected = unlist(options[c("mcDonaldItem", "alphaItem", "guttman2Item", "glbItem",
                                           "itemRestCor", "meanItem", "sdItem")]),

    namesEstimators     = list(
      tables = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Greatest Lower Bound", 
                 "McDonald's \u03C9", "Average interitem correlation", "mean", "sd"),
      tables_item = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Greatest Lower Bound", 
                 "McDonald's \u03C9", "Item-rest correlation", "mean", "sd"),
      coefficients = c("Cronbach's \u03B1", "Guttman's \u03BB2", "Greatest Lower Bound", 
                       "McDonald's \u03C9", "Item-rest correlation"),
      plots = list(expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]), 
                   "Greatest Lower Bound", expression("McDonald's"~omega))
    )

  )
  # order to show in JASP
  derivedOptions[["order"]] <- match(c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2",
                                       "Greatest Lower Bound", "Average interitem correlation", "mean", "sd"), 
                                     derivedOptions[["namesEstimators"]][["tables"]])
  derivedOptions[["order_item"]] <- match(c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2",
                                       "Greatest Lower Bound", "Item-rest correlation", "mean", "sd"), 
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
        nvar <- length(variables)
        key <- rep(1, nvar)
        key[match(.v(unlist(options[["reverseScaledItems"]])), nvar)] <- -1
        dataset <- dataset %*% diag(key, nvar, nvar)
      }

      model[["footnote"]] <- .BayesianReliabilityCheckLoadings(dataset, variables)
      relyFit <- try(Bayesrel::strel(x = dataset, estimates=c("alpha", "lambda2", "glb", "omega"), 
                                     n.iter = options[["noSamples"]], freq = F,
                                     item.dropped = TRUE))
      
      # add the scale info
      corsamp <- apply(relyFit$bay$covsamp, 1, cov2cor)
      relyFit$bay$samp$avg_cor <- coda::mcmc(apply(corsamp, 2, function(x) mean(x[x!=1])))
      relyFit$bay$est$avg_cor <- median(relyFit$bay$samp$avg_cor)
      
      relyFit$bay$samp$mean <- coda::mcmc(c(0, 0))
      relyFit$bay$est$mean <- mean(dataset)
      relyFit$bay$samp$sd <- coda::mcmc(c(0, 0))
      relyFit$bay$est$sd <- sd(apply(dataset, 2, mean))
      
      # now the item statistics
      ###### how to do this? I dont know, the item-rest correlation needs some thinking, mean and sd are straightforward 
      relyFit$bay$ifitem$samp$ircor <- .reliabilityItemRestCor(dataset, options[["noSamples"]])
      relyFit$bay$ifitem$est$ircor <- apply(relyFit$bay$ifitem$samp$ircor, 1, median)
      
      relyFit$bay$ifitem$est$mean <- apply(dataset, 2, mean)
      relyFit$bay$ifitem$est$sd <- apply(dataset, 2, sd)  
      relyFit$bay$ifitem$samp$mean <- coda::mcmc(matrix(0, ncol(dataset), 2))
      relyFit$bay$ifitem$samp$sd <- coda::mcmc(matrix(0, ncol(dataset), 2))

      
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
        stateObj$dependOn(options = c("variables", "reverseScaledItems", "noSamples"))
        jaspResults[["modelObj"]] <- stateObj
      }
    }
  } else {
    print("model from state")
  }

	if (is.null(model[["error"]])) {
	  criState <- jaspResults[["criObj"]]$object
	  if (is.null(criState) && !is.null(relyFit)) {

	    scaleCri <- .BayesianReliabilityCalcCri(relyFit[["bay"]][["samp"]],             
	                                            options[["CredibleIntervalValue"]])
	    itemCri  <- .BayesianReliabilityCalcCri(relyFit[["bay"]][["ifitem"]][["samp"]], 
	                                            options[["CredibleIntervalValue"]])

	    
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

.BayesianReliabilityCalcCri <- function(samps, criValue) {

  cri <- vector("list", length(samps))
  names(cri) <- names(samps)
  
  for (nm in names(samps)) {
    cri[[nm]] <- coda::HPDinterval(samps[[nm]], prob = criValue)
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


# tables ----
.BayesianReliabilityScaleTable <- function(jaspResults, model, options) {

  if (!is.null(jaspResults[["scaleTable"]])) {
    print("jaspResults[['scaleTable']] from state")
    return()
  }

  scaleTable <- createJaspTable("Bayesian Scale Reliability Statistics")
  scaleTable$dependOn(options = c("variables", "mcDonaldScale", "alphaScale", "guttman2Scale", 
                                  "glbScale", "reverseScaledItems", "CredibleIntervalValue", "noSamples", 
                                  "averageInterItemCor", "meanScale", "sdScale"))

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


.BayesianReliabilityItemTable <- function(jaspResults, model, options) {

	if (!is.null(jaspResults[["itemTable"]]) || !any(model[["derivedOptions"]][["itemDroppedSelected"]])) {
	  print("jaspResults[['hasItemTable']] from state or not wanted")
		return()
	}

  derivedOptions <- model[["derivedOptions"]]
  itemDroppedSelected <- derivedOptions[["itemDroppedSelected"]]
  order <- derivedOptions[["order_item"]]
  overTitles <- format(derivedOptions[["namesEstimators"]][["tables_item"]][order], digits = 3, drop0trailing = T)
  
  cred <- format(100*options[["CredibleIntervalValue"]], digits = 3, drop0trailing = TRUE)
  itemTable <- createJaspTable("Bayesian If Item Dropped Scale Reliability Statistics")
  itemTable$dependOn(options = c("variables",
                                 "mcDonaldScale", "alphaScale", "guttman2Scale", "glbScale", 
                                 "averageInterItemCor", "meanScale", "sdScale",
                                 "mcDonaldItem",  "alphaItem",  "guttman2Item", "glbItem",
                                 "reverseScaledItems", "CredibleIntervalValue", 
                                 "itemRestCor", "meanItem", "sdItem"))
  itemTable$addColumnInfo(name = "variable", title = "Item", type = "string")

  idxSelected <- which(itemDroppedSelected)
  estimators <- derivedOptions[["namesEstimators"]][["tables_item"]][order]
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]
  
  for (i in idxSelected) {
    if (estimators[i] %in% coefficients) {
      itemTable$addColumnInfo(name = paste0("postMean", i), title = "Posterior Median", type = "number", 
                              overtitle = overTitles[i])
      itemTable$addColumnInfo(name = paste0("lower", i), title = paste0("Lower ", cred, "%"), type = "number", 
                              overtitle = overTitles[i])
      itemTable$addColumnInfo(name = paste0("upper", i), title = paste0("Upper ", cred, "%"), type = "number", 
                              overtitle = overTitles[i])
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
      if (idx %in% c(1:5)) {
        newtb <- cbind(postMean = relyFit$bay$ifitem$est[[idx]], cris[[idx]])
      } else {
        newtb <- cbind(postMean = relyFit$bay$ifitem$est[[idx]])
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
  probTable$dependOn(options = c("variables", "mcDonaldScale", "alphaScale", "guttman2Scale",
                                 "glbScale", "reverseScaledItems", "probTableValueLow", "probTable",
                                 "probTableValueHigh"))
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
  
	n.item <- dim(relyFit$bay$covsamp)[2]
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
      probsPost[i] <- mean(relyFit[["bay"]][["samp"]][[idx[i]]] > options[["probTableValueLow"]]) -
                      mean(relyFit[["bay"]][["samp"]][[idx[i]]] > options[["probTableValueHigh"]])
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

# plots ----
.BayesianReliabilityPosteriorPlot <- function(jaspResults, model, options) {

	if (!options[["plotPosterior"]])
		return()

  plotContainer <- jaspResults[["plotContainer"]]
  if (is.null(plotContainer)) {
    print("Plotcontainer remade")
    plotContainer <- createJaspContainer("Posteriors Plots")
    plotContainer$dependOn(options = c("variables", "reverseScaledItems", "plotPosterior", "shadePlots",
                                       "probTable", "probTableValueLow", "probTableValueHigh", "fixXRange", 
                                       "dispPrior", "noSamples"))
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
	prior <- priors[[as.character(n.item)]][c(4, 1, 2, 3)] ##### change this when more estiamtors are included!!!


	if (options[["shadePlots"]] && options[["probTable"]]) {
	  plotContainer$dependOn("probTableValueLow")
	  plotContainer$dependOn("probTableValueHigh")
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
	  plotContainer[["Posterior Distribution"]] <- createJaspPlot()
	}

	return()
}

.BayesianReliabilityMakeSinglePosteriorPlot <- function(relyFit, scaleCri, i, nms, fixXRange, 
                                                        shade = NULL, priorTrue, priorSample) {

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
	  datFilter <- datDens[datDens[["x"]] >= shade[1] & datDens[["x"]] <= shade[2], ]
	  g <- g + ggplot2::geom_ribbon(data = datFilter, mapping = ggplot2::aes(ymin = 0,ymax = y), 
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

.reliabilityItemRestCor <- function(dataset, n.iter) {
  help_dat <- array(0, c(ncol(dataset), nrow(dataset), 2))
  
  for (i in 1:ncol(dataset)) {
    idx <- seq(1, ncol(dataset))
    idx <- idx[idx!=i]
    help_dat[i, , ] <- cbind(dataset[, i], apply(dataset[, idx], 1, mean))
  }
  
  ircor_samp <- apply(help_dat, c(1), .WishartCorTransform, n.iter = n.iter, n.burnin = 50)
  return(t(ircor_samp))
}


.WishartCorTransform <- function(x, n.iter, n.burnin = 50) {
  tmp_cov <- Bayesrel:::covSamp(x, n.iter, n.burnin)
  tmp_cor <- apply(tmp_cov, 1, cov2cor)
  out <- apply(tmp_cor, 2, function(x) mean(x[x!=1]))
  return(out)
}
