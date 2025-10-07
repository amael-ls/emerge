
#### Aim of prog: Gathering tool functions, note that the required packages are not loaded here
## Table of Contents
#	- getParams: get fixed values parameters
#	- lazyTrace: Bayesplot is having troubles on my mac (Arial font not always found), so I create my own traces plot
#	- lazyPosterior: plot the prior and posterior of a parameter
#	- gamma_sum_approx: approximates the sum of gamma distributions according to https://doi.org/10.1007/bf02481123
#
## Comments
# This R file contains tool functions only that help me to analyse the results and do some check-up. Note that
#	some functions are quite similar to bayesplot, but I use base plot rather than ggplot. Moreover, Bayesplot
#	is having troubles on my mac (Arial font not always found), so I created my own traces and posterior plots.

#### Tool functions
## Get fixed values parameters (will not work for draws with third dimension > 1)
getParams = function(model_cmdstan, params_names, type = "mean", ...)
{
	if (!(type %in% c("all", "chain-iter", "mean", "median", "quantile")))
		stop("Unknown type. Please choose all, iter-chain, mean, median, or quantile")

	args = list(...)

	if (type %in% c("chain-iter", "mean", "median"))
	{
		vals = numeric(length(params_names))
		names(vals) = params_names

		if (type == "chain-iter")
		{
			if (!all(c("iter", "chain") %in% names(args)))
				stop("You must provide iter and chain when using the option iter-chain")
			iter = args[["iter"]]
			chain = args[["chain"]]

			if (iter > model_cmdstan$metadata()$iter_sampling)
				stop("iter cannot be larger than iter_sampling")

			if (chain > model_cmdstan$num_chains())
				stop("chain cannot be larger than num_chains")

			for (i in seq_along(params_names))
			{
				draws = model_cmdstan$draws(params_names[i])
				if (dim(draws)[3] != 1)
					stop("This function is not yet designed to handle a multi-params")
				vals[i] = draws[iter, chain, 1]
			}
			return(vals)
		} else {
			for (i in seq_along(params_names))
			{
				vals[i] = ifelse(type == "mean",
					mean(model_cmdstan$draws(params_names[i])),
					median(model_cmdstan$draws(params_names[i])))
			}
			return(vals)
		}
	} else if (type == "quantile") {
		if (!("probs" %in% names(args)))
		{
			probs = c(0, 0.05, 0.5, 0.95, 1)
			vals = data.table(parameter = params_names, min = 0, q05 = 0, med = 0, avg = 0, q95 = 0, max = 0)
		}
		if ("probs" %in% names(args))
		{
			warning("I have not yet coded the general output, so far I use probs = c(0, 0.05, 0.5, 0.95, 1)")
			probs = c(0, 0.05, 0.5, 0.95, 1)
			vals = data.table(parameter = params_names, min = 0, q05 = 0, med = 0, avg = 0, q95 = 0, max = 0)
		}
		for (i in seq_along(params_names))
		{
			vals[i, c("min", "q05", "med", "q95", "max") := as.list(quantile(model_cmdstan$draws(params_names[i]), c(0, 0.05, 0.5, 0.95, 1)))]
			vals[i, avg := mean(model_cmdstan$draws(params_names[i]))]
		}
		return(vals)
	}

	vals = model_cmdstan$draws(params_names)

	return(vals)
}

## Bayesplot is having troubles on my mac (Arial font not always found), so I create my own traces plot
lazyTrace = function(draws, filename = NULL, ...)
{
	if (!is.array(draws) && !all(class(draws) %in% c("draws_array", "draws", "array")))
		stop("The class of draws should be either array, or compatible with cmdstanr (draws_array, draws, array)")

	n_chains = dim(draws)[2]
	n_iter = dim(draws)[1]
	colours = MetBrewer::met.brewer("Hokusai3", n_chains)
	colours_str = grDevices::colorRampPalette(colours)(n_chains)

	min_val = min(draws)
	max_val = max(draws)

	providedArgs = list(...)
	nbArgs = length(providedArgs)

	ls_names = names(providedArgs)

	val_ind = stri_detect(str = ls_names, regex = "val[[:digit:]]")
	xlab_ind = (ls_names == "xlab")
	ylab_ind = (ls_names == "ylab")
	main_ind = (ls_names == "main")
	label_ind = stri_detect(str = ls_names, regex = "label")
	iter_ind = stri_detect(str = ls_names, regex = "iter[[:digit:]]")

	scaling_ind = (ls_names == "scaling")
	if (any(scaling_ind)) scaling = providedArgs[["scaling"]] else scaling = 1

	if (any(label_ind))
		par(mar = c(5, 4, 4, 4))

	# Plot
	if (!is.null(filename))
	{
		pdf(paste0(filename, ".pdf"))
		print(paste0("Figure saved under the name: ", filename, ".pdf"))
	}

	plot(0, pch = "", xlim = c(0, n_iter), ylim = scaling*c(min_val, max_val), axes = TRUE, bg = "transparent",
		xlab = ifelse(any(xlab_ind), providedArgs[["xlab"]], ""),
		ylab = ifelse(any(ylab_ind), providedArgs[["ylab"]], ""),
		main = ifelse(any(main_ind), providedArgs[["main"]], ""))

	for (chain in 1:n_chains)
	{
		if (all(class(draws) %in% c("draws_array", "draws", "array")))
			lines(1:n_iter, scaling*draws[, chain, ], type = "l", col = colours_str[chain])
		if (is.array(draws) && !all(class(draws) %in% c("draws_array", "draws", "array")))
			lines(1:n_iter, scaling*draws[, chain], type = "l", col = colours_str[chain])
	}

	if (any(val_ind))
	{
		for (val in ls_names[val_ind])
			abline(h = scaling*providedArgs[[val]], col = "#CD212A", lwd = 4)

		if (any(label_ind))
		{
			num_vals = stri_sub(str = ls_names[val_ind], from = stri_locate(str = ls_names[val_ind], regex = "val")[, "end"] + 1)
			for (label in ls_names[label_ind])
			{
				num_label = stri_sub(str = label, from = stri_locate(str = label, regex = "label")[, "end"] + 1)
				corresponding_val = (ls_names[val_ind])[num_vals == num_label]
				axis(4, at = scaling*providedArgs[[corresponding_val]], providedArgs[[label]], las = 1)
			}
		}
	}

	if (any(iter_ind))
	{
		for (iter in ls_names[iter_ind])
			abline(v = providedArgs[[iter]], col = "#66666644", lwd = 0.2)
	}

	if (!is.null(filename))
		dev.off()
}

## Function to plot the prior and posterior of a parameter
lazyPosterior = function(draws, fun = NULL, expand_bounds = FALSE, filename = NULL, multi = FALSE, ls_nfi = NULL, ...)
{
	# Dealing with draws as a list
	if (is.list(draws))
	{
		# --- Check-up
		nDraws = length(draws)
		len = data.table(nIter = integer(nDraws), nChains = integer(nDraws), nVariables = integer(nDraws), variableNames = character(nDraws))
		for (i in seq_along(draws))
		{
			if (!all(class(draws[[i]]) == c("draws_array", "draws", "array")))
				stop(paste0("draws[", i, "] is not an array extracted from a CmdStanMCMC object"))

			len[i, c("nIter", "nChains", "nVariables", "variableNames") := append(as.list(dim(draws[[i]])), summary(draws[[i]])[["variable"]])]
		}
		if (unique(len)[, .N] != 1)
			stop("Dimensions or variables' name mismatch within the provided draws")

		rm(len)
		# --- Format
		draws_unlist = posterior::bind_draws(draws[[1]], along = "iteration")
		for (j in 2:length(draws))
			draws_unlist = posterior::bind_draws(draws_unlist, draws[[i]], along = "iteration")
		draws = draws_unlist
		rm(draws_unlist)
	}

	# Check-up
	if (!all(class(draws) == c("draws_array", "draws", "array")))
		stop("Draws should be an array (a list) extracted from a CmdStanMCMC object (CmdStanMCMC objects)")

	if (!is.null(fun))
	{
		if (!isTRUE(all.equal(fun, dunif)) && # Decreasing alphabetical order
			!isTRUE(all.equal(fun, dnorm)) &&
			!isTRUE(all.equal(fun, dlnorm)) &&
			!isTRUE(all.equal(fun, dgamma)) &&
			!isTRUE(all.equal(fun, dexp)) &&
			!isTRUE(all.equal(fun, dbeta))) # isFALSE will not work here, hence !isTRUE
		{
			stop("This function only accepts dnorm, dlnorm, dexp, dgamma, or dbeta as priors")
		}
	}

	# Get list of arguments
	providedArgs = list(...)
	ls_names = names(providedArgs)
	nbArgs = length(providedArgs)

	# Get the argument for density if provided
	n = 512
	n_ind = (ls_names == "n")
	if (any(n_ind))
	{
		n = providedArgs[["n"]]
		print(paste0("Using n = ", n, " for the density plot"))
	}

	# Get the parameter's name if provided
	params = ""
	params_ind = (ls_names == "params")
	if (any(params_ind))
		params = providedArgs[["params"]]

	# Get the index of the x-axis label
	xlab_ind = (ls_names == "xlab")

	# Get the indices of the x-axis limits
	min_x_ind = (ls_names == "min_x")
	max_x_ind = (ls_names == "max_x")

	# Get values indices
	val_ind = stri_detect(str = ls_names, regex = "val[[:digit:]]") 

	# Get the scaling on the x-axis if provided
	scaling_ind = (ls_names == "scaling")
	if (any(scaling_ind)) scaling = providedArgs[["scaling"]] else scaling = 1

	# Get parameters for prior
	if (isTRUE(all.equal(fun, dunif)))
	{
		if ((!all(c("min", "max") %in% names(providedArgs))) && (!all(c("arg1", "arg2") %in% names(providedArgs))))
			stop("You must provide min and max for dunif")

		if (all(c("min", "max") %in% names(providedArgs)))
		{
			arg1 = providedArgs[["min"]]
			arg2 = scaling*providedArgs[["max"]]
		} else {
			print("args 1 and 2 provided; it is assumed that they are min and max, respectively")
			arg1 = providedArgs[["arg1"]]
			arg2 = scaling*providedArgs[["arg2"]]
		}
	}

	if (isTRUE(all.equal(fun, dnorm)))
	{
		if ((!all(c("mean", "sd") %in% names(providedArgs))) && (!all(c("arg1", "arg2") %in% names(providedArgs))))
			stop("You must provide mean and sd for dnorm")

		if (all(c("mean", "sd") %in% names(providedArgs)))
		{
			arg1 = providedArgs[["mean"]]
			arg2 = scaling*providedArgs[["sd"]]
		} else {
			print("args 1 and 2 provided; it is assumed that they are mean and sd, respectively")
			arg1 = providedArgs[["arg1"]]
			arg2 = scaling*providedArgs[["arg2"]]
		}
	}

	if (isTRUE(all.equal(fun, dlnorm)))
	{
		if ((!all(c("mean", "sd") %in% names(providedArgs))) && (!all(c("arg1", "arg2") %in% names(providedArgs))) && 
			(!all(c("meanlog", "sdlog") %in% names(providedArgs))))
			stop("You must provide mean and sd or meanlog and sdlog for dlnorm")

		if (all(c("mean", "sd") %in% names(providedArgs)))
		{
			dlnorm_mean = providedArgs[["mean"]]
			dlnorm_sd = providedArgs[["sd"]]

			arg1 = log(dlnorm_mean^2/sqrt(dlnorm_sd^2 + dlnorm_mean^2)) + log(scaling)
			arg2 = sqrt(log(dlnorm_sd^2/dlnorm_mean^2 + 1))
		} else if (all(c("meanlog", "sdlog") %in% names(providedArgs))) {
			arg1 = providedArgs[["meanlog"]] + log(scaling)
			arg2 = providedArgs[["sdlog"]]
		} else {
			print("args 1 and 2 provided; it is assumed that they are meanlog and sdlog, respectively")
			arg1 = providedArgs[["arg1"]] + log(scaling)
			arg2 = providedArgs[["arg2"]]
		}
	}

	if (isTRUE(all.equal(fun, dgamma)))
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) && (!all(c("shape", "rate") %in% names(providedArgs)))
			&& (!all(c("arg1", "arg2") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape and rate for dgamma")

		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]] # Squared because in this case it is a variance, not a std. dev.

			arg1 = temp1^2/temp2 # shape
			arg2 = temp1/(temp2*scaling) # rate
		}

		if (all(c("shape", "rate") %in% names(providedArgs)))
		{
			arg1 = providedArgs[["shape"]]
			arg2 = providedArgs[["rate"]]/scaling
		}

		if (all(c("arg1", "arg2") %in% names(providedArgs)))
		{
			print("args 1 and 2 provided; it is assumed that they are shape and rate, respectively")
			arg1 = providedArgs[["arg1"]]
			arg2 = providedArgs[["arg2"]]/scaling
		}
	}

	if (isTRUE(all.equal(fun, dexp)))
	{
		if (!("rate" %in% names(providedArgs) && !("arg1" %in% names(providedArgs))))
			stop("You must provide rate for dunif")

		arg1 = ifelse("rate" %in% names(providedArgs), providedArgs[["rate"]], providedArgs[["arg1"]])
		if (arg1 %in% names(providedArgs))
			print("arg 1 provided; it is assumed that it is the rate")
		arg2 = NULL
	}

	if (isTRUE(all.equal(fun, dbeta)))
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) && (!all(c("shape1", "shape2") %in% names(providedArgs)))
			&& (!all(c("arg1", "arg2") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape1 and shape2 for dbeta")

		if (scaling != 1)
			warning("I have not coded the scaling for the beta distribution. Your plot might be out of the window")

		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]]

			arg1 = ((1 - temp1)/temp2 - 1/temp1)*temp1^2 # shape 1
			arg2 = arg1*(1/temp1 - 1) # shape 2
		}

		if (all(c("shape1", "shape2") %in% names(providedArgs)))
		{
			arg1 = providedArgs[["shape1"]]
			arg2 = providedArgs[["shape2"]]
		}

		if (all(c("arg1", "arg2") %in% names(providedArgs)))
		{
			print("args 1 and 2 provided; it is assumed that they are shape1 and shape2, respectively")
			arg1 = providedArgs[["arg1"]]
			arg2 = providedArgs[["arg2"]]
		}

		max_y_prior = optimise(f = fun, interval = c(0, 1), maximum = TRUE, shape1 = arg1, shape2 = arg2)[["objective"]]
	}

	# Get posterior
	if (multi)
	{
		info = summary(draws)
		setDT(info)
		length_params = info[, .N]
		density_from_draws = vector(mode = "list", length = length_params)
		x = vector(mode = "list", length = length_params)
		y = vector(mode = "list", length = length_params)
		names(density_from_draws) = info[, variable]
		names(x) = info[, variable]
		names(y) = info[, variable]
		for (varName in info[, variable])
		{
			density_from_draws[[varName]] = density(scaling*draws[, , varName], n = n)
			x[[varName]] = density_from_draws[[varName]]$x
			y[[varName]] = density_from_draws[[varName]]$y
		}
		min_x = min(sapply(x, min))
		max_x = max(sapply(x, max))
		max_y = max(sapply(y, max))

	} else {
		density_from_draws = density(draws, n = n)
		x = density_from_draws$x
		y = density_from_draws$y
		min_x = min(x)
		max_x = max(x)
		max_y = max(y)
	}

	min_x = ifelse(min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
	max_x = ifelse(max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x

	if (isTRUE(all.equal(fun, dunif)))
		max_y_prior = arg2 # Can be computed analytically  this one...

	if (isTRUE(all.equal(fun, dnorm)))
	{
		max_y_prior = optimise(f = fun, interval = c(min_x, max_x), maximum = TRUE, mean = arg1, sd = arg2)[["objective"]]
		if (expand_bounds)
		{
			check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x,
				mean = arg1, sd = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_min_bound$value > 0.1)
			{
				min_x = ifelse(min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
				check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x,
					mean = arg1, sd = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}

			check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x),
				mean = arg1, sd = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_max_bound$value > 0.1)
			{
				max_x = ifelse(max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x
				check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x),
					mean = arg1, sd = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}
		}
	}

	if (isTRUE(all.equal(fun, dlnorm)))
	{
		max_y_prior = optimise(f = fun, interval = c(min_x, max_x), maximum = TRUE,
			meanlog = arg1, sdlog = arg2)[["objective"]]
		if (expand_bounds)
		{
			check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x, meanlog = arg1,
				sdlog = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_min_bound$value > 0.1)
			{
				min_x = ifelse(min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
				check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x, meanlog = arg1,
					sdlog = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}

			check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x), meanlog = arg1,
				sdlog = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_max_bound$value > 0.1)
			{
				max_x = ifelse(max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x
				check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x), meanlog = arg1,
					sdlog = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}
		}
	}

	if (isTRUE(all.equal(fun, dgamma)))
	{
		max_y_prior = optimise(f = fun, interval = c(min_x, max_x), maximum = TRUE, shape = arg1, rate = arg2)[["objective"]]
		if (expand_bounds)
		{
			check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x,
				shape = arg1, rate = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_min_bound$value > 0.1)
			{
				min_x = ifelse(min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
				check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x,
					shape = arg1, rate = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}

			check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x),
				shape = arg1, rate = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_max_bound$value > 0.1)
			{
				max_x = ifelse(max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x
				check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x),
					shape = arg1, rate = arg2, subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}
		}
	}

	if (isTRUE(all.equal(fun, dexp)))
	{
		max_y_prior = optimise(f = fun, interval = c(min_x, max_x), maximum = TRUE, rate = arg1)[["objective"]]
		if (expand_bounds)
		{
			check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x, rate = arg1,
				subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_min_bound$value > 0.1)
			{
				min_x = ifelse(min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
				check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x, rate = arg1,
					subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}

			check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x), rate = arg1,
				subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_max_bound$value > 0.1)
			{
				max_x = ifelse(max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x
				check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x), rate = arg1,
					subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}
		}
	}

	if (isTRUE(all.equal(fun, dbeta)))
	{
		max_y_prior = optimise(f = fun, interval = c(min_x, max_x), maximum = TRUE,
			shape1 = arg1, shape2 = arg2)[["objective"]]
		if (expand_bounds)
		{
			check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x, rate = arg1,
				subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_min_bound$value > 0.1)
			{
				min_x = ifelse(min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
				check_min_bound = integrate(fun, lower = ifelse(min_x < 0, 10*min_x, -10*min_x), upper = min_x, rate = arg1,
					subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}

			check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x), rate = arg1,
				subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			while (check_max_bound$value > 0.1)
			{
				max_x = ifelse(max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x
				check_max_bound = integrate(fun, lower = max_x, upper = ifelse(max_x < 0, -10*max_x, 10*max_x), rate = arg1,
					subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
			}
		}
	}
	
	if (!is.null(fun))
		max_y = max(max_y, max_y_prior)
	max_y = ifelse(max_y < 0, 0.9*max_y, 1.1*max_y) # To extend 10% from max_y

	# Plot
	if (!is.null(filename))
	{
		filename = paste0(filename, ".pdf")
		pdf(filename)
		print(paste0("Figure saved under the name: ", filename))
	}

	if (any(min_x_ind))
		min_x = providedArgs[["min_x"]]
	if (any(max_x_ind))
		max_x = providedArgs[["max_x"]]

	# Plot posterior
	if (multi)
	{
		colours = MetBrewer::met.brewer("Hokusai3", length_params)
		colours_str = grDevices::colorRampPalette(colours)(length_params)
		colours_str_pol = paste0(colours_str, "66")
		plot(0, type = "n", xlim = c(min_x, max_x), ylim = c(0, max_y), axes = FALSE,
			xlab = ifelse(any(xlab_ind), providedArgs[["xlab"]], ""), ylab = "density", main = "")
		for (i in 1:length_params)
		{
			lines(x = density_from_draws[[i]]$x, y = density_from_draws[[i]]$y, col = colours_str[i], lwd = 2)
			polygon(density_from_draws[[i]], col = colours_str_pol[i])
		}
	} else {
		plot(density_from_draws, xlim = c(min_x, max_x), col = "#295384", lwd = 2,
			xlab = ifelse(any(xlab_ind), providedArgs[["xlab"]], ""), ylab = "density", main = "", axes = FALSE)
		polygon(density_from_draws, col = "#29538466")
	}
	axis(1)
	axis(2, las = 1)

	# Plot prior
	if (!is.null(fun))
	{
		curve(fun(x, arg1, arg2), add = TRUE, lwd = 2, col = "#E9851D")
		if (!is.null(arg2))
			DescTools::Shade(fun(x, arg1, arg2), breaks = c(min_x, max_x), col = "#E9851D66", density = NA)
		if (is.null(arg2))
			DescTools::Shade(fun(x, arg1), breaks = c(min_x, max_x), col = "#E9851D66", density = NA)
	}

	if (any(val_ind))
	{
		for (val in ls_names[val_ind])
			abline(v = providedArgs[[val]], col = "#CD212A", lwd = 4)
	}

	# Add legend
	if (multi && !is.null(ls_nfi))
		if (length(ls_nfi) != length_params)
			warning("Dimension mismatch between ls_nfi and length_params! The legend might not be correctly printed")

	if (!multi && !is.null(ls_nfi))
		if (length(ls_nfi) != 1)
			warning("Too many NFI provided in ls_nfi! The legend might not be correctly printed")

	if (multi)
	{
		legend_text = ifelse(is.null(fun), paste("Posterior", if (!is.null(ls_nfi)) ls_nfi else 1:length_params),
			c("Prior", paste("Posterior", if (!is.null(ls_nfi)) ls_nfi else 1:length_params)))
		legend_colours = ifelse(is.null(fun), colours_str, c("#E9851D", colours_str))
	} else {
		legend_text = if(is.null(fun)) "Posterior" else c("Prior", "Posterior")
		legend_colours = if(is.null(fun)) "#295384" else c("#E9851D", "#295384")
	}
	legend(x = "topright", legend = legend_text, fill = legend_colours, box.lwd = 0)

	if (!is.null(filename))
		dev.off()

	if (is.null(fun))
		return(list(arg1 = NA, arg2 = NA, min_x = min_x, max_x = max_x, max_y = max_y, max_y_prior = NA, filename = filename))

	return(list(arg1 = arg1, arg2 = arg2, min_x = min_x, max_x = max_x, max_y = max_y, max_y_prior = max_y_prior, filename = filename))
}

# Function to approximate the sum of gamma distributions
gamma_sum_approx = function(x, order, shape_vec, scale_vec, tol = sqrt(.Machine$double.eps), rate = FALSE)
{
	delta_func = function(lim, shape_vec, scale_vec) # lim is up to where compute delta!
	{
		is.wholenumber = function(x, tol = .Machine$double.eps^0.5)
			abs(x - round(x)) < tol

		if (lim < 0 || !is.wholenumber(lim))
			stop(paste0("lim = ", lim, " must be a positive integer!!"))

		delta = numeric(lim + 1) #! +1 because R starts its index at 1 insted of 0
		delta[1] = 1 # delta 0

		gamma_vec = numeric(lim)

		for (k in seq_len(lim))
		{
			gamma_vec[k] = 1/k * sum(shape_vec * (1 - scale_vec[1]/scale_vec)^k)
			delta[k + 1] = 1/k * sum(1:k * gamma_vec[1:k] * delta[k:1]) #! Watch out, R starts its index at 1 instead of 0, hence k + 1
		}

		return(delta)
	}

	if (rate)
		scale_vec = 1/scale_vec # The formula was developped for shape and scale, with scale = 1/rate

	if (isTRUE(all.equal(min(scale_vec), max(scale_vec), tolerance = tol)))
	{
		ans_url = "https://math.stackexchange.com/questions/250059/sum-of-independent-gamma-distributions-is-a-gamma-distribution"
		warning(paste("When all the scale parameters (or rate) are equal, there exists an analytical formula:", ans_url))
	}

	beta1 = scale_vec[1]
	if (min(scale_vec) != beta1)
	{
		warning("Without loss of generality, beta_1 was assumed to be the min in https://doi.org/10.1007/bf02481123. Permutating 1st with min")
		if (rate)
			warning("You are using rate, but the formula was developped with scale = 1/rate. Maybe this is why beta_1 is not the minimum")

		# Permutation within scale_vec
		min_loc = which.min(scale_vec)
		scale_vec[1] = scale_vec[min_loc]
		scale_vec[min_loc] = beta1
		beta1 = scale_vec[1]

		# Permutating shape_vec accordingly
		temporary = shape_vec[1]
		shape_vec[1] = shape_vec[min_loc]
		shape_vec[min_loc] = temporary
	}

	rho = sum(shape_vec)
	C = prod((beta1/scale_vec)^shape_vec)

	delta_vec = delta_func(lim = order, shape_vec = shape_vec, scale_vec = scale_vec)

	results = numeric(length(x))

	for (i in seq_along(x))
		results[i] = C * exp(-x[i]/beta1) * sum(delta_vec * x[i]^(rho + 0:order - 1)/(gamma(rho + 0:order) * beta1^(rho + 0:order)))

	return(results)
}

## Function to print pretty summary CmdStanR model for quarto
pretty_summary = function(fit, params)
{
	draws = fit$draws(params)
	var_names = dimnames(draws)$variable

	if (length(var_names) != length(params))
		warning("Watch out, some params were vectors in stan model")

	# See Bayesian Data Analysis, 3rd ed., p. 283-285 for explanation
	# B is the between sequence variance
	# W is the within sequence variance
	n_iter = fit$metadata()$iter_sampling
	n = floor(n_iter/2)
	m = 2*fit$num_chains()

	info_dt = data.table(params_name = var_names, mean_params = apply(X = draws, MARGIN = 3, FUN = mean),
		sd_params = apply(X = draws, MARGIN = 3, FUN = sd), r_hat_params = 0.0)
	setkey(info_dt, params_name)
	for (current_var in var_names)
	{
		split_chains = vector(mode = "list", length = m)
		for (cc in seq_len(fit$num_chains()))
		{
			split_chains[[2*cc - 1]] = draws[1:n, cc, current_var]
			split_chains[[2*cc]] = draws[(n + 1):n_iter, cc, current_var]
		}

		mean_within_seq = sapply(X = split_chains, FUN = mean)
		mean_btw_seq = mean(mean_within_seq)
		B = n/(m - 1) * sum((mean_within_seq - mean_btw_seq))

		s_j = numeric(length = m)
		for (j in seq_len(m))
			s_j[j] = 1/(n - 1) * sum((split_chains[[j]] - mean_within_seq[j])^2)

		W = mean(s_j)

		# Estimate marginal posterior variance with weighted sum
		var_plus = (n - 1)/n*W + 1/n*B
		info_dt[(current_var), r_hat_params := sqrt(var_plus/W)]
	}

	return(info_dt)
}
