## Get fixed values parameters (will not work for draws with third dimension > 1)
getParams = function(model_cmdstan, params_names, type = "mean", ...)
{
	if (!(type %in% c("all", "chain", "chain-iter", "mean", "median", "quantile")))
		stop("Unknown type. Please choose all, iter-chain, mean, median, or quantile")

	args = list(...)

	if (type %in% c("chain", "chain-iter", "mean", "median"))
	{
		vals = numeric(length(params_names))
		names(vals) = params_names

		if (type == "chain")
		{
			if (!("chain" %in% names(args)))
				stop("You must provide chain number when using the option chain")

			chain = args[["chain"]]
			if (chain > model_cmdstan$num_chains())
				stop("chain cannot be larger than num_chains")

			for (i in seq_along(params_names))
			{
				draws = model_cmdstan$draws(params_names[i])
				if (dim(draws)[3] != 1)
					stop("This function is not yet designed to handle a multi-params")
				vals[i] = mean(draws[, chain, 1])
			}
			return(vals)
		} else if (type == "chain-iter") {
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
			vals = data.table::data.table(parameter = params_names, min = 0, q05 = 0, med = 0, avg = 0, q95 = 0, max = 0)
		}
		if ("probs" %in% names(args))
		{
			warning("I have not yet coded the general output, so far I use probs = c(0, 0.05, 0.5, 0.95, 1)")
			probs = c(0, 0.05, 0.5, 0.95, 1)
			vals = data.table::data.table(parameter = params_names, min = 0, q05 = 0, med = 0, avg = 0, q95 = 0, max = 0)
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
	colours = MetBrewer::met.brewer("Egypt", n_chains)
	colours_str = grDevices::colorRampPalette(colours)(n_chains)

	min_val = min(draws)
	max_val = max(draws)

	providedArgs = list(...)
	nbArgs = length(providedArgs)

	ls_names = names(providedArgs)

	val_ind = stringi::stri_detect(str = ls_names, regex = "val[[:digit:]]")
	xlab_ind = (ls_names == "xlab")
	ylab_ind = (ls_names == "ylab")
	main_ind = (ls_names == "main")
	label_ind = stringi::stri_detect(str = ls_names, regex = "label")
	iter_ind = stringi::stri_detect(str = ls_names, regex = "iter[[:digit:]]")
	chain_ind = stringi::stri_detect(str = ls_names, regex = "chain_id")

	scaling_ind = (ls_names == "scaling")
	if (any(scaling_ind)) scaling = providedArgs[["scaling"]] else scaling = 1

	if (any(label_ind))
		par(mar = c(5, 4, 4, 4))

	if (any(chain_ind)) chain_id = providedArgs[["chain_id"]] else chain_id = seq_len(n_chains)

	if (!all(seq_len(n_chains) %in% chain_id))
	{
		if (!all(chain_id %in% seq_len(n_chains)))
		{
			warning("Some provided chain_id are beyond n_chains. Trimmed")
			chain_id = chain_id[chain_id %in% seq_len(n_chains)]
		}

		draws = draws[, chain_id, ]
		n_chains = dim(draws)[2]
	}

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

	for (chain in seq_len(n_chains))
	{
		if (all(class(draws) %in% c("draws_array", "draws", "array")))
			lines(seq_len(n_iter), scaling*draws[, chain, ], type = "l", col = colours_str[chain])
		if (is.array(draws) && !all(class(draws) %in% c("draws_array", "draws", "array")))
			lines(seq_len(n_iter), scaling*draws[, chain], type = "l", col = colours_str[chain])
	}

	if (any(val_ind))
	{
		for (val in ls_names[val_ind])
			abline(h = scaling*providedArgs[[val]], col = "#CD212A", lwd = 4)

		if (any(label_ind))
		{
			num_vals = stringi::stri_sub(str = ls_names[val_ind], from = stringi::stri_locate(str = ls_names[val_ind], regex = "val")[, "end"] + 1)
			for (label in ls_names[label_ind])
			{
				num_label = stringi::stri_sub(str = label, from = stringi::stri_locate(str = label, regex = "label")[, "end"] + 1)
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

	legend("topright", legend = paste("Chain", seq_len(n_chains)), fill = colours)

	if (!is.null(filename))
		dev.off()
}

## Function to plot divergences (when any during Bayesian run)
plot_divergences = function(fit, div, simplif = FALSE)
{
	params = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"), type = "all")

	div_c = posterior::subset_draws(params[, , "c"], draw = which(div == 1))
	div_j = posterior::subset_draws(params[, , "j"], draw = which(div == 1))
	div_m = posterior::subset_draws(params[, , "m"], draw = which(div == 1))
	div_k = posterior::subset_draws(params[, , "k"], draw = which(div == 1))
	div_s = posterior::subset_draws(params[, , "s"], draw = which(div == 1))

	# Plot divergence for c and j
	plot(params[, , "c"], params[, , "j"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "j", main = "c and j")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# Plot divergence for c and k
	plot(params[, , "c"], params[, , "k"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "k", main = "c and k")
	points(div_c, div_k, pch = 18, col = "#0F7BA2")

	# Plot divergence for c and m
	plot(params[, , "c"], params[, , "m"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "m", main = "c and m")
	points(div_c, div_m, pch = 18, col = "#0F7BA2")

	# Plot divergence for j and k
	plot(params[, , "j"], params[, , "k"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "j", ylab = "k", main = "j and k")
	points(div_j, div_k, pch = 18, col = "#0F7BA2")

	# Plot divergence for j and m
	plot(params[, , "j"], params[, , "m"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "j", ylab = "m", main = "j and m")
	points(div_j, div_m, pch = 18, col = "#0F7BA2")

	# Plot divergence for k and m
	plot(params[, , "k"], params[, , "m"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "k", ylab = "m", main = "k and m")
	points(div_k, div_m, pch = 18, col = "#0F7BA2")
	
	# Bonus plots of susceptible correlated parameters: k and s
	plot(params[, , "k"], params[, , "s"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "k", ylab = "s", main = "k and s")
	points(div_k, div_s, pch = 18, col = "#0F7BA2")
}

## Function to plot a species fit (simplif = TRUE for submodel)
plot_sp = function(fit, sp, forest = tree_dt[.(sp)], pred = TRUE, simplif = FALSE,
	n_bins = 4, pal = "Hiroshige", selected_variable = NULL, print_plot = TRUE)
{
	# Get parameters
	n_sampling = fit$metadata()$iter_sampling
	n_chains = fit$num_chains()

	if (n_chains != nrow(fit$metadata()$time))
	{
		warning("Some chains are missing probably because they could not start! Using the number of successful chains")
		n_chains = nrow(fit$metadata()$time)
	}

	if (n_chains != 4)
		warning("This function will bug if plotting chains (rhat > 1.05) and n_chains != 4")

	# Add colour to the measured points according to selected variable
	if (!is.null(selected_variable))
	{
		if (!(selected_variable %in% colnames(forest)))
			stop(paste0("'", selected_variable, "' is not a column of forest"))
		
		colours = MetBrewer::met.brewer(pal, n_bins)[1:n_bins]
		forest[, colour_ind := as.numeric(cut(.SD[[selected_variable]], breaks = n_bins))] # Indices per category
		forest[, colour := colours[colour_ind]] # Map selected variable to colour using indices
	} else {
		forest[, colour := "#3355AA33"] # Default colour
	}

	# Define fitted functions (depends if simplified model 2 stages vs 3 stages)
	if (!simplif)
	{
		paramsVec = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"), type = "mean")
	} else {
		paramsVec_simplif = getParams(model_cmdstan = fit,
			params_names = c("alpha", "beta_", "gamma", "delta"), type = "mean")
		paramsVec = c(
			c = unname(paramsVec_simplif["alpha"]),
			j = 1,
			k = unname(paramsVec_simplif["beta_"]),
			m = unname(exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"]),
			n = unname(paramsVec_simplif["delta"] + paramsVec_simplif["alpha"]),
			s = unname(paramsVec_simplif["beta_"])
		)
	}

	r_func = function(x, pars)
	{
		return ((pars["m"] - pars["c"]) * exp(pars["j"] - pars["k"]*x) * (pars["k"]*x/pars["j"])^pars["j"] +
			pars["c"] - (pars["c"] - pars["n"])*exp(-pars["s"]*x));
	}
	
	psi = function(x, pars)
	{
		r = r_func(x, pars)

		return ( 1/pars["c"] * x^( 1 - (log(r) - log(pars["c"])) / log(x) ) )
	}

	# Check rhats
	pars_names = names(paramsVec)
	if (simplif)
		pars_names = c("alpha", "beta_", "gamma", "delta")
	rhats = bayesplot::rhat(object = fit, pars = pars_names)

	# Divergence check up
	params = getParams(model_cmdstan = fit, params_names = pars_names, type = "all")
	div = posterior::subset_draws(fit$sampler_diagnostics(), variable = "divergent__")

	if (!print_plot & !pred)
		return (list(any_div = any(div != 0), loc = div, rhats = rhats))

	# Plot figure r vs bole and fitted curve
	if (print_plot)
	{
		plot(forest[, bole_volume_m3], forest[, r], pch = 19, col = forest[, colour], axes = FALSE,
			xlab = "Bole volume", ylab = "Ratio bole/tot volumes", lwd = 0, cex = 0.75,
			main = ifelse(is.null(selected_variable), sp, paste(sp, "(", selected_variable, ")")))
		if (!is.null(selected_variable))
			legend("topright", fill = colours, title = paste(selected_variable, "classes"),
				legend = levels(forest[, cut(.SD[[selected_variable]], breaks = n_bins, dig.lab = 1)]))
		curve(r_func(x, paramsVec), add = TRUE, lwd = 3)
		axis(1)
		axis(2, las = 1)
	}

	# Plot predictions
	if (pred)
	{
		sim = apply(X = fit$draws("v_gen_mean"), MARGIN = 3, FUN = mean)

		if (print_plot)
		{
			plot(sim, forest[, total_volume_m3], pch = 19, col = forest[, colour], axes = FALSE,
				xlab = "Predicted total volume", ylab = "Observed total volume",
				main = ifelse(is.null(selected_variable), sp, paste(sp, "(", selected_variable, ")")))
			if (!is.null(selected_variable))
				legend("topright", fill = colours, title = paste(selected_variable, "classes"),
					legend = levels(forest[, cut(.SD[[selected_variable]], breaks = n_bins, dig.lab = 1)]))
			abline(a = 0, b = 1, lty = "dashed", lwd = 1.5, col = "#2E2E2E")
			axis(1)
			axis(2, las = 1)

			plot(forest[, bole_volume_m3], forest[, total_volume_m3], pch = 19, col = forest[, colour],
				axes = FALSE, xlab = "Observed bole volume", ylab = "Observed total volume",
				main = ifelse(is.null(selected_variable), sp, paste(sp, "(", selected_variable, ")")))
			if (!is.null(selected_variable))
				legend("topright", fill = colours, title = paste(selected_variable, "classes"),
					legend = levels(forest[, cut(.SD[[selected_variable]], breaks = n_bins, dig.lab = 1)]))
			abline(a = 0, b = 1, lty = "dashed", lwd = 0.85, col = "#2E2E2E")
			curve(psi(x, paramsVec), add = TRUE, lwd = 3)
			axis(1)
			axis(2, las = 1)
		}
	}

	if (!print_plot & pred)
	{
		return (list(any_div = any(div != 0), loc = div, rhats = rhats, pred_tot = sim,
			obs_bole = forest[, bole_volume_m3], obs_tot = forest[, total_volume_m3], params = paramsVec))
	}

	# Plot chains if any rhat problem
	if (any(rhats > 1.05))
	{
		if (!simplif)
		{
			lazyTrace(fit$draws("c"), main = paste0(sp, ", c"))
			lazyTrace(fit$draws("j"), main = paste0(sp, ", j"))
			lazyTrace(fit$draws("k"), main = paste0(sp, ", k"))
			lazyTrace(fit$draws("m"), main = paste0(sp, ", m"))
			lazyTrace(fit$draws("n"), main = paste0(sp, ", n"))
			lazyTrace(fit$draws("s"), main = paste0(sp, ", s"))

			paramsVec1 = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"),
				type = "chain", chain = 1)
			paramsVec2 = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"),
				type = "chain", chain = 2)
			paramsVec3 = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"),
				type = "chain", chain = 3)
			paramsVec4 = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"),
				type = "chain", chain = 4)
		} else {
			lazyTrace(fit$draws("alpha"), main = "alpha")
			lazyTrace(fit$draws("beta_"), main = "beta_")
			lazyTrace(fit$draws("gamma"), main = "gamma")
			lazyTrace(fit$draws("delta"), main = "delta")

			paramsVec_simplif = getParams(model_cmdstan = fit, params_names = c("alpha", "beta_", "gamma", "delta"),
				type = "chain", chain = 1)
			paramsVec1 = c(
				c = paramsVec_simplif["alpha"],
				j = 1,
				k = paramsVec_simplif["beta_"],
				m = exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"],
				n = paramsVec_simplif["delta"] + paramsVec_simplif["alpha"],
				s = paramsVec_simplif["beta_"]
			)
			paramsVec_simplif = getParams(model_cmdstan = fit, params_names = c("alpha", "beta_", "gamma", "delta"),
				type = "chain", chain = 2)
			paramsVec2 = c(
				c = paramsVec_simplif["alpha"],
				j = 1,
				k = paramsVec_simplif["beta_"],
				m = exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"],
				n = paramsVec_simplif["delta"] + paramsVec_simplif["alpha"],
				s = paramsVec_simplif["beta_"]
			)
			paramsVec_simplif = getParams(model_cmdstan = fit, params_names = c("alpha", "beta_", "gamma", "delta"),
				type = "chain", chain = 3)
			paramsVec3 = c(
				c = paramsVec_simplif["alpha"],
				j = 1,
				k = paramsVec_simplif["beta_"],
				m = exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"],
				n = paramsVec_simplif["delta"] + paramsVec_simplif["alpha"],
				s = paramsVec_simplif["beta_"]
			)
			paramsVec_simplif = getParams(model_cmdstan = fit, params_names = c("alpha", "beta_", "gamma", "delta"),
				type = "chain", chain = 4)
			paramsVec4 = c(
				c = paramsVec_simplif["alpha"],
				j = 1,
				k = paramsVec_simplif["beta_"],
				m = exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"],
				n = paramsVec_simplif["delta"] + paramsVec_simplif["alpha"],
				s = paramsVec_simplif["beta_"]
			)
		}
		plot(forest[, bole_volume_m3], forest[, r], pch = 19, axes = FALSE, lwd = 0,
			xlab = "Bole volume", ylab = "Ratio bole/tot volumes", main = sp)
		curve(r_func(x, paramsVec), add = TRUE, lwd = 3)
		curve(r_func(x, paramsVec1), lwd = 3, add = TRUE, col = MetBrewer::met.brewer("Egypt")[1])
		curve(r_func(x, paramsVec2), lwd = 3, add = TRUE, col = MetBrewer::met.brewer("Egypt")[2])
		curve(r_func(x, paramsVec3), lwd = 3, add = TRUE, col = MetBrewer::met.brewer("Egypt")[3])
		curve(r_func(x, paramsVec4), lwd = 3, add = TRUE, col = MetBrewer::met.brewer("Egypt")[4])
		axis(1)
		axis(2, las = 1)
	}

	return (list(any_div = any(div != 0), loc = div, rhats = rhats))
}

## Same function as above but for groups
plot_gr = function(fit, gr, forest = tree_dt[.(gr)], pred = TRUE, simplif = FALSE,
	n_bins = 4, pal = "Hiroshige", selected_variable = NULL, print_plot = TRUE)
{
	# Get parameters
	n_sampling = fit$metadata()$iter_sampling
	n_chains = fit$num_chains()

	if (n_chains != nrow(fit$metadata()$time))
	{
		warning("Some chains are missing probably because they could not start! Using the number of successful chains")
		n_chains = nrow(fit$metadata()$time)
	}

	# Add colour to the measured points according to selected variable
	if (!is.null(selected_variable))
	{
		if (!(selected_variable %in% colnames(forest)))
			stop(paste0("'", selected_variable, "' is not a column of forest"))
		
		colours = MetBrewer::met.brewer(pal, n_bins)[1:n_bins]
		forest[, colour_ind := as.numeric(cut(.SD[[selected_variable]], breaks = n_bins))] # Indices per category
		forest[, colour := colours[colour_ind]] # Map selected variable to colour using indices
	} else {
		forest[, colour := "#3355AA33"] # Default colour
	}

	# Define fitted functions (depends if simplified model 2 stages vs 3 stages)
	if (!simplif)
	{
		paramsVec = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"), type = "mean")
	} else {
		paramsVec_simplif = getParams(model_cmdstan = fit,
			params_names = c("alpha", "beta_", "gamma", "delta"), type = "mean")
		paramsVec = c(
			c = unname(paramsVec_simplif["alpha"]),
			j = 1,
			k = unname(paramsVec_simplif["beta_"]),
			m = unname(exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"]),
			n = unname(paramsVec_simplif["delta"] + paramsVec_simplif["alpha"]),
			s = unname(paramsVec_simplif["beta_"])
		)
	}

	r_func = function(x, pars)
	{
		return ((pars["m"] - pars["c"]) * exp(pars["j"] - pars["k"]*x) * (pars["k"]*x/pars["j"])^pars["j"] +
			pars["c"] - (pars["c"] - pars["n"])*exp(-pars["s"]*x));
	}
	
	psi = function(x, pars)
	{
		r = r_func(x, pars)

		return ( 1/pars["c"] * x^( 1 - (log(r) - log(pars["c"])) / log(x) ) )
	}

	# Check rhats
	pars_names = names(paramsVec)
	if (simplif)
		pars_names = c("alpha", "beta_", "gamma", "delta")
	rhats = bayesplot::rhat(object = fit, pars = pars_names)

	# Divergence check up
	params = getParams(model_cmdstan = fit, params_names = pars_names, type = "all")
	div = posterior::subset_draws(fit$sampler_diagnostics(), variable = "divergent__")

	if (!print_plot & !pred)
		return (list(any_div = any(div != 0), loc = div, rhats = rhats))

	# Plot figure r vs bole and fitted curve
	if (print_plot)
	{
		plot(forest[, bole_volume_m3], forest[, r], pch = 19, col = forest[, colour], axes = FALSE,
			xlab = "Bole volume", ylab = "Ratio bole/tot volumes", lwd = 0, cex = 0.75,
			main = ifelse(is.null(selected_variable), gr, paste(gr, "(", selected_variable, ")")))
		if (!is.null(selected_variable))
			legend("topright", fill = colours, title = paste(selected_variable, "classes"),
				legend = levels(forest[, cut(.SD[[selected_variable]], breaks = n_bins, dig.lab = 1)]))
		curve(r_func(x, paramsVec), add = TRUE, lwd = 3)
		axis(1)
		axis(2, las = 1)
	}

	# Plot predictions
	if (pred)
	{
		sim = apply(X = fit$draws("v_gen_mean"), MARGIN = 3, FUN = mean)

		if (print_plot)
		{
			plot(sim, forest[, total_volume_m3], pch = 19, col = forest[, colour], axes = FALSE,
				xlab = "Predicted total volume", ylab = "Observed total volume",
				main = ifelse(is.null(selected_variable), gr, paste(gr, "(", selected_variable, ")")))
			if (!is.null(selected_variable))
				legend("topright", fill = colours, title = paste(selected_variable, "classes"),
					legend = levels(forest[, cut(.SD[[selected_variable]], breaks = n_bins, dig.lab = 1)]))
			abline(a = 0, b = 1, lty = "dashed", lwd = 1.5, col = "#2E2E2E")
			axis(1)
			axis(2, las = 1)

			plot(forest[, bole_volume_m3], forest[, total_volume_m3], pch = 19, col = forest[, colour],
				axes = FALSE, xlab = "Observed bole volume", ylab = "Observed total volume",
				main = ifelse(is.null(selected_variable), gr, paste(gr, "(", selected_variable, ")")))
			if (!is.null(selected_variable))
				legend("topright", fill = colours, title = paste(selected_variable, "classes"),
					legend = levels(forest[, cut(.SD[[selected_variable]], breaks = n_bins, dig.lab = 1)]))
			abline(a = 0, b = 1, lty = "dashed", lwd = 0.85, col = "#2E2E2E")
			curve(psi(x, paramsVec), add = TRUE, lwd = 3)
			axis(1)
			axis(2, las = 1)
		}
	}

	if (!print_plot & pred)
	{
		return (list(any_div = any(div != 0), loc = div, rhats = rhats, pred_tot = sim,
			obs_bole = forest[, bole_volume_m3], obs_tot = forest[, total_volume_m3], params = paramsVec))
	}

	print("Hello 3")

	# Plot chains if any rhat problem
	if (any(rhats > 1.05))
	{
		if (!simplif)
		{
			lazyTrace(fit$draws("c"), main = paste0(gr, ", c"))
			lazyTrace(fit$draws("j"), main = paste0(gr, ", j"))
			lazyTrace(fit$draws("k"), main = paste0(gr, ", k"))
			lazyTrace(fit$draws("m"), main = paste0(gr, ", m"))
			lazyTrace(fit$draws("n"), main = paste0(gr, ", n"))
			lazyTrace(fit$draws("s"), main = paste0(gr, ", s"))

			paramsVec1 = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"),
				type = "chain", chain = 1)
			paramsVec2 = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"),
				type = "chain", chain = 2)
			paramsVec3 = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"),
				type = "chain", chain = 3)
			paramsVec4 = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s"),
				type = "chain", chain = 4)
		} else {
			lazyTrace(fit$draws("alpha"), main = "alpha")
			lazyTrace(fit$draws("beta_"), main = "beta_")
			lazyTrace(fit$draws("gamma"), main = "gamma")
			lazyTrace(fit$draws("delta"), main = "delta")

			paramsVec_simplif = getParams(model_cmdstan = fit, params_names = c("alpha", "beta_", "gamma", "delta"),
				type = "chain", chain = 1)
			paramsVec1 = c(
				c = paramsVec_simplif["alpha"],
				j = 1,
				k = paramsVec_simplif["beta_"],
				m = exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"],
				n = paramsVec_simplif["delta"] + paramsVec_simplif["alpha"],
				s = paramsVec_simplif["beta_"]
			)
			paramsVec_simplif = getParams(model_cmdstan = fit, params_names = c("alpha", "beta_", "gamma", "delta"),
				type = "chain", chain = 2)
			paramsVec2 = c(
				c = paramsVec_simplif["alpha"],
				j = 1,
				k = paramsVec_simplif["beta_"],
				m = exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"],
				n = paramsVec_simplif["delta"] + paramsVec_simplif["alpha"],
				s = paramsVec_simplif["beta_"]
			)
			paramsVec_simplif = getParams(model_cmdstan = fit, params_names = c("alpha", "beta_", "gamma", "delta"),
				type = "chain", chain = 3)
			paramsVec3 = c(
				c = paramsVec_simplif["alpha"],
				j = 1,
				k = paramsVec_simplif["beta_"],
				m = exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"],
				n = paramsVec_simplif["delta"] + paramsVec_simplif["alpha"],
				s = paramsVec_simplif["beta_"]
			)
			paramsVec_simplif = getParams(model_cmdstan = fit, params_names = c("alpha", "beta_", "gamma", "delta"),
				type = "chain", chain = 4)
			paramsVec4 = c(
				c = paramsVec_simplif["alpha"],
				j = 1,
				k = paramsVec_simplif["beta_"],
				m = exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] + paramsVec_simplif["alpha"],
				n = paramsVec_simplif["delta"] + paramsVec_simplif["alpha"],
				s = paramsVec_simplif["beta_"]
			)
		}
		plot(forest[, bole_volume_m3], forest[, r], pch = 19, axes = FALSE, lwd = 0,
			xlab = "Bole volume", ylab = "Ratio bole/tot volumes", main = gr)
		curve(r_func(x, paramsVec), add = TRUE, lwd = 3)
		curve(r_func(x, paramsVec1), lwd = 3, add = TRUE, col = met.brewer("Egypt")[1])
		curve(r_func(x, paramsVec2), lwd = 3, add = TRUE, col = met.brewer("Egypt")[2])
		curve(r_func(x, paramsVec3), lwd = 3, add = TRUE, col = met.brewer("Egypt")[3])
		curve(r_func(x, paramsVec4), lwd = 3, add = TRUE, col = met.brewer("Egypt")[4])
		axis(1)
		axis(2, las = 1)
	}

	return (list(any_div = any(div != 0), loc = div, rhats = rhats))
}

## Function to rebuild comparison between models
rebuild_comp = function(save_ls)
{
	ls_species = names(save_ls)
	
	comp = data.table(species = ls_species, best = "", elpd_diff = -Inf, se_diff = -Inf,
		warning = FALSE, key = "species")

	weights_dt = data.table(species = ls_species, best = "", W_full = -Inf, W_sub = -Inf, key = "species")

	rhat_dt = data.table(species = ls_species, Rhat_full = -Inf, Rhat_sub = -Inf, key = "species")

	R2D2 = data.table(species = ls_species,
		R2_full = -Inf, R2_sub = -Inf, R2_Vtot_full = -Inf, R2_Vtot_sub = -Inf,
		R2_loo_full = -Inf, R2_loo_sub = -Inf, R2_loo_Vtot_full = -Inf, R2_loo_Vtot_sub = -Inf,
		key = "species")

	for (sp in ls_species)
	{
		weights_dt[sp, c("W_full", "W_sub") := as.list(save_ls[[sp]]$weights)]
	  
	  comp[.(sp), c("best", "elpd_diff", "se_diff", "warning") :=
			.(save_ls[[sp]]$comploo[1, "model"], save_ls[[sp]]$comploo[2, "elpd_diff"],
			save_ls[[sp]]$comploo[2, "se_diff"], save_ls[[sp]]$warning)]

		R2D2[.(sp), c("R2_full", "R2_sub") := .(
			median(save_ls[[sp]][["rsq_distrib"]][["full"]]),
			median(save_ls[[sp]][["rsq_distrib"]][["sub"]])
		)]

		R2D2[.(sp), c("R2_Vtot_full", "R2_Vtot_sub") := .(
			median(save_ls[[sp]][["rsq_vtot_distrib"]][["full"]][["rsq_vtot"]]),
			median(save_ls[[sp]][["rsq_vtot_distrib"]][["sub"]][["rsq_vtot"]])
		)]
		
		R2D2[.(sp), c("R2_loo_full", "R2_loo_sub") := .(
			median(save_ls[[sp]][["rsq_loo_distrib_r"]][["full"]]),
			median(save_ls[[sp]][["rsq_loo_distrib_r"]][["sub"]])
		)]

		R2D2[.(sp), c("R2_loo_Vtot_full", "R2_loo_Vtot_sub") := .(
			median(save_ls[[sp]][["rsq_loo_distrib_v"]][["full"]]),
			median(save_ls[[sp]][["rsq_loo_distrib_v"]][["sub"]])
		)]

		rhat_dt[sp, c("Rhat_full", "Rhat_sub") := .(save_ls[[sp]][["rhat_full"]], save_ls[[sp]][["rhat_sub"]])]
	}

	return(list(comp = comp, weights_dt = weights_dt, R2D2 = R2D2, rhat = rhat_dt))
}

## Function to compare model new version with submodel, based on PSIS-LOO
comparison_full_sub = function(sp, tree_dt, path_models, path_output,
	woodstock_seed = 1969 - 08 - 18, printPlot = FALSE)
{
	#### Generate quantitities
	## Compile generator models
	genQ_full = cmdstanr::cmdstan_model(paste0(path_models, "fullmodel-genQ.stan"))
	genQ_sub = cmdstanr::cmdstan_model(paste0(path_models, "submodel-genQ.stan"))

	## Prepare data
	stanData = list(
		N = tree_dt[.(sp)][, .N],
		N_new = tree_dt[.(sp)][, .N],
		bole_volume_m3 = tree_dt[.(sp)][, bole_volume_m3],
		bole_volume_m3_new = tree_dt[.(sp)][, bole_volume_m3],
		total_volume_m3 = tree_dt[.(sp)][, total_volume_m3],
		total_volume_m3_new = tree_dt[.(sp)][, total_volume_m3]
	)

	## Load fitted model
	full_file = paste0(path_output, stringi::stri_replace(str = sp, regex = " ", replacement = "-"), "_fullmodel")
	sub_file = paste0(path_output, stringi::stri_replace(str = sp, regex = " ", replacement = "-"), "_submodel")

	if (sp == "Quercus sp.")
		full_file = stringi::stri_replace(str = full_file, regex = "._", replacement = "_")

	full = readRDS(paste0(full_file, ".rds"))
	sub = readRDS(paste0(sub_file, ".rds"))

	n_chains = full$num_chains()
	n_draws = n_chains*full$metadata()$iter_sampling

	## Compute max rhat for full model and submodel
	rhat_full = posterior::summarise_draws(full$draws(), posterior::rhat)
	rhat_full_m = max(rhat_full[, "posterior::rhat"])
	rhat_sub = posterior::summarise_draws(sub$draws(), posterior::rhat)
	rhat_sub_m = max(rhat_sub[, "posterior::rhat"])

	## Predict on new data
	sim_full = genQ_full$generate_quantities(fitted_params = full, data = stanData,
		seed = woodstock_seed, parallel_chains = min(n_chains, 4))
	
	sim_sub = genQ_sub$generate_quantities(fitted_params = sub, data = stanData,
		seed = woodstock_seed, parallel_chains = min(n_chains, 4))

	## Save simulations
	v_gen_mean_full = apply(X = sim_full$draws("v_gen_mean"), MARGIN = 3, FUN = mean)
	v_gen_mean_sub = apply(X = sim_sub$draws("v_gen_mean"), MARGIN = 3, FUN = mean)

	if (printPlot)
	{	
		plot(x = v_gen_mean_sub, y = v_gen_mean_full, pch = 19, xlab = "Sub", ylab = "New", axes = FALSE)
		abline(a = 0, b = 1, col = "#CD212A")
		axis(1)
		axis(2, las = 1)
	}

	## Compute PSIS-LOO...
	# ... for full model
	r_eff_full = loo::relative_eff(exp(sim_full$draws("log_lik")), cores = 8)
	loo_full = loo::loo(x = sim_full$draws("log_lik"), r_eff = r_eff_full, cores = 8, save_psis = TRUE)
	
	# ... for submodel
	r_eff_sub = loo::relative_eff(exp(sim_sub$draws("log_lik")), cores = 8)
	loo_sub = loo::loo(x = sim_sub$draws("log_lik"), r_eff = r_eff_sub, cores = 8, save_psis = TRUE)

	warning_loo_full = FALSE
	if (any(loo_full$diagnostics$pareto_k >= 0.7))
		warning_loo_full = TRUE
	
	warning_loo_sub = FALSE
	if (any(loo_sub$diagnostics$pareto_k >= 0.7))
		warning_loo_sub = TRUE

	warn = warning_loo_full | warning_loo_sub

	## Compare both models based on psisloo and weights
	comp = loo::loo_compare(list(full = loo_full, sub = loo_sub))
	best = comp[1, "model"]

	weights = head(loo::loo_model_weights(x = list(full = loo_full, sub = loo_sub),
		method = "stacking")) # I use head to return a named vector

	## R squared based on Gelman 2019...
	# ... for the ratio
	var_fit_full = apply(X = posterior::as_draws_matrix(sim_full$draws("r_gen_mean")),
		MARGIN = 1, FUN = var) # The var contains the correction 1/(n - 1) already!
	var_res_full = apply(X = posterior::as_draws_matrix(sim_full$draws("sigma_var")),
		MARGIN = 1, FUN = mean)
	
	var_fit_sub = apply(X = posterior::as_draws_matrix(sim_sub$draws("r_gen_mean")),
		MARGIN = 1, FUN = var)
	var_res_sub = apply(X = posterior::as_draws_matrix(sim_sub$draws("sigma_var")),
		MARGIN = 1, FUN = mean)

	rsq_full = var_fit_full / (var_fit_full + var_res_full)
	rsq_sub = var_fit_sub / (var_fit_sub + var_res_sub)

	# ... for the volume. This requires rescaling (see notebook 4, entry 29.07.2026)
	rsq_vtot_fct = function(sim, Vbole = stanData[["bole_volume_m3_new"]])
	{
		# Get parameters of r ~ BetaDistribution(a, b)
		a_mat = posterior::as_draws_matrix(sim$draws("shape1_new")) # N_draws x N_obs
		b_mat = posterior::as_draws_matrix(sim$draws("shape2_new")) # N_draws x N_obs
	
		warning_a_mat = FALSE
		if (min(a_mat) <= 2)
		{
			warning_a_mat = TRUE
			warning("The analytical calculus does not work!")
		}
	
		# Variance of 1/r, dim: N_draws x N_obs, each column being the draws for one individual
		var_1_over_r = b_mat * (a_mat + b_mat - 1) / ((a_mat - 1)^2 * (a_mat - 2))

		# Now multiply by Vbole^2: 1st column with Vbole_1^2, ..., ith column with Vbole_1^2
		sigma_vtot_mat = sweep(var_1_over_r, MARGIN = 2, STATS = Vbole^2, FUN = "*", check.margin = TRUE)

		var_fit_vtot = apply(posterior::as_draws_matrix(sim$draws("v_gen_mean")), 1, var)
		var_res_vtot = apply(sigma_vtot_mat, 1, mean)

		rsq_vtot = var_fit_vtot / (var_fit_vtot + var_res_vtot)

		return(list(rsq_vtot = rsq_vtot, warning = warning_a_mat))
	}

	rsq_vtot_full = rsq_vtot_fct(sim_full)
	rsq_vtot_sub = rsq_vtot_fct(sim_sub)

	## R2 based on LOO, adapted from https://avehtari.github.io/bayes_R2/bayes_R2.html
	loo_R2 = function(sim, psis_object, n_draws, ratio = NULL, vtot = NULL, n = tree_dt[.(sp), .N])
	{
		# LOO weighted predictive expectation for each observation, using posterior draws + PSIS weights
		if ((is.null(ratio)) && !(is.null(vtot)))
		{
			corresponding_variable = "v_gen_mean"
			y = vtot
			if (any(vtot < 0))
				stop("Presence of negative total volumes!")
		}			
		
		if (!(is.null(ratio)) && (is.null(vtot)))
		{
			corresponding_variable = "r_gen_mean"
			y = ratio
			if (any(ratio < 0) || any(ratio > 1))
				stop("Ratio not between 0 and 1!")
		}
		
		if ((!(is.null(ratio)) && !(is.null(vtot))))
			stop("Either ratio or vtot should be left NULL")
		
		mu_loo = loo::E_loo(posterior::as_draws_matrix(sim$draws(corresponding_variable)),
			psis_object = psis_object, type = "mean",
			log_ratios = -posterior::as_draws_matrix(sim$draws("log_lik")))$value
		# Residuals LOO: Observation - weighted predictive expectation (computed above)
		e_loo = y - mu_loo

		rd = bayesboot::rudirichlet(n = n_draws, d = n)
		# Explanation:
		# Rubin's Bayesian bootstrap. The Dirichlet draws represent the posterior
		# uncertainty over what the "true" weight of each observed data point should be in
		# computing a population statistics (like a mean, variance). Each row of rd can be
		# thought of as a "resampling" of the dataset. rowSums(rd) should give a vector of ones.
		# In other words, when I sum the weight of one individual (column), then it is one,
		# as in the whole (unknown) population. However, in each draw (row of the matrix),
		# the weight of individuals does not sum to 1, as some individuals have a heavier
		# weight in the sample, i.e., they represent more individuals of the whole population.

		vary = n/(n - 1) * ( # small-sample bias correction
			rowSums(sweep(rd, 2, y^2, FUN = "*", check.margin = TRUE)) - # weighted mean of y^2 
			rowSums(sweep(rd, 2, y, FUN = "*", check.margin = TRUE))^2 # weighted mean of y
		) # Gives the weighted variance, where rd assign each observation a continuous random weight.

		vareloo = n/(n - 1) * ( # small-sample bias correction
			rowSums(sweep(rd, 2, e_loo^2, FUN = "*", check.margin = TRUE)) - # weighted mean by individual
			rowSums(sweep(rd, 2, e_loo, FUN = "*", check.margin = TRUE))^2 # weighted mean of LOO residuals
		)

		looR2 = 1 - vareloo/vary

		looR2[looR2 < -1] = -1
		looR2[looR2 > 1] = 1

		return(looR2)
	}

	r2_loo_full_r = loo_R2(sim = sim_full, psis_object = loo_full$psis_object, n_draws = n_draws,
		ratio = stanData[["bole_volume_m3_new"]]/stanData[["total_volume_m3_new"]])
	r2_loo_sub_r = loo_R2(sim = sim_sub, psis_object = loo_sub$psis_object, n_draws = n_draws,
		ratio = stanData[["bole_volume_m3_new"]]/stanData[["total_volume_m3_new"]])
	
	r2_loo_full_vtot = loo_R2(sim = sim_full, psis_object = loo_full$psis_object, n_draws = n_draws,
		vtot = stanData[["total_volume_m3_new"]])
	r2_loo_sub_vtot = loo_R2(sim = sim_sub, psis_object = loo_sub$psis_object, n_draws = n_draws,
		vtot = stanData[["total_volume_m3_new"]])

	## Return results
	return (list(
		best = best, comploo = comp, weights = weights,
		warning_loo_full = warning_loo_full, warning_loo_sub = warning_loo_sub, warning = warn,
		rsq_distrib = list(full = rsq_full, sub = rsq_sub),
		rsq_vtot_distrib = list(full = rsq_vtot_full, sub = rsq_vtot_sub),
		rsq_loo_distrib_r = list(full = r2_loo_full_r, sub = r2_loo_sub_r),
		rsq_loo_distrib_v = list(full = r2_loo_full_vtot, sub = r2_loo_sub_vtot),
		rhat_full = rhat_full_m, rhat_sub = rhat_sub_m))
}