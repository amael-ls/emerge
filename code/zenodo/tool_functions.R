#### Aim of script: Contains miscellaneous functions
## Comments
# Line 5: form_vallet, predicts total volume from Vallet 2006
# Line 23: getParams, extract/sum-up draws from Bayesian fit

#' Predict total tree volume using the Vallet et al. (2006) model
#'
#' Predicts total aboveground stem volume from stem circumference and tree
#' height using the model of Vallet et al. (2006). Different model variants
#' are selected through the number of parameters, with a specific formulation
#' available for Douglas-fir.
#'
#' @param circumference_cm Numeric. Stem circumference in centimetres.
#' @param height Numeric. Total tree height in metres.
#' @param params Named numeric vector containing the model parameters
#'   `alpha`, `beta`, and, depending on the model, `gamma` and `delta`.
#' @param n_params Integer. Number of model parameters (2, 3, or 4).
#' @param is_douglas Logical or integer. Whether to use the Douglas-fir
#'   formulation.
#'
#' @return Numeric. Predicted total aboveground volume in cubic metres.
#'
#' @references
#' Vallet, P., Dhôte, J.-F., Moguédec, G. L., Ravart, M., & Pignard, G.
#' (2006). Development of total aboveground volume equations for seven
#' important forest tree species in France. *Forest Ecology and Management*,
#' *229*(1–3), 98–110. \doi{10.1016/j.foreco.2006.03.013}
#'
#' @export

form_vallet = function(circumference_cm, height, params, n_params, is_douglas)
{
	res = params["alpha"] + params["beta"]*circumference_cm

	if (is_douglas == 1)
		return(res*(1 + params["delta"]/circumference_cm^2))
	
	if (n_params >= 3)
		res = res + params["gamma"]*sqrt(circumference_cm)/height

	if (n_params == 4)
		res = res*(1 + params["delta"]/circumference_cm^2)
	
	res = res/(40000*pi) * circumference_cm^2*height # Conversion to volume
	return(res)
}

## Get fixed values parameters (will not work for draws with third dimension > 1
#' Extract parameter values from a CmdStan model
#'
#' Extracts parameter values or summaries from a CmdStan model. Parameters can
#' be returned as all draws, by chain and iteration, as means or medians, or
#' as a set of standard quantiles.
#'
#' @param model_cmdstan A fitted CmdStan model object.
#' @param params_names Character vector of parameter names to extract.
#' @param type Character. Type of values to return: `"all"` for all draws,
#'   `"chain"` for the mean of a specified chain, `"chain-iter"` for a
#'   specific iteration and chain, `"mean"` for the posterior mean,
#'   `"median"` for the posterior median, or `"quantile"` for posterior
#'   quantiles.
#' @param ... Additional arguments required by some `type` values. `chain`
#'   specifies the chain number for `"chain"` and `"chain-iter"`; `iter`
#'   specifies the iteration for `"chain-iter"`. `probs` is currently ignored
#'   for `"quantile"`.
#'
#' @return A numeric vector for `"chain"`, `"chain-iter"`, `"mean"`, and
#'   `"median"`; a data.table of summary statistics for `"quantile"`; or the
#'   CmdStan draws object for `"all"`.
#'
#' @note
#' The function currently only supports scalar parameters (i.e. parameters
#' whose draws have a third dimension of 1).
#'
#' @export

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

#' Plot MCMC trace plots
#'
#' Creates trace plots for MCMC draws, with a separate colour for each chain.
#' The function provides a lightweight alternative to `bayesplot` and can
#' optionally save the plot as a PDF.
#'
#' @param draws An MCMC draws array, such as a `draws_array` object from
#'   `cmdstanr`, with dimensions iteration, chain, and parameter.
#' @param filename Character. Optional filename (without extension) for saving
#'   the plot as a PDF.
#' @param ... Optional plotting arguments. `xlab`, `ylab`, and `main` specify
#'   axis labels and the title; `scaling` scales the values plotted; `chain_id`
#'   selects chains to display; `val1`, `val2`, etc. add horizontal reference
#'   lines; `label1`, `label2`, etc. label the corresponding reference lines;
#'   and `iter1`, `iter2`, etc. add vertical lines at specified iterations.
#'
#' @return No value is returned. A trace plot is produced in the current
#'   graphics device and, if `filename` is provided, saved as a PDF.
#'
#' @export

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

#' Plot parameter pairs associated with divergent transitions
#'
#' Produces pairwise diagnostic plots of selected model parameters, highlighting
#' posterior draws associated with divergent transitions during Bayesian
#' sampling.
#'
#' @param fit A fitted CmdStan model object.
#' @param div An MCMC draws array, such as a `draws_array` object from
#'   `cmdstanr`, with dimensions iteration, chain. Values equal to `1` identify
#'   divergent draws.
#'
#' @return No value is returned. A series of pairwise scatter plots is
#'   produced in the current graphics device.
#'
#' @details
#' The diagnostic plots show pairwise relationships among the parameters
#' `c`, `j`, `k`, `m`, `n`, `s`, and `tau`, with divergent draws highlighted
#' separately from the remaining posterior draws.
#'
#' @export

plot_divergences = function(fit, div)
{
	params = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s", "tau"), type = "all")

	div_c = posterior::subset_draws(params[, , "c"], draw = which(div == 1))
	div_j = posterior::subset_draws(params[, , "j"], draw = which(div == 1))
	div_k = posterior::subset_draws(params[, , "k"], draw = which(div == 1))
	div_m = posterior::subset_draws(params[, , "m"], draw = which(div == 1))
	div_n = posterior::subset_draws(params[, , "n"], draw = which(div == 1))
	div_s = posterior::subset_draws(params[, , "s"], draw = which(div == 1))
	div_tau = posterior::subset_draws(params[, , "tau"], draw = which(div == 1))

	# Plot divergence for m and...
	# ... c
	plot(params[, , "m"], params[, , "c"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "c", main = "m and c")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and tau
	plot(params[, , "m"], params[, , "tau"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "tau", main = "m and tau")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and j
	plot(params[, , "m"], params[, , "j"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "j", main = "m and j")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and n
	plot(params[, , "m"], params[, , "n"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "n", main = "m and n")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and s
	plot(params[, , "m"], params[, , "s"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "s", main = "m and s")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")



	# Plot divergence for c and...
	# ... and tau
	plot(params[, , "c"], params[, , "tau"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "tau", main = "c and tau")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and j
	plot(params[, , "c"], params[, , "j"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "j", main = "c and j")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and n
	plot(params[, , "c"], params[, , "n"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "n", main = "c and n")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and s
	plot(params[, , "c"], params[, , "s"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "c", ylab = "s", main = "c and s")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")



	# Plot divergence for tau and...
	# ... and j
	plot(params[, , "tau"], params[, , "j"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "tau", ylab = "j", main = "tau and j")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and n
	plot(params[, , "tau"], params[, , "n"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "tau", ylab = "n", main = "tau and n")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and s
	plot(params[, , "tau"], params[, , "s"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "tau", ylab = "s", main = "tau and s")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")



	# Plot divergence for j and...
	# ... and n
	plot(params[, , "j"], params[, , "n"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "j", ylab = "n", main = "j and n")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")

	# ... and s
	plot(params[, , "j"], params[, , "s"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "j", ylab = "s", main = "j and s")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")



	# Plot divergence for n and...
	# ... and s
	plot(params[, , "n"], params[, , "s"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "n", ylab = "s", main = "n and s")
	points(div_c, div_j, pch = 18, col = "#0F7BA2")
}

## Function to plot joint posterior of two parameters
plot_joint = function(fit, params, div)
{
	if (length(params) > 6)
		stop("Too many parameters, choose 6 params at most")
	
}

## Function to plot joint posterior of two parameters
plot_joint_dirty = function(fit, div)
{
	params = getParams(model_cmdstan = fit, params_names = c("alpha", "beta_", "gamma", "delta"), type = "all")
	div_a = posterior::subset_draws(params[, , "alpha"], draw = which(div == 1))
	div_b = posterior::subset_draws(params[, , "beta_"], draw = which(div == 1))
	div_c = posterior::subset_draws(params[, , "gamma"], draw = which(div == 1))
	div_d = posterior::subset_draws(params[, , "delta"], draw = which(div == 1))

	# Plot divergence for alpha and...
	# ... and beta_
	plot(params[, , "alpha"], params[, , "beta_"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "alpha", ylab = "beta_", main = "alpha and beta_")
	points(div_a, div_b, pch = 18, col = "#0F7BA2")

	# ... and gamma
	plot(params[, , "alpha"], params[, , "gamma"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "alpha", ylab = "gamma", main = "alpha and gamma")
	points(div_a, div_c, pch = 18, col = "#0F7BA2")

	# ... and delta
	plot(params[, , "alpha"], params[, , "delta"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "alpha", ylab = "delta", main = "alpha and delta")
	points(div_a, div_d, pch = 18, col = "#0F7BA2")



	# Plot divergence for beta_ and...
	# ... and gamma
	plot(params[, , "beta_"], params[, , "gamma"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "beta_", ylab = "gamma", main = "beta_ and gamma")
	points(div_b, div_c, pch = 18, col = "#0F7BA2")

	# ... and delta
	plot(params[, , "beta_"], params[, , "delta"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "beta_", ylab = "delta", main = "beta_ and delta")
	points(div_b, div_d, pch = 18, col = "#0F7BA2")



	# Plot divergence for gamma and...
	# ... and delta
	plot(params[, , "gamma"], params[, , "delta"], pch = 19, cex = 0.65, col = "#FAB255",
		xlab = "gamma", ylab = "delta", main = "gamma and delta")
	points(div_c, div_d, pch = 18, col = "#0F7BA2")
}

## Function to plot a species fit (simplif = TRUE for submodel)
#' Plot species-specific model fit and diagnostics
#'
#' Produces diagnostic plots for a fitted species-specific Bayesian model,
#' including the observed bole-to-total volume relationship and, optionally,
#' predicted versus observed total volume. Divergent draws and problematic
#' parameter convergence are also reported and visualised.
#'
#' @param fit A fitted CmdStan model object.
#' @param sp Character. Species identifier used to select the corresponding
#'   forest data.
#' @param forest Data table containing tree observations. By default, the
#'   observations for `sp` are selected from `tree_dt`.
#' @param pred Logical. Whether to plot and return model predictions of total
#'   volume.
#' @param simplif Logical. Whether the simplified model formulation is used.
#' @param n_bins Integer. Number of classes used to colour observations when
#'   `selected_variable` is provided.
#' @param pal Character. Name of the `MetBrewer` palette used for colouring
#'   observations.
#' @param selected_variable Character or `NULL`. Optional column of `forest`
#'   used to classify and colour observations.
#' @param print_plot Logical. Whether to produce the diagnostic plots. If
#'   `FALSE`, only diagnostic information and, when requested, predictions are
#'   returned.
#'
#' @return A list containing `any_div`, indicating whether divergent
#'   transitions occurred, `loc`, containing the divergence indicators, and
#'   `rhats`, containing parameter R-hat values. When `print_plot = FALSE`
#'   and `pred = TRUE`, the list also contains predicted and observed volumes
#'   and the fitted parameter values.
#'
#' @export

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
#' Plot group-specific model fit and diagnostics
#'
#' Produces diagnostic plots for a fitted group-specific Bayesian model,
#' including the observed bole-to-total volume relationship and, optionally,
#' predicted versus observed total volume. Divergent draws and problematic
#' parameter convergence are also reported and visualised.
#'
#' @param fit A fitted CmdStan model object.
#' @param gr Character. Group identifier used to select the corresponding
#'   forest data.
#' @param forest Data table containing tree observations. By default, the
#'   observations for `gr` are selected from `tree_dt`.
#' @param pred Logical. Whether to plot and return model predictions of total
#'   volume.
#' @param simplif Logical. Whether the simplified model formulation is used.
#' @param n_bins Integer. Number of classes used to colour observations when
#'   `selected_variable` is provided.
#' @param pal Character. Name of the `MetBrewer` palette used for colouring
#'   observations.
#' @param selected_variable Character or `NULL`. Optional column of `forest`
#'   used to classify and colour observations.
#' @param print_plot Logical. Whether to produce the diagnostic plots. If
#'   `FALSE`, only diagnostic information and, when requested, predictions are
#'   returned.
#'
#' @return A list containing `any_div`, indicating whether divergent
#'   transitions occurred, `loc`, containing the divergence indicators, and
#'   `rhats`, containing parameter R-hat values. When `print_plot = FALSE`
#'   and `pred = TRUE`, the list also contains predicted and observed volumes
#'   and the fitted parameter values.
#'
#' @export

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
#' Rebuild model comparison results
#'
#' Reconstructs model comparison tables from saved results for each species.
#' Summarises model weights, LOO comparison statistics, R-squared measures,
#' and R-hat diagnostics.
#'
#' @param save_ls Named list containing the saved model results for each
#'   species. Each element is expected to contain the model weights, LOO
#'   comparison results, R-squared distributions, and R-hat diagnostics used
#'   by the function.
#'
#' @return A list containing four data tables: `comp` with LOO model
#'   comparisons, `weights_dt` with model weights, `R2D2` with R-squared
#'   summaries, and `rhat` with R-hat diagnostics.
#'
#' @export

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

## Function to compare full model with submodel, based on PSIS-LOO
#' Compare full and submodels using PSIS-LOO
#'
#' Compares the full and submodel using PSIS-LOO, model weights, R-hat
#' diagnostics, and several Bayesian R-squared measures. Posterior predictive
#' quantities are generated for both models and optional diagnostic plots can
#' be produced.
#'
#' @param sp Character. Species identifier used to select the corresponding
#'   observations and fitted models.
#' @param tree_dt Data table containing tree observations, including
#'   `bole_volume_m3` and `total_volume_m3`.
#' @param path_models Character. Path to the Stan model files.
#' @param path_output Character. Path to the saved fitted models.
#' @param woodstock_seed Numeric. Random seed used for posterior predictive
#'   simulation.
#' @param printPlot Logical. Whether to plot predicted values from the two
#'   models against each other.
#'
#' @return A list containing the best model, PSIS-LOO comparison and model
#'   weights, LOO warnings, Bayesian R-squared distributions for the ratio and
#'   total volume, LOO-based R-squared distributions, and maximum R-hat values
#'   for the full and submodels.
#'
#' @export

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
	full_file = paste0(path_output, stringi::stri_replace(str = sp, regex = " ", replacement = "-"), "_fullmodel_theta")
	sub_file = paste0(path_output, stringi::stri_replace(str = sp, regex = " ", replacement = "-"), "_submodel")

	if (sp == "Quercus sp.")
	{
		full_file = stringi::stri_replace(str = full_file, regex = "._", replacement = "_")
		sub_file = stringi::stri_replace(str = sub_file, regex = "._", replacement = "_")
	}

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

## Function to compute "traditional" RMSE
#' Compute RMSE and MAPE from posterior predictions
#'
#' Computes the traditional root mean squared error (RMSE) and mean absolute
#' percentage error (MAPE) between posterior predictions of total tree volume
#' and observed total volume.
#'
#' @param sp Character. Species identifier used to select the observations and
#'   fitted model.
#' @param tree_dt Data table containing tree observations, including
#'   `bole_volume_m3` and `total_volume_m3`.
#' @param path_output Character. Path to the saved fitted model.
#' @param path_models Character. Path to the Stan model files used to generate
#'   posterior predictions.
#' @param is_simplif Logical. Whether to use the simplified (submodel)
#'   formulation instead of the full model.
#' @param woodstock_seed Numeric. Random seed used for posterior predictive
#'   simulation.
#'
#' @return A list containing posterior distributions of `rmse` and `mape`,
#'   with one value per posterior draw.
#'
#' @export

RMSE_bayes = function(sp, tree_dt, path_output, path_models, is_simplif = FALSE,
	woodstock_seed = 1969 - 08 - 18)
{
	#### Generate quantitities
	## Compile generator models
	genQ = cmdstanr::cmdstan_model(paste0(path_models, "fullmodel-genQ.stan"))
	filename = paste0(path_output, stringi::stri_replace(str = sp, regex = " ", replacement = "-"),
		"_fullmodel_theta.rds")

	if (is_simplif)
	{
		genQ = cmdstanr::cmdstan_model(paste0(path_models, "submodel-genQ.stan"))
		filename = paste0(path_output, stringi::stri_replace(str = sp, regex = " ", replacement = "-"),
			"_submodel.rds")
	}

	if (sp == "Quercus sp.")
	{
		filename = stri_replace_first(str = filename, regex = "\\._fullmodel",
			replacement = "_fullmodel")
		filename = stri_replace_first(str = filename, regex = "\\._submodel",
			replacement = "_submodel")
	}

	## Prepare data
	stanData = list(
		N = tree_dt[.(sp)][, .N],
		N_new = tree_dt[.(sp)][, .N],
		bole_volume_m3 = tree_dt[.(sp)][, bole_volume_m3],
		bole_volume_m3_new = tree_dt[.(sp)][, bole_volume_m3],
		total_volume_m3 = tree_dt[.(sp)][, total_volume_m3],
		total_volume_m3_new = tree_dt[.(sp)][, total_volume_m3]
	)

	## Compute rmse with respect to the mean
	fit = readRDS(filename)
	n_chains = fit$num_chains()
	sim = genQ$generate_quantities(fitted_params = fit, data = stanData,
		seed = woodstock_seed, parallel_chains = min(n_chains, 4))

	sim = posterior::as_draws_matrix(sim$draws("v_gen_mean"))
	rmse = numeric(dim(sim)[1])
	mape = numeric(dim(sim)[1])

	n_indiv = tree_dt[.(sp)][, .N]

	for (i in seq_along(rmse))
	{
		rmse[i] = sqrt(1/n_indiv * sum((sim[i,] - tree_dt[.(sp)][, total_volume_m3])^2))
		mape[i] = 100/n_indiv * sum(abs((sim[i,] - tree_dt[.(sp)][, total_volume_m3])/tree_dt[.(sp)][, total_volume_m3]))
	}

	return(list(rmse = rmse, mape = mape))
}

## Function to compute the average ratio without params/residual uncertainty
#' Predict the average bole-to-total volume ratio
#'
#' Computes the average bole-to-total volume ratio from model parameters,
#' without accounting for parameter or residual uncertainty.
#'
#' @param x Numeric. Predictor value, typically tree size or bole volume.
#' @param pars Named numeric vector or one-row data table containing the
#'   parameters `c`, `j`, `k`, `m`, `n`, and `s`.
#' @param expansion_factor Logical. If `TRUE`, returns the expansion factor
#'   (the inverse of the predicted ratio) instead of the ratio.
#'
#' @return Numeric. Predicted average bole-to-total volume ratio, or its inverse
#'   when `expansion_factor = TRUE`.
#'
#' @export

pred_ratio = function(x, pars, expansion_factor = FALSE)
{
	if (is.data.table(pars))
	{
		if (pars[, .N] != 1)
			stop("pars shoud correspond to one species only")
		pars = c(c = pars[, c], j = pars[, j], k = pars[, k], m = pars[, m], n = pars[, n], s = pars[, s])
	}
	ratio = (pars["m"] - pars["c"]) * exp(pars["j"] - pars["k"]*x) * (pars["k"]*x/pars["j"])^pars["j"] +
		pars["c"] - (pars["c"] - pars["n"])*exp(-pars["s"]*x)

	if (expansion_factor)
		return (1/ratio)
	
	return(ratio);
}

## Function to compute the average above-ground pred volume without params/residual uncertainty
#' Predict average above-ground volume
#'
#' Computes the average above-ground tree volume from the predictor and model
#' parameters, without accounting for parameter or residual uncertainty.
#'
#' @param x Numeric. Bole volume used to predict total volume.
#' @param pars Named numeric vector or one-row data table containing the
#'   parameters required by `pred_ratio`.
#'
#' @return Numeric. Predicted average above-ground tree volume.
#'
#' @export

pred_vol = function(x, pars)
	return(x/pred_ratio(x, pars))

## Logit and inv_logit function
logit = function(x)
	return(log(x/(1 - x)))

inv_logit = function(x)
	return(1/(1 + exp(-x)))
