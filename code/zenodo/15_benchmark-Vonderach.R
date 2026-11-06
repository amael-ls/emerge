#### Aim of script: Compare my predictions with Vonderach 2018
## Comments
# Vonderach et al. 2018,
#	Consistent set of additive biomass functions for eight tree species in Germany fit
#	by nonlinear seemingly unrelated regression.
#	Ann For Sci 75:49. https://doi.org/10.1007/s13595-018-0728-4
# published allometric equations for biomass. Combining with Cuny 2025 et al.
#	https://doi.org/10.6084/m9.figshare.c.7515396
# I can deduce a volume
#

#### Load packages
options(max.print = 500)
source("./global_variables.R")

library(data.table)

#### Tool functions
## Function to predict above-ground biomass and volume from Vonderach2018 and Cuny2025
pred_vonderach = function(dbh_cm, height_m, sp, params_agb, wood_params)
{
	coarse_wood_including_bark = params_agb[.(sp, "cww", "intercept"), value] *
		dbh_cm^params_agb[.(sp, "cww", "dbh"), value] *
		height_m^params_agb[.(sp, "cww", "height"), value] + # Coarse wood without bark
		params_agb[.(sp, "cwb", "intercept"), value] *
		dbh_cm^params_agb[.(sp, "cwb", "dbh"), value] *
		height_m^params_agb[.(sp, "cwb", "height"), value] # Coarse wood bark

	small_wood = params_agb[.(sp, "sw", "base"), value] +
		params_agb[.(sp, "sw", "intercept"), value] *
		dbh_cm^params_agb[.(sp, "sw", "dbh"), value] *
		height_m^params_agb[.(sp, "sw", "height"), value] # Small wood

	agb = coarse_wood_including_bark + small_wood
	return(list(agb = agb, vtot = agb/wood_params[sp, dens]))
}

#### Parameters
## From Vonderach2018
sp_params = data.table(species = rep("Picea abies", 10),
	compartment = c(rep("cww", 3), rep("cwb", 3), rep("sw", 4)),
	parameter = c(rep(c("intercept", "dbh", "height"), 3), "base"),
	value = c(0.0142, 1.7414, 1.2401,
		0.0038, 1.6076, 1.0588,
		0.0243, 2.9671, -0.8183, 1.8472),
	key = c("species", "compartment", "parameter"))

## From Cuny2025
sp_wood = data.table(species = "Picea abies", dens = 393.8, key = "species")

#### Predictions on the above-ground volume dataset
## Load data
# Fitting dataset of my model, with truth
tree_dt = readRDS(paste0(path_data, "tree_dt_14species.rds"))[speciesName_sci == "Picea abies",
	.(speciesName_sci, plot_id, tree_id, vtot_obs = total_volume_m3,
		dbh_cm = 100*circumference_m/pi, height_m = height)]

# My predictions for Picea abies
my_pred = fread(paste0(path_pgfplotsfig, "vallet-mine_Picea-abies.csv"))

## Predictions
tree_dt[, vonderach := pred_vonderach(dbh_cm, height_m, unique(speciesName_sci),
	params_agb = sp_params, wood_params = sp_wood)[["vtot"]], by = speciesName_sci]

pdf("./vonderach_Picea-abies.pdf", height = 6, width = 6)
plot(my_pred[, obs], my_pred[, mine], axes = FALSE, pch = 19, cex = 0.85, col = "#0F7BA266",
	xlab = "Observations", ylab = "Predictions", main = "Picea abies")
points(tree_dt[, vtot_obs], tree_dt[, vonderach], pch = 19, cex = 0.65, col = "#FAB25544")
axis(1)
axis(2, las = 1)
abline(a = 0, b = 1, lwd = 0.8, lty = "dashed", col = "#5A5A5A")
dev.off()

#### Predictions on the French NFI
## Load data
vtot_alamod = readRDS("./U_V0_ALAMOD-base-test.rds") |>
	setnames(new = c("plot_id", "tree_id", "vtot_alamod"))

db = inventR::connect_db(base = "exploitation", serveur = "inv-exp") # Connect to test server

nfi_data = inventR::exec_req(conn = db, req = paste0("
	SELECT
		tree.NPP AS plot_id, tree.A AS tree_id,
		tree.c13 * 100 / PI() AS dbh_cm, tree.htot AS height_m,
		tree.V AS bole_volume, u_tree.U_V0 AS vtot_nfi

	FROM
		inv_exp_nm.g3arbre AS tree
	
	LEFT JOIN
		inv_exp_nm.u_g3arbre AS u_tree
	ON
		tree.NPP = u_tree.NPP AND tree.A = u_tree.A

	WHERE
		tree.VEGET = '0' AND -- Living standing trees
		tree.ess = '62' AND -- Picea abies
		tree.htot IS NOT NULL AND -- Height is measured
		tree.V > 0 AND -- I found some that are 0! That's because of culls set to 100%
		tree.INCREF IN ('15', '16', '17', '18', '19') -- Campaign from 2020 to 2024

	ORDER BY
		plot_id, tree_id;"), DT = TRUE)

inventR::disconnect_db(db)

## Predict Vonderach
nfi_data = merge.data.table(nfi_data, vtot_alamod, by = c("plot_id", "tree_id"), all.x = TRUE)
nfi_data[, any(is.na(vtot_alamod))]

nfi_data[, vtot_vonderach := pred_vonderach(dbh_cm, height_m, "Picea abies",
	params_agb = sp_params, wood_params = sp_wood)[["vtot"]]]

y_rg = c(0, max(nfi_data[, .(vtot_nfi, vtot_alamod, vtot_vonderach)]))

pdf("./nfi-vonderach_Picea-abies.pdf", height = 6, width = 6)
plot(nfi_data[, bole_volume], nfi_data[, vtot_alamod], axes = FALSE, pch = 19, cex = 1, col = "#0F7BA266",
	xlab = "Bole volume", ylab = "Predictions", main = "Picea abies", ylim = y_rg)
points(nfi_data[, bole_volume], nfi_data[, vtot_vonderach], pch = 19, cex = 0.75, col = "#FAB25544")
points(nfi_data[, bole_volume], nfi_data[, vtot_nfi], pch = 19, cex = 0.15, col = "#DD5129")
axis(1)
axis(2, las = 1)
abline(a = 0, b = 1, lwd = 0.8, lty = "dashed", col = "#5A5A5A")
dev.off()

