#### Aim of script: Compute averaged VEF from my own values
## Comments:
# One of the problem, in computing individual VEFs, then averaging by plot, to finally average within
#	biogeographical subregions (SER86 variable in the French NFI) is that each plot get the same weight,
#	which is not necessarily the case, especially in terms of resources (an almost empty plot should weight
#	less than a full plot!). This is where the function `ocre_nm` is useful.
#	This is the output I use in the article.
#

#### Load packages
library(data.table)
library(inventR)
library(stringi)
library(terra)

#### Prepare data
## Get path
source("./global_variables.R")

## Get map
ser_map_ocre = vect(inventR::ser_simp)

## Connect to db
db = connect_db(base = "exploitation", serveur = "test-inv-exp")

## Get species names from database
speciesName_nfi = inventR::ListUnitMod("ESPAR1")[, .(sp_code = mode, speciesName_sci = definition)]
setkey(speciesName_nfi, sp_code)

## Modify few names, removing subsp. and var.
speciesName_nfi[stri_detect(str = speciesName_sci, regex = "var."),
	speciesName_sci := stri_sub(str = speciesName_sci,
		to = stri_locate(str = speciesName_sci, regex = "var.")[, "start"] - 2)]

speciesName_nfi[stri_detect(str = speciesName_sci, regex = "subsp."),
	speciesName_sci := stri_sub(str = speciesName_sci,
		to = stri_locate(str = speciesName_sci, regex = "subsp.")[, "start"] - 2)]

## Manually set for 29AF (broadleaf sp.) and 68CE (conifer sp.)
speciesName_nfi["29AF", speciesName_sci := "Broadleaf"]
speciesName_nfi["68CE", speciesName_sci := "Conifer"]

## Manually set for 53CA (Pinus nigra var. calabrica) and 53CO (Pinus nigra var. corsicana)
speciesName_nfi["53CA", speciesName_sci := "Pinus laricio"]
speciesName_nfi["53CO", speciesName_sci := "Pinus laricio"]

## Manually set for 58 (Pinus mugo subsp. uncinata)
speciesName_nfi["58", speciesName_sci := "Pinus uncinata"]

#### Compute expansion factors
## Average expansion factors per ser
ser_expF_ocre = ocre_nm(vars = c("U_V0_ALAMOD", "V"), vent_p = "ser_86", ans = 2020:2024,
	test = TRUE, us_nm = "1")[visite == 1, .(ser_86 = ser_86.ser_86, vars, tot)] |>
	dcast(formula = ser_86 ~ vars, value.var = "tot")

ser_expF_ocre[, avg_exp_factor := U_V0_ALAMOD/V]

ser_map_ocre = terra::merge(ser_map_ocre, ser_expF_ocre[, .(ser_86, avg_exp_factor)],
	by.x = "code_ser", by.y = "ser_86", all.x = TRUE)

## Average expansion factors per tree functional type (broadleaf/conifer)
broadConif_expF_ocre = ocre_nm(vars = c("U_V0_ALAMOD", "V"), vent_a = "ess:essd1", ans = 2020:2024,
	test = TRUE, us_nm = "1")[visite == 1, .(fct_type = ess.essd1, vars, tot)] |>
	dcast(formula = fct_type ~ vars, value.var = "tot")

broadConif_expF_ocre[, avg_exp_factor := U_V0_ALAMOD/V]

## Average expansion factors per ser and tree functional type (broadleaf/conifer)
serBC_expF_ocre = ocre_nm(vars = c("U_V0_ALAMOD", "V"), vent_p = "ser_86", vent_a = "ess:essd1", ans = 2020:2024,
	test = TRUE, us_nm = "1")[visite == 1, .(ser_86 = ser_86.ser_86, fct_type = ess.essd1, vars, tot)] |>
	dcast(formula = ser_86 + fct_type ~ vars, value.var = "tot")

serBC_expF_ocre[, avg_exp_factor := U_V0_ALAMOD/V]

dt_BC_expF = serBC_expF_ocre[, c(as.list(quantile(avg_exp_factor, probs = c(0, 0.025, 0.5, 0.975, 1)),
	mean = mean(avg_exp_factor))),
	by = fct_type]

dt_BC_expF = serBC_expF_ocre[, {
	q = quantile(avg_exp_factor, probs = c(0, 0.025, 0.5, 0.975, 1))
	names(q) = c("min", "q025", "med", "q975", "max")
	c(as.list(q), mean = mean(avg_exp_factor))
}, by = fct_type]

filename = paste0(path_pgfplotstable, "expF_broad-conif.csv")
if (!file.exists(filename))
	fwrite(dt_BC_expF, filename)

## Average expansion factors per species
spSpecific_expF_ocre = ocre_nm(vars = c("U_V0_ALAMOD", "V"), vent_a = "ess:ess", ans = 2020:2024,
	test = TRUE, us_nm = "1")[visite == 1, .(ess = ess.ess, vars, tot)] |>
	dcast(formula = ess ~ vars, value.var = "tot")

spSpecific_expF_ocre[, avg_exp_factor := U_V0_ALAMOD/V]

spSpecific_expF_ocre = spSpecific_expF_ocre[, .(ess, avg_exp_factor)]

# nfi_data[, n_indiv := .N, by = speciesName_sci]
# select_sp_code = nfi_data[n_indiv > 1000, unique(sp_code)]

spSpecific_expF_ocre[ess < 50, round(mean(avg_exp_factor), 2)] # 1.35
# spSpecific_expF_ocre[select_sp_code][ess < 50][which.min(avg_exp_factor)] # 05, Quercus pubescens, 1.19
# spSpecific_expF_ocre[select_sp_code][ess < 50][which.max(avg_exp_factor)] # 40 Arbutus unedo, 1.73
spSpecific_expF_ocre[ess > 50, round(mean(avg_exp_factor), 2)] # 1.19
# spSpecific_expF_ocre[select_sp_code][ess > 50][which.min(avg_exp_factor)] # 63 Larix decidua, 1.09
# spSpecific_expF_ocre[select_sp_code][ess > 50][which.max(avg_exp_factor)] # 51 Pinus pinaster, 1.25

spSpecific_expF_ocre[
	ess %in% c("02", "03", "09", "17", "51", "52", "53", "54", "56", "58", "61", "62", "63", "64")][
	ess < 50, round(mean(avg_exp_factor), 2)] # 1.24

spSpecific_expF_ocre[
	ess %in% c("02", "03", "09", "17", "51", "52", "53", "54", "56", "58", "61", "62", "63", "64")][
	ess > 50, round(mean(avg_exp_factor), 2)] # 1.16

## For Abies
sp_ind = unique(speciesName_nfi[stri_detect(str = speciesName_sci, regex = "^Abies"),
	stri_sub(str = sp_code, to = 2)]) # to = 2 in order to remove the sub-species letters and keep rough groups
spSpecific_expF_ocre[
	ess %in% sp_ind, round(mean(avg_exp_factor), 2)] # 1.19

## For Fagus
sp_ind = unique(speciesName_nfi[stri_detect(str = speciesName_sci, regex = "Fagus"),
	stri_sub(str = sp_code, to = 2)]) # to = 2 in order to remove the sub-species letters and keep rough groups
spSpecific_expF_ocre[
	ess %in% sp_ind, round(mean(avg_exp_factor), 2)] # 1.29

## For Fraxinus
sp_ind = unique(speciesName_nfi[stri_detect(str = speciesName_sci, regex = "Fraxinus"),
	stri_sub(str = sp_code, to = 2)]) # to = 2 in order to remove the sub-species letters and keep rough groups
spSpecific_expF_ocre[
	ess %in% sp_ind, round(mean(avg_exp_factor), 2)] # 1.29

## For Larix
sp_ind = unique(speciesName_nfi[stri_detect(str = speciesName_sci, regex = "Larix"),
	stri_sub(str = sp_code, to = 2)]) # to = 2 in order to remove the sub-species letters and keep rough groups
spSpecific_expF_ocre[
	ess %in% sp_ind, round(mean(avg_exp_factor), 2)] # 1.09

## For Picea
sp_ind = unique(speciesName_nfi[stri_detect(str = speciesName_sci, regex = "Picea"),
	stri_sub(str = sp_code, to = 2)]) # to = 2 in order to remove the sub-species letters and keep rough groups
spSpecific_expF_ocre[
	ess %in% sp_ind, round(mean(avg_exp_factor), 2)] # 1.17

## For Pinus
sp_ind = unique(speciesName_nfi[stri_detect(str = speciesName_sci, regex = "Pinus"),
	stri_sub(str = sp_code, to = 2)]) # to = 2 in order to remove the sub-species letters and keep rough groups
spSpecific_expF_ocre[
	ess %in% sp_ind, round(mean(avg_exp_factor), 2)] # 1.16

## For Pseudotsuga
sp_ind = unique(speciesName_nfi[stri_detect(str = speciesName_sci, regex = "Pseudotsuga"),
	stri_sub(str = sp_code, to = 2)]) # to = 2 in order to remove the sub-species letters and keep rough groups
spSpecific_expF_ocre[
	ess %in% sp_ind, round(mean(avg_exp_factor), 2)] # 1.21

## For Quercus
sp_ind = unique(speciesName_nfi[stri_detect(str = speciesName_sci, regex = "Quercus"),
	stri_sub(str = sp_code, to = 2)]) # to = 2 in order to remove the sub-species letters and keep rough groups
spSpecific_expF_ocre[
	ess %in% sp_ind, round(mean(avg_exp_factor), 2)] # 1.23

#? -----------------------------------------------------------------------------------------------------
#! -----------------------------------------------------------------------------------------------------
## ---------------------------------------    CRASH TEST ZONE    ---------------------------------------
#! -----------------------------------------------------------------------------------------------------
#? -----------------------------------------------------------------------------------------------------

# test = ocre_nm(vars = c("U_V0_ALAMOD", "V"), vent_a = c("ess:ess", "clad:u_dim_g3"), ans = 2020:2024,
# 	test = TRUE, us_nm = "1")[visite == 1, .(ess = ess.ess, clad = clad.u_dim_g3, vars, tot)] |>
# 	dcast(formula = ess + clad ~ vars, value.var = "tot")

# test[, vef := U_V0_ALAMOD/V]
# setkey(test, ess, clad)

# #### French NFI data query
# db = connect_db(base = "exploitation", serveur = "test-inv-exp") # Connect to test server

# nfi_data = exec_req(conn = db, req = paste0("
# 	SELECT
# 		tree.NPP AS plot_id, tree.A AS tree_id, tree.ESPAR AS sp_code,
# 		CASE
# 			WHEN tree.clad::numeric < 12.5 THEN
# 				'PB'
# 			WHEN 12.5 <= tree.clad::numeric AND tree.clad::numeric < 22.5 THEN
# 				'PB'
# 			WHEN 22.5 <= tree.clad::numeric AND tree.clad::numeric < 37.5 THEN
# 				'GB'
# 			WHEN 37.5 <= tree.clad::numeric AND tree.clad::numeric < 52.5 THEN
# 				'GB'
# 			ELSE
# 				'GB'
# 		END AS class,
# 		tree.ESS AS sp_code_simp, tree.W, tree.c13 AS circumference_m, tree.htot AS height,
# 		tree.V AS bole_volume,
# 		CASE
# 			WHEN tree.orig = '0' THEN
# 				'seed'
# 			WHEN tree.orig = '1' THEN
# 				'sprout'
# 		END AS origin,
# 		CASE
# 			WHEN tree.ESS < '50' THEN
# 				'broadleaf'
# 			ELSE
# 				'conifer'
# 		END AS type,
# 		u_tree.u_v0_alamod AS total_volume,
# 		forest.ESSPRE AS dom_sp, forest.SVER AS structure, forest.GEST AS management,
# 		coords1.XL AS x_EPSG_27572, coords1.YL AS y_EPSG_27572,
# 		coords2.SER_86, coords2.GRECO

# 	FROM
# 		inv_exp_nm.g3arbre AS tree

# 		LEFT JOIN inv_exp_nm.u_g3arbre AS u_tree
# 		ON tree.NPP = u_tree.NPP AND tree.A = u_tree.A
		
# 		LEFT JOIN inv_exp_nm.g3foret AS forest
# 		ON tree.NPP = forest.NPP
		
# 		LEFT JOIN inv_exp_nm.e1point AS coords1
# 		ON tree.NPP = coords1.NPP

# 		LEFT JOIN inv_exp_nm.e2point AS coords2
# 		ON tree.NPP = coords2.NPP

# 	WHERE
# 		tree.VEGET = '0' AND -- Living standing trees
# 		tree.V > 0 AND -- I found some that are 0!
# 		tree.INCREF IN ('15', '16', '17', '18', '19') -- Campaign from 2020 to 2024

# 	ORDER BY
# 		plot_id, tree_id;"), DT = TRUE)

# disconnect_db(db)


# plot(nfi_data[(sp_code_simp == "02") & bole_volume < 1, bole_volume],
# 	nfi_data[(sp_code_simp == "02") & bole_volume < 1, circumference_m],
# 	axes = FALSE, xlab = "Volume BFT", ylab = "Circonférence", pch = 19, cex = 0.75, main = "Ess = 02")
# axis(1)
# axis(2, las = 1)
# abline(h = 22.5*pi/100, lwd = 2, col = "#CD212A", lty = "dashed")
# abline(v = 0.55, lwd = 2, col = "#CD212A", lty = "dashed")


# nfi_data[, vef := total_volume/bole_volume]
# setkey(nfi_data, sp_code_simp, class)




# dd1 = density(nfi_data[.("02", "PB")][, sum(w*vef)/sum(w), by = .(plot_id)][, V1])
# dd2 = density(nfi_data[.("02", "GB")][, sum(w*vef)/sum(w), by = .(plot_id)][, V1])

# plot(dd1$x, dd1$y, lwd = 4, col = "#FAB255", type = "l", ylim = c(0, 10))
# lines(dd2$x, dd2$y, lwd = 4, col = "#0F7BA2")
