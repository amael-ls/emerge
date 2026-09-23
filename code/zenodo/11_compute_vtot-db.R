#### Aim of script: Compute the total volume for trees in database (campaigns 2020--2024)
## Comments
# This file will not run outside of the French NFI
# This script predicts the above-ground volume of trees in the French NFI database.
# Note that they do not support posteriors, so predictions are only done with
# the averaged coefficients. I checked and found that average of posterior predictions
# matches medians and also predictions done with averaged coefficients.

#### Load packages
library(data.table)
library(inventR)
library(stringi)

#### Load data
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

## French NFI data query
db = connect_db(base = "exploitation", serveur = "inv-exp") # Connect to test server

nfi_data = exec_req(conn = db, req = paste0("
	SELECT
		tree.NPP AS plot_id, tree.A AS tree_id, tree.ESPAR AS sp_code,
		tree.ESS AS sp_code_simp, tree.V AS bole_volume,
		CASE
			WHEN tree.ESS < '50' THEN
				'broadleaf'
			ELSE
				'conifer'
		END AS type

	FROM
		inv_exp_nm.g3arbre AS tree

	WHERE
		tree.VEGET = '0' AND -- Living standing trees
		tree.V > 0 AND -- I found some that are 0! That's because of culls set to 100%
		tree.INCREF IN ('15', '16', '17', '18', '19') -- Campaign from 2020 to 2024

	ORDER BY
		plot_id, tree_id;"), DT = TRUE)

disconnect_db(db)

## Get species names from database
species_group = readRDS(paste0(path_data, "groups_NFI.rds"))

## Modify few names, removing subsp. and var.
species_group[stri_detect(str = speciesName_sci, regex = "var."),
	speciesName_sci := stri_sub(str = speciesName_sci,
		to = stri_locate(str = speciesName_sci, regex = "var.")[, "start"] - 2)]

species_group[stri_detect(str = speciesName_sci, regex = "subsp."),
	speciesName_sci := stri_sub(str = speciesName_sci,
		to = stri_locate(str = speciesName_sci, regex = "subsp.")[, "start"] - 2)]

## Manually set for 29AF (broadleaf sp.) and 68CE (conifer sp.)
setkey(species_group, sp_code)
species_group["29AF", speciesName_sci := "Broadleaf"]
species_group["68CE", speciesName_sci := "Conifer"]

## Manually set for 53CA (Pinus nigra var. calabrica) and 53CO (Pinus nigra var. corsicana)
species_group["53CA", speciesName_sci := "Pinus laricio"]
species_group["53CO", speciesName_sci := "Pinus laricio"]

## Manually set for 58 (Pinus mugo subsp. uncinata)
species_group["58", speciesName_sci := "Pinus uncinata"]

## Remove French name column
species_group[, nom_espar := NULL]

## Merge with nfi_data
nfi_data = merge.data.table(nfi_data, species_group[, .(speciesName_sci, sp_code, group)], by = "sp_code")

## Load params
params_dt_species = readRDS(paste0(path_output, "avg_params.rds"))
ls_species = params_dt_species[, speciesName_sci]
params_dt_group = readRDS(paste0(path_output, "avg_params_full.rds"))[!ls_species, .(group, c, j, k, m, n, s, tau)] |>
	unique() |>
	setkey(group)

#### Compute total above-ground volume...
nfi_data[, U_V0_ALAMOD := NA_real_]

## ...for parametrised species
setkey(nfi_data, speciesName_sci)

nfi_data[, U_V0_ALAMOD := pred_vol(bole_volume, params_dt_species[.(.BY$speciesName_sci)]),
	by = speciesName_sci]

## ...for groups AND NOT parametrised species
nfi_data[!(ls_species), U_V0_ALAMOD := pred_vol(bole_volume, params_dt_group[.(.BY$group)]),
	by = group]

## ...for Pinus uncinata which species-specific model is unreliable (lack of data in the tail), apply group
nfi_data["Pinus uncinata", U_V0_ALAMOD := pred_vol(bole_volume, params_dt_group[.(.BY$group)]),
	by = group]

## ...for the group A2 that does not exist
nfi_data[is.na(U_V0_ALAMOD), unique(group)]
nfi_data[is.na(U_V0_ALAMOD), unique(speciesName_sci)] # Should be conifers only, Cupressus and Juniperus
nfi_data[is.na(U_V0_ALAMOD), U_V0_ALAMOD := pred_vol(bole_volume, params_dt_group["conifer-generic"])]

if (nfi_data[, any(is.na(U_V0_ALAMOD))])
	stop("There should not be any NA left")

saveRDS(nfi_data[, .(npp = plot_id, a = tree_id, U_V0_ALAMOD)], "U_V0_ALAMOD-base-test.rds")
