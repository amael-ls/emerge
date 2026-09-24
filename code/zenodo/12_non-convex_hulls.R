#### Aim of script: Compute a non-convex hull in the c-h space and check overlaps NFI vs fitting data
## Comments
# I expect the bole volumes of both datasets (NFI and training_dt) to overlap quite well,
#	while it might not be the case for the circumference -- height space range
#

#### Load packages
library(data.table)
library(concaveman)
library(terra)
# library(sf)

#### Load common data and tool functions
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

tree_dt = readRDS(paste0(path_data, "tree_dt_14species.rds"))
ls_species = tree_dt[, unique(speciesName_sci)]

## French NFI data query, 2020 -- 2024 campaigns
nfi_codes = readRDS(paste0(path_data, "nfi_codes.rds"))
sp_code_nfi = paste0("'", paste(nfi_codes[, species_code_nfi], collapse = "','"), "'")

db = inventR::connect_db(base = "exploitation", serveur = "inv-exp")

nfi_data = inventR::exec_req(conn = db, req = paste0("
	SELECT
		tree.ESS AS species_code_nfi, tree.C13 AS circumference_m,
		tree.HTOT AS height, tree.V AS bole_volume_m3

	FROM
		inv_exp_nm.g3arbre AS tree

	WHERE
		tree.ESS IN (", sp_code_nfi, ") AND -- Selected tree species
		tree.VEGET = '0' AND -- Living standing trees
		tree.V > 0 AND -- I found some that are 0! That's because of culls set to 100%
		tree.INCREF IN ('15', '16', '17', '18', '19') -- Campaign from 2020 to 2024

	ORDER BY
		species_code_nfi;"), DT = TRUE)

inventR::disconnect_db(db)

## Merge training dataset and NFI data
# Add sp_code
tree_dt = merge.data.table(tree_dt, nfi_codes, by = "speciesName_sci")
tree_dt = tree_dt[, .SD, .SDcols = colnames(nfi_data)]

# Row binding
tree_dt = rbindlist(l = list(training = tree_dt, nfi = nfi_data), idcol = "origin")
rm(nfi_data)
