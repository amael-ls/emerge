#### Aim of script: Compute a non-convex hull in the c-h space and check overlaps NFI vs fitting data
## Comments
# I expect the bole volumes of both datasets (NFI and training_dt) to overlap quite well,
#	while it might not be the case for the circumference -- height space range
#

#### Load packages
library(data.table)
library(terra)
library(ks)

#### Load common data and tool functions
## Tool functions
source("./tool_functions.R")

## Function to build species-specific concave hull
sp_hull = function(dt, param = 0.2, type = "concave_ratio", trim_nfi = TRUE, level = 0.95, allowHoles = FALSE)
{
	if (param < 0 || param > 1)
		stop("param must be between 0 and 1")

	if (type != "concave_ratio" && type != "concave_length")
		stop("I want a concave hull!")

	pts_nfi = as.matrix(unique(dt[origin == "nfi", .(circumference_m, height)]))
	if (trim_nfi) {
		kd_nfi = kde(x = pts_nfi)
		dens_nfi = predict(kd_nfi, x = pts_nfi)
		thresh_nfi = contourLevels(kd_nfi, cont = 100 * level)
		pts_nfi = pts_nfi[dens_nfi >= thresh_nfi, , drop = FALSE]
	}

	# Hull NFI
	dt_vect_nfi = vect(dt[origin == "nfi", .(circumference_m, height)])
	dt_vect_nfi_trimmed = vect(pts_nfi, type = "points")
	nfi = hull(dt_vect_nfi, type = type, param = param, allowHoles = allowHoles)
	nfi_trimmed = hull(dt_vect_nfi_trimmed, type = type, param = param, allowHoles = allowHoles)

	# Hull training
	dt_vect_train = vect(unique(dt[origin == "training", .(circumference_m, height)]))
	training = hull(dt_vect_train, type = type, param = param, allowHoles = allowHoles)

	# Compute overlap
	overlap_full = expanse(intersect(nfi, training), unit = "m")
	overlap_trimmed = expanse(intersect(nfi_trimmed, training), unit = "m")
	nfi_area = expanse(nfi, unit = "m")
	train_area = expanse(training, unit = "m")

	overlap_2d = data.table(target = c("training", "NFI"),
		area_target = c(train_area, nfi_area),
		area_intersect_full = overlap_full,
		area_intersect_trimmed = overlap_trimmed,
		overlap_percent_full = 100 * c(overlap_full/train_area, overlap_full/nfi_area),
		overlap_percent_trimmed = 100 * c(overlap_trimmed/train_area, overlap_trimmed/nfi_area)
	)

	return (list(nfi_hull = nfi, nfi_pt = dt_vect_nfi,
		training_hull = training, training_pt = dt_vect_train, overlap = overlap_2d))
}

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
		tree.SIMPLIF = '0' AND -- Trees with height measured
		tree.INCREF IN ('15', '16', '17', '18', '19') -- Campaign from 2020 to 2024

	ORDER BY
		species_code_nfi;"), DT = TRUE)

inventR::disconnect_db(db)

## Merge training dataset and NFI data
# Add sp_code
tree_dt = merge.data.table(tree_dt, nfi_codes, by = "speciesName_sci")
tree_dt = tree_dt[, .SD, .SDcols = colnames(nfi_data)]

# Row binding
tree_dt = rbindlist(l = list(training = tree_dt, nfi = nfi_data), idcol = "origin")[
	!is.na(circumference_m)][!is.na(height)]
setkey(tree_dt, species_code_nfi)

rm(nfi_data)

#### Overlap circumference -- height space
## Build hulls, compute overlaps
ls_hulls = vector(mode = "list", length = length(ls_species))
names(ls_hulls) = ls_species

temp = vector(mode = "list", length = length(ls_species))
names(temp) = ls_species

tree_dt[, .N, by = species_code_nfi][order(N)]
sp = "Pinus strobus"

for (sp in ls_species)
{
	print(paste("Running", sp))

	param = 0.2 # Seems ok visually
	if (sp %in% c("Pinus laricio", "Pinus strobus"))
		param = 0.3 # Adjusted visually

	ls_hulls[[sp]] = sp_hull(dt = tree_dt[nfi_codes[.(sp), species_code_nfi]],
		param = param, level = 0.975, trim_nfi = FALSE)

	# Extract overlap
	temp[[sp]] = ls_hulls[[sp]][["overlap"]]

	# Plot
	ext_nfi = ext(ls_hulls[[sp]][["nfi_hull"]])
	ext_train = ext(ls_hulls[[sp]][["training_hull"]])

	xlim = c(min(ext_nfi[1], ext_train[1]), max(ext_nfi[2], ext_train[2]))
	ylim = c(min(ext_nfi[3], ext_train[3]), max(ext_nfi[4], ext_train[4]))

	plot(ls_hulls[[sp]][["nfi_hull"]], col = "#0F7BA255", main = sp, xlim = xlim, ylim = ylim,
		axes = FALSE, xlab = "Circumference", ylab = "Height", clip = TRUE)
	plot(ls_hulls[[sp]][["training_hull"]], col = "#FAB25588", add = TRUE)

	# points(ls_hulls[[sp]][["nfi_pt"]], pch = 19, cex = 0.5, col = "#5A5A5AAA")

	axis(1)
	axis(2, las = 1)
}

overlap = rbindlist(temp, idcol = "speciesName_sci")

#### Overlap bole volume space
## Span of measurements
range_dt = tree_dt[, .(min_bole_v = min(bole_volume_m3), max_bole_v = max(bole_volume_m3)),
	by = .(species_code_nfi, origin)] |>
	dcast(species_code_nfi ~ origin, value.var = c("min_bole_v", "max_bole_v")) # NA for 01 as not in NFI 2020--2024

## Overlap (i.e., the smallest max - the largest min)
range_dt[, overlap := min(max_bole_v_training, max_bole_v_nfi) -
	max(min_bole_v_training, min_bole_v_nfi), by = species_code_nfi]
range_dt[, overlap_percent := 100*overlap/(max_bole_v_nfi - min_bole_v_nfi)]

#### Merge 2D space and bole volume overlaps
## Keep only NFI target and columns of interest
overlap = overlap[target == "NFI", .(speciesName_sci, circum_height = overlap_percent_full)]

## Merging
range_dt = merge.data.table(range_dt, nfi_codes, by = "species_code_nfi")
range_dt = range_dt[, .(speciesName_sci, bole_volume = overlap_percent)]

overlap = merge.data.table(overlap, range_dt)

filename = paste0(path_pgfplotstable, "overlap.csv")
if (!file.exists(filename))
	fwrite(overlap, filename)
