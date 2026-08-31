#### Aim of script: Make directories, check data
## Packages
library(data.table)

## Check/make directories
path_data = "./data/"
path_output = "./results/"

if (!dir.exists(path_data))
	stop(paste0("The folder <", path_data, ">is missing. Please check your working directory"))

if (!dir.exists(path_output))
	dir.create(path_output)

## Check data
if (!file.exists("tree_dt.rds"))
{
	# Loading Emerge data
	inra = readRDS(paste0(path_data, "inra.rds"))
	inra = inra[, .(speciesName_sci, tree_id = unique_id, plot_id, fct_type, year, circumference_m, height,
		taper_height = taper_height_flo, bole_volume_m3, total_volume_m3)]
	setkey(inra, speciesName_sci)
	
	# Loading Swiss data
	swiss = readRDS(paste0(mnt_point, "data/switzerland.rds"))
	swiss = unique(swiss[, .(speciesName_sci, tree_id, plot_id, fct_type, year, circumference_m, height,
		taper_height = hdec, bole_volume_m3, total_volume_m3)])
	setkey(swiss, speciesName_sci)
	
	## Modern Emerge data
	emerge = readRDS(paste0(mnt_point, "data/emerge_2009-2010.rds"))
	emerge[, year := as.integer(stringi::stri_replace_all(str = dataset, replacement = "", regex = "emerge_"))]
	emerge = emerge[, .(speciesName_sci, tree_id, fct_type, year, circumference_m, height, taper_height,
		bole_volume_m3 = bole_volume_conic_m3, total_volume_m3)]
	
	## Bind everything
	tree_dt = rbindlist(list(inra = inra, swiss = swiss, emerge = emerge), idcol = "origin", fill = TRUE)

}