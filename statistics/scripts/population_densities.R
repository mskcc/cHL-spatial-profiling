suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))
log_threshold(DEBUG)

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript fov_counts_and_fractions.R 

          [REQUIRED (may be defined on command line OR in manifest file)]
            --cell_data_dir                path to annotated cells RDA files that each 
                                           contain single table where rows are cells
                                           and columns are all data for that single 
                                           cell (one file per sample)
            --statistics_conditions_index  XLSX file containing all conditions to be analyzed,
                                           previously indexed
            --statistics_conditions_file   XLSX file containing all conditions to be analyzed, unindexed
            --meta_dir                     path to meta files in XLSX format, required IF 
                                           meta_data_file is NULL
            --fov_metrics_dir              root directory where all FOV metrics are saved 
            --annotation_config_file       YAML file containing cell annotation settings (see docs)

          [OPTIONAL]
            --force_recalc        logical, default = FALSE; when TRUE, recount and recalc even if 
                                  files already contain data
            --area_unit           ['px'|'um'|mm'], default = 'mm'; unit in which area should be returned
            --cell_dive_id        a single Cell Dive ID for which cell conditions should be counted
                                  and fractions calculated
            --calc_unit           unit of one count and fraction; default = FOV_ID
            --manifest            YAML file containing one or more parameter; NOTE: arguments on command
                                  line override manifest arguments!!!
            --meta_files          comma-delimited list of meta files; can be supplied instead of meta_dir
            --number_threads      number of threads to use for parallel processes
        \n"
    )
}

## set up required args & defaults 
minReq   <- list(c("meta_dir","meta_files"), 
                 "cell_data_dir", 
                 "statistics_conditions_file",
                 "statistics_conditions_index",
                 "fov_metrics_dir")

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), paramDflts, minReq, usage)
} else {
    args <- read_yaml("input/config/study_config.yaml")
    args$cell_dive_id <- "HL_8"
    args <- processCMD(args, paramDflts, minReq, usage)
}
logParams(args, names(args))


metaFiles <- getCurrentMetaFiles(metaDir = args$meta_dir)
sampAnn <- loadStudyAnnotations(metaFiles = metaFiles)$flat
markers <- getAllMarkers(markerFile = getMetaFile("Markers", metaFiles)) 
cdids   <- getSampleIDs(metaFiles = metaFiles, cdid = args$cell_dive_id)
conds   <- getConditionsIndex(args$statistics_conditions_index,
                              args$statistics_conditions_file,
                              read_yaml(args$annotation_config_file)$arrange_annotation) %>%
           filter(grepl("dens", CalcType))
populations <- conds %>% select(Subpopulation, Population) %>% unlist %>% unique
populations <- populations[!is.na(populations)]

threads     <- min(args$number_threads, length(cdids))
cl <- startCluster(threads, exports = c("args", "sampAnn", "conds", "markers", "populations"))
allDens <- parLapply(cl, seq(cdids), function(x){

    cdid  <- cdids[x]
    afile <- dir(args$cell_data_dir, full.names = T, pattern = paste0("^", cdid, "_"))
    rfile <- file.path(args$fov_metrics_dir, paste0(cdid, "_total_area_per_", args$calc_unit ,".rda"))
    cfile <- file.path(args$fov_metrics_dir, paste0(cdid, "_population_counts_per_", args$calc_unit, ".rda"))
    dfile <- file.path(args$fov_metrics_dir, paste0(cdid, "_population_densities_per_", args$calc_unit, ".rda"))

    tryCatch({
        getPopulationCounts(afile, populations, markers,
                            calcUnit = args$calc_unit,
                            numThreads = threads,
                            outFile = cfile,
                            forceRecalc = FALSE) %>%
        getPopulationDensities(readRDS(rfile) %>% left_join(sampAnn %>% select(CellDive_ID, FOV_ID, FOV_number)),
                               calcUnit = args$calc_unit,
                               outFile = dfile,
                               forceRecalc = args$force_recalc) %>%
        mutate(CellDive_ID = cdid) %>%
        select(CellDive_ID, FOV_ID, Population, Count, TotalCounts, TotalArea, Density)
      }, error = function(e){
        print(e)
        log_error("Counts and fractions FAILED for ", cdid)
    })

})
stopCluster(cl)

log_info("Done.")
