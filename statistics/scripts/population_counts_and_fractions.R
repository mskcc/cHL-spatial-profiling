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
            --cell_dive_id                 CellDive_ID of sample to process
            --statistics_conditions_index  XLSX file containing all conditions to be analyzed,
                                           previously indexed
            --statistics_conditions_file   XLSX file containing all conditions to be analyzed (pre-indexing)
            --meta_dir                     path to meta files in XLSX format, required IF 
                                           meta_data_file is NULL
            --fov_metrics_dir              directory where all FOV counts and fractions should be saved
            --annotation_config_file       YAML file containing cell annotation settings (see docs)

          [OPTIONAL]
            --cell_region         calculate metrics only for cells that fall within a certain neighborhood;
                                  a neighborhood is defined as the 30-micron radius around a cell of a certain
                                  type or subtype (e.g., 'HRS' or 'T_All')
            --force_recalc        logical, default = FALSE; when TRUE, recount and recalc even if 
                                  files already contain data
            --calc_unit           unit of one count and fraction; default = FOV_ID
            --manifest            YAML file containing one or more parameter; NOTE: arguments on command
                                  line override manifest arguments!!!
            --meta_files          comma-delimited list of meta files; can be supplied instead of meta_dir
            --neighbor_dir        directory containing cell-to-cell distance files for each sample
            --neighborhood_radius radius in microns defining a cell neighborhood; default = 30
            --number_threads      number of threads to use for parallel processes
        \n"
    )
}

## set up required args & defaults 
minReq   <- list(c("meta_dir","meta_files"), 
                 "cell_data_dir",
                 "cell_dive_id", 
                 "statistics_conditions_file",
                 "statistics_conditions_index",
                 "fov_metrics_dir")

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), paramDflts, minReq, usage)
} else {
    args <- read_yaml("input/config/study_config.yaml")
    args$cell_dive_id <- "HL_10"
    args <- processCMD(args, paramDflts, minReq, usage)
}
logParams(args, names(args))

## load study data
stDat <- loadStudyData(args,
                       analyses = T,
                       annotatedCells = T,
             #          cellsInTumorNeighborhood = T,
                       tmeCellStatus = F)

populations <- stDat$conds %>% select(Subpopulation, Population) %>% unlist %>% unique 
populations <- populations[!is.na(populations)]
cdids       <- getSampleIDs(metaFiles = stDat$metaFiles, cdid = args$cell_dive_id)
threads     <- min(args$number_threads, length(cdids))

cl <- startCluster(threads, exports = c("args", "stDat", "populations"))
allFracs <- parLapply(cl, seq(cdids), function(x){

    cdid  <- cdids[x]
    afile <- dir(args$cell_data_dir, full.names = T, pattern = paste0("^", cdid, "_"))
    cfile <- file.path(args$fov_metrics_dir, paste0(cdid, "_population_counts_per_", args$calc_unit, ".rda"))
    ffile <- file.path(args$fov_metrics_dir, paste0(cdid, "_population_fractions_per_", args$calc_unit, ".rda"))

    tryCatch({
        #annDat <- loadAnnotatedCells(stDat$sampAnn, annFile = afile) 
        annDat <- stDat$annCells

        if(!is.null(args$cell_region)){
            annDat <- annDat %>% 
                      filter(!!!cellStateFilter(args$cell_region))

            #if(args$cell_region != "HRS_nbhd"){
            #    log_error(paste0("Unrecognized cell_region: '", args$cell_region,
            #                     "'. Currently only supporting 'HRS_nbhd'."))
            #    return(NULL)
            #}
            #if(is.null(stDat$tumorNbhdCells) | length(stDat$tumorNbhdCells) == 0){
            #    log_warn(paste0("Zero cells are in a tumor neighborhood."))
            #    return(NULL)
            #}
            #log_info(paste0("Filtering data for cells within ", 
            #                args$neighborhood_radius, " microns of an HRS cell"))
            #annDat <- annDat %>% filter(UUID %in% stDat$tumorNbhdCells)
        }

        missing_markers <- stDat$sampAnn %>% 
                           filter(CellDive_ID == cdid) %>%
                           pull(Final_missing_markers_functional) %>%
                           strsplit(',') %>%
                           unlist() %>% unique

        ### REMOVE POPULATIONS THAT INCLUDE MISSING MARKERS IN THIS SAMPLE
        na_pops <- na_population(populations, missing_markers) 
        pops <- populations[!na_pops]

        annDat %>%
        getPopulationCounts(pops, 
                            stDat$markers, 
                            calcUnit = args$calc_unit,
                            numThreads = threads, 
                            outFile = cfile, 
                            forceRecalc = args$force_recalc) %>%
        getPopulationFractions(stDat$conds %>% filter(grepl("fraction", CalcType)), 
                               calcUnit = args$calc_unit, 
                               outFile = ffile, 
                               forceRecalc = args$force_recalc)
      }, error = function(e){
        print(e)
        log_error("Counts and fractions FAILED for ", cdid)
    })

})
stopCluster(cl)


