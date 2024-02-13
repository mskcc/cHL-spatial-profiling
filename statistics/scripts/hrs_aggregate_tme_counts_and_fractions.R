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
            --meta_dir                     path to meta files in XLSX format, required IF 
                                           meta_data_file is NULL
            --microenvironment_dir         path to RDA file containing microenvironment assignments
            --microenvironment_metrics_dir directory where all ME counts and fractions should be saved
            --neighbor_cell_class          Category, Cell_type or Subtype of neighbor cells
                                           of interest
            --statistics_conditions_index  XLSX file containing all conditions to be analyzed,
                                           previously indexed
            --statistics_conditions_file   XLSX file containing all conditions to be analyzed (pre-indexing)


          [OPTIONAL]
            --cell_dive_id                a single Cell Dive ID for which cell conditions should be counted
                                          and fractions calculated
            --calc_unit                   unit of one count and fraction; default = FOV_ID
            --force_recalc                logical, default = FALSE; when TRUE, recount and recalc even if 
                                          files already contain data
            --manifest                    YAML file containing one or more parameter; NOTE: 
                                          arguments on command line override manifest arguments!!!
            --meta_files                  comma-delimited list of meta files; can be 
                                          supplied instead of meta_dir
            --min_microenv_cells                                                
            --number_threads              number of threads to use for parallel processes
            --neighbor_cell_marker_combo  comma-delimited marker combination that together 
                                          with neighbor_cell_class define neighborhood of 
                                          interest
            --postive_cutoff 
        \n"
    )
}

## set up required args & defaults 
minReq   <- list(c("meta_dir","meta_files"), 
                 "cell_data_dir", 
                 "statistics_conditions_file",
                 "statistics_conditions_index",
                 "microenvironment_dir",
                 "microenvironment_metrics_dir",
                 "neighbor_cell_class")

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), paramDflts, minReq, usage)
} else {
    args <- read_yaml("input/config/study_config.yaml")
    args$cell_dive_id <- "HL_20"
    args$neighbor_cell_class <-  "HRS"
    args$neighbor_cell_marker_combo <- "aggregated"
    args <- processCMD(args, paramDflts, minReq, usage)
}

logParams(args, sort(names(args)))

log_info("Loading study data...")
## load study data
stDat <- loadStudyData(args,
                       analyses = T,
                       annotatedCells = F)
 
populations <- stDat$conds %>% select(Subpopulation, Population) %>% unlist %>% unique
populations <- populations[!is.na(populations)]
cdids       <- getSampleIDs(metaFiles = stDat$metaFiles, cdid = args$cell_dive_id)
threads     <- min(args$number_threads, length(cdids))
me          <- paste(c(args$neighbor_cell_class, args$neighbor_cell_marker_combo), collapse = ",")
cu          <- paste0(args$calc_unit, "___", me)

cl <- startCluster(threads, exports = c("args", "stDat", "populations", "me", "cu"))
allFracs <- parLapply(cl, seq(cdids), function(x){

    cdid  <- cdids[x]
    afile <- dir(args$cell_data_dir, full.names = T, pattern = annotated_cell_file_pattern(cdid = cdid))
    mfile <- dir(args$microenvironment_dir, full.names = T, 
                 pattern = microenv_file_pattern(cdid = cdid, 
                                                 neighbor = args$neighbor_cell_class,
                                                 state = args$neighbor_cell_marker_combo,
                                                 frac_pos = args$positive_cutoff, 
                                                 min_cells = args$min_microenv_cells))
 
    mdat <- readRDS(mfile) %>% 
            ungroup() %>%
            select(UUID = C.UUID, FOV_number, `HRS,aggregated_microEnv`)

    cfile <- file.path(args$microenvironment_metrics_dir, paste0(cdid, "_population_counts_per__FOV_ID___HRS__aggregated.rda"))
    ffile <- file.path(args$microenvironment_metrics_dir, paste0(cdid, "_population_fractions_per__FOV_ID___HRS__aggregated.rda"))

    missing_markers <- stDat$sampAnn %>%                                    
                       filter(CellDive_ID == cdid) %>%                      
                       pull(Final_missing_markers_functional) %>%           
                       strsplit(',') %>%                                    
                       unlist() %>% unique                                  
                                                                                
    ### REMOVE POPULATIONS THAT INCLUDE MISSING MARKERS                     
    na_pops <- na_population(populations, missing_markers)                  
    pops <- populations[!na_pops]                                           
     
    dat <- loadAnnotatedCells(stDat$sampAnn, annFile = afile) %>%
           left_join(mdat) %>%
           rename(`cell_HRS,aggregated` = `HRS,aggregated_microEnv`) %>%
           unite(!!as.name(cu), args$calc_unit, paste0("cell_",me), sep="___", remove = F)

    tryCatch({
        dat %>%
        getPopulationCounts(pops, stDat$markers, 
                            calcUnit = cu,
                            numThreads = threads, 
                            outFile = cfile, 
                            forceRecalc = args$force_recalc) %>%
        getPopulationFractions(stDat$conds %>% filter(grepl("fraction", CalcType)), 
                               calcUnit = cu, 
                               outFile = ffile, 
                               forceRecalc = args$force_recalc)
      }, error = function(e){
        print(e)
        log_error("Counts and fractions FAILED for ", cdid)
    })

})
stopCluster(cl)


