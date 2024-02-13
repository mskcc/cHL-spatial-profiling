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
            --microenv_pos_cutoff         minimum fraction of neighbor cells with 'neighbor_cell_marker_combo'
                                          needed in order for a center cell to have 'Pos.Env' (used to filter
                                          files in 'microenvironment_dir')
            --min_cell_count              minimum number of neighbor cells of type 'neighbor_cell_class'
                                          needed in order for a cell to have 'Pos.Env' (used to filter
                                          files in 'microenvironment_dir')      
            --calc_unit                   unit of one count and fraction; default = FOV_ID
            --force_recalc                logical, default = FALSE; when TRUE, recount and recalc even if 
                                          files already contain data
            --manifest                    YAML file containing one or more parameter; NOTE: 
                                          arguments on command line override manifest arguments!!!
            --meta_files                  comma-delimited list of meta files; can be 
                                          supplied instead of meta_dir
            --number_threads              number of threads to use for parallel processes
            --neighbor_cell_marker_combo  comma-delimited marker combination that together 
                                          with neighbor_cell_class define neighborhood of 
                                          interest
 
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
    #args$neighbor_cell_marker_combo <- "aggregated"
    args$neighbor_cell_marker_combo <- "PDL1"
    args <- processCMD(args, paramDflts, minReq, usage)
}

logParams(args, sort(names(args)))

log_info("Loading study data...")
## load study data
stDat <- loadStudyData(args,
                       analyses = T,
                       annotatedCells = F)
 
populations <- stDat$conds %>% select(Subpopulation, Population) %>% unlist %>% unique
cdids       <- getSampleIDs(metaFiles = stDat$metaFiles, cdid = args$cell_dive_id)
threads     <- min(args$number_threads, length(cdids))

me <- args$neighbor_cell_class
if(!gtools::invalid(args$neighbor_cell_marker_combo)){
    me <- paste0(args$neighbor_cell_class, ",", args$neighbor_cell_marker_combo)
}
cu <- paste0(args$calc_unit, "___", me)


cl <- startCluster(threads, exports = c("args", "stDat", "populations", "me", "cu"))
allFracs <- parLapply(cl, seq(cdids), function(x){

    cdid  <- cdids[x]
    afile <- dir(args$cell_data_dir, full.names = T, pattern = annotated_cell_file_pattern(cdid))

    ## output files
    cfile <- file.path(args$microenvironment_metrics_dir, paste0(cdid, "_population_counts_per__", cu, ".rda"))
    ffile <- file.path(args$microenvironment_metrics_dir, paste0(cdid, "_population_fractions_per__", cu, ".rda"))

### addTME here adds a columns for ALL TME's for which it finds data in the microenv dir; 
### need to change that but leave it for now
    dat <- loadAnnotatedCells(stDat$sampAnn, annFile = afile) %>%
           addTME(args$microenvironment_dir, 
                  filePattern = microenv_file_pattern(cdid = cdid, 
                                                      frac_pos = args$positive_cutoff, 
                                                      min_cells = args$min_microenv_cells)) %>%
           unite(!!as.name(cu), args$calc_unit, paste0("cell_",me), sep="___", remove = F)

    missing_markers <- dat %>%                                    
                       pull(Final_missing_markers_functional) %>%           
                       strsplit(',') %>%                                    
                       unlist() %>% unique                                  
                                                                                
    ### REMOVE POPULATIONS THAT INCLUDE MISSING MARKERS                     
    na_pops <- na_population(populations, missing_markers)                  
    pops <- populations[!na_pops]                                           

    tryCatch({
        dat %>%
        getPopulationCounts(pops, 
                            stDat$markers, 
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


