library(funr)
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_threshold(DEBUG)
log_debug(paste0("loading source files from: ",sdir))

#####################################
#####       GET USER INPUT      #####
#####################################

usage <- function(){

 cat("\nUsage:  Rscript report_statistics.R 
            
    [REQUIRED (may be defined on command line OR in manifest file)] 
      --annotation_config_file       YAML file describing how conditions are to be 
                                     arranged and indexed
      --cell_data_dir                path to RDA files, each containing a single 
                                     table where rows are cells and columns are 
                                     all data for that single cell
      --fov_area_dir                 path to RDA files, each containing table of FOVs and 
                                     total FOV area for all FOVs in a single sample 
      --fov_metrics_dir              root directory for FOV area, fractions and 
                                     densities
      --meta_dir                     path to meta files in XLSX format
      --statistics_conditions_file   XLSX file listing all cell states/conditions to 
                                     compare between two sample groups
      --statistics_config_file       YAML file containing stats config such as filters and
                                     calculations to run stats on (see docs for details)
      --statistics_conditions_index  XLSX file with pre-indexed cell states/conditions
      --statistics_tables_dir        output directory where XLSX files of results should 
                                     be written

    [OPTIONAL]
      --question             a single QuestionNumber from statistics_questions_file to run stats on
      --manifest             YAML file containing one or more parameter; NOTE: arguments 
                             on command line override manifest arguments!!!        
      --microenvironment_dir directory containing files of tumor microenvironment assignments by cell, 
                             required if one or more question to be run involves grouping or filtering
                             on these assignments
      --microenvironment_metrics_dir directory for FOV-level microenvironment fractions
                                     and densities                              
      --neighbor_dir         directory of cell to cell distances, specifically distances 
                             between tumor cells and immune cells
      --number_threads       number of threads to use for parallel processes

  \n"
 )

}

## names of required args
minReq <- list("annotation_config_file", 
               "cell_data_dir", "fov_area_dir",
               c("meta_dir","meta_files","meta_data_file"),
               "fov_metrics_dir",
               "statistics_conditions_file", 
               "statistics_conditions_index",
               "statistics_config_file",
               "statistics_tables_dir")

defaults <- list(number_threads = 1, calc_unit = "FOV_ID")

used <- c(unlist(minReq), "question", "manifest", "number_threads", 
          "neighbor_dir", "microenvironment_dir", "calc_unit")

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

#####################################
###   CONFIGURATION & DATA INIT   ###
#####################################

cfg <- resolveConfig(args, 
                     read_yaml(args$annotation_config_file), 
                     read_yaml(args$statistics_config_file))
 
logParams(cfg, used) 
mkdir(cfg$statistics_tables_dir)

log_debug("Loading study data...")
stDat <- loadStudyData(cfg, 
                       questions = T,
                       analyses = T, 
                       annotatedCells = F) 

log_debug("Gathering all SAMPLE annotation...")
allAnnot <- getAllDataAnnotation(dir(cfg$cell_data_dir, full.names=T), 
                                 stDat$sampAnn, 
                                 getAllSampleFilters(stDat$allQuestions), 
                                 calcUnit = cfg$calc_unit, 
                                 threads = cfg$number_threads,
                                 filter_stats_exclusions = TRUE)

######################################
###        ANSWER  QUESTIONS       ###
######################################
log_info(paste0("Running stats on ", cfg$calc_unit, " level..."))

for(q in names(stDat$allQuestions)){
    log_info(paste0("QUESTION: ",q))

    question    <- stDat$allQuestions[[q]]
    allAnalyses <- stDat$analysisList
    metDir      <- cfg$fov_metrics_dir

    ## for tumor neighborhood analyses, only run fractions analysis 
    if(any(grepl("Neighborhood", names(unlist(stDat$allQuestions[[q]]))))){  
        allAnalyses <- allAnalyses[names(allAnalyses) == "fractions"]
    }

    res <- compareSampleGroups(question, allAnalyses, allAnnot, 
                               cfg$results_filters[[cfg$use_filter]], 
                               metDir,
                               calcUnit = cfg$calc_unit,
                               filterID = cfg$use_filter)

    ### Save results
    outFile <- file.path(cfg$statistics_tables_dir, paste0(q,".xlsx"))
    log_debug("writing XLSX file...")
    writeStatsQuestionXLSX(res, outFile)
    log_info("done.")
}



