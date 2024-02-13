suppressMessages(library(R.utils))
srcDir <- dirname(dirname(commandArgs(asValue=TRUE)$file))
source(file.path(srcDir, "source_all.R"))

## get git revision
#gitrev <- system(paste("git rev-parse HEAD", srcDir), intern=TRUE)
# can't get from /juno to /home?? 
# not important right now


log_threshold(DEBUG)

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript configure_study.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --study_name           character string identifying the study; must match prefix
                                   included in all meta data files (e.g., 'HodgkinsLymphoma')

          [OPTIONAL]
            --project_data_root     path to project directory containing the standard subdirectory
                                    structure: child folders are named according to date of last change;
                                    dated subfolders contain all input data folders including 
                                    drift, haloCoor, haloCSV, haloImages, haloRDA, info, meta (required unless
                                    passing individual data directories)
            --source_meta_dir       full path to project meta data directory; required when --project_data_root 
                                    not provided
            --source_halo_data_dir  full path to final Halo RDA files; required when --project_data_root 
                                    not provided
            --source_halo_csv_dir   full path to Halo CSV files; required when --project_data_root 
                                    not provided
            --source_halo_image_dir full path to Halo image files; required when --project_data_root 
                                    not provided
            --source_halo_coord_dir full path to Halo coordinate (XML) files; required when --project_data_root 
                                    not provided
            --source_drift_summary_dir full path to drift directories drift_mask, drift_summary, drift_pdf; 
                                    required when --project_data_root not provided
            --source_neighbor_dir   full path to cell neighbor RDA files; required when --project_data_root
                                    not proided
            --study_config          YAML file containing custom config
            --plot_config_file      YAML file containing plot colors to be applied to all analyses
        \n"
    )
}

### SET DEFAULTS:
defaults <- list(
                 meta_data_file = "flat_meta_data.xlsx",

                 control_marker = "DAPI",
                 pad = 20,
                 drift_threshold = 0.1,
                 max_g = 5.0,
                 number_threads = 6.0,
                 image_identifier_string = 'mono_dapi_reg_pyr16',
                 neighborhood_radius = 30,
                 min_microenv_cells = 1,
                 positive_cutoff = 0.95,  ## positive microenvironment
                 flex_neg_req = TRUE, ## allow removal of required negatives for missing markers
                 debug = "yes",

                 plot_config_file = "input/config/global_plot_colors.yaml",
                 detail_figure_config_file = "input/config/stats_detail_figure_config.yaml",
                 statistics_config_file = "input/config/statistics_config.yaml",
                 annotation_config_file = "input/config/annotation_config.yaml",
                 condition_detail_config = "input/config/condition_detail_config.yaml",

                 statistics_conditions_file = "input/config/cell_state_conditions.xlsx",
                 detail_figure_conditions_file = NULL,

                 config_dir = "input/config",
                 meta_dir = "input/meta",
                 halo_csv_dir = "input/halo_csv",
                 halo_image_dir = "input/halo_images",
                 qc_dir = "qc",
                 log_dir = "logs",

                 study_figure_dir = "results/figures",
                 drift_summary_dir = "preprocessing/drift/drift_summary",
                 drift_mask_dir = "preprocessing/drift/drift_mask",
                 drift_pdf_dir = "preprocessing/drift/drift_pdf",
                 data_dir = "preprocessing/exclusions",
                 cell_data_dir = "preprocessing/annotated",
                 neighbor_dir = "preprocessing/neighbors",
                 microenvironment_dir = "preprocessing/microenvironments",
                 statistics_conditions_index = "results/statistics/cell_state_conditions_index.xlsx",
                 metrics_dir = "processed/metrics",
                 fov_metrics_dir = "processed/metrics/fovs",
                 fov_area_dir = "processed/metrics/fovs",
                 microenvironment_metrics_dir = "processed/metrics/fovs/microenvironments",

                 counts_dir = "results/counts",
                 cell_type_counts_file = "results/counts/marker_combination_cell_type_counts.xlsx",
                 statistics_dir = "results/statistics",
                 statistics_general_tables_dir = "results/statistics/general/tables",
                 statistics_spatial_tables_dir = "results/statistics/spatial/tables",
                 statistics_overview_dir = "results/statistics/figures/overviews",
                 statistics_general_detail_dir = "results/statistics/general/figures/detail",
                 statistics_spatial_detail_dir = "results/statistics/spatial/figures/detail",
                 statistics_misc_dir = "results/statistics/figures/misc"

            )

minReq <- list("study_name", 
               c("project_data_root", "source_meta_dir", "source_drift_summary_dir")) ## really it's project_data_root or all source_*_dir files, 
                                                          ## but processCMD doesn't support that kind of thing at the moment
used <- unlist(c(minReq,  names(defaults), 
               "source_meta_dir", "source_halo_image_dir", "source_halo_csv_dir", "source_halo_data_dir",
               "source_halo_coord_dir", "source_drift_summary_dir", "source_neighbor_dir")) %>% unique

if(interactive()){
    args <- list(project_data_root = "/juno/res/bic/shared/Multiomyx/Projects/Hodgkins",
                 study_name = "HodgkinLymphoma")
    cfg <- processCMD(args, unlist(defaults), minReq, usage)
} else {
    cfg <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
}
cfg <- cfg[!names(cfg) %in% c("", "no-echo", "no-restore", "file", "args")]


## validate command line input
srcDirs <- paste0("source_", c("meta_dir", 
                               "halo_csv_dir", 
                               "halo_image_dir",
                               "halo_coord_dir", 
                               "halo_data_dir",
                               "drift_summary_dir",
                               "neighbor_dir"))
given <- cfg[names(cfg) %in% srcDirs]

if(is.null(cfg$project_data_root) && length(given) < length(srcDirs)){
    stop("Must provide either project_data_root or ALL of the following: ", 
         paste(srcDirs, collapse = ","))
} else if(!is.null(cfg$project_data_root) && length(given) > 0 && length(given) < length(srcDirs)) {
    stop("Must provide EITHER project_data_root or ALL of the following: ", 
         paste(srcDirs, collapse = ","), ". Please do not provide both.")
}


if(!is.null(cfg$project_data_root)){
    current <- getProjectDataDir(cfg$project_data_root)
    cfg$source_halo_data_dir     <- file.path(current, "haloRDA/preNorm/exclusionsMarked")
    cfg$source_meta_dir          <- file.path(current, "meta")
    cfg$source_halo_image_dir    <- file.path(current, "haloImages")
    cfg$source_halo_csv_dir      <- file.path(current, "haloCSV")
    cfg$source_halo_coord_dir    <- file.path(current, "haloCoor")
    cfg$source_drift_summary_dir <- file.path(current, "drift/drift_summary")
    cfg$source_neighbor_dir      <- file.path(current, "haloRDA/neighbors/exclusionsRemoved")
}

if(is.null(cfg$cell_states_file)){
    cfg$cell_states_file <- file.path(cfg$meta_dir, paste0(cfg$study_name, "_CellStates.xlsx"))
}

logParams(cfg, used)

cat("\nDirectory structure will be created according to the configuration above.  Type 'stop' to write config to file for editing. Press return to proceed: ")

proceed <- readLines("stdin", n=1)
if(tolower(proceed) == "stop"){
    cat("
         Writing template config to file
        
               TEMPLATE_CONFIG.yaml

         Please modify as necessary, rename and rerun configuration with '--study_config [config_file]'\n\n")
    write_yaml(cfg, "TEMPLATE_CONFIG.yaml")
    q()
}

if(cfg$debug == FALSE){ log_threshold(INFO) }

cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
log_info("Creating directories...")
dir_setup(cfg)

cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
cfg_file <- file.path(cfg$config_dir, "study_config.yaml")
log_info(paste0("Saving all study configuration to: ", cfg_file))
write_yaml(cfg, cfg_file)



cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
log_info("Validating study config...")
tmp <- file.symlink(dir(cfg$source_meta_dir, full.names = T), cfg$meta_dir)  ## this must be done before validation
scValid <- validateStudyConfig(cfg)


cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
log_info("Validating all meta data...")
validateInput(cfg)

cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
log_info("Compiling all meta data...")
meta <- loadStudyAnnotations(metaFiles = getCurrentMetaFiles(metaFiles = cfg$meta_files, metaDir = cfg$meta_dir)) 
if(!is.null(cfg$meta_data_file)){
    write.xlsx(meta, gsub(".rda$", ".xlsx", cfg$meta_data_file), check.names = F)
    saveRDS(meta, gsub(".xlsx$", ".rda", cfg$meta_data_file))
}

cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
log_info("Linking all source files locally...")
# link images, CSV files and RDA files
inputPaths <- list(drift  = list(src = cfg$source_drift_summary_dir, local = cfg$drift_summary_dir, pattern = ".txt"),
                   csv    = list(src = cfg$source_halo_csv_dir, local = cfg$halo_csv_dir, pattern = ".csv"),
                   #coord  = list(src = cfg$source_coord_dir, local = cfg$halo_coord_dir),  ## maybe needed later
                   data   = list(src = cfg$source_halo_data_dir, local = cfg$data_dir, pattern = ".rda"),
                   neighbors = list(src = cfg$source_neighbor_dir, local = cfg$neighbor_dir, pattern = ".rda"))

log_info("Linking images...")
## images are a little complicated to find and link, so they need to be done separately;
## also they are not really needed after drift is done, which is currently upstream
## of this script so don't bother linking them; 

#linkCount <- link_halo_images(meta$flat, sourceDir = cfg$source_halo_image_dir, 
#                              localDir = cfg$halo_image_dir, 
#                              imageIdentifier = cfg$image_identifier_string)
#log_info("Linked [", linkCount, "] image files.")

log_info("Linking all other input files...\n") 
tmp <- link_existing_files(inputPaths)


## linking pipeline scripts
cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
log_info("Linking launch script(s) in working directory.")
tmp <- file.symlink(file.path(srcDir, "pipeline/run.sh"), getwd())
tmp <- file.symlink(file.path(srcDir, "pipeline/run_drift.sh"), getwd())  ## this soon will not be necessary once it is combined with run.sh

log_info("Done!")

