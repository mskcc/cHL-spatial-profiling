suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))

#####################################
#####        SOURCE CODE        #####
#####################################

usage <- function(){

    cat("\nUsage:  Rscript calculate_area.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)]
            --data_dir         directory containing all halo data RDA files for which area should
                               be calculated; alternatively, provide a vector of specific files to
                               --data_files argument
            --fov_area_dir     output directory for FOV area files (one per sample)
            --meta_dir         path to meta files in XLSX format, required IF meta_data_file is NULL           

          [OPTIONAL]
            --area_unit        ['px'|'um'|mm'], default = 'mm'; unit in which area should be returned
            --cell_dive_id     CellDive_ID of sample to be processed; to be used only if a 
                               single data file is provided
            --data_files       halo data RDA file; must be specified if data_dir is not
            --force_recalc     recalculate area even if output file already exists; default = FALSE
            --manifest         YAML file containing one or more parameter;
                               NOTE: arguments on command line override manifest arguments!!!
            --max_g            default = 5
            --meta_data_file   RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                               meta_files are NULL
            --number_threads   number of threads to use for parallel processes
        \n"
    )
}

### process command line
minReq <- list("fov_area_dir",
               c("data_files", "data_dir"),
               c("meta_dir","meta_files","meta_data_file"))
defaults <- list(max_g = 5, number_threads = 1, cell_dive_id = NULL, area_unit = 'mm',
                 force_recalc = FALSE, calc_unit = "FOV_ID")
args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(args, names(args))


### load all annotation data 
annot <- loadStudyAnnotations(metaFiles = getCurrentMetaFiles(metaDir = args$meta_dir),
                              metaDataFile = args$meta_data_file)$flat
dataFiles <- getFiles(path = args$data_dir, files = args$data_files) 
threads   <- min(args$number_threads, length(dataFiles))


###
### Calculate and save FOV areas
###
cl <- makeCluster(threads, type="FORK", outfile="")
clusterExport(cl, c("annot", "dataFiles", "args"), envir=environment())
x <- parLapply(cl, unique(annot$CellDive_ID), function(cdid){

        haloFile <- dataFiles[grepl(paste0("^", cdid, "_"), basename(dataFiles))]
        if(length(haloFile) == 0){ return(NULL) }
        log_info("Calculating area for FOVs in file ", haloFile)
        fAreas <- getAllSampleFOVAreas(haloFile, 
                                       unit = args$area_unit, 
                                       maxG = args$max_g,
                                       trimFOVs = TRUE, 
                                       outFile = file.path(args$fov_area_dir, paste0(cdid, "_total_area_per_FOV_ID.rda")),
                                       forceRecalc = args$force_recalc)
        fAreas
     })
stopCluster(cl)    

