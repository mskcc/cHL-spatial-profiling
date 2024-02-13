suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))

log_threshold(DEBUG)

usage <- function(){

    cat("\nUsage:  Rscript calculate_area.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --meta_dir            path to meta files in XLSX format, required if meta_data_file is NULL
            --meta_data_file      RDA or XLSX file to which newly flattened data should be saved

          [OPTIONAL]
            --meta_files          comma-delimited character string with each element containing path
                                  to one meta file; use this when there are multiples of one or more meta
                                  file in a directory
            --manifest            YAML file containing one or more required parameter; NOTE: arguments  
                                  on command line override manifest arguments!!!         
        \n"
    )
}

dflt <- list()
minReq <- list(c("meta_dir","meta_files"), "meta_data_file")

###
### process user input
###
if(!interactive()){
    ## if manifest of parameters is provided, any additional parameters given on command line
    ## will override those provided in manifest
    suppressMessages(library(R.utils))
    args <- processCMD(commandArgs(asValue=TRUE), dflt, minReq, usage)
} else {
    args <- list(manifest       = "input/config/study_config.yaml",
                 meta_data_file = "all_meta_data.rda")
    args <- resolveConfig(dflt, read_yaml(args$manifest), args)
}


###
### validate user input 
###
if(!any(is.null(args$meta_dir), is.null(args$meta_files))){
    err <- "Values were found for both 'meta_dir' and 'meta_files'. Please provide one or the other."
    log_error(err)
    stop(err)
}

metaFiles <- getCurrentMetaFiles(metaDir = args$meta_dir, metaFiles = args$meta_files)
log_debug("Meta files:")
tmp <- sapply(metaFiles, function(x) log_debug(paste(" ", x))) 

###
### process & save meta data
###
meta <- loadStudyAnnotations(metaFiles = metaFiles, pathMap = pathConv)
write.xlsx(meta, gsub(".rda$", ".xlsx", args$meta_data_file), check.names = F)
saveRDS(meta, gsub(".xlsx$", ".rda", args$meta_data_file))
log_debug("Compiled & flattened meta data saved in files:")
log_debug(" ", gsub(".rda$", ".xlsx", args$meta_data_file))
log_debug("   and")
log_debug(" ", gsub(".xlsx$", ".rda", args$meta_data_file))

