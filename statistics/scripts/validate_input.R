suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))

log_threshold(DEBUG)

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript validate_input.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --meta_dir     directory containing all required XLSX meta files

          [OPTIONAL]
            --manifest     YAML file containing custom config
        \n"
    )
}

minReq   <- c("meta_dir")
used     <- c(minReq) 
defaults <- list(debug = TRUE)
args     <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
if(args$debug == FALSE){ log_threshold(INFO) }
logParams(args, used)

cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
log_info("\nValidating all meta data...\n")
validateInput(args)


