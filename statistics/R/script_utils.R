#' Get latest project directory
#'
#' Given the root directory for a Halo project, where the latest
#' data will always be in a subfolder named with the most recent
#' data in the form YYYY-DD-MM, return the full path to that
#' subfolder
#'
#' @param  rootDir  full path to root project directory
#' @return full path to latest subdirectory under the root
#' @export
getProjectDataDir <- function(rootDir){
    dir(rootDir, pattern = "\\d{4}(\\-\\d{2}){2}", full.names = T) %>%
    sort %>% 
    tail(1)
}


#' Set up and start cluster
#' 
#' Wrapper around makeCluster() for cleaner output
#'
#' @param threads
#' 
#' @return cluster object
startCluster <- function(threads, exports = NULL){
    cl <- makeCluster(threads, type="FORK", outfile="")

    if(!is.null(exports) && length(exports) > 0){
        clusterExport(cl, exports, envir = environment())
    }

    cl
}

#' Set log threshold
#' 
#' Set log threshold
#'
#' @param config  script configuration in list format, possibly with value set for key 'debug'; default = NULL
#' @param debug   logical or character; [T|F|'yes'|'no'|'on'|'off'] indicating whether debug messages
#'                should be printed; default = FALSE
#' @return nothing
setLogThreshold <- function(config=NULL, debug=NULL){
    log_threshold(INFO)
    on = c(TRUE, 'yes', 'on')
    off = c(FALSE, 'no', 'off')

    vals = c(config$debug, debug)
    if(!(all(vals %in% on) || all(vals %in% off))){
        log_warn(paste0("Conflicting or invalid debug argument(s): config$debug = ", config$debug, " and debug=", debug))
        log_warn(paste0("Setting debug mode to OFF"))
    } else {
        log_threshold(DEBUG)
    }
}


#' Log study parameters
#'
#' Log study parameters
#' 
#' @param args    named list of arguments provided via manifest or command line
#' @param used    named list of arguments to be used in analyses
#' @return nothing 
logParams <- function(args, used){
    cat("\n")
    cat("#~~~~~~~~~~~~~~~~~~~~~  CONFIGURATION  ~~~~~~~~~~~~~~~~~~~~~#\n")
    log_info("Working directory: ", getwd())
    for(param in used){
        val <- args[[param]]
        if(is.null(val)){
            val <-  "~"
        } else if(length(val) > 1){
            val <- paste0("[", paste(val, collapse = ","), "]")
        }
        log_info(paste0(param,": ",val,"\n"))
    }
    cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
    cat("\n")
}


#' Resolve multiple sources of configuration by priority
#'
#' Each script may be passed configuration parameters from either a manifest
#' file, at the command line or both. If both are provided, any parameter
#' passed from the command line overrides parameter found in manifest file.
#' Any required parameter that is not specified in either a manifest or at
#' the command line will be set to a default
#'
#' @param cArgs       command line arguments list
#' @param defaults    defaults set 'defaults' directory under code directory
#' @param required    vector of required parameters
#' @param sourceFile  R file containing script to source defaults and pipeline code
#' @param usage       usage string of script
#' 
#' @return configuration in list form
processCMD <- function(cArgs, defaults, required, usage, sourceFile = NULL){

    translate <- list("~" = NULL, "TRUE" = TRUE, "FALSE" = FALSE, "yes" = TRUE, "no" = FALSE)

    ## if no args, print help
    if(!interactive() && length(cArgs) == 4){ usage(); quit(save="no", status=0) }

    if(!is.null(sourceFile)){
        suppressMessages(source(sourceFile))
    }

    ## manifest args override defaults and command line args 
    ## override manifest args (if even provided)
    if("manifest" %in% names(cArgs) && !is.null(cArgs$manifest)){
        man <- read_yaml(cArgs$manifest)
        args <- resolveConfig(defaults, man, cArgs)
    } else {
        args <- resolveConfig(defaults, cArgs)
    }
    args <- args[which(!is.na(names(args)) & !names(args) == "")]
    args <- args[!names(args) %in% c("no-echo", "no-restore", "args", "file")]

    ## make sure all arguments are being interpreted correctly 
    args <- lapply(args, function(x){
                if(is.logical(x)){ return(x) }
                y <- ifelse(x %in% names(translate), translate[[x]], x)
                y <- ifelse(!is.logical(y) && !is.na(suppressWarnings(as.numeric(y))), as.numeric(y), y) 
                y
            })

    checkRequiredInput(args, required, usage())

    args
}

#' Resolve multiple configuration lists with overlapping parameters
#' 
#' Given two or more lists containing key-value pairs of study configuration,
#' build one list with parameters overwritten 
#' 
#' @param ...   vector of all lists of configuration parameters, sorted from
#'              least important to most important; duplicate paramters found in 
#'              list x will override those in list x - 1
#' @return complete list of parameters from all lists provided
resolveConfig <- function(...){

    configs <- list(...)
    finCfg <- configs[[1]]
    for(x in 2:length(configs)){
        for(param in names(configs[[x]])){
            p <- configs[[x]][[param]]
            if(!is.null(p) && length(p) > 0){
                finCfg[[param]] <- p
            }
        }
    }
    finCfg
}

#' Check that all required input exists
#' 
#' Report any missing parameters and if there are any, exit.
#' 
#' @param args    full list of named arguments
#' @param minReq  vector of all parameters required to run script
#' @param usage   function that prints usage information for script
#'
#' @return nothing; if script doesn't exit, minimum requirements have been met 
#'         and it is OK to continue
checkRequiredInput <- function(args, minReq, usage=NULL){
    missing <- c()
    for(mr in minReq){
        if(!any(mr %in% names(args))){
            missing <- c(missing, paste0("[", paste(mr, collapse = "|"), "]"))
        }
    }

    if(length(missing) > 0){
        cat(paste0("\nERROR: The following required parameters are missing from command line and/or manifest: ",missing),"\n")
        if(!is.null(usage)){
            usage()
        }

        missing <- paste(missing, collapse = ", ")
        cat(paste0("\nERROR: The following required parameters are missing from command line and/or manifest: ",missing),"\n")
        quit(save="no", status=1)
    }
}

