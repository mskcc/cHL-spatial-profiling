#' Generate a single tibble from several files
#' 
#' Provided a directory path or a vector of files, read all RDA files
#' and create a single tibble by binding the rows of all individual tibbles
#'
#' @param dataDir   full path to directory of RDA files; all files found in directory
#'                  will be included (default = NULL)
#' @param files     vector of paths to RDA files to be included
#' @param keepCols  vector of column names to keep; all others will be excluded from
#'                  final tibble
#' @return a single tibble 
tibbleFromMultipleFiles <- function(dataDir = NULL, files = NULL, keepCols=NULL){

    if(is.null(files) && !is.null(dataDir)){
        files <- file.path(dataDir, dir(dataDir))
    }
    tbl <- tibble()
    for(fl in files){
        flTbl <- readRDS(fl)
        if(!is.null(keepCols) && nrow(flTbl) > 0){
            flTbl <- flTbl %>% select_at(keepCols)
        }
        tbl <- tbl %>% bind_rows(flTbl %>% unique())
    }
    tbl %>% unique()
}


#' Convert a string containing underscores to separate words, all capitalized
#' 
#' Convert a string containing underscores to separate words, all capitalized
#' 
#' @param str    character string containing underscores (e.g., 'Column_title')
#'
#' @return capitalized string (e.g., 'Column Title')
underscoreToUpper <- function(str){
    mtch <- str_match_all(str, "_.{1}")[[1]][,1]
    for(m in mtch){
        str <- gsub(m, toupper(gsub("_"," ",m)), str)
    }
    str
}


#' Get a list of study meta files
#' 
#' Given either a directory containing xlsx meta files or study config,
#' get a full vector of meta file paths. When giving metaDir, the vector
#' will be paths to all XLSX files in the directory. When giving config,
#' if parameter 'meta_files' is specified, that value will be returned. 
#' Otherwise all files in directory specified by parameter 'meta_dir' 
#' will be returned.
#' 
#' @param metaDir    directory containing meta files
#' @param config     configuration with either specific files given in 'meta_files'
#'                   or directory given in 'meta_dir'
#'
#' @return vector of XLSX meta data files
getMetaFiles <- function(metaDir=NULL, config=NULL){
    if(!is.null(metaDir)){
        dir(metaDir, full.names = T, pattern = ".xlsx") 
    } else if(!is.null(config)){
        if(!is.null(config$meta_files)){
            config$meta_files      
        } else if(!is.null(config$meta_dir)){
            metaDir <- config$meta_dir
            file.path(metaDir, dir(metaDir)[grepl(".xlsx",dir(metaDir))])
        } else {
            stop("No meta data given/found.")
        }
    }
}

getCurrentMetaFiles <- function(metaFiles = NULL, metaDir = NULL, config = NULL){
    if(is.null(metaFiles)){
        metaFiles <- getMetaFiles(metaDir = metaDir, config = config)
    }
    patterns <- c("CellStates", "CellTypes", "FOVs", "Markers", 
                  "Paths", "Questions", "Samples", "StudyOverview")

    lapply(patterns, function(x){
        getMetaFile(x, metaFiles = metaFiles, metaDir = metaDir, config = config)
    }) %>% unlist
}

#' Get any one of the meta files
#'
#' Utility function to extract a specific meta file using that file's unique
#' identifier (e.g., 'Markers', 'Samples', etc.). If there are multiple files
#' matching the identifier, assume they include version numbers following pattern
#' '__V*.xlsx' and return the one with the largest number. Must provide at least one
#' argument from options 'metaFiles', 'metaDir', and 'config'. If more than one of
#' these is given, metaFiles takes priority over metaDir, which takes priority over
#' config.
#'
#' @param fileID      unique character string identifying the file of interest
#' @param metaFiles   vector of all meta files; required if metaDir AND config are NULL
#' @param metaDir     path to directory containing all meta files; required if metaFiles
#'                    AND config are NULL
#' @param config      study configuration in list form including key/value pairs for either
#'                    'metaFiles' or 'metaDir'; use when both metaFiles and metaDir are NULL
#' @return file matching pattern [fileID]
getMetaFile <- function(fileID, metaFiles = NULL, metaDir = NULL, config = NULL){
    if(is.null(metaFiles)){
        metaFiles <- getMetaFiles(metaDir = metaDir, config = config)
    }
    matches <- metaFiles[grep(fileID, metaFiles)]
    if(length(matches) == 1){ return(matches) }
    currentV <- gsub(".*__V", "", gsub(".xlsx", "", matches)) %>% as.integer %>% max
    matches[grep(paste0("_V", currentV, ".xlsx"), matches)]
}


#' Add backticks around a string
#' 
#' Mostly for column names with special characters, add a backtick
#' before and after string 
#' 
#' @param str   character string
#'
#' @return str with backticks added to beginning and end
add_ticks <- function(str){
    paste0("`",str,"`")
}

#' Remove leading or trailing backticks from string
#'
#' Remove leading or trailing backticks from string
#'
#' @param str   character string containing backticks
#'
#' @return str with backticks removed
remove_ticks <- function(str){
    gsub("^`","",gsub("`$","",str))
}


#' Get standard star notation of a numeric value
#' 
#' Get standard star notation of a numeric value, generally a p.value; 
#' val < 0.001 = "***", val < 0.01 = "**", val < 0.05 = "*"
#'
#' @param val   numeric value
#' @return character string of stars, the length of which signifies significance
signifStars <- function(val){
    ifelse(val < 0.001, "***",
           ifelse(val < 0.01, "**",
                  ifelse(val < 0.05, "*", "")))
}


#' Insert newlines into a long character string
#'
#' Wrap a long character string by inserting newline characters wherever appropriate
#' 
#' @param longStr    character string to be wrapped
#' @param maxChar    integer; the maximum number of characters allowed on one line
#'
#' @return wrapped character sring
wrapTitle <- function(longStr, maxChar=60){
    if(nchar(longStr) <= maxChar || grepl("\n",longStr)){ return(longStr) }

    finalStr <- ""
    for(x in 1:(ceiling(nchar(longStr)/maxChar))){
        part <- substr(longStr, maxChar*(x-1)+1, maxChar*x)
        if(!nchar(part) < maxChar){
            spcs  <- grep(" |/",unlist(strsplit(part, "")))
            lnBrk <- spcs[which(abs(spcs - maxChar) == min(abs(spcs - maxChar)))]
            part  <- paste0(replace(unlist(strsplit(part,"")), lnBrk, c("\n")),collapse="")
        }
        finalStr <- paste0(finalStr, part)
    }
    finalStr
}



