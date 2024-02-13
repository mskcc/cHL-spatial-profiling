read.xlsx <- function(xl_file, sheet = 1, check.names = F){
    tbl <- xlsx::read.xlsx(xl_file, sheet, check.names = check.names) 
    tbl[, names(tbl) != 'NA'] %>%
    mutate_if(is.character, stringr::str_trim)
}




#' Extract valid file extension from path
#'
#' Separate file extension from file path basename. If vector of 
#' valid/expected file extensions is provided and extension is
#' not equal to one of them, an error is thrown.
#'
#' @param path      path containing basename from which file extension 
#'                  will be extracted
#' @param expected  optional character vector of expected file extensions
#'
#' @return if path contains an expected file type, the file extension is returned
validFileType <- function(path, expected = NULL){
    ext <- tools::file_ext(path)
    if(!is.null(expected) && !ext %in% expected){
        msg <- paste0("Unexpected file type: [", ext, "]")
        log_error(msg)
        stop(msg)
    }
    ext
}

#' If given PDF file is not NULL, open it
#' 
#' If given PDF file is not NULL, open it
#' 
#' @param pdfFile     path to file
#' @param pdfHeight   height of PDF in inches; default = 8.5
#' @param pdfWidth    width of PDF in inches; default = 11
#' @param singlePage  logical; when TRUE, a blank page will appear before anything is 
#'                    printed to file (?????)
#'
#' @return nothing
openPDF <- function(pdfFile, pdfHeight=8.5, pdfWidth=11, singlePage=FALSE){
    if(!is.null(pdfFile)){
        pdf(pdfFile, height=pdfHeight, width=pdfWidth, onefile = !singlePage)
    }
}

#' If given PDF is not null, close it
#' 
#' If given PDF is not null, close it
#' 
#' @param pdfFile  path to file to be closed
#'
#' @return nothing
closeOpenPDF <- function(pdfFile){
    if(!is.null(pdfFile)){
        dev.off()
    }
}

#' Get subset of files either from a directory or from a list of files
#' that have names matching a specific pattern
#' 
#' Given either a directory path or a vector of file names, return either
#' all files in the file list or, if a pattern is provided, the subset
#' of files matching that pattern
#'
#' @param path    full path to directory containing files of interest; if 'pattern' is 
#'                provided, return only the files within the directory that match; default: NULL
#' @param files   vector of files; if 'pattern' is provided, return only the files that
#'                match; default = NULL
#' @param pattern character pattern of interest; all files that match pattern will be returned
#'
#' @return  vector of the full set or a subset of given file list
getFiles <- function(path=NULL, files=NULL, pattern=NULL){

    if(is.null(path) && is.null(files)){
        stop("Need either a directory or a vector of files and a pattern to get files.")
    }
    fls <- c()
    if(!is.null(files)){
        fls <- files
    } else if(!is.null(path)){
        fls <- file.path(path, dir(path))
    } else {
        stop("Must provide either path or files argument.")
    }
    if(!is.null(pattern)){
        fls <- fls[grepl(pattern, basename(fls))]
    }
    unique(fls)

}

#' Determine whether a file exists and is not empty
#'
#' Determine whether a file exists and is not empty
#' 
#' @param fileName  file name

#' @return logical
fileDone <- function(fileName){
    !is.na(fileName) && !is.null(fileName) && length(fileName) > 0 && file.exists(fileName) && file.size(fileName) > 0
}


#' Determine whether a file in a configuration list needs to/should be written
#' 
#' @param lst  list of study parameters
#' @param key  key of file to test

#' @return logical
needToWrite <- function(lst, key){
    is.null(lst[[key]]) || length(lst[[key]]) == 0 || !file.exists(lst[[key]])
}


#' Create a directory recursively without warnings
#' 
#' Create a directory recursively without warnings
#' 
#' @param path    full path of directory to create
#'
#' @return nothing
mkdir <- function(path){
    dir.create(path, recursive = T, showWarnings = F)
}


#' Make sure a vector of files is not empty
#'
#' Given a vector of files, throw an error if vector is empty
#' or if any one of the files does not exist
#' 
#' @param files  vector of files to check
#' @param path   directory containing files of interest
#' @param desc   a character string describing the expected file type
#'               (for logging purposes, can be anything)
#' 
#' @return nothing 
checkFilesFound <- function(files, path, desc, error = TRUE){
    msg = NULL
    if(is.null(files) || length(files) == 0){
        msg <- paste("No", desc, "files found in directory:", path)
    } else if(!all(file.exists(files))){
        msg <- paste("One or more file does not exist: ", 
                     paste(files[which(!file.exists(files))], collapse = ", "))
    }
    if(!is.null(msg)){
        if(error){
            log_error(msg)
            stop()
        }
        log_warn(msg)
    }
}

#' Convert file paths from one location to another
#' 
#' Convert a table of file paths from one mount to another; in our case,
#' we're converting paths from shared drives to server paths. Paths are
#' mapped in a list named by the paths existing in the paths table, with
#' values being the paths to which they should be converted in meta data
#'
#' @param  paths    tibble of file paths to be converted
#' @param  pathMap  list of new file paths, named by the ones existing in paths tibble
convertServerPaths <- function(paths, pathMap){
    lapply(paths, function(x){
        for(pth in names(pathMap)){
            x <- gsub(pth, pathMap[[pth]], x)
        }
        x
    }) %>% as_tibble()
}


formatted_cycles <- function(cycle1, cycleN){
    paste0('S',
           formatC(seq(gsub("S0*", "", cycle1), gsub("S0*", "", cycleN)),
                   width = 3, format = 'd', flag = '0'))
}

link_halo_csv <- function(flatMeta, localCSVdir){

    for(cdid in unique(flatMeta$CellDive_ID)){

        sampDat <- flatMeta %>% filter(CellDive_ID == cdid)

        csvPath <- sampDat$path_to_halo_csv_AI %>% unique
        csvPath <- gsub("/$", "", csvPath) # remove trailing slashes
        csvFile <- dir(file.path(csvPath), full.names=T, pattern = paste0(cdid, ".csv")) # this should find the file whether it has *.gz ext or not
        if(length(csvFile) > 1){
            log_error(paste0("More than one CSV file matching pattern '", paste0(cdid, ".csv"), "'. Don't know which to link. Skipping."))
            next 
        } else if(length(csvFile) == 0){
            log_error(paste0("No CSV file matching pattern '", paste0(cdid, ".csv"), "' found. Skipping."))
            next
        }
        file.symlink(csvFile, localCSVdir) 
    }

}


link_halo_images <- function(flatMeta, sourceDir, localDir, imageIdentifier){

    itMap <- list(FOV = 'spot', WS = 'region')

    linkCount <- 0
    for(cdid in unique(flatMeta$CellDive_ID)){

        sampDat <- flatMeta %>% filter(CellDive_ID == cdid)

        ## link only images to be kept in analysis
        cycle1    <- sampDat$DAPI_first %>% unique
        cycleN    <- sampDat$DAPI_last %>% unique
        fovs      <- sampDat$FOV_number %>% unique
        imageType <- itMap[[sampDat$Image_type %>% unique]]

        ## link scanplan
        origPath <- file.path(sourceDir, cdid, 'scanplan', 'scanplan.tif')
        localPath <- file.path(localDir, cdid, 'scanplan.tif')
        if(file.exists(origPath)){ 
            file.symlink(origPath, localPath) 
            linkCount <- linkCount + 1
        }

        ## link moves
        origPath <- file.path(sourceDir, cdid, 'moves.dat')
        localPath <- file.path(localDir, cdid, 'moves.dat')
        if(file.exists(origPath)){ 
            file.symlink(origPath, localPath) 
            linkCount <- linkCount + 1
        }
 
        for(cycleF in formatted_cycles(cycle1, cycleN)){
            origPath  <- file.path(sourceDir, cdid, 'RegisteredImages', cycleF)
            localPath <- file.path(localDir, cdid, cycleF)
            mkdir(localPath)

            for(fov in fovs){
                fovF <- formatC(fov, width = 3, format = 'd', flag = '0')
                fileName = paste(cycleF, imageIdentifier, imageType, paste0(fovF, '.tif'), sep = "_")
                if(file.exists(file.path(origPath, fileName))){
                    file.symlink(file.path(origPath, fileName), localPath)            
                    linkCount <- linkCount + 1 
                } 
            }
        }

    }
    linkCount 
}

link_halo_files <- function(flatMeta, localCSVdir, localImageDir, sourceImageIdentifier){

    link_halo_csv(flatMeta, localCSVdir)
    link_halo_images(flatMeta, sourceImageIdentifier, localImageDir)

}

link_existing_files <- function(locationList){

    lapply(locationList, function(x){
        srcFiles <- dir(x$src, pattern = x$pattern)
        log_info("Linking to [ ", length(srcFiles), " ] files in ", x$src)
        if(length(srcFiles) > 0){
            file.symlink(file.path(x$src, srcFiles), file.path(x$local, srcFiles))
        }
    })

}

