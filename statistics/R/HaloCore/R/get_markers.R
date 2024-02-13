#' Read marker meta file
#' 
#' Read XLSX meta file containing all marker data
#'
#' @param marker_dat  tibble of all marker data; default = NULL 
#' @param marker_xlsx  XLSX file containing marker meta data; when
#'                    both marker_dat and marker_xlsx are given, 
#'                    marker_xlsx will be ignored; default = NULL
#'
#' @return tibble of all marker info
get_marker_data <- function(marker_xlsx = NULL, marker_dat = NULL){
    if(!is.null(marker_dat)){ return(marker_dat) }
    if(is.null(marker_xlsx)){ 
        stop("Must provide either marker_dat or marker_xlsx")
    }
    read.xlsx(marker_xlsx, 1, check.names = F) %>% as_tibble
}

#' Get vector of all markers used in study
#'
#' Pull names of all markers used in study listed in Marker_name column
#'
#' @param marker_dat  table read from Markers meta file
#' @param marker_xlsx  XLSX file containing marker meta data; when
#'                    both marker_dat and marker_xlsx are given, 
#'                    marker_xlsx will be ignored; default = NULL
#'
#' @return vector of all marker names 
get_all_markers <- function(marker_dat = NULL, marker_xlsx = NULL){
    get_marker_data(marker_dat = marker_dat, marker_xlsx = marker_xlsx) %>% 
    pull(Marker_name)
}


#' Get vector of all 'identity' markers used in study
#'
#' Get vector of all 'identity' markers used in study, those marked with
#' 'Y' in column Identity.
#'
#' @param marker_dat  tibble of all marker data; default = NULL 
#' @param marker_xlsx  XLSX file containing marker meta data; when
#'                    both marker_dat and marker_xlsx are given, 
#'                    marker_xlsx will be ignored; default = NULL
#'
#' @return vector of identity marker names 
get_identity_markers <- function(marker_dat = NULL, marker_xlsx = NULL){
    get_marker_data(marker_dat = marker_dat, marker_xlsx = marker_xlsx) %>% 
    filter(Identity == "Y") %>% 
    pull(Marker_name)
}


#' Get vector of all functional markers used in study
#' 
#' Get vector of all 'functional' markers used in study, those marked with
#' 'Y' in column Functional.
#'
#' @param marker_dat  tibble of all marker data; default = NULL 
#' @param marker_xlsx  XLSX file containing marker meta data; when
#'                    both marker_dat and marker_xlsx are given, 
#'                    marker_xlsx will be ignored; default = NULL
#'
#' @return vector of functional marker names
get_functional_markers <- function(marker_dat = NULL, marker_xlsx = NULL){
    get_marker_data(marker_dat = marker_dat, marker_xlsx = marker_xlsx) %>% 
    filter(Functional == "Y") %>% 
    pull(Marker_name)
}



