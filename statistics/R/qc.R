#' Determine whether marker exists in annotated cell data
#'
#' If marker name occurs at least one time in EITHER the 'PositiveMarkers'
#' or 'NegativeMarkers' column, return TRUE. Otherwise, return FALSE.
#' 
#' @param dat     annotated cell table
#' @param marker  a single marker name
#' 
#' @return logical
markerExistsInData <- function(dat, marker){
    if(all(c("PositiveMarkers", "NegativeMarkers") %in% names(dat))){
        if(!any(grepl(getClassifierPattern(marker, delim=","), dat$PositiveMarkers)) &&
           !any(grepl(getClassifierPattern(marker, delim=","), dat$NegativeMarkers))){
            return(FALSE)
        }
    }
    if("Marker" %in% names(dat) && !marker %in% dat$Marker){
        return(FALSE)
    }
    TRUE
}

#' Check for markers missing in each sample
#'
#' For each sample, find out which markers are not in the data at all (i.e.,
#' no cell has marker listed under 'PositiveMarkers' OR 'NegativeMarkers' in
#' annotated cell table. Return a table of all samples with a column of known
#' missing markers (as listed in annotation), a column of actual missing markers,
#' and an 'ERROR' column indicating if the columns do not match.
#'
#' @param dataFiles  vector of all annotated cell files to check
#' @param sampAnn    flattened sample annotation table
#' @param markers    vector of all markers to check
#' @param threads
#' 
#' @return tibble with columns 'CellDive_ID', 'Expected_missing', 'Observed_missing',
#'         and 'ERROR'
checkMissingMarkers <- function(dataFiles, sampAnn, markers, threads = 1){

    cl <- makeCluster(threads, type="FORK", outfile="")
    clusterExport(cl, c("dataFiles", "sampAnn", "markers"), envir=environment())

    missing <- parLapply(cl, unique(sampAnn$CellDive_ID), function(cdid){
                   dataFile <- dataFiles[grepl(paste0("^", cdid, "_"), basename(dataFiles))]
                   if(length(dataFile) == 0){ return(NULL) }
                   log_debug(cdid, ": checking for missing markers")
                   mDat <- readRDS(dataFile)
                   if("marker.data" %in% names(mDat)){
                       mDat <- mDat$marker.data
                   }
                   oMissing <- lapply(markers, function(marker){
                                  if(!markerExistsInData(mDat, marker)){
                                      marker
                                  }
                               }) %>% unlist
                   eMissing <- getMissingMarkersInSample(sampAnn, cdid)
                   tibble(CellDive_ID = cdid,
                          Expected_missing = paste0(sort(eMissing), collapse = ","),
                          Observed_missing = paste0(sort(oMissing), collapse = ","),
                          ERROR = !identical(sort(eMissing), sort(oMissing)))
               }) %>% 
               bind_rows
    stopCluster(cl)

    missing

}
