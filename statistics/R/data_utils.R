#' Get list of CellDive_IDs
#'
#' Get the unique list of study CellDive_IDs from meta data OR
#' a single sample ID to use in downstream calculations
#' 
#' @param metaFiles  list of all meta data files including *Samples.xlsx
#' @param cdid       a single cell dive id; when not NULL, this value
#'                   overrides the full list pulled from meta data
getSampleIDs <- function(metaFiles, cdid = NULL){
    if(!is.null(cdid)){
        return(cdid)
    }
    loadStudyAnnotations(metaFiles = metaFiles)$flat %>%
    pull(CellDive_ID) %>%
    unique
}

#' Replace Min/Max X and Y coordinates from Halo with cell midpoints
#'
#' Given a tibble containing columns XMax, Xmin, YMax, YMin, replace
#' those columns with 'X' and 'Y' columns which contain the midpoint
#' between those min and max values
#' 
#' @param ds  tibble containing columns XMax, Xmin, YMax, YMin, Sample,
#'            SPOT, UUID, Marker, Value
#' @return  tibble containing only those columns mentioned above
convertCellMinMaxToMidpoints <- function(ds){
    ds %>%
        mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2)
}

#' Subset table of center/neighborhood cell pairs based on specific marker combinations
#' within a single center cell type
#'
#' Neighborhood counts are stored for MAIN center cell types/subtypes (e.g., Macrophage, Tconv8, etc.).  Fitler data based on a specific marker combinations within CENTER cells only. 
#'  
#' @param neighborhoodCounts  table of center/neighborhood cell pairs including at least columns 
#'                            'CenterCellType' and 'C.PositiveMarkers' 
#' @param subCenterCellType   comma-delimited character string where the first element is the main 
#'                            cell type or subtype of a cell and all subsequent items are markers. 
#'                            NOTE: here is is ASSUMED that the additional parts are markers because
#'                            all cell types and subtypes are already included in neighborhood counts
#'                            table as center cells. 
#' @return filtered neighborhood counts table with 'CenterCellType' renamed to subCenterCellType
getCenterCellTypeSubset <- function(neighborhoodCounts, subCenterCellType){

    allN <- unique(neighborhoodCounts$NeighborhoodCellType)
    allF <- unique(neighborhoodCounts$FOV_ID)
    if(subCenterCellType %in% neighborhoodCounts$CenterCellType){
        return(neighborhoodCounts %>% filter(CenterCellType == subCenterCellType))
    }

    mainAndMarkers <- unlist(strsplit(subCenterCellType,","))
    main <- mainAndMarkers[1]
    markers <- mainAndMarkers[-1]
    subDat <- neighborhoodCounts %>% filter(CenterCellType == main)
    if(nrow(subDat) == 0){
        wrn <- paste0("Center cell type [",main,"] not found in data")
        log_warn(wrn)
        warning(wrn)
        return(subDat)
    }

    if(!is.null(markers) && length(markers) > 0){
        for(m in markers){
            if(grepl("-", m)){
                subDat <- subDat %>% filter(!grepl(m, C.PositiveMarkers))
            } else {
                subDat <- subDat %>% filter(grepl(m, C.PositiveMarkers))
            }
        }
    }

    if(nrow(subDat) == 0){
        wrn <- paste0("Neighborhood data for center cell type [", subCenterCellType, "] DOES NOT EXIST within data given.")
        log_warn(wrn)
        warning(wrn)
    } else {
        log_debug(paste0("Found ",nrow(subDat)," rows of neighborhood counts for center cell type: ",
                   subCenterCellType))
        subDat$CenterCellType <- subCenterCellType
    }
    subDat
}

subDivideLongSegments<-function(boundary,maxSegLength) {

    if(nrow(boundary) == 0){
        return()
    }

    bline=boundary %>% dplyr::select(X,Y) %>% mutate(SegNo=seq(nrow(.)))
    newPts=list()
    for(is in seq(nrow(bline)-1)) {

        pt1=bline[is,]
        pt2=bline[is+1,]
        subDivide=ceiling(distance(pt1,pt2)/maxSegLength)
        if(subDivide>1) {
            ##cat(is,"Dist=",distance(pt1,pt2),"subDivide =",subDivide,"\n")
            dx=pt2$X-pt1$X
            dy=pt2$Y-pt1$Y

            newPts[[is]]=data.frame(
                    X=pt1$X+dx*seq(subDivide-1)/subDivide,
                    Y=pt1$Y+dy*seq(subDivide-1)/subDivide,
                    SegNo=is+seq(subDivide-1)/subDivide
                    )
        }
    }

    bind_rows(bline,bind_rows(newPts)) %>% arrange(SegNo)

}

#' Validate and normalize character string representing a specific cell region
#'
#' Check that input cell region string is valid and can be matched to a normalized
#' character string to be used downstream
#'
#' @param cr   character string representing a specific cell region 
#'
#' @return lowercase character string that can be mapped downstream to a single
#'         cell region
getCellRegion <- function(cr = NULL){
    if(is.null(cr)){ return('fov') }
    validCR <- list('fov' = c('all', 'fov'),
                    'interface' = 'interface',
                    'interface inside' = c('interface inside', 'inside', 'inside interface'),
                    'interface outside' = c('interface outside', 'outside', 'outside interface'),
                    'neighborhood' = c('neighborhood', 'tumor neighborhood', 'tme', 'nbhd'))
    vcr <- names(validCR)[which(lapply(validCR, function(x) tolower(cr) %in% x) %>% unlist())]
    if(is.null(vcr) || length(vcr) == 0){
        msg <- paste0("Unrecognized cell region: ",cr)
        log_error(msg)
        stop(msg)
    }
    vcr
}

#' Get character string representing a single calculation unit
#' 
#' Given one or more column name, paste those names together and
#' return a single character string representing a unique unit
#' for which calculations should be made (e.g., 'FOV_ID' and 'Band')
#' 
#' @param indiv  names of columns in data to combine 
#' @param sep    character string to be used as delimiter between multiple
#'               column names; default: "__"
#' @return a single character string to be used as column name of new
#'         combined column that will contain unique calculation units
getCalcUnit <- function(indiv, sep = "__"){
    paste(indiv, collapse = sep)
}


#' Separate previously combined calculation unit character strings
#' 
#' Undo getCalcUnit()
#' 
#' @param combinedUnit  character string that is the result of separate column
#'                      names being joined to create unique calculation units
#' @param sep           character string used as delimiter between multiple
#'                      column names; default: "__"
#' @return vector of character strings
separateCalcUnit <- function(combinedUnit, sep="__"){
    unlist(strsplit(combinedUnit, sep))
}


#' Combine columns in data to create a new calculation unit column
#' 
#' Combine columns in data to create a new calculation unit column
#' 
#' @param dat   tibble of data containing columns to combine
#' @param cols  vector of column names of those to be combined
#' @param sep   character string; delimiter with which to join column names
#' @return  tibble with a new additional column that is a combination
#'          of multiple columns; original separate columns will NOT be removed 
transformCalcUnit <- function(dat, cols, sep="__"){
    cu <- getCalcUnit(cols)
    dat %>%
    unite(!!as.name(cu), cols, sep = sep, remove = FALSE)
}


