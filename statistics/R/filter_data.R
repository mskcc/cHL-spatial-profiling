#' Filter annotated cell data for cells in a specific neighborhod
#' 
#' Get all cells that fall in a specific neighborhood. For example
#' extract cells that are within 30 microns of at least one tumor cell.
#' The neighborhood criteria is provided as one or more quosure of expressions. 
#' e.g., quos(Cell_type == 'Tumor', Dij.micron < 30, grepl('MHCII', PositiveMarkers))
#' 
#' @param annDat   tibble of annotated cell data to be filtered
#' @param nbhrDat  tibble of all cell-cell pairs and the distance between them in microns
#' @param neighborhood  quosure or list of quosures containing expressions that define 
#'                 a neighboring cell of interest  
filterForCellsInNeighborhood <- function(annDat, nbhrDat, neighborhood, minNeighbors = 1){

    centers <- annDat %>%
               filter(!!!cellStateFilter(neighborhood)) %>% 
               pull(UUID)

    nbrs <- nbhrDat %>% 
            filter(C.UUID %in% centers) %>%
            group_by(N.UUID) %>%
            summarize(num = n()) %>%
            filter(num >= minNeighbors)

    annDat %>% 
    filter(UUID %in% nbrs$N.UUID) 

}


cellStateFilter <- function(condStr){
    pts <- unlist(strsplit(condStr, ","))
    ct_st <- pts[1]
    pts <- pts[-1]
    pos <- pts[-grep("\\-", pts)]
    neg <- gsub("\\-", "", pts[grep("\\-", pts)])

    pos_ptrn <- paste0("(^|,)", paste(pos, collapse = "(,.*|$)"), "(,|$)")
    neg_ptrn <- paste(paste0("(^|,)", neg, "(,|$)"), collapse = "|")

    quos((Category == ct_st | Cell_type == ct_st | Subtype == ct_st),
         grepl(pos_ptrn, PositiveMarkers),
         !grepl(neg_ptrn, PositiveMarkers))
}


#' Filter annotated cell data for a single condition
#' 
#' Get all cells that meet all criteria outlined in condition string
#' 
#' @param dat           table of annotated cell data, where one row represents a single cell
#'                      and, at minimum, columns for Classifiers and PositiveMarkers 
#' @param conditionStr  comma-delimited string containing a combination of cell types (classifiers) 
#'                      and/or markers
#' @param markers       vector of individual markers used to distinguish markers from classifiers
#'
#' @param cellTypeCols  vector of column names for those columns containing cell identifiers; default:
#'                       c("Category", "Cell_type", "Subtype")
#'
#' @return subset of dat including only cells that match ALL criteria included in conditionStr
filterDataForCondition <- function(dat, conditionStr, markers, cellTypeCols = c("Category", "Cell_type", "Subtype")){

    if(conditionStr == "DAPI"){ return(dat) }

    pts <- unlist(strsplit(conditionStr,","))
    pts <- pts[pts != "DAPI"] 

    class  <- pts[!gsub("\\-", "", pts) %in% markers]
    if(length(class) > 1){
        wrn <- paste0("Unable to parse condition '", conditionStr, "' -- Skipping.")
        log_warn(wrn)
        warning(wrn)
        return(NULL)
    }
    mCombo <- paste(sort(pts[pts != class]), collapse=",")

    dat %>%
    filter_at(all_of(cellTypeCols), any_vars(. == class)) %>%
    filterForMarkerCombo(mCombo)

}

#' Filter data tibble for a selection of calculation units
#'
#' Filter data tibble for a selection of calculation units
#' 
#' @param dat       data tibble to be filtered
#' @param calcUnit  name of column corresponding to the calculation unit of 
#'                  interest (default: FOV_ID)
#' @param include   vector of calculation units to keep
#' @param ignoreWarnings logical; when TRUE, warnings about missing units will NOT 
#'                       be printed
#'
#' @return subset of dat consisting of only calcuation units in 'include' vector
filterForCalculationUnits <- function(dat, calcUnit = "FOV_ID", include = NULL, ignoreWarnings = FALSE){
    if(!calcUnit %in% names(dat)){
        err <- paste0("Invalid calculation unit: ", calcUnit)
        log_error(err)
        stop(err)
    }
    if(!is.null(include)){
       if(!all(include %in% dat[[calcUnit]]) && !ignoreWarnings){
           wrn <- paste0("Calculation units missing from data: ",
                           paste(include[!include %in% dat[[calcUnit]] ], collapse=", "))
           log_warn(wrn)
           warning(wrn)
       }
       log_debug(paste0("Number of FOVs INCLUDED in analysis: [", length(include[include %in% dat[[calcUnit]]]), "]"))
       return( dat %>% filter_at(all_of(calcUnit), all_vars(. %in% include)) )
    }
    log_debug(paste0("Number of FOVs INCLUDED in analysis: [", length(unique(dat[[calcUnit]])), "]"))
    dat
}


#' Filter data for cells that are part of a single comparison group
#'
#' Given a list that defines a single group of cells, filter data for only those cells
#'
#' @param dat             data table to be filtered
#' @param group           list describing data to be included in group; list names are
#'                        either column names from dat or 'Cell Region' and values are
#'                        the values to keep
#' @param tumorNbhdCells  vector of UUIDs for cells that fall within the neighborhood of at least
#'                        one tumor cell
#'
#' @return filtered data table
filterForComparisonGroup <- function(dat, group, nbhdDir = NULL){ #, tumorNbhdCells = NULL){
    grp <- dat

    ## filter for neighbor cells first if necessary
    if(isNeighborhoodQuestion(group)){
        log_debug("Filtering for neighborhood")
        neighborCells <- loadCellNeighborhood(nbhdDir, 
                             group$`Center Cell Type`, 
                             centerMarkerCombo = group$`Center Marker Combo`,
                             centerClass = group$`Center Cell Class`,
                             neighborCellType = group$`Neighbor Cell Type`,
                             neighborMarkerCombo = group$`Neighbor Marker Combo`,
                             neighborClass = group$`Neighbor Cell Class`,
                             pullColumn = "neighborUUIDs")
        grp <- grp %>% filter(UUID %in% neighborCells)
    } 

    otherFilters <- names(group)[!grepl("center cell|neighbor cell", tolower(names(group)))]
    for(filt in otherFilters){
        if(filt == "Cell Region"){
            log_debug(paste0("Filtering for cells in region: ", group[[filt]]))
            cr <- getCellRegion(tolower(group$`Cell Region`))
            grp <- grp %>% filterDataForCellRegion(cr)
        } else if(filt %in% names(dat)){
            log_debug(paste0("Filtering ",filt," for ",paste(group[[filt]], collapse = ",")))
            grp <- grp %>%
                   filter_at(all_of(filt), all_vars(. %in% group[[filt]]))
        } else {
            warning(paste0("Can not filter on column that does not exist: ",filt,"...skipping."))
        }
    }
    grp
}

#' Filter for cells that are positive for certain markers and negative for others
#' 
#' Given a combination of positive and negative markers, filter data for cells
#' that match that match that combination
#'
#' @param dat           tibble of cell data to filter
#' @param combo         comma-delimited string of positve markers (indicated by 
#                       lack of '-' character) and negative markers (indicated 
#                       by '-' appended to end of marker name)
#' 
#' @return tibble of all cells that match the marker combination given
filterForMarkerCombo <- function(dat, combo){

    if(is.null(combo) || length(combo) == 0 || combo == ""){
        return(dat)
    }

    pts <- unlist(strsplit(combo,","))
    pts <- pts[pts != "DAPI"]

    neg     <- gsub("-", "", pts[grep("-",pts)])
    negPtrn <- sapply(neg, function(x){
                       paste(getClassifierPattern(x,delim=","))
                 }) %>%
               unlist %>% paste(collapse = "|")

    pos     <- sort(pts[!grepl("\\-", pts)])
    posPtrn <- paste0("(^|,)", paste0(pos, collapse = ",.*"), "(,|$)")

    dat %>%
    { if(length(neg) > 0) filter(., !grepl(negPtrn, PositiveMarkers)) else . } %>%
    { if(length(pos) > 0) filter(., grepl(posPtrn, PositiveMarkers, perl=TRUE)) else . }

}

#' Filter for cells that are in a specific cell region
#' 
#' Remove from data any cells that fall outside the region of interest
#' 
#' @param dat             any halo data table; in order to filter for cells within
#'                        the interface, table must have a 'Band' column
#' @param cellRegion      the region to filter for
#'                          - fov: region is total FOV, return entire dat table
#'                          - interface: keep any row that is assigned to a Band
#'                          - interface inside: keep any row that is in a negative distance Band
#'                          - interface outside: keep any row that is in a positive distance Band
#' @param tumorNbhdCells  character vector of cell UUIDs for cells that are in the neighborhood of
#'                        at least one tumor cell; to be used for filtering data for 'neighborhood'
#'                        cell region
#'
#' @return  filtered data table
filterDataForCellRegion <- function(dat, cellRegion){ 
    valid <- c("fov","interface","interface inside", "interface outside", "stroma") 
    if(!tolower(cellRegion) %in% valid){
        stop(paste0("Unrecognized cell region: ",cellRegion))
    }
    dat <- switch(tolower(cellRegion),
                  "interface" = dat %>% filter(!is.na(Band)),
                  "interface inside" = dat %>% filter(!is.na(Band), grepl("\\-", Band)),
                  "interface outside" = dat %>% filter(!is.na(Band), !grepl("\\-", Band)),
                  "stroma" = dat %>% filter(Distance > 0),
                  dat)
    dat
}


#' Extract area information from annotated cell data for each FOV
#' for a specific cell region
#'
#' Pull out only area information from annotated cell data for each FOV
#' depending on the cell region of interest
#' 
#' @param dat           annotated cell data that has been joined with cell region areas
#'                      for cell region of interest
#' @param cellRegion    region of interest [fov|interface|interface inside|interface outside]
#' @param groupBy       data column name(s) containing values for which a single area value
#'                      should be calculated (default: FOV_ID)
#'
#' @return  a table containing columns [groupBy], Area 
filterForAreas <- function(dat, cellRegion = "fov", groupBy = "FOV_ID"){
    groupCols  <- unique(c("Sample_ID","FOV_ID","Band",groupBy)) ## cols required to identify unique values;
                                                           ## must keep them so that two rows
                                                           ## that happen to be identical without any
                                                           ## one of them are all kept (e.g., a sample 
                                                           ## in which Band1 in FOV1 has the same area
                                                           ## as Band1 in FOV2
    areaCol    <- "FOVBandArea"

    if(tolower(cellRegion) == "fov"){
        groupCols <- groupBy
        areaCol <- "FOVArea"
    }

    dat %>%
    select_at(c(groupCols, areaCol)) %>%
    unique() %>%
    rename_at(vars(areaCol), list(~ (. = "Area"))) %>%
    filterDataForCellRegion(cellRegion) %>%
    group_by_at(groupBy) %>%
    summarize(Area = sum(Area, na.rm = T))
}


#' Filter conditions index for all conditions of a single analysis type
#' 
#' Given a table of all conditions including column 'AnalysisType', filter
#' for conditions of a single type (e.g., spatial or general)
#' 
#' @param condIdx       table of conditions including column 'AnalysisType'
#' @param analysisType  character string matching a single value in column
#'                      'AnalysisType'
#'
#' @return  subset of condIdx
filterConditionsByAnalysisType <- function(condIdx, analysisType){
    condIdx %>%
    filter(AnalysisType == analysisType) %>%
    select_if(function(x){ !all(is.na(x)) })
}

#' Filter table of neighborhood data for neighborhood cells
#' that satisfy all filter criteria including classifiers and/or positive markers
#' 
#' @param dat   data tibble of neighborhood data including columns for Center cells
#'              and columns for Neighborhood cells
#' @param nbhd  comma-delimited character string describing neighborhood cells to keep
#'              (e.g., "Tconv8,PD1,LAG3,TIM3-")
#' @param markers character vector containing all study markers, used to distinguish
#'                markers from cell types in {nbhd} string
#' 
#' @return  filtered data tibble
filterForNeighborhood <- function(dat, nbhd, markers){
    n <- unlist(strsplit(nbhd,","))
    classes <- n[!n %in% c(markers, paste0(markers, "-"))]
    mrkrs <- n[n %in% c(markers, paste0(markers,"-"))]

    if(length(classes) > 0){
        for(cls in classes){
            if(cls == "EXM0"){
                dat <- dat %>%
                       filter(!grepl("EXM",N.Classifiers))
            } else {
                dat <- dat %>%
                       filter(grepl(getClassifierPattern(cls,delim=";"),N.Classifiers))
            }
        }
    }
    if(length(mrkrs) > 0){
        for(mrk in mrkrs[!grepl("\\-", mrkrs)]){ dat <- dat %>%
                                                  filter(grepl(getClassifierPattern(mrk,delim=","),
                                                             N.FullPosMarkerStr)) }
        for(mrk in mrkrs[grepl("\\-", mrkrs)]){ dat <- dat %>%
                                                filter(!grepl(getClassifierPattern(gsub("-","",mrk), delim=","),
                                                              N.FullPosMarkerStr)) }
    }
    return(dat)
}


