#' Generate vector of positive marker combinations
#' 
#' Given a vector of markers, expand into vector of all combinations of those
#' markers, with the option of including the empty combination (no positive markers)
#' 
#' @param  markers     vector of individual markers
#' @param  into        character string; type of vector 'into' which markers should be expanded; currently function
#'                     only supports expanding into 'allCombinations'
#' @param  emptyCombo  logical; TRUE indicates the empty combo should be included in final vector; default=FALSE
#' @return  list of expanded cell types
#' @export
expandMarkerCombos <- function(markers, into="allCombinations", emptyCombo=FALSE){
    markers <- unlist(strsplit(markers,","))

    combos <- switch(tolower(into),
                     "allcombinations" = getAllCombos(markers, superNeg = emptyCombo),
                     "singlemarkers" = markers,
                     NULL)

    if(is.null(combos)){
        log_error("Code currently can only expand markers into 'allCombinations' or 'singleMarkers'")
        stop(paste0("Unsupported 'into' argument: ",into))    
    }

    #if(emptyCombo){
    #    combos <- c(combos, "")
    #}

    return(combos)
}


#' Parse cell types XLSX file                                                   
#'                                                                              
#' Parse cell types XLSX file, expanding cell types and assigning labels when necessary
#'                                                                              
#' @param cellTypesXLSX   xlsx file in cell types format (see docs)             
#'                                                                              
#' @return table containing cell type definitions to be used for cell annotation
#' @export                                                                      
getCellTypes <- function(cellTypesXLSX){   
    read.xlsx(cellTypesXLSX, 1, check.names = F) %>%
    filter(!is.na(Category)) %>%
    mutate(Cell_type_figures_subscript = ifelse(is.na(Cell_type_figures_subscript), 
                                                "", Cell_type_figures_subscript),
           Subtype_figures_subscript = ifelse(is.na(Subtype_figures_subscript), 
                                              "", Subtype_figures_subscript)) 
} 

#' Parse cell types XLSX file
#' 
#' Parse cell types XLSX file, expanding cell types and assigning labels when necessary
#' 
#' @param cellTypesXLSX   xlsx file in cell types format (see docs)
#' 
#' @return table containing parsed/expanded cell types to be used for cell annotation
#' @export
getCellTypes_OLD <- function(cellTypesXLSX){

    if(!file.exists(cellTypesXLSX) || file.size(cellTypesXLSX) == 0){
        stop(paste0("Cell Types file ", cellTypesXLSX, " is missing or empty."))
    }

    cellTypes   <- as_tibble(read.xlsx(cellTypesXLSX,1))
    toExpand    <- cellTypes %>% filter(Pos_required != 'all')
    allExpanded <- cellTypes %>% filter(Pos_required == 'all')

    for(te in 1:nrow(toExpand)){
        curRow <- toExpand[te,]
        markers <- curRow %>% pull(Pos_markers)
        into <- NULL
        into <- ifelse(curRow$Pos_required == "+", 
                       "allCombinations", 
                       "singleMarkers")
        combos <- expandMarkerCombos(markers, into=into, emptyCombo=FALSE)
        neg <- curRow$Neg_markers[!is.na(curRow$Neg_markers)]
        expanded <- combos %>%
                    lapply(., function(x){
                                   x <- sort(unlist(strsplit(x,",")))
                                   curRow %>%
                                   mutate(Pos_markers = paste(x[!grepl("-",x)],
                                                              collapse=","),
                                          Neg_markers = paste(c(neg, 
                                                                gsub("\\-", "", x[grepl("-",x)])),
                                                                collapse=","))
                              }
                    ) %>%
                    bind_rows()
        empt <- which(is.null(expanded$Neg_markers) | expanded$Neg_markers == "")
        expanded$Neg_markers[!empt] <- paste(expanded$Neg_markers[!empt], curRow$Neg_markers, sep=",")
        expanded$Neg_markers[empt] <- curRow$Neg_markers

        allExpanded <- bind_rows(allExpanded, expanded)

    }

    allExpanded[allExpanded == "MΦ+"] <- "MΦ"
    ########################

    return(allExpanded)
}

#' Format cell types to show positive/negative markers explicitly
reformatCellTypes <- function(cellTypes){
    cts <- lapply(1:nrow(cellTypes), function(ctRow){
             ct = cellTypes[ctRow,]
             res <- ct %>% select(!dplyr::matches("Pos|Neg|arkers"))
             pos <- unlist(strsplit(ct$Pos_markers, ","))
             neg <- unlist(strsplit(ct$Neg_markers, ","))
             if(ct$Pos_required == "all"){
                 mrkrs <- c(rep(1, length(pos)), rep(0, length(neg)))
                 names(mrkrs) <- c(pos, neg) 
                 return(bind_cols(res, as_tibble(as.list(mrkrs))))
             } else if(ct$Pos_required == "+"){ ## one or more required
                 tmp <- tibble()
                 for(x in 1:length(pos)){
                     mrkrs <- c(1, rep(0, length(neg)))
                     names(mrkrs) <- c(pos[x], neg)
                     tmp <- tmp %>% bind_rows(as_tibble(as.list(mrkrs)))
                 }
                 return(bind_cols(res, tmp))
             }
           }) %>% 
           bind_rows 

    allMrkrs <- cts %>% select_if(is.numeric) %>% names
    superNeg <- rep(0, length(allMrkrs))
    names(superNeg) <- allMrkrs

    cts %>% 
    bind_rows(as_tibble(as.list(superNeg)) %>% 
              mutate(Category = "Other", 
                     Cell_type = "superNeg", 
                     Subtype = "superNeg"))
}


#' Add class(es) to pre-existing classes in Classifier column of annotated cell data
#' 
#' Given a vector of character strings where each element is a delimited collection of
#' classifiers for a single cell, add to every element a delimited string of new classes
#'
#' @param  existingClasses   a character vector where each element is a delimited (default delim=";")
#'                           string of classes to which a single cell belongs
#' @param  newClasses        a delimited character string of new classes to add to every item in 
#'                           {existingClasses}; NOTE: setting unique at this point will NOT 
#' @param  delim             delimiter separating individual classes in both {existingClasses} and {newClasses}
#'
#' @return  vector of expanded class lists
appendClassifiers <- function(newClasses, existingClasses = NULL, unique=FALSE, delim=";"){

    exC <- existingClasses
    if(!is.null(exC) && !is.na(exC)){
        exC <- unlist(strsplit(existingClasses, delim))
        exC <- exC[!is.na(exC)]
    } else {
        exC <- NULL
    } 
    if(unique){
        newClasses <- unlist(strsplit(newClasses, delim))
        newClasses <- newClasses[!is.na(newClasses)]
        nwC <- paste(unique(c(exC, newClasses)), collapse=delim)
    } else {
        nwC <- paste(c(exC, newClasses), collapse=delim)
    }
    return(nwC)
}

#' Get vector of all markers used in study
#'
#' Get vector of all markers used in study
#'
#' @param markerDesc  table read from Markers meta file
#' @param markerFile  XLSX file containing marker meta data; when
#'                    both markerDesc and markerFile are given, 
#'                    markerFile will be ignored
#'
#' @return vector of all markers under Marker_name column of markers file
getAllMarkers <- function(markerDesc = NULL, markerFile = NULL){
    if(is.null(markerDesc) & is.null(markerFile)){
        stop("Must provide either markerDesc or markerFile")
    }
    if(is.null(markerDesc)){
        markerDesc <- read.xlsx(markerFile, 1, check.names = F) 
    }
    markerDesc %>% pull(Marker_name)
}


#' Get vector of all 'identity' markers used in study
#'
#' Get vector of all 'identity' markers used in study
#'
#' @param markerDesc  read from Markers meta file
#'
#' @return vector of markers classified as 'Identity' markers by the word 'Identity'
#'         in the Description column of markers file 
getIdentityMarkers <- function(markerDesc){
    markerDesc %>% filter(Identity == "Y") %>% pull(Marker_name)
}


#' Get vector of all functional markers used in study
#' 
#' Get vector of all 'identity' markers used in study
#'
#' @param markerDesc  read from Markers meta file
#'
#' @return vector of markers classified as 'Identity' markers by the word 'Identity'
#'         in the Description column of markers file 
getFunctionalMarkers <- function(markerDesc){
    markerDesc %>% filter(Functional == "Y") %>% pull(Marker_name)
}



#' Add columns to indicate all positive markers and positive cell identity markers
#'
#' Given a tibble where each row is a unique cell and columns for each marker have
#' a value of 1 (marker is positive in that cell) or 0 (marker is negative in that cell),
#' add columns 'PositiveMarkers' and 'CTMarkers'. 'PositiveMarkers' is a comma-delimited
#' string of ALL positive markers in a cell and 'CTMarkers' is a comma-delimited string
#' of ONLY cell identity markers
#'
#' @param dat         tibble where each row is a unique cell and columns for each marker have
#'                    a value of 1 (marker is positive in that cell) or 0 (marker is negative 
#'                    in that cell)
#' @param markerDesc  tibble with at minimum columns Marker_name, Cell_type; markers where Cell_type == 'X'
#'                    will be included in values under CTMarkers
#'
#' @return  original data table with PositiveMarkers and CTMarkers columns added
addCellPositiveMarkers <- function(dat, markerDesc){

    markers   <- getAllMarkers(markerDesc = markerDesc)
    idMarkers <- getIdentityMarkers(markerDesc)

    if(!any(markers %in% names(dat)) && all(c("Marker", "Value") %in% names(dat))){
        dat <- dat %>% spread(Marker, Value)
    }

    if(!"PositiveMarkers" %in% names(dat)){
        dat$PositiveMarkers <- apply(dat[,which(names(dat) %in% markers)], 1,
                                      function(x){
                                         paste(sort(names(x)[x==1]), collapse=",")
                                      })

        dat$IDMarkers <- unlist(lapply(dat$PositiveMarkers,
                                         function(x){
                                             x <- unlist(strsplit(x,","))
                                             paste0(x[x %in% idMarkers], collapse=",")
                                       }))
    }

    dat

}

addCellNegativeMarkers <- function(dat, markerDesc){

    markers   <- getAllMarkers(markerDesc)
    idMarkers <- getIdentityMarkers(markerDesc)

    if(!any(markers %in% names(dat)) && all(c("Marker", "Value") %in% names(dat))){
        dat <- dat %>% spread(Marker, Value)
    }

    if(!"NegativeMarkers" %in% names(dat)){
        dat$NegativeMarkers <- apply(dat[,which(names(dat) %in% markers)], 1,
                                      function(x){
                                         paste(sort(names(x)[x==0]), collapse=",")
                                      })

        dat$NegativeIDMarkers <- unlist(lapply(1:nrow(dat), 
                                         function(x){
                                             neg <- unlist(strsplit(dat$NegativeMarkers[x],","))
                                             msng <- unlist(strsplit(dat$MissingMarkers[x],","))               
                                             paste0(neg[which(neg %in% idMarkers & !neg %in% msng)], collapse=",")
                                       }))
    }
    dat
}

#' Throw error if attempting to assign one or more classifications to 
#' cells that already contain conflicting classes
#'
#' Given a vector of existing cell class assignments (and/or NAs), if any classes
#' exist and are different from the ones to be added, throw an error 
#'
#' @param dat      vector of cell classifications to be checked for conflicts
#' @param classes  a vector of classes named by their classification type (Category, 
#'                 Cell_type, Subtype, etc.) to be assigned to cells
#'
#' @return nothing
cellClassConflicts <- function(dat, classes){
    for(cl in names(classes)){  ## Category Cell_type Subtype Tag
        if(any(!is.na(dat[[cl]]))){
            vals <- unique(dat[[cl]])
            if(vals[!is.na(vals)] != classes[cl]){
                log_error(paste("Conflict betwen", cl, vals[!is.na(vals)], "and", classes[cl]))
                stop(paste("Conflict betwen", cl, vals[!is.na(vals)], "and", classes[cl]))
            }
        }
    }
}


#' Set all class columns in data table to NA
#' 
#' If class columns exist in data, make sure all data
#' in those columns are removed. If they do not exist, initialize
#' them to NA
#'
#' @param  dat         cell-level data table
#' @param  classTypes  column names of all classification types to be assigned
#'                     (generally Category, Cell_type, Subtype, Tag, Classifiers)
#'
#' @return data table with empty columns for all class types
resetDataClassification <- function(dat, classTypes){
    for(type in classTypes){
        dat[[type]] <- as.character(NA)
    }
    dat
}

#' Get all marker combinations that are not sufficient to identify a cell type
#' 
#' Cells with a postive marker combination that contains no cell type markers and
#' no identity markers are classified as "Unknown super negative" cells. Pull out
#' the unique set of these combinations from halo data
#' 
#' @param dat       cell level halo data including columns CTMarkers and PositiveMarkers
#' @param idMarkers vector of Identity markers pulled from marker descriptions xlsx file
#' 
#' @return vector of unique marker combinations existing in data that are not sufficient
#'         to identify the cell as a certain type
getAllUnknownSupernegCombos <- function(dat, idMarkers){
    dat %>%
    filter(IDMarkers == "", !grepl(getClassifierPattern(idMarkers, delim=","), PositiveMarkers)) %>%
    select(PositiveMarkers) %>% 
    unique() %>% 
    pull()
}


addCellTypes <- function(dat, cellTypes){
    cts <- lapply(1:nrow(cellTypes), function(x){ 
             log_debug(paste0(cellTypes[x,c("Category", "Cell_type", "Subtype")], collapse=";"))
             mrkrs <- cellTypes[x,] %>% 
                      select_if(~sum(!is.na(.)) > 0) %>%
                      select(!dplyr::matches("Category|Cell_type|Subtype")) %>%
                      names
             if(!all(mrkrs %in% names(dat))){
                 ## for superNeg cell type, assign if all identity markers existing in data are negative
                 ## (ignore missing markers in this case) 
                 if(!cellTypes[x,"Cell_type"] == "superNeg"){
                     log_debug("Skipping cell type [", x,"]! No data for required marker: ", setdiff(mrkrs, names(dat)))
                     return(NULL) 
                 }
             } 
             cells <- cellTypes[x,] %>%
                      select_if(~sum(!is.na(.)) > 0) %>%
                      left_join(dat, by = intersect(names(.), names(dat))) 
           }) %>%
           bind_rows %>%
           unique()       ## need to unique here because of cell types with multiple definitions
    if(nrow(cts) == 0){
        log_error("NO cell types assigned!!!")
        stop()   
    }
    unk <- dat %>% 
           filter(!UUID %in% cts$UUID) %>% 
           mutate(Category = "Other", Cell_type = "UNKNOWN", Subtype = "UNKNOWN")
    bind_rows(cts, unk)
}


#' Check that all group names exist in a set of column names 
#' 
#' Throw an error if not all column names to be used as grouping
#' variables exist in set of column names pulled from data
#'
#' @param cols     vector of tibble column names
#' @param groupBy  groupBy vector of group names
#'
#' @return nothing; throw error if not all group names are in col names
checkGrouping <- function(cols, groupBy){
    if(!all(groupBy %in% cols)){
        stop(paste0("Invalid groups: ", groupBy))
    }
}


#' Get vector of UUIDs for cells that fall into a certain class
#' 
#' Get vector of UUIDs for cells that fall into a certain class; this is
#' especially helpful for ensuring that there is no overlap/double counting when
#' comparing two groups of cells
#'
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param class           class/cell type to find
#'
#' @return vector of UUIDs of cells in class
getClassUUIDs <- function(dat, class){
    dat %>% 
    filter(grepl(getClassifierPattern(class), Classifiers)) %>%
    pull(UUID)
}

#' Build regex pattern needed to properly pull cell types out of Classifiers column
#' 
#' Build regex pattern needed to properly pull cell types out of Classifiers column; this
#' is to ensure special characters in the class name are escaped and that the class and 
#' ONLY the class are being counted
#' 
#' @param c   character string of class being searched for
#'
#' @return  a character string to be used when grepping for class in Classifiers column of data
getClassifierPattern <- function(c, delim=";"){

    pat <- paste(c(paste0("^",c,delim),paste0(delim,c,delim),paste0(delim,c,"$"),paste0("^",c,"$")), collapse="|")
    pat <- gsub("\\+","\\\\+",pat)
    pat <- gsub("\\(","\\\\(",pat)
    pat <- gsub("\\)","\\\\)",pat)

    return(pat)
}

#' Check that there are no overlaps in UUIDs of cells in multiple groups
#' 
#' In order to ensure there is no double counting of cells, check for overlaps
#' in lists of UUIDs in each group being counted
#' 
#' @param UUIDlist    list containing vectors of UUIDs, each representing a cell group
#'                    that should have no intersections with any other vectors in the list
#' 
#' @return list where each element is a pair of vectors that contains one or more overlapping cells
cellGroupOverlaps <- function(UUIDlist){

    allOverlaps <- list()

    for(x in 1:(length(UUIDlist)-1)){ 
        s1 <- UUIDlist[[x]]
        tmp <- UUIDlist[-x]
        for(y in 1:length(tmp)){
            s2 <- tmp[[y]]
            overlap <- intersect(s1,s2)
            if(!is.null(overlap) && length(overlap) > 0){
                if(!is.null(names(UUIDlist))){
                    allOverlaps[[length(allOverlaps) + 1]] <- c(names(UUIDlist)[x], names(tmp)[y])
                } else {
                    idx2 <- which(unlist(lapply(UUIDlist, function(lst){ identical(lst, s2) })))
                    allOverlaps[[length(allOverlaps) + 1]] <- c(x, idx2) 
                }
            }
        }
    }

    allOverlaps
}


#' Get vector of UUIDs that match one or more class/cell type
#' 
#' Get vector of UUIDs that match one or more class/cell type
#' 
#' @param dat        table of halo data where a row represents a single cell and contains,
#'                   at minimum, columns UUID and Classifiers
#' @param cellTypes  vector of one or more cell types to search for in Classifiers column
#'
#' @return  vector of all UUIDs that match one or more value in {cellTypes}
uuidsAnyClass <- function(dat, cellTypes){
    fullPat <- getClassifierPattern(cellTypes[1])
    if(length(cellTypes) > 1){
        for(ct in cellTypes[2:length(cellTypes)]){
            fullPat <- paste(fullPat, getClassifierPattern(ct), sep="|")       
        }
    }    
    dat %>% filter(grepl(fullPat, Classifiers)) %>% pull(UUID)
} 

#' Get vector of UUIDs that match ALL classes/cell types
#' 
#' Get vector of UUIDs that match ALL classes/cell types
#' 
#' @param dat        table of halo data where a row represents a single cell and contains,
#'                   at minimum, columns UUID and Classifiers
#' @param cellTypes  vector of one or more cell types to search for in Classifiers column
#'
#' @return  vector of all UUIDs that match ALL values in {cellTypes}
uuidsAllClasses <- function(dat, cellTypes){
    tmp <- dat    
    for(ct in cellTypes){
        tmp <- tmp %>% filter(grepl(getClassifierPattern(ct), Classifiers))
    }
    if(nrow(tmp) > 0){
        return(tmp$UUID)
    } else {
        return(NULL)
    }
}

#' Separate Classifiers column into individual category columns according to *CellTypes.xlsx file
#' 
#' Separate Classifiers column in annotated cells file into columns described
#' in cell types table
#' 
#' @param annDat    table of annotated cell data where each row represents one unique cell and 
#'                  contains a Classifiers column containing all classifiers assigned to a cell, 
#'                  delimited by semi-colons
#' @param cellTypes parsed and expanded table version of *CellTypes.xlsx
#' @param classCols vector of columns from cellTypes table; default=c("Category", "Cell_type", "Subtype", "Tag")
#'
#' @return  annDat table with columns for Category, CellType, Subtype and Tag
separateClassifiers <- function(annDat, cellTypes, classCols=NULL){

    classCols <- c("Category","Cell_type","Subtype")

    for(col in classCols){
        annDat[[col]] <- NA
        cts <- unique(cellTypes[[col]])
        for(ct in cts){
            uuids <- uuidsAllClasses(annDat, ct)
            annDat[[col]][annDat$UUID %in% uuids] <- ct
        }
    }
    annDat

}

#' Get all possible combinations of markers
#'
#' Given a vector of single markers, return a list
#' of all possible combinations of those markers,
#' including both positive and negative variations
#' 
#' @param markers  vector of single markers
#' @param superNeg include combination of all negative markers; default = TRUE
#'
#' @return  list of all combinations
getAllCombos <- function(markers, superNeg = TRUE){

    all.combos <- list()

    ## get all combinations of markers
    for(x in 1:length(markers)){
        ## get all combos of x
        combos = combn(markers,x,simplify=FALSE)
        all.combos = c(all.combos, combos)
    }

    for(c in 1:length(all.combos)){
        combo = all.combos[[c]]
        if(length(combo) < length(markers)){
            full <- c(combo,paste0(setdiff(markers,combo),"-"))
            srtd <- c()
            for(m in 1:length(markers)){
                srtd <- c(srtd, full[full %in% c(markers[m], paste0(markers[m], "-"))])
            }
            all.combos[[c]] <- srtd
        }
        all.combos[[c]] <- paste0(all.combos[[c]],collapse=",")
    }

    if(superNeg){
        ## add all neg combo 
        return(c(list(paste0(markers, "-", collapse = ",")), all.combos))
    }
    all.combos
}

#' Remove any markers determined to be problematic after analysis
#' 
#' Most exclusions are determined prior to rethresholding, but this function
#' allows users to exclude markers determined to be problematic during
#' QC or any other analysis step.
#'
#' @param dat    data tibble loaded from halo object analysis RDA files, joined
#'               with all sample and FOV IDs
#' @param annot  flat annotation data returned from loadStudyAnnotations(),
#'               including columns FOV_exclusion (see docs for 
#'               details on format)
#' @param return data tibble with any additional or markers removed 
postAnalysisMarkerExclusions <- function(dat, annot){
   
    if(all(is.na(annot$Marker_exclusion) | annot$Marker_exclusion == "")){ return(dat) }

    ## filter out markers
    mx <- lapply(annot$Marker_exclusion[!is.na(annot$Marker_exclusion) & annot$Marker_exclusion != ""], 
                 function(x) {
                   length(unlist(strsplit(x, ",")))
                 }
                ) %>% unlist %>% max
    sepCols <- paste0("Marker.", seq(1:mx))

    mExcl <- annot %>%
             filter(!is.na(Marker_exclusion), Marker_exclusion != "") %>%
             select(CellDive_ID, FOV_number, Marker_exclusion) %>%
             separate(Marker_exclusion, sepCols, sep=",") %>%
             gather(all_of(sepCols), key = "tmp", value = Marker) %>%
             select(-tmp) %>%
             filter(!is.na(Marker)) %>%
             unique %>%
             mutate(REMOVE = "YES", Marker = gsub(" ", "", Marker))

    tmp  <- dat %>% left_join(mExcl, by = intersect(names(.), names(mExcl)))

    ## log what is to be removed
    excl <- tmp %>%
            filter(REMOVE == "YES") %>%
            group_by(CellDive_ID, FOV_number, Marker) %>%
            summarize(Count = n()) %>%
            mutate(LOG = paste("Removed marker [", Marker, "] from", Count, 
                               "cells in FOV", FOV_number, "of sample", CellDive_ID))
    sapply(1:nrow(excl), function(x) log_debug(excl$LOG[x]))

    tmp %>%
    filter(is.na(REMOVE)) %>%
    select(-REMOVE)
}

#' Remove any FOVs determined to be problematic after analysis
#' 
#' Most exclusions are determined prior to rethresholding, but this function
#' allows users to exclude FOVs determined to be problematic during
#' QC or any other analysis step.
#'
#' @param dat    data tibble loaded from halo object analysis RDA files
#' @param annot  flat annotation data returned from loadStudyAnnotations(),
#'               including columns FOV_exclusion (see docs for 
#'               details on format)
#' @param return data tibble with any additional FOVs removed 
postAnalysisFOVexclusions <- function(dat, annot){

    totalCellsRemoved <- 0
    excl <- dat

    ## filter out FOVs
    metaExcl <- annot %>%
                filter(!is.na(FOV_exclusion) | FOV_exclusion != "") %>%
                unique() %>% 
                select(CellDive_ID, FOV_number)

    if(nrow(metaExcl) > 0 ){
        log_info("FOV(s) marked for exclusion in meta data: [ ", paste(unique(metaExcl$FOV_number), collapse = ","), " ]. EXCLUDING!")
    }

    fExcl <- dat %>% 
             select(CellDive_ID, FOV_number) %>%
             left_join(annot %>% select(CellDive_ID, FOV_number, FOV_ID, FOV_exclusion)) %>%
             filter(is.na(FOV_ID)) %>%
             unique() %>%
             select(CellDive_ID, FOV_number)

    if(nrow(fExcl) > 0){
        log_warn("FOV(s) still in data but MISSING from meta data: [ ", paste(unique(fExcl$FOV_number), collapse = ","), " ]. EXCLUDING!") 
    }

    fExcl <- bind_rows(fExcl, metaExcl)

    if(nrow(fExcl) != 0){
        for(x in 1:nrow(fExcl)){
            rmv <- excl %>% 
                   filter(CellDive_ID == fExcl$CellDive_ID[x],
                          FOV_number == fExcl$FOV_number[x]) %>%
                   pull(UUID)
   
            log_debug(paste("Removing [", length(rmv), "] cells",
                            "from FOV", fExcl$FOV_number[x], 
                             "of sample", fExcl$CellDive_ID[x]))  
        
            totalCellsRemoved <- totalCellsRemoved + length(rmv)

            excl <- excl %>% filter(!UUID %in% rmv)
        }
    }
    log_debug(paste0("Removed ", nrow(fExcl), " entire FOVs for a total of ", totalCellsRemoved, " cells"))
    
    excl
}

#' Get markers missing in a particular sample
#' 
#' Get vector of marker names that are known to be missing from a sample
#'
#' @param sampAnn    sample annotation table containing columns
#'                   'Final_missing_markers_functional' and 'Final_missing_markers_identity'
#' @param cellDiveID Celldive ID of sample to check for missing markers
#'
#' @return vector of all markers missing in sample
getMissingMarkersInSample <- function(sampAnn, cellDiveID){
    msng <- sampAnn %>%
            filter(CellDive_ID == cellDiveID) %>%
            select(Final_missing_markers_functional, Final_missing_markers_identity) %>%
            unique %>%
            select_if(~sum(!is.na(.)) > 0)
    if(length(msng) == 0 ){
        return(NULL)
    }
    msng %>% unlist %>% strsplit(",") %>% unlist(use.names=F)
}


getMissingMarkers <- function(sampAnn){
    sampAnn %>% 
    select(CellDive_ID, dplyr::matches("missing_markers")) %>%
    unique() %>%
    gather(2:ncol(.), key = "X", val = "MissingMarkers") %>% 
    select(-X) %>% 
    group_by(CellDive_ID) %>% 
    filter(!is.na(MissingMarkers)) %>% 
    summarize(MissingMarkers = paste0(MissingMarkers, collapse = ","))
}

getExcludedMarkers <- function(sampAnn){
    maxExcl <- lapply(sampAnn$Marker_exclusion, function(x){
                   unlist(strsplit(x, ",")) %>% length
               }) %>% unlist %>% max
    suppressWarnings(sampAnn %>% 
                     select(CellDive_ID, FOV_ID, Marker_exclusion) %>%
                     unique %>%
                     separate(Marker_exclusion, paste0("mExcl.", seq(1:maxExcl)), sep=",", remove=F))
}

actualMissingIDMarkers <- function(vals, idMarkers){
    paste(idMarkers[which(is.na(vals))], collapse = ",")
}

unexpectedMissingMarkers <- function(actualMissing, expectedMissing){
    lapply(1:length(actualMissing), function(x){
        paste(setdiff(unlist(strsplit(actualMissing[x], ",")), 
                      unlist(strsplit(expectedMissing[x], ","))), 
              collapse = ",")
    }) %>%
    unlist
}

annotationComplete <- function(dat, idMarkers){

    ## get actual missing markers (all NA or not in data at all?)
    dat$ActualMissing <- dat %>% 
                         select(all_of(idMarkers)) %>% 
                         apply(., 1, actualMissingIDMarkers, idMarkers) 
    dat <- dat %>%
           mutate(UnexpectedMissing = unexpectedMissingMarkers(ActualMissing, MissingMarkers))

    if(any(is.na(dat$UnexpectedMissing)) || any(dat$UnexpectedMissing != "")){
        log_error("Found unexpected missing marker(s)!!! [ ", unique(dat$UnexpectedMissing), " ]")
        return(FALSE) 
    }
    TRUE
}


#' Add to cell-level tibble columns for positive markers and cell classifications
#' 
#' Given a tibble complete table of marker positivity for each cell, add columns
#' consolidating list of positive markers in each cell, positive cell identity markers
#' for each cell, and assign cell classes/types based on those positive markers
#' 
#' @param annotatedCellsFile   file that either already contains a cell-level tibble of
#'                             Halo data or the file to which said table will be saved
#' @param dataDir              directory of Halo files to be annotated; when NOT NULL, 
#'                             ALL RDA files in this directory will be loaded and annotated
#' @param dataFiles            vector of RDA files to be annotated; use EITHER dataDir or dataFiles
#' @param metaFiles            vector of meta data XLSX files, including sample annotation
#'                             marker information and cell type definitions files (default = NULL)
#' @param metaDataFile         RDA file of pre-compiled meta data (default = NULL)
#' @param numThreads           integer; number of threads
#' @param filterExclusions     logical indicating whether to remove cells marked with any text
#'                             in 'EXCLUDE' column; default = FALSE
#' @param controlMarker        marker whose negativity indicated the cell is not usable; these
#'                             cells are removed
#'
#' @return all annotated data
annotateCells <- function(annotatedCellsFile, dataDir = NULL, dataFiles = NULL, 
                          metaFiles = NULL, numThreads = 1, 
                          filterExclusions = TRUE, controlMarker = "DAPI", 
                          forceReannotation = FALSE){

    log_debug("Reading cell type and marker files")
    cellTypes  <- getMetaFile("CellTypes", metaFiles) %>%
                  getCellTypes() %>% 
                  reformatCellTypes()
    markerDesc <- getMetaFile("Markers", metaFiles) %>% 
                  read.xlsx(1, check.names = F) %>% as_tibble()
    idMarkers  <- getIdentityMarkers(markerDesc)

    ###
    ### load data
    ###
    if(fileDone(annotatedCellsFile) && !forceReannotation){
        log_info(annotatedCellsFile, "already exists. Returning existing data.")
        return(readRDS(annotatedCellsFile))   
    }
    annot <- loadStudyAnnotations(metaFiles = metaFiles)
    missingMarkers <- getMissingMarkers(annot$flat)

    log_debug("Loading object analysis data...")
    dataFiles <- getFiles(path = dataDir, files = dataFiles, pattern=".rda")
    nThreads <- min(length(dataFiles), numThreads) 

    allDat <- loadAllHaloData(dataFiles = dataFiles,
                              nThreads = nThreads,
                              filterExclusions = filterExclusions, 
                              controlMarker = controlMarker) %>%
              joinIDs(annot$IDs) %>%
              left_join(missingMarkers, by = "CellDive_ID") %>%
              spread(Marker, Value)

    log_debug(paste0("Starting with [ ", formatC(nrow(allDat), big.mark = ","), " ] ", controlMarker, "+ cells"))

    ### 
    ###  add annotation columns including pos/neg markers and cell type assignments
    ###
    annDat <- allDat %>%
              addCellTypes(cellTypes) %>%
              filter(!is.na(UUID)) %>%
              addCellPositiveMarkers(markerDesc) %>%
              addCellNegativeMarkers(markerDesc) %>%
              unite("Classifiers", c("Category", "Cell_type", "Subtype"), sep=";", remove = F)

    log_debug("Verifying annotation is complete.")
    if(!annotationComplete(annDat, idMarkers)){
        return(NULL)
    }

    ###
    ### organize and select needed columns
    ###
    annDat <- annDat %>%
              select(UUID, CellDive_ID, Patient_ID, Sample_ID, FOV_ID, 
                     XMin, XMax, YMin, YMax,
                     PositiveMarkers, NegativeMarkers, MissingMarkers, 
                     IDMarkers, NegativeIDMarkers, 
                     Category, Cell_type, Subtype, Classifiers, everything()) %>%
              select(!dplyr::matches("Exclude"))

    if(nrow(annDat) != nrow(allDat)){ 
        log_error("Resulting annotated data has different number of cells than original data!")
        log_error("Original data: ", nrow(allDat), "cells;  Annotated data: ", nrow(annDat), "cells.")
        stop()
    }

    log_info(paste0("Saving newly annotated cells in file: ", annotatedCellsFile))
    ## save data
    saveRDS(annDat, annotatedCellsFile)

    annDat

}

#' Determine whether a cell meets the definition of a specific phenotype
#' 
#' Given a comma-delimited string representing a cell phenotype, in the form 
#' [Cell Type/Subtype/Tag],[Marker1],[Marker2],[Marker3], and the cells' cell type
#' classifiers and positive markers, return a 1 if the cell fits the phenotype and
#' a 0 if it does not.
#' 
#' @param phenStr        comma-delimited string where the first element contains a cell
#'                       type classifier (Cell_type, Subtype, Category or Tag) and the 
#'                       remaining elements are any combination of positive and negative 
#'                       markers
#' @param classifierStr  semi-colon delimited string containing all classifiers that 
#'                       apply to a cell
#' @param posMarkerStr   comma-delimited string of positive and negative markers
#' @param nameMap        list of out-dated classifiers named by the current ones used
#'                       in phenStr
cellFitsPhenotype <- function(phenStr, classifierStr, posMarkerStr, nameMap = NULL){

    if(length(classifierStr) != length(posMarkerStr)){
        stop("Vectors of classifier strings and positive marker strings differ in length.")
    }

    class <- gsub(",.*", "", phenStr)
    if(class %in% names(nameMap)){ class <- nameMap[[class]] }

    mrkrs <- unlist(strsplit(phenStr, ","))[-1]
    pos <- mrkrs[!grepl("\\-", mrkrs)]
    neg <- mrkrs[grepl("\\-", mrkrs)]

    lapply(1:length(classifierStr), function(x){


        if(!grepl(getClassifierPattern(class), classifierStr[x])){
            return(0)
        }

        if(length(pos) > 0){
            for(pm in pos){
                if(!grepl(getClassifierPattern(pos, delim = ","), posMarkerStr[x])){
                    return(0)
                }
            }
        }

        if(length(neg) > 0){
            if(grepl(getClassifierPattern(neg, delim = ","), posMarkerStr[x])){
                return(0)
            }
        }

        1

    }) %>% unlist()
}


