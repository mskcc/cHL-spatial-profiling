#' Summarize counts of cell type marker combinations and their classifiers (cell types)
#' 
#' Count cell type marker combinations in a number of ways, by Sample and FOV for example. 
#' Output will include classification based on Classifiers column in dat and cell type class columns 
#' specified to be output explicitly in their own columns (e.g., Cell_type, Tag)
#' 
#' 
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param cellTypes       tibble of cell type assignments
#' @param classCols       column names from cellTypes tibble indicating which categories should be 
#'                        summarized; default=c("Category","Cell_type","Subtype","Tag")
#' @param summarizeBy     vector of ways in which data should be grouped; values must be equal to
#'                        column names in dat; default=c("Sample","FOV")
#' @param xlsxFile        name of XLSX file to which counts summaries should be written; if NULL, no
#'                        XLSX file will be written
#'
#' @return list of summary tables
getMarkerComboCounts <- function(dat, cellTypes, classCols=NULL, summarizeBy=NULL, xlsxFile=NULL){

    if(is.null(classCols)){ classCols <- c("Category", "Cell_type", "Subtype") }
    if(is.null(summarizeBy)){ summarizeBy <- c("CellDive_ID", "Sample_ID", "FOV_ID") }

    ctClasses <- cellTypes %>%
                 filter(Category != "Functional") %>%
                 select_at(classCols) %>%
                 unlist() %>%
                 unique()
    ctClasses <- ctClasses[!is.na(ctClasses)]

    ## remove any non-cell type classifiers from Classifiers column
    dat <- dat %>% unite("Cell type(s)", Category, Cell_type, Subtype, sep=";", remove = F)
    dat$IDMarkers[which(dat$IDMarkers == "")] <- "superNeg"

    #### get all unique cell type combinations with counts and classifiers
    ctsTotal <- dat %>%
                group_by(IDMarkers, `Cell type(s)`) %>%
                summarize(Total = n()) %>%
                arrange(desc(Total))

    sheetLst <- list()

    for(sb in summarizeBy){
        ## get counts for each markerPosTag+ctClassifier combo by sample
        sheetLst[[paste0("combos_by_",sb)]] <- ctsTotal %>%
                                               left_join(dat %>%
                                                         group_by_at(c(sb, "IDMarkers", "NegativeIDMarkers", "Cell type(s)")) %>%
                                                         summarize(Count = n()) %>%
                                                         spread(sb, Count),
                                                     by = c("IDMarkers", "Cell type(s)"))
    }
    return(sheetLst)
}

#' Count cells in each class/cell type
#'
#' Summarize cell counts by cell types and data level (e.g., all, Sample, FOV)
#' 
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param cellTypes       table of cell type definitions
#' @param classCols       column names from cellTypes tibble indicating which categories should be 
#'                        summarized; default=c("Category","Cell_type","Subtype","Tag")
#' @param summarizeBy     vector of ways in which data should be grouped; values must be equal to
#'                        column names in dat; default=c("Sample","FOV")
#'
#' @return list of tables, each containing a different summary of counts
getClassCountSummaries <- function(dat, cellTypes, classCols=c("Category", "Cell_type", "Subtype"), 
                                    summarizeBy=c("Sample_ID", "FOV_ID")){


    tbls <- list()
    for(sb in summarizeBy){
        for(cl in classCols){
            ## summarize standard classes (those found in individual class columns in data)
            sheet <- summarizeCellCounts(dat, summarizeBy = c(cl, sb)) %>%
                     spread(sb,CellCount) %>%
                     ungroup() %>%
                     mutate(Total = rowSums(.[,-1], na.rm=T)) %>%
                     select(all_of(cl), Total, everything()) %>%
                     filter_at(all_of(cl), all_vars(!is.na(.)))

            ## summarize additional misc classes that are only in Classifiers column
            addlVals <- cellTypes %>% 
                        filter_at(vars(all_of(cl)), all_vars(!. %in% sheet[[cl]])) %>%
                        pull(cl) %>%
                        unique()
            sheet <- sheet %>%
                     bind_rows(lapply(addlVals, function(x){
                                  summarizeClassCounts(dat, x, cl, summarizeBy = sb) 
                              })
                     ) 
            tbls[[paste0(cl,"_by_",sb)]] <- sheet 
        }
    }

    return(tbls)
}

#' Create a cell count summary
#'
#' Create a summary of cell counts for each group formed by a specified set of
#' columns
#' 
#' @param dat          cell level tibble (each row represents a single unique cell)
#'                     to be summarized
#' @param summarizeBy  vector of column names from dat on which data should be grouped
#'
#' @return a tibble of cell/row counts with one row for each combination of column 
#'         values
summarizeCellCounts <- function(dat, summarizeBy=NULL){
    
    checkGrouping(names(dat), summarizeBy)

    dat %>%
    group_by_at(summarizeBy) %>%
    summarize(CellCount = n()) 

}


#' Create a count summary for a single cell class
#' 
#' Filter data for a specific cell class, group by any one column, get count
#' for each value in that column & calculate the total cell count
#
#' @param dat          cell level tibble (each row represents a single unique cell)
#'                     to be summarized
#' @param class        cell class to summarize
#' @param classType    type of class (Cell_type|Subtype|etc)
#' @param summarizeBy  vector of column names from dat on which data should be grouped
#
#' @return a tibble with columns: Class, Total, and one column per value in 'summarizeBy'
summarizeClassCounts <- function(dat, class, classType, summarizeBy=NULL){

    checkGrouping(names(dat), summarizeBy)

    dat <- dat %>%
           getClassifierCounts(class, summarizeBy = summarizeBy) 

    if(!is.null(summarizeBy)){
        dat <- dat %>% spread(summarizeBy, Count) 
    }

    dat <- dat %>%
           mutate(Total = rowSums(.), Class = class) %>%
           select(Class, Total, everything()) %>%
           rename(!!classType := Class)
    dat
}

#' Get count(s) for a single marker
#' 
#' Get number of cells that are positive for a single marker
#' 
#' @param dat      table of cell level data where a row represents a single unique cell 
#'                 and at minimum, columns include any values in {groupBy} vector and PositiveMarkers
#' @param marker   marker to count
#' @param groupBy  vector of column names to group data by before counting; default=NULL; set to
#'                 NULL to get total count
#'
#' @return tibble of count(s) with one row per group
countMarkerPositiveCells <- function(dat, marker, groupBy = "CellDive_ID"){
    dat %>% 
    select_at(c("UUID", "PositiveMarkers", groupBy)) %>%
    mutate(Marker = marker, 
           MarkerPos = ifelse(grepl(getClassifierPattern(marker, delim=","), PositiveMarkers), 1, 0)) %>%
    group_by_at(groupBy) %>%
    mutate(`Total Cell Count` = n()) %>%
    group_by_at(c(groupBy, "Marker", "Total Cell Count")) %>%
    summarize(`Marker Positive Cells` = sum(MarkerPos)) %>%
    mutate(PCT = `Marker Positive Cells`/`Total Cell Count`)
}


#' Get table of individual marker counts and percentages
#'  
#' Given a table of annotated cell data, count total cells in each
#' group, the number of cells positive for each marker in each group, and
#' calculate the pertages of cells positive for each marker in each group
#'
#' @param dat      table of annotated cell data
#' @param markers  vector of marker names to count in column PositiveMarkers
#' @param groupBy  column(s) to use to identify groups for which cells should
#'                 be counted
#'
#' @return table of total cell counts, marker positive cell counts and marker
#'         positive cell percentages in each group
getMarkerCounts <- function(dat, markers, groupBy = "CellDive_ID"){
 print(dat$CellDive_ID[1])
    lapply(markers, function(m){
        if(!any(grepl(getClassifierPattern(m, delim=","), dat$PositiveMarkers)) &&
           !any(grepl(getClassifierPattern(m, delim=","), dat$NegativeMarkers))){
           msg <- paste0("ZERO data found for marker ", m, " in sample(s) ",  
                           paste(unique(dat$CellDive_ID), collapse=","), "!!!")
           log_warn(msg)
           return(dat %>% group_by_at(groupBy) %>% summarize(`Total Cell Count` = n()) %>%  
                  mutate(Marker = m, `Marker Positive Cells` = NA))
        }
        countMarkerPositiveCells(dat, m, groupBy = groupBy)
    }) %>%
    bind_rows()
}

getAllMarkerCounts <- function(samples, datFiles, markers, datFilters = NULL, 
                               knownMissing = NULL, threads = 1, groupBy = "CellDive_ID"){

    cl <- makeCluster(threads, type="FORK", outfile="")
    clusterExport(cl, c("datFiles", "markers", "knownMissing", "datFilters", "groupBy"), envir=environment())
    
    mc <- parLapply(cl, samples, function(smp){

              df <- datFiles[grepl(paste0("^", smp, "_"), basename(datFiles))]
              dat <- readRDS(df)

              oMissing <- NULL
              if(!is.null(knownMissing)){ 
                  oMissing <- knownMissing %>%
                              filter(CellDive_ID == smp) %>%
                              pull(Observed_missing) %>%
                              strsplit(",") %>% unlist
              }
             
              if(!is.null(datFilters)){
                  for(filt in names(datFilters)){
                      dat <- dat %>% 
                             filter_at(vars(all_of(filt)), all_vars(. %in% datFilters[[filt]]))
                  }
              }
              getMarkerCounts(dat, setdiff(markers, oMissing), groupBy = groupBy)
          }) %>%
          bind_rows
    stopCluster(cl)

    mc
}

#' Get table of individual marker counts and percentages for
#' a specific cell type or state
#'
#' Given a full data set, get positive marker counts and percentages
#' for a subset of the data
#'
#' @param dat           complete table of annotated cell data
#' @param conditionStr  comma-delimited condition string specifying a specific cell
#'                      type and/or state for which to count positive markers
#' @param markers       vector of marker names to count in column PositiveMarkers
#' @param groupBy       column(s) to use to identify groups for which cells should
#'                      be counted
#'
#' @return table of total cell counts, marker positive cell counts and marker
#'         positive cell percentages in each group of condition data
getConditionMarkerCounts <- function(dat, conditionStr, markers, groupBy = "CellDive_ID"){
    condDat <- dat %>% filterDataForCondition(conditionStr, markers) 
    if(nrow(condDat) == 0){ 
        grps <- dat %>% select(all_of(groupBy)) %>% unique
        grps <- lapply(seq(markers), function(x){ grps }) %>% bind_rows

        counts <- cbind(grps) %>%
                  tibble(Marker = rep(markers, nrow(unique(grps)))) %>%
                  mutate(TotalCells = 0, Count = NA, PCT = NA)
    } else {
        counts <- condDat %>%
                  getMarkerCounts(markers, groupBy = groupBy)
    }
    counts
}



#' Count and summarize cell types by looking at positive marker combinations
#' 
#' Determine the markers expressed in each cell and count and summarize the unique combinations of
#' expressed markers, cell types/subtypes, and cells that touch or do not touch macrophages
#'
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param cellTypes       tibble of cell type assignments
#' @param classCols       column names from cellTypes tibble indicating which categories should be 
#'                        summarized; default=c("Category","Cell_type","Subtype","Tag")
#' @param summarizeBy     vector of ways in which data should be grouped; values must be equal to
#'                        column names in dat; default=c("Sample","FOV")
#' @param xlsxFile        name of XLSX file to which counts summaries should be written; if NULL, no
#'                        XLSX file will be written
#' @param markerCombos    logical; when TRUE, all cell type marker combinations will be counted
#' @param classes         logical; when TRUE, counts will be summarized by classes
#' @param touchingMacros  logical; when TRUE, class counts will be repeated on the TM cells only 
#'
#' @return list of tibbles, each a different summary of counts
getCountSummaries <- function(dat, cellTypes, classCols=NULL, summarizeBy=NULL, xlsxFile=NULL,
                              markerCombos=TRUE, classes=TRUE, touchingMacros=FALSE){

    if(is.null(classCols)){ classCols <- c("Category", "Cell_type", "Subtype") }
    if(is.null(summarizeBy)){ summarizeBy <- c("CellDive_ID", "Sample_ID", "FOV_ID") }

    allSheets <- list()
   
    if(markerCombos){
        comboSheets <- getMarkerComboCounts(dat, cellTypes, classCols=classCols, summarizeBy=summarizeBy)
        allSheets <- c(allSheets, comboSheets)
    }
    if(classes){
        classSheets <- getClassCountSummaries(dat, cellTypes, classCols=classCols, summarizeBy=summarizeBy)
        allSheets <- c(allSheets, classSheets)
    
        if(touchingMacros){
            ## break TM category by Cell_type and by Subtype
            TMdat <- dat %>% filter(grepl(getClassifierPattern("TM"), Classifiers))
            TMclassSheets <- getClassCountSummaries(TMdat, cellTypes, classCols=c("Cell_type", "Subtype"), 
                                                    summarizeBy=summarizeBy)
            names(TMclassSheets) <- paste0("TM_",names(TMclassSheets))
            allSheets <- c(allSheets, TMclassSheets)
        }
    }
    if(!is.null(xlsxFile)){
        write.xlsx(allSheets, file=xlsxFile, check.names = F)
    }

    return(allSheets)

}

#' Summarize counts for a single cell type included in Classifiers column
#' 
#' Get counts for a single cell type for entire data set, or summarized by one or
#' more columns
#'
#' @param dat         table where column headers include Classifiers, and each row
#'                    represents a single cell
#' @param cellType    character string to pull from Classifiers column
#' @param summarizeBy name of column in dat by which counts should be grouped; default
#'                    is NULL, so counts returned will be for complete data set
#'
#' @return table of counts 
getClassifierCounts <- function(dat, cellType, summarizeBy=NULL){

    checkGrouping(names(dat), summarizeBy)

    dat %>% 
    filter(grepl(getClassifierPattern(cellType), Classifiers)) %>%
    group_by_at(summarizeBy) %>%
    summarize(Count=n())
}


