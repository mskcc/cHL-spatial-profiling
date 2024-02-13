suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))

log_threshold(DEBUG)

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript annotate_cells.R 

          [REQUIRED (may be defined on command line OR in manifest file)]
            --cell_data_dir           path to annotated cells RDA files that each 
                                      contain single table where rows are cells
                                      and columns are all data for that single 
                                      cell (one file per sample)
            --cell_type_counts_file   path to output file - XLSX file to contain counts summaries
            --meta_dir                path to meta files in XLSX format, required IF 
                                      meta_data_file is NULL

          [OPTIONAL]
            --flex_neg_req        logical; when TRUE, negative requirement of missing markers
                                  will be ignored; default = FALSE
            --manifest            YAML file containing one or more parameter; NOTE: arguments on command
                                  line override manifest arguments!!!
            --meta_files          comma-delimited list of meta files; can be supplied instead of meta_dir
            --number_threads      number of threads to use for parallel processes
        \n"
    )
}

## set up required args & defaults 
minReq   <- list(c("meta_dir","meta_files"), "cell_data_dir", "cell_type_counts_file")
defaults <- list(number_threads = 1, flex_neg_req = FALSE)
used     <- c(unlist(minReq), names(defaults), "manifest", "meta_files")
args     <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(args, used)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

summarizeSampleCellTypeCounts <- function(cdid, sampAnn, annDir, cellTypes){

    annCells <- loadAnnotatedCells(sampAnn, annDir = annDir, cellDiveID = cdid)

    if(is.null(annCells) | nrow(annCells) == 0){ 
        log_warn("No annotated cells found for sample ", cdid, ". Excluding sample from summaries.")
        return(NULL) 
    }

    getCountSummaries(annCells, cellTypes, summarizeBy = c("CellDive_ID", "FOV_ID"))
}

combineMarkerComboCounts <- function(allSumms, tblName = "combos_by_CellDive_ID", 
                                      missingMarkers = NULL){

    fin <- lapply(names(allSumms), function(x){
               if(is.null(allSumms[[x]])){ return(NULL) }
               allSumms[[x]][[tblName]] %>%
               mutate(`Cell type(s)` = ifelse(IDMarkers %in% c("", "superNeg"), 
                                               "superNeg", `Cell type(s)`)) %>%
               select(-Total, -NegativeIDMarkers) %>%
               gather(all_of(x), key = 'CellDive_ID', value = 'Count') %>%
               group_by(IDMarkers, `Cell type(s)`, CellDive_ID) %>%
               summarize(Count = sum(Count)) %>%
               ungroup
           }) %>%
           bind_rows %>%
           spread(CellDive_ID, Count, fill = 0) 
print(1)
    ## back fill NA combinations
    for(x in 1:nrow(missingMarkers)){
        cdid <- missingMarkers$CellDive_ID[x]
        if(!cdid %in% names(fin)){ next }
        mm   <- unlist(strsplit(missingMarkers$missing[x], ","))
        ptrn <- paste(paste0(mm,"(,|$)"), collapse = "|")
        idx  <- grep(ptrn, fin$IDMarkers)
        fin[[cdid]][idx] <- NA
    }
print(2)
    fin %>%
    mutate(Total = rowSums(select_if(., is.numeric), na.rm = T)) %>%            
    mutate(DF = rowSums(!is.na(select_if(., is.numeric))) - 2) %>%  ## subtract an extra 1 for 'Total' column
    arrange(desc(Total)) %>%                                                    
    select(any_of(c("IDMarkers", "Cell type(s)", "Cell_type", "Subtype", "Total", "DF")), everything())

}


getSampleNACellTypes <- function(cellTypes, markersMissingInSample, ctLevel = "Cell_type",
                                  flex_neg_req = FALSE){

    cellTypes %>%
    filter_at(vars(markersMissingInSample), any_vars(. == 1)) %>% #any_vars(!is.na(.))) %>%
    pull(ctLevel) %>%
    unique

}

fillInCellTypeNAs <- function(countsTbl, ctLevel, missingIDmarkers = NULL, cellTypes = NULL,
                                flex_neg_req = FALSE){

    if(is.null(missingIDmarkers) | is.null(cellTypes)){
        log_debug("No missing markers or cell types provided. No NA counts to be filled in.") 
        return(allSumms)
    }

    for(cdid in names(countsTbl)[-1]){
        mm <- missingIDmarkers %>% filter(CellDive_ID == cdid) %>% pull(missing)
        if(length(mm) == 0){ next }
        naCT <- getSampleNACellTypes(cellTypes, mm, ctLevel = ctLevel, flex_neg_req = flex_neg_req)
        idxs <- which(countsTbl[[ctLevel]] %in% naCT)
        countsTbl[[cdid]][idxs] <- NA
    }
    return(countsTbl)
}



combineCountSummaries <- function(allSumms, tblNames = NULL, missingIDmarkers = NULL, cellTypes = NULL,
                                    flex_neg_req = FALSE){

    if(is.null(tblNames)){
        tblNames <- lapply(allSumms, function(x){ names(x) }) %>% unlist %>% unique
    }

    fin <- list()
    for(x in tblNames){
        for(cdid in names(allSumms)){
            if(is.null(allSumms[[cdid]])){ next }
            tbl <- allSumms[[cdid]][[x]] %>% select(-Total)
            if(is.null(fin[[x]])){
                fin[[x]] <- tbl
            } else {
                fin[[x]] <- fin[[x]] %>% 
                            full_join(tbl, by = intersect(names(.), names(tbl)))  
            }
        }

        fin[[x]][is.na(fin[[x]])] <- 0 

        ctLevel = gsub("_by_.*", "", x)
        if(ctLevel %in% c("Cell_type", "Subtype")){
            fin[[x]] <- fillInCellTypeNAs(fin[[x]], ctLevel, 
                                          missingIDmarkers = missingIDmarkers, 
                                          cellTypes = cellTypes, 
                                          flex_neg_req = flex_neg_req)
        }

        fin[[x]] <- fin[[x]] %>%
                    mutate(Total = rowSums(select_if(., is.numeric), na.rm=TRUE),
                           DF = rowSums(!is.na(select_if(., is.numeric))) - 1, 
                           DF = ifelse(DF < 0, 0, DF)) %>%
                    arrange(desc(Total)) %>% 
                    select(any_of(c("IDMarkers", "Cell type(s)", "Cell_type", "Subtype", "Total", "DF")), everything())
    }

    fin
}

stDat   <- loadStudyData(args, annotatedCells = FALSE)
threads <- min(args$number_threads, detectCores() - 2)
cdids   <- as.list(unique(stDat$sampAnn$CellDive_ID))
names(cdids) <- cdids

cl <- makeCluster(threads, type="FORK", outfile="")
clusterExport(cl, c("args", "stDat"), envir=environment())
allSumms <- parLapply(cl, cdids, function(cdid){
              sampInData <- ifelse(length(grep(paste0("^", cdid, "_"), dir(args$cell_data_dir))) > 0, TRUE, FALSE)
              if(!sampInData){ return(NULL) } 
              tryCatch({
                  summarizeSampleCellTypeCounts(cdid, stDat$sampAnn, args$cell_data_dir, stDat$cellTypes)
                }, error = function(e){
                    log_error(e)
                    log_error("Summaries do NOT include sample ", cdid)
              })
            })
stopCluster(cl)

tbls <- c("combos_by_CellDive_ID", "Category_by_CellDive_ID",
          "Cell_type_by_CellDive_ID", "Subtype_by_CellDive_ID")

missingID <- stDat$sampAnn %>%
             filter(!is.na(Final_missing_markers_identity)) %>%
             select(CellDive_ID, Final_missing_markers_identity) %>%
             unique() %>%
             rename(missing = Final_missing_markers_identity)

final <- list()

final[["combos_by_CellDive_ID"]] <- combineMarkerComboCounts(allSumms, 
                                                             tblName = "combos_by_CellDive_ID", 
                                                             missingMarkers = missingID)

ctSumms <- combineCountSummaries(allSumms, 
                                 tblNames = tbls[-1], 
                                 missingIDmarkers = missingID, 
                                 cellTypes = stDat$cellTypes,
                                 flex_neg_req = args$flex_neg_req)

fin <- c(final, ctSumms)

write.xlsx(fin, args$cell_type_counts_file, check.names = F)

