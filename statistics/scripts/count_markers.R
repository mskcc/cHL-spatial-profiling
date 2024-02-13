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
            --counts_dir              output directory 
            --meta_dir                path to meta files in XLSX format, required IF
                                      meta_data_file is NULL

          [OPTIONAL]
            --manifest            YAML file containing one or more parameter; NOTE: arguments on command
                                  line override manifest arguments!!!
            --meta_files          comma-delimited list of meta files; can be supplied instead of meta_dir
            --number_threads      number of threads to use for parallel processes
        \n"
    )
}

## set up required args & defaults 
minReq   <- list(c("meta_dir","meta_files"), "cell_data_dir", "counts_dir")
defaults <- list(number_threads = 1)

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(args, names(args))
threads   <- min(args$number_threads, detectCores() - 2)

#########

stDat       <- loadStudyData(args, annotatedCells = FALSE)
annFiles    <- dir(args$cell_data_dir, full.names = T)
funcMarkers <- read.xlsx(getMetaFile("Markers", metaFiles = stDat$metaFiles), 1, check.names = F) %>%
               getFunctionalMarkers()

missing <- checkMissingMarkers(annFiles, stDat$sampAnn, stDat$markers, threads = threads) 
if(nrow(missing) > 0){
    print(missing %>% filter(ERROR))
    #stop()
}

formatMarkerCounts <- function(mCounts, groupBy = "CellDive_ID"){

    grpN <- paste0("N.", paste(groupBy, collapse = "__"))
    mrkrN <- "N.Markers"
    samplesPerMarker <- mCounts %>% group_by(Marker) %>% summarize(!!as.name(grpN) := n())
    markersPerSample <- mCounts %>% group_by_at(groupBy) %>% summarize(!!as.name(mrkrN) := n())

    mc <- mCounts %>%
          ungroup() %>%
          gather(c("Marker Positive Cells", "PCT"), key = "ValType", value= "Val") %>% 
          mutate(ValType = ifelse(ValType == "Marker Positive Cells", "N", "PCT"), 
                 Col = paste(Marker, ValType, sep=".")) %>%
          select_at(c(groupBy, "Total Cell Count", "Col", "Val")) %>%
          spread(Col, Val) %>%
          left_join(markersPerSample) %>%
          arrange(desc(`Total Cell Count`)) %>%
          select(groupBy, `Total Cell Count`, mrkrN, everything())

    res <- list(mc, samplesPerMarker)
    names(res) <- c(paste0("marker_counts_by_", paste(groupBy[length(groupBy)], collapse = "__")), 
                    paste0("marker_summary_by_", paste(groupBy[length(groupBy)], collapse = "__")))
    res
}


### get counts for ALL markers out of all cells
all <- getAllMarkerCounts(unique(stDat$sampAnn$CellDive_ID), annFiles, stDat$markers, 
                          knownMissing = missing, threads = args$number_threads, groupBy = "CellDive_ID") %>%
       formatMarkerCounts(groupBy = "CellDive_ID")

write.xlsx(all, file.path(args$counts_dir, "all_marker_counts_per_sample.xlsx"), check.names = F)


### ALL MARKER COUNTS IN TUMOR CELLS
filt <- list('Cell_type' = 'HRS')

## by sample
tum <- getAllMarkerCounts(unique(stDat$sampAnn$CellDive_ID), annFiles, funcMarkers, 
                          datFilters = filt, knownMissing = missing, 
                          threads = args$number_threads, groupBy = "CellDive_ID") %>%
       formatMarkerCounts(groupBy = "CellDive_ID")
tum[[1]] <- tum[[1]] %>% dplyr::rename(`HRS Cells` = `Total Cell Count`)

## by FOV
tumFOV <- getAllMarkerCounts(unique(stDat$sampAnn$CellDive_ID), annFiles, funcMarkers,       
                          datFilters = filt, knownMissing = missing,
                          threads = args$number_threads, groupBy = c("CellDive_ID", "FOV_ID")) %>%
          formatMarkerCounts(groupBy = c("CellDive_ID", "FOV_ID"))
tumFOV[[1]] <- tumFOV[[1]] %>% dplyr::rename(`HRS Cells` = `Total Cell Count`)

tumRes <- c(tum, tumFOV)
write.xlsx(tumRes, file.path(args$counts_dir, "tumor_functional_marker_counts_per_sample.xlsx"), check.names = F)

