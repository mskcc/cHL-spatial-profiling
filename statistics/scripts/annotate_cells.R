suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript annotate_cells.R 
            
          [REQUIRED] 
            --data_file             path to processed, exclusion-marked RDA file of 
                                    formatted halo object data (Nick's RDAs)
            --meta_dir              path to meta files in XLSX format

          [OPTIONAL]
            --annotated_cells_file output RDA file; default = 'annotated_cells.rda' 
            --meta_files           comma-delimited list of meta files; default = NULL
            --number_threads       number of threads to use for parallel processes; default = 4
            --flex_neg_req         logical; when TRUE, if a definition includes the negative
                                   requirement of a missing marker, that requirement will be 
                                   ignored; if a missing marker is a required positive,
                                   no cells will be assigned to that cell type; default = FALSE
        \n"
    )
}

minReq   <- list("data_file", "meta_dir")
defaults <- list(annotated_cells_file = "annotated_cells.rda",
                 flex_neg_req = FALSE, 
                 number_threads = 4, 
                 debug = TRUE)
allUsed  <- unique(c(unlist(minReq), names(defaults), "manifest"))

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

setLogThreshold(debug = args$debug)
logParams(args, allUsed) 


#### prep input data
metaFiles <- getCurrentMetaFiles(metaDir = args$meta_dir, metaFiles = args$meta_files)
ctFile <- getMetaFile("CellTypes", metaFiles = metaFiles)
mFile <- getMetaFile("Markers", metaFiles = metaFiles)

markerDesc <- read.xlsx(mFile, 1, check.names = F) %>% 
              as_tibble 
idMarkers <- getIdentityMarkers(markerDesc)
cellTypes <- read.xlsx(ctFile, 1, check.names = F) %>% 
             as_tibble %>%
             filter(!is.na(Cell_type))

study <- loadStudyAnnotations(metaFiles = metaFiles) 
missingMarkers <- getMissingMarkers(study$flat)


### TODO: put everything below in a function

#### start building annotated_cells file in the current format ####
aDat <- loadAllHaloData(dataFiles = args$data_file,                            
                        nThreads = args$number_threads,                              
                        filterExclusions = TRUE, 
                        controlMarker = "DAPI") %>%                
        joinIDs(study$IDs) %>%                                            
        left_join(missingMarkers, by = "CellDive_ID") %>%                 
        spread(Marker, Value) %>%
        addCellPositiveMarkers(markerDesc) %>%                               
        addCellNegativeMarkers(markerDesc)          


#### ASSIGNMENTS ###
dat <- readRDS(args$data_file) 
ctDat <- dat$geom.data %>% 
         filter(!Exclude) %>% 
         select(UUID, Sample, FOV = SPOT) %>%
         left_join(dat$marker.data %>%
                   filter(Marker %in% idMarkers), by = "UUID") %>%
         assign_cell_types_by_fov(cellTypes, idMarkers, flexible = args$flex_neg_req)

markerCols <- getAllMarkers(markerDesc = markerDesc)
annotated <- aDat %>% 
             left_join(ctDat, by = "UUID") %>%
             unite("Classifiers", c("Category", "Cell_type", "Subtype"), sep=";", remove = F) %>%
             select(-all_of(markerCols)) %>%
             select(UUID, CellDive_ID, Patient_ID, Sample_ID, FOV_ID,          
                     XMin, XMax, YMin, YMax,                                    
                     PositiveMarkers, NegativeMarkers, MissingMarkers,          
                     IDMarkers, NegativeIDMarkers,                              
                     Category, Cell_type, Subtype, Classifiers, everything()) %>%
              select(!dplyr::matches("Exclude"))    

expCells <- length(unique(dat$geom.data %>% filter(!Exclude) %>% pull(UUID)))
if(nrow(annotated) != expCells){                                           
    log_error("Resulting annotated data has different number of cells than original data!")
    log_error("Original data: ", nrow(allDat), "cells;  Annotated data: ", nrow(annDat), "cells.")

    quit(save = "no", status = 1, runLast = FALSE)
}

log_info(paste0("Saving annotated cells to ", args$annotated_cells_file))
saveRDS(annotated, args$annotated_cells_file)

log_debug("Annotated:")
log_debug(paste("  UUIDs:\t", nrow(annotated)))
log_debug(paste("  FOV_IDs:\t", length(unique(annotated$FOV_ID))))
log_debug(paste("  Sample_IDs:\t", length(unique(annotated$Sample_ID))))
log_info("Done.")
