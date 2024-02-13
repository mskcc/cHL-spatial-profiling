#' Set cell type definition table
#'
#' Read and parse cell types XLSX file and save table in tibble.
#' In case of error, catch and print error message, and warn that 
#' variable is not set.
#'
#' @param ctFile   XLSX file containing cell type definitions
#' 
#' @return nothing
setCellTypes <- function(ctFile){
    log_info("Setting cell type definition table...") 
    tryCatch({
        getCellTypes(ctFile)
      }, error = function(e){
        log_error(e)
        log_warn("Cell type definition table NOT set.")
    })
}

#' Set vector of all study markers 
#'
#' Read marker description file and save marker names in a vector. 
#' In case of error, catch and print error message, and warn that 
#' variable is not set.
#'
#' @param mFile   XLSX file containing all marker information, at least Marker_name
#'
#' @return nothing 
setMarkers <- function(mFile){
    log_info("Setting vector of all study markers...") 
    tryCatch({
        read.xlsx(mFile, 1, check.names = F) %>%
        pull(Marker_name)
      }, error = function(e){
        log_error(e)
        log_warn("Vector of study markers NOT set.")
     })
}


#' Set flat table of all sample annotation
#' 
#' Read, parse and flatten all sample meta files and save in  
#' table 'sampAnn'. If parameter 'id' is provided, filter annotation
#' on CellDive_ID for given value. In case of error, catch and print error 
#' message, and warn that variable is not set.
#'
#' @param   metaFiles   vector of all sample meta data files
#' @param   id          CellDive_ID of sample to keep; default = NULL (all samples kept)
#'
#' @return
setSampleAnnotation <- function(metaFiles, id = NULL){
    log_info("Setting flattened table of sample annotation...")
    tryCatch({
        sa <- loadStudyAnnotations(metaFiles)$flat %>%
              filter(is.na(FOV_exclusion))
        if(!is.null(id) && tolower(id) != "all"){
             sa <- sa %>% filter(CellDive_ID == id)
        }
        sa
      }, error = function(e){
        log_error(e)
        log_warn("Flattened table of sample annotation NOT set.")
    })
}

#' Set list of study questions
#'
#' Read, parse and reformat statistics questions. In case of error, 
#' catch and print error message, and warn that variable is not set.
#' 
#' @param  qFile       XLSX file describing all study questions
#' @param  question    QuestionNumber to return; default = NULL (all questions)
#' 
#' @return nothing
setQuestions <- function(qFile, question = NULL){
    log_info("Setting list of study questions...")
    tryCatch({
        qs <- parseStatsQuestions(qFile)
        if(!is.null(question) && tolower(question) != 'all'){ qs <- qs[names(qs) == question] }
        qs
    }, error = function(e){
        log_error(e)
        log_warn("Study questions not set.")
    })
}


#' Set variable 'conds'
#' 
#' Generate and index or read all previously indexed cell states/conditions
#'
#' @param cFile              XLSX file containing human-generated cell states/conditions
#' @param cIndex             XLSX file that does or will contain indexed version of cFile
#' @param arrangeAnnotation  list describing how conditions should be sorted
#'
#' @return nothing
setConditions <- function(cFile, cIndex, arrangeAnnotation){
    log_info("Setting table of indexed conditions...")
    tryCatch({
        getConditionsIndex(cIndex, cFile, arrangeAnnotation)
      }, error = function(e){
        log_error(e)
      log_warn("Table of indexed conditions NOT set.")
    })
}


#' Set variable 'analysisList'
#'
#' Divide indexed conditions into list according to calculation type
#'
#' @param cFile   XLSX file containing human-generated cell states/conditions
#' @param cTbl    table containing all indexed conditions
#'
#' @return nothing
setAnalyses <- function(cFile, cTbl){
    log_info("Setting list of analyses by calculation type...") 
    tryCatch({
        getStatsAnalyses(cFile, cTbl)
      }, error = function(e){
        log_warn("List of all analyses NOT set.")
    })
}


#' Set variable 'annCells'
#' 
#' Load all annotated cells tables and save in tibble 'annCells'. 
#' In case of error, catch and print error message and warn that variable is
#' not set.
#' 
#' @param cellDir     directory containing all annotated cells tables
#' @param fovAreaDir  directory containing FOV area files
#' @param studyAnn    flattened table of all study annotation
#' @param id          CellDive_ID of sample of interest; default = NULL, all samples
#'                    loaded
#'
#' @return
setAnnCells <- function(cellDir, fovAreaDir, studyAnn, id = NULL){
    log_info("Setting study annotated cells table...")
    tryCatch({
        loadAnnotatedCells(studyAnn, annDir = cellDir,
                           fovAreaDir = fovAreaDir,
                           cellDiveID = id)
      }, error = function(e){
        log_error(e)
        log_warn("Study table of annotated cells not set.")
    })
}


#' Set variable 'nbhdCounts'
#'
#' Read all macrophage neighborhood files and save counts in table
#' 'nbhdCounts'. In case of error, catch and print error message and warn
#' that variable is not set.
#'
#' @param nbhdDir        directory containing subfolders "C2" and "C3", each
#'                       containing files of cell to cell distances for 
#'                       cell pairs within 30 microns of each other
#' @param nbhdCountsDir  directory that already contains or will contain
#'                       neighborhood cell type counts for each center cell
#'                       of a specific type (e.g., counts of T cells within
#'                       30 microns of MHCIIpos_macro cells, per FOV)
#' @param cells          annotated cells table
#' @param id             CellDive_ID of sample to keep
#' @param analyses       list of all analyses
#' 
#' @return nothing
setNbhdCounts <- function(nbhdDir, nbhdCountsDir, cells, analyses,
                          id = NULL, threads = 1){

    log_info("Setting table of macrophage neighborhood counts...")
    tryCatch({
        nbhdDirs <- file.path(nbhdDir, c("C2", "C3"))
        loadNeighborhoodCounts(nbhdDirs, 
                               analyses, 
                               nbhdCountsDir, 
                               cells,
                               cellDiveID = id, 
                               numThreads = threads)
      }, error = function(e){
        log_error(e)
        log_warn("Table of neighborhood counts NOT set.")
    })
}

#' Set study variables 'tumorNbhd' and 'tumorNbhdCells'
#'
#' Read tumor neighborhood file and save all data in table 'tumorNbhd'. 
#' For convenience, save unique list of neighborhood cell UUIDs in  
#' vector 'tumorNbhdCells'. In case of error, catch and print error 
#' message and warn that variable is not set.
#' 
#' @param cellAnn    annotated cell table; must contain, at minimum, columns 'UUID'
#'                   and 'Cell_type'
#' @param nbhdFile   RDA file of cell to cell distances for 
#'                   cell pairs within {maxRadius} microns of each other
#' @param questions  study questions list, used to determine whether it is 
#'                   necessary to load tumor neighborhood data
#' @param cellType   name of tumor in Cell_type annotation column of {cellAnn}; 
#'                   default = 'Tumor'
#'
#' @return
setTumorNeighborhood <- function(cellAnn, nbhdFile, questions = NULL,
                                  cellType = "Tumor", maxRadius = 30){

    if(!is.null(questions)){
        if(!any(tolower(unlist(questions)) %in% c("nbhd", "neighborhood"))){
            log_warn("No questions are restricted to tumor neighborhood cell regions. Study variables 'tumorNbhd' and 'tumorNbhdCells' not set.")
            return()
        }
    }

    log_info("Setting study variables: tumorNbhd and tumorNbhdCells")
    tumorNbhd <- NULL
    tumorNbhdCells <- NULL

    tryCatch({

        log_info(paste0("Loading cells in tumor neighborhood from file ", nbhdFile))
        neighbors <- loadNeighbors(nbhdFile, maxRadius = maxRadius)
  print("loaded neighbors")      
        tumorCells <- cellAnn %>%
                      filter(Cell_type == cellType) %>%
                      pull(UUID)
print("got tumor cells")
        if(length(tumorCells) == 0){
            log_warn(paste0("No tumor cells found in sample(s) ", 
                            paste(unique(cellAnn$CellDive_ID), collapse = ", "), 
                            "!"))
        }

        tumorNbhd <- neighbors %>%
                     filter(N.UUID %in% tumorCells)
print("filtered neighbor table for tumor neighbors")
        tumorNbhdCells <- tumorNbhd %>% pull(C.UUID) %>% unique()
print("pulled out cells that are in tumor neighborhood")

      }, error = function(e){
        log_error(e)
        log_warn("Tumor neighborhood table and vector of cells in tumor neighborhood not set.")
    })

    list('tumorNbhd' = tumorNbhd, 'tumorNbhdCells' = tumorNbhdCells)        

}


#' Add tumor microenvironment assignments to annotated cells table
#' 
#' Read each TME file, and join to study variable 'annCells'. Depenging on
#' assignment level specified, prepend level to column name (e.g., cell_MHCI). 
#' In case of error, catch and print error message and warn that variable 
#' is not set.
#'
#' @param tmeDir           directory containing tumor microenvironment assignment
#'                         files, each containing one column named 'microEnv.*". 
#'                         Column in final table will be named according to assignment
#'                         level and microenvironment (e.g., cell_MHCI)
#' @param annCells         tibble of annotated cells to which TME column should be added
#' @param assignmentLevel  'sample' or 'cell', the level at which TME assignments were made
#' @param questions        study list of questions
#' @param filePattern      regex to use when identifying files to read in tmeDir
#'
#' @return nothing
addTME <- function(annCells, tmeDir, assignmentLevel = "cell", questions = NULL, filePattern = NULL){

    if(!is.null(questions)){
        colPtrn <- paste0(assignmentLevel, "_")
        if(!any(grepl(colPtrn, names(unlist(questions))))){
            log_warn("No questions involve TME on ", assignmentLevel, " level. No TME columns added to 'annCells'.")
            return(annCells)
        }
    }

    tryCatch({

        log_info(paste0("Adding ", assignmentLevel, " level TME assignments to study variable 'annCells'"))

        files <- dir(tmeDir, full.names = T, pattern = filePattern)
        if(gtools::invalid(files)){
            log_warn(paste0("No microenvironment files found for [", id, "]"))
            return(annCells)
        }

        for(fl in files){
            log_debug("Adding info from file: ", fl)
            tme <- readRDS(fl) %>% select(UUID = C.UUID, everything())
            meCol <- names(tme)[grepl("_microEnv", names(tme))]
            me <- gsub("_microEnv","",meCol)
            log_debug("  Adding column ", paste0(assignmentLevel, "_", me))

            annCells <- annCells %>%
                        left_join(tme %>%
                                  select(all_of(c("UUID", meCol))) %>%
                                  rename_at(all_of(meCol),
                                            list(~(. = paste0(assignmentLevel, "_",me)))),
                                  by = "UUID")
        }

        annCells

      }, error = function(e){
        log_error(e)
        log_warn(paste0("TME assignments on ", assignmentLevel, " level not added to 'annCells'"))
    })

}


#' Load all study data 
#'
#' Given just the study configuration in list form, load into a list all 
#' study data including parsed cell types, marker names, condition analysis list, 
#' sample annotation, annotated cell data, parsed study questions, neighborhood counts
#' and neighborhood distances. NOTE: this assumes all meta data have been processed,
#' cells have been annotated and metrics have been calculated prior to running
#' this function.
#'
#' @param config                  study configuration in list form (see docs for details)
#' @param cellTypes               logical; when TRUE, parse and expand cell types XLSX 
#'                                file and store in variable 'cellTypes'; default: TRUE
#' @param markerList              logical; when TRUE, pull Marker_name column from *Markers.xlsx 
#'                                file and store in variable 'markers'; default: TRUE  
#' @param analyses                logical; when TRUE, read stats conditions XLSX file into variable
#'                                'analysisList'; default: TRUE
#' @param conditions              logical; when TRUE, load conditions index into variable 'conds'; default: TRUE
#' @param sampleAnnotation        logical; when TRUE, load all sample annotations into variable 'sampAnn';
#'                                default: TRUE
#' @param annotatedCells          logical; when TRUE, load annotated cell data into variable 'annCells';
#'                                default: TRUE
#' @param questions               logical; when TRUE, parse stats questions XLSX file into list form and 
#'                                store in variable 'allQuestions'; default=FALSE
#' @param neighborhoodCounts      logical; when TRUE, read neighborhood counts file into variable 'nbhdCounts';
#'                                default: FALSE
#' @param neighborhoodDistances   logical; when TRUE, read neighborhood distances file into variable 
#'                                'nbhdDistances'; default: FALSE
#' @param cellsInTumorNeighborhood logical; when TRUE, read tumor neighborhood file; save table in 'tumorNbhd'
#                                  and save UUIDs of all neighborhood cells in vector 'tumorNbhdCells'
#' @param tmeSampleStatus         logical; when TRUE, add to 'annCells' table columns for all categorized tumor
#'                                microenvironments (e.g., MHCI+/- TME) available in directory referenced by
#'                                'tme_by_sample_dir' key in config
#' @param tmeCellStatus           logical; when TRUE, add to 'annCells' table columns for all categorized tumor
#'                                microenvironments (e.g., MHCI+/- TME) available in directory referenced by
#'                                'microenvironment_dir' key in config
#'
#' @return nothing
loadStudyData <- function(config,
                          all = FALSE,
                          cellTypeDefinitions = TRUE,
                          markerList = TRUE,
                          sampleAnnotation = TRUE,
                          annotatedCells = TRUE,
                          analyses = FALSE,
                          conditions = FALSE,
                          questions = FALSE,
                          neighborhoodCounts = FALSE,
                          cellsInTumorNeighborhood = FALSE,
                          tmeSampleStatus = FALSE,
                          tmeCellStatus = FALSE){

    if(all){
        cellTypeDefinitions <- markerList <- analyses <- conditions <- sampleAnnotation <- TRUE
        annotatedCells <- questions <- neighborhoodCounts <- cellsInTumorNeighborhood <- TRUE
        #tmeSampleStatus <- TRUE 
        tmeCellStatus <- TRUE
    }

    stDat <- list()

    metaFiles <- getCurrentMetaFiles(metaDir = config$meta_dir, metaFiles=config$meta_files)
    stDat$metaFiles <- metaFiles

    if(is.null(config$cell_dive_id)){
        config$cell_dive_id <- "all"
    }

    if(!is.null(config$debug) && (config$debug == "yes" || config$debug)){
        log_threshold(DEBUG)
    } else {
        log_threshold(INFO)
    }

    if(cellTypeDefinitions){
        ctFile <- getMetaFile("CellTypes", metaFiles)
        stDat$cellTypes <- setCellTypes(ctFile)
    }

    if(markerList){
        stDat$markers <- setMarkers(getMetaFile("Markers", metaFiles))
    }

    if(questions){
        stDat$allQuestions <- setQuestions(getMetaFile("Questions", metaFiles),
                                           question = config$question)
    }

    if(conditions || analyses){
        stDat$conds <- setConditions(config$statistics_conditions_file,
                                     config$statistics_conditions_index,
                                     config$arrange_annotation)
    }

    if(analyses){
        stDat$analysisList <- setAnalyses(config$statistics_conditions_file, stDat$conds)
    }

    if(sampleAnnotation || annotatedCells || neighborhoodCounts || tmeCellStatus){
        stDat$sampAnn <- setSampleAnnotation(metaFiles, id = config$cell_dive_id)
    }

#### EVERYTHING BELOW (INSIDE THIS FUNCTION) IS OLD AND IS NO LONGER USED
    if(annotatedCells || cellsInTumorNeighborhood || neighborhoodCounts || 
      tmeSampleStatus || tmeCellStatus){
        stDat$annCells <- setAnnCells(config$cell_data_dir,
                                      config$fov_area_dir,
                                      stDat$sampAnn,
                                      id = config$cell_dive_id)
        if(tmeCellStatus){
            stDat$annCells <- stDat$annCells %>%
                              addTME(config$microenvironment_dir,
                                     assignmentLevel = "cell",
                                     questions = stDat$allQuestions,
                                     pattern = paste0("^", config$cell_dive_id, "_"))
        }
    }

#    if(neighborhoodCounts){
#        nAnalyses <- stDat$analysisList[c("navgcounts", "nfracs")]
#        stDat$nbhdCounts <- setNbhdCounts(config$neighborhood_dir, 
#                                          config$neighborhood_counts_dir, 
#                                          cells = stDat$annCells,
#                                          nAnalyses,
#                                          id = config$cell_dive_id, 
#                                          threads = config$number_threads)
#    }

    if(cellsInTumorNeighborhood || tmeSampleStatus || tmeCellStatus){
        nbhdFile <- dir(config$neighbor_dir, full.names = T,
                        pattern = paste0("^", config$cell_dive_id, "___"))
        if(length(nbhdFile) == 0){
            log_warn(paste0("Could not file neighborhood file for sample ", config$cell_dive_id))
        } else {
            stDat <- c(stDat, 
                       setTumorNeighborhood(stDat$annCells,
                                            nbhdFile, 
                                            questions = stDat$allQuestions,
                                            cellType = "HRS",    ## move this to config eventually
                                            maxRadius = config$neighborhood_radius))
        }
    }

#    if(tmeSampleStatus){
#        stDat$annCells <- stDat$annCells %>%
#                          addTME(config$tme_by_sample_dir,
#                                 assignmentLevel = "sample",
#                                 questions = stDat$allQuestions)
#    }

    stDat
}

#' Load cell level data with all available cell and study annotations
#' 
#' Read pre-processed cell level data from one or more RDA files and join
#' additional information including sample annotation, band assignments
#' and areas
#' 
#' @param annFile      full path to a single RDA file of annotated cell data; 
#'                     when annFile is not NULL, annDir will be ignored; if
#'                     loading multiple files, use param 'annDir' instead
#' @param annDir       directory path containing one or more cell-level RDA files
#' @param fovAreaDir   directory path containing one or more RDA file of total area 
#'                     of each FOV
#' @param sampAnn      tibble of all sample annotation data to be joined to cell data
#' @param cellDiveID   cellDiveID for which data should be loaded; default: All 
#' 
#' @return single tibble of all cell and study annotations
loadAnnotatedCells <- function(sampAnn, annFile = NULL, annDir = NULL, fovAreaDir = NULL, 
                               cellDiveID = "all"){ 

    if(is.null(annFile) && is.null(annDir)){
        stop("Both params 'annFile' and 'annDir' are NULL. Must provide a value for one of these.")
    }

    faFiles <- NULL

    annFiles <- annFile
    if(is.null(annFiles)){
        annFiles <- file.path(annDir, dir(annDir))
    }
    #checkFilesFound(annFiles, annDir, "annotated cells", error = FALSE)
    #if(length(annFiles) == 0){
    #    return(NULL)
    #}

    if(!is.null(fovAreaDir)){
        faFiles  <- file.path(fovAreaDir, dir(fovAreaDir))
        tryCatch({
            checkFilesFound(faFiles, fovAreaDir, "FOV area")
         }, error = function(e){
            log_warn("No FOV areas added to annotated cells table.")
            faFiles <- NULL
        })
    }

    if(!is.null(cellDiveID) && tolower(cellDiveID[1]) != "all"){
        cdidPat  <- paste0("^", cellDiveID, "_")
        annFiles <- annFiles[grepl(cdidPat, basename(annFiles))]
        checkFilesFound(annFiles, annDir, "annotated cells", error = FALSE)
        if(is.null(annFiles) || length(annFiles) == 0){
            return(NULL)
        }
    }

    dat <- tibble()
    for(af in annFiles){
        log_debug(paste0("loading cells from file ",af))

        sampCells <- readRDS(af) %>%
                     left_join(sampAnn, by = intersect(names(.), names(sampAnn))) %>%
                     select(-dplyr::matches("MSK"), -dplyr::matches("detailed"),
                            -dplyr::matches("exclusion"))

        ## get band data for sample, if it exists
        cdids <- sampCells %>% pull(CellDive_ID) %>% unique()

        for(id in cdids){
            faFile <- NULL
            ptrn <- paste0("^", id, "_")
            if(!is.null(faFiles)){ faFile <- faFiles[grep(ptrn, basename(faFiles))] }
            if(!is.null(faFile) && length(faFile) == 1 && file.exists(faFile)){
                log_debug(paste0("joining FOV areas from file: ", faFile))
                fa <- readRDS(faFile)
                sampCells <- sampCells %>%
                             left_join(fa, by = intersect(names(sampCells), names(fa)))
            }
        }

        dat <- dat %>% bind_rows(sampCells)
    }

    dat
}

#' Load or generate index of all fraction and density conditions to be analyzed
#'
#' If an XLSX index file already exists, read and return the index; if not, generate
#' an index based on the Fraction and Density sheets in a manually-created conditions
#' XLSX file and condition configuration
#' 
#' @param indexFile         XLSX file containing indexed conditions or file to which newly
#'                          generated index will be writeen
#' @param conditionsFile    XLSX file containing a sheet for Fraction conditions and one for
#'                          Density conditions (see docs for details)
#' @param conditionConfig   config in nested list form (key 'arrange_annotation' in study 
#'                          overview config file); see docs for details
#'
#' @return tibble of indexed conditons including a Cell State ID
getConditionsIndex <- function(indexFile, conditionsFile, conditionConfig){

    if(!fileDone(indexFile)){
        condIdx <- indexConditions(conditionsFile, conditionConfig)
        write.xlsx(condIdx, indexFile, row.names=F, check.names = F)
    } else {
        log_debug(paste0("Loading previously indexed conditions from file ",indexFile))
        condIdx <- read.xlsx(indexFile, 1, check.names=F) %>% as_tibble()
    }
    condIdx
}

#' Load a single Halo megatable from RDA file
#' 
#' Load all halo data from a single RDA file. Only necessary columns are kept
#' in order to minimize memory. Only positive marker data is kept and any cells NOT
#' positive for DAPI are removed. Optionally, all data marked EXCLUDE are also
#' removed. Finally, column names are adjusted for consistency throughout pipeline, 
#' specifically 'Sample' is changed to 'CellDive_ID' and 'SPOT' is changed to 'FOV_number'.
#' 
#' @param file               full path to RDA file
#' @param sampAnn            flattened sample and FOV meta data
#' @param filterExclusions   logical; when set to TRUE, any row with text in column EXCLUDE
#'                           will be filtered out
#' @param controlMarker      name of marker used to identify individual, usable cells
#'
#' @return  filtered & formatted Halo data
loadHaloDataFile <- function(file, sampAnn, filterExclusions = FALSE, controlMarker = "DAPI"){

    dat <- readRDS(file)
    if(filterExclusions){
        pre <- nrow(dat$geom.data)
        dat$geom.data <- dat$geom.data %>% filter(Exclude == FALSE)
        log_debug(paste0("Excluded [", formatC(pre - nrow(dat$geom.data), big.mark = ","), "] cells."))
    }
 
    ## check for dapi neg cells
    dapiNeg <- dat$marker.data %>% filter(Marker == "DAPI", Positive_Classification != 1)
    if(nrow(dapiNeg) > 0){
        log_warn(paste0(nrow(dapiNeg), " DAPI negative cells being filtered out."))
        negCells <- unique(dapiNeg$UUID)
        dat$marker.data <- dat$marker.data %>% filter(!UUID %in% negCells)
        dat$geom.data <- dat$geom.data %>% filter(!UUID %in% negCells) 
    }
 
    dat$geom.data %>%
    left_join(dat$marker.data %>% select(UUID, Marker, Positive_Classification), ## leave out intensity for now
              by = "UUID") %>%
    select(CellDive_ID = Sample, 
           FOV_number = SPOT, 
           UUID, 
           dplyr::matches("X|Y"), 
           PixelScaling, 
           x0, y0, 
           Marker,
           Value = Positive_Classification)

}

#' Load all halo data files into a single table
#' 
#' Load all halo data files into a single table
#'
#' @param dataDir            directory containing ONLY halo RDA files, all of which will be loaded
#' @param dataFiles          vector of specific file paths, possibly a subset of files in a directory
#' @param nThreads           number of threads
#' @param filterExclusions   logical; when set to TRUE, any row with text in column EXCLUDE
#'                           will be filtered out
#' @param controlMarker      name of marker used to identify individual, usable cells
#'
#' @return a single table of filtered & formatted Halo data
loadAllHaloData <- function(dataDir = NULL, dataFiles = NULL, nThreads = 1, 
                            filterExclusions = TRUE, controlMarker = "DAPI"){
    allDat <- NULL
    log_info("Loading all pre-processed halo data...")

    if(is.null(dataFiles) && !is.null(dataDir)){
        dataFiles <- file.path(dataDir, dir(dataDir))
    }

    cl <- makeCluster(nThreads, type = "FORK", outfile = "")
    clusterExport(cl, c('filterExclusions', 'controlMarker'), envir = environment())

    allDat <- parLapply(cl, dataFiles, function(df){
                    log_info(paste0("Reading file ",df))
                    loadHaloDataFile(df, 
                                     filterExclusions = filterExclusions, 
                                     controlMarker = controlMarker)
              }) %>%
              bind_rows()

    stopCluster(cl)
    return(allDat)
}

#' Read, format and assemble all patient/sample/fov meta data
#' 
#' Read XLSX files for Patients, Samples and FOVs. Generate unique identifiers
#' for each sample in the form [Patient_ID]_[Sample_number] and for each FOV
#' in the form [Patient_ID]_[Sample_number]_[FOV_number]. Join and flatten all
#' data into a single complete table and create a separate table containing a 
#' map of all identifiers.
#'
#' @param metaFiles     vector of XLSX meta files, including at minimum *_Patients.xlsx,
#'                      *_Samples.xlsx and *_FOVs.xlsx; see docs for details on the formats
#'                      of these files
#' @param metaDataFile  RDA file that already contains or will contain the list of tables
#'                      returned from this function; if file already exists, data in the
#'                      file are what will be returned; otherwise, data will be compiled
#'                      and saved in this file
#' @param pathMap       list of file paths named by the corresponding paths existing in 
#'                      *_Paths.xlsx meta data file
#' @return list of five tables: Patients, Samples, FOVs, flat, IDs
loadStudyAnnotations <- function(metaFiles = NULL, metaDataFile = NULL, pathMap = NULL){

    if(fileDone(metaDataFile)){
        log_info(paste0("Reading pre-processed meta data from file ",metaDataFile))
        return(readRDS(metaDataFile))
    }

    sampleAnnFile <- getMetaFile("Samples", metaFiles)
    fovAnnFile    <- getMetaFile("FOVs", metaFiles) 
    pathFile      <- getMetaFile("Paths", metaFiles) 

    log_debug("Reading Sample annotations...")
    if(is.null(sampleAnnFile) || length(sampleAnnFile) != 1){
        log_error("No *_Samples.xlsx file found.")
    }
    samples <- as_tibble(read.xlsx(sampleAnnFile, 1, check.names = FALSE)) %>% 
               filter(!grepl("_WS$", CellDive_ID), !is.na(CellDive_ID))

    log_debug("Reading FOV annotations...")
    if(is.null(fovAnnFile) || length(fovAnnFile) != 1){
        log_error("No *_FOVs.xlsx file found.")
    }
    fovs <- as_tibble(read.xlsx(fovAnnFile, 1, check.names = FALSE)) %>%
            filter(!grepl("_WS$", CellDive_ID)) %>%
            mutate(FOV_number = ifelse(FOV_number == "Whole_slide", "1", FOV_number)) %>%
            mutate(FOV_number = as.numeric(FOV_number))

    ##
    ### FORMAT IDENTIFIERS
    ##
    ptd <- paste0("%0",max(nchar(as.character(samples$Patient_ID))),"d")
    fnd <- paste0("%0",max(nchar(as.character(fovs$FOV_number))),"d")

    ##
    ### JOIN ALL DATA INTO A SINGLE FLATTENED TABLE
    ##
    log_debug("Flattening meta data...")
    flat <- samples %>%
            left_join(fovs, by = "CellDive_ID") %>%
            mutate(Sample_ID = Patient_ID,
                   FOV_ID = paste(Sample_ID,
                                  sprintf(fnd, FOV_number),
                                  sep = "_")) %>%
            filter(is.na(FOV_exclusion))

    ##
    ## CREATE A TABLE OF ONLY IDENTIFIERS
    ##
    log_debug("Creating map of all IDs...")
    idMap <- flat %>%
             select(Patient_ID, CellDive_ID, Sample_ID, FOV_ID, FOV_number) %>%
             unique()

    ##
    ### RETURN A LIST OF ALL TABLES
    ##
    annot <- list(Samples = samples, 
                  FOVs = fovs, 
                #  Paths = paths,
                  flat = flat, 
                  IDs = idMap)

    return(annot)
}



#' Compile all analyses to be run
#' 
#' Read 'conditions' XLSX file, in which each tab contains a list of conditions to 
#' be analyzed using a specific calculation type. These types include general fractions
#' and densities in addition to spatial analyses 'neighborhood' fractions and 'neightobhood'
#' average counts
#' 
#' @param condFile    XLSX 'conditions' file described in the docs
#' @param condIdx     table of indexed conditions including Cell State IDs
#'
#' @return list of indexed and formatted conditions where each element contains
#'         conditions pertaining to one calculation type
getStatsAnalyses <- function(condFile, condIdx){ #, cellTypes, funcMarkers, funcCombos){

    fracConds  <- read.xlsx(condFile, "fractions", check.names=F) %>%
                  as_tibble() %>%
                  left_join(condIdx %>%
                            filterConditionsByAnalysisType("general"),
                            by=intersect(names(.), names(condIdx)))

    densConds <- getDensityConditions(condIdx) %>%
                 left_join(condIdx %>% select(-Population, -Subpopulation) %>% unique()) %>%
                 left_join(read.xlsx(condFile, "densities", check.names=F) %>%
                           as_tibble())

    list(fractions  = fracConds,
         densities  = densConds)
}


#' Convert XLSX table of study 'question' info to a list
#' 
#' Read XLSX file where each row contains info describing a question
#' and how to subset data in order to answer the question. Convert all
#' rows into a nested list
#'
#' @param  questionXLSXfile  XLSX file
#'
#' @return  question data in list form
parseStatsQuestions <- function(questionsXLSXfile){

    qList <- list()
    questions <- read.xlsx(questionsXLSXfile, 1, check.names=FALSE) %>% as_tibble %>%
                 filter(!is.na(QuestionNumber))

    for(x in 1:nrow(questions)){
        grp  <- as.list(questions[x,])
        if(!grp$QuestionNumber %in% names(qList)){
            q <- list(question  = grp$Question,
                      groupVar  = grp$ComparisonVariable,
                      qNum      = grp$QuestionNumber,
                      calcUnit  = grp$CalculationUnit,
                      `Group 1` = list(),
                      `Group 2` = list())
            qList[[grp$QuestionNumber]] <- q
        }
        filt <- grp[5:length(grp)]
#        filt <- lapply(filt[!is.na(filt)], function(x) trimws(unlist(strsplit(as.character(x), ";;"))) )
filt <- lapply(filt[!is.na(filt)], function(x) trimws(unlist(strsplit(as.character(x), ", "))) )
        qList[[grp$QuestionNumber]][[paste0("Group ",grp$Group)]] <- filt
    }
    qList
}



loadNeighborData <- function(dataFile, cellAnn = NULL){

    cdid <- gsub("___.*", "", basename(dataFile))
    nbhrs <- readRDS(dataFile)$neighbors %>%
             mutate(CellDive_ID = cdid, FOV_Number = SPOT)

    if(!is.null(cellAnn)){
        cellAnn <- cellAnn %>% select(UUID, Category, Cell_type, Subtype, PositiveMarkers)
        cAnn <- cellAnn
        names(cAnn) <- paste0("C.", names(cAnn))
        nAnn <- cellAnn
        names(nAnn) <- paste0("N.", names(nAnn))
        nbhrs <- nbhrs %>%
                 left_join(cAnn, by = c("C.UUID")) %>%
                 left_join(nAnn, by = c("N.UUID"))
    }

}



