## required columns (constants)
req_cols <- list(patients = c("Patient_ID"), # "Patient_response"),
                 samples  = c("CellDive_ID", "Patient_ID", "Microscope", "DAPI_first", "DAPI_last"),
                 fovs     = c("CellDive_ID", "FOV_number", "FOV_exclusion", #"FOV_exclusion_stage", 
                              "Marker_exclusion", "Num_manual_exclusions"), #"Num_epidermis_exclusions", 
                              #"Num_interface_areas"),
                 cell_types = c("Category", "Cell_type", "Subtype"), 
                                #"Pos_markers", "Pos_required", "Neg_markers"),
                 markers = c("Marker_name", "Identity", "Functional", "Threshold_compartment"), 
                 paths = c("CellDive_ID", "path_to_halo_csv_nonAI", "path_to_images"),
                 annotated_cells = c("Patient_ID", "Sample_ID", "FOV_ID"),
                 sample_annotation = c("Sample", "FOV"),
                 questions = c("Question", "QuestionNumber", "ComparisonVariable", "Group")) 


#' Check that a table includes all columns required for that dataset
#' 
#' @param cols   vector of cols in the table to be checked
#' @param dataset  name of dataset being checked; values allowed:
#'                      ['patients'|'samples'|'fovs'|'cell_types'|'markers'|
#'                        'annotated_cells'|'sample_annotation']
#' @return vector of error messages (NULL if none found)
checkRequiredCols <- function(cols, dataset){
    errs <- c()
    req <- req_cols[[dataset]]
    missing <- req[!req %in% cols]
    if(length(missing) > 0){
        for(m in missing){ errs <- c(errs, paste0("Column missing from ",dataset," table: ",m)) }
    }
    return(as.list(errs))
}

#' Check for existence of all required meta files
#' 
#' Given entire study configuration, ensure that either {meta_files} or {meta_dir} is defined.
#' If {meta_dir} given and not {meta_files}, populate {meta_files} with all xlsx files in 
#' directory defined by {meta_dir}. Minimum required files are: 
#'    *_CellTypes*.xlsx
#'    *_Markers*.xlsx
#'    *_Patients*.xlsx
#'    *_Samples*.xlsx
#'    *_FOVs*.xlsx
#'
#' @param sCfg    study config in list format
#' @return nested list of errors, warnings and messages
validateMetaFileList <- function(sCfg){

    log_debug("  validating meta file list..\n")
    errors = list()
    warnings = list()
    messages = list()

    ## meta_dir/files
    metaFiles <- c()
    if(is.null(sCfg$meta_files)){
        if(is.null(sCfg$meta_dir)){
            errors[[length(errors) + 1]] <- "Must specify meta_dir OR/AND meta_files"
        } else if(!dir.exists(sCfg$meta_dir) || length(dir(sCfg$meta_dir)) == 0) {
            errors[[length(errors) + 1]] <- "meta_dir must exists and must not be empty"
        } else {
            metaFiles <- getCurrentMetaFiles(config = sCfg) 
        }
    } else {
        metaFiles <- sCfg$meta_files
    }

    req <- c("_CellTypes.*\\.xlsx",
             "_CellStates.*\\.xlsx", 
             "_Markers.*\\.xlsx",
             "_Samples.*\\.xlsx",
             "_FOVs.*\\.xlsx",
             "_Paths.*\\.xlsx")

    for(fs in req){
        if(!any(grepl(fs, metaFiles))){
            errors[[length(errors) + 1]] <- paste0("Can not find required file *",fs," in list of meta files.")
        }
    }

    return(list(errs=errors, wrns=warnings, msgs=messages))
}


#' Check for access to halo boundary information
#' 
#' Given study config if list format, check for either a pre-processed
#' *.rda file or a directory of XML files directly from Halo; if a directory is specified,
#' it must contain one or more subdirectories, each containing XML files for one type of boundary.
#' 
#' @param sCfg    study config in list format
#' @return nested list of errors, warnings and messages
validateHaloBoundaryFiles <- function(sCfg){
    log_debug("  validating halo boundary files..\n")
    errors = list()
    warnings = list()
    messages = list()

    hbf <- sCfg$halo_boundaries_file
    hbd <- sCfg$halo_boundaries_xml_dir

    ## must have either halo boundary annot file OR annotation dirs
    if(is.null(hbf) | !file.exists(hbf)){
        ## validate annotations_dirs
        if(is.null(hbd) || !dir.exists(hbd)){
            errors[[length(errors) + 1]] <- "Must provide either halo_boundaries_file or halo_boundaries_xml_dir"
        } else {
            if(file.size(hbd) == 0 || length(dir(hbd)) == 0){
                errors[[length(errors) + 1]] <- "Param [halo_boundaries_xml_dir] provided but is empty"
            }
            if(!any(sapply(file.path(hbd, dir(hbd)), function(x){ dir.exists(x) }))){
                errors[[length(errors) + 1]] <- "Param [halo_boundaries_xml_dir] provided but does not contain any subdirectories. This directory should contain a subdirectory for each type of Halo boundary (Interface, Exclusion, Epidermis, Glass)"
            }
        }
    } else if(!is.null(hbf)) {
        if(!file.exists(hbf) || file.size(hbf) == 0){
            errors[[length(errors) + 1]] <- paste0("Param [halo_boundaries_file] is given but file is empty")    
        }
    }

    return(list(errs=errors, wrns=warnings, msgs=messages))
}


#' List files that already exist and will be used as input and files that are to be 
#' written by pipeline
#'
#' Given entire study config, for each *_file paramater print info on whether the file
#' exists and will be used as input or does not yet exist (or is empty) and will be
#' generated during analysis.
#' 
#' @param   sCfg    study config in list format
#' @return nested list of errors, warnings and messages
validatePipelineFiles <- function(sCfg){
    log_debug("  validating provided files...\n")
    warnings <- list() 
    existing <- list()
    non      <- list()
    fileFields <- names(sCfg)[grepl("_file$",names(sCfg))]
    for(ff in fileFields){
        if(fileDone(sCfg[[ff]])){ 
            existing[[ff]] <- sCfg[[ff]] 
        } else {
            non[[ff]] <- sCfg[[ff]]
        }
    }
    if(length(existing) > 0){
        warnings <- as.list("The following files already exist and will be used (NOT OVERWRITTEN) in associated analyses: ")
        warnings <- c(warnings, paste0("\t", names(existing), ": ", existing)) 
        warnings[[length(warnings) + 1]] <- "\n"
    }
    if(length(non) > 0){
        warnings[[length(warnings) + 1]] <- "The following files do NOT exist or are empty and contents will be generated from scratch:"
        warnings <- c(warnings, paste0("\t", names(non), ": ",non)) 
        warnings[[length(warnings) + 1]] <- "\n"
    }

    return(list(errs=list(), wrns=warnings, msgs=list()))  
}


#' Indicate duplicate values within a vector, usually a column from a dataframe
#' 
#' Given a vector of values, print info on which values are found multiple times.
#' 
#' @param vals  vector of values
#' @param label  label (usually column name) to include in message
#' @param rtrn   ['message'|'values'] what should be returned
#' @return  character string with message of duplicates or vector of
#'          duplicated values (NULL if none found)
duplicatesWithinCol <- function(vals, label, rtrn='message'){
    e <- NULL
    if(length(vals[duplicated(vals)]) > 0){
        e <- switch(rtrn,
                    'message' = paste0("Duplicated ",label,"(s): ",
                        paste(vals[duplicated(vals)], collapse=", ")),
                    'values' = vals[duplicated(vals)],
                    NULL)
    }
    e
}


#' Indicate duplicates across two columns in a data frame
#' 
#' Given a data frame of two columns, print info on values that occur in both
#' columns.
#' 
#' @param df  data frame of two columns
#' @return  character string with message of values in both columns (NULL if none found)
duplicatesAcrossCols <- function(df){
    e <- NULL
    if(any(df[,1] %in% df[,2])){
        e <- paste0(paste(df[,1][df[,1] %in% df[,2]], collapse=", "),
                    " in both ",names(df)[1]," and ",names(df)[2]," columns.")
    }
    e
}

#' Check for conflicts between marker combinations of duplicate classes
#' 
#' @param df  data frame of cell types with a single identifier and at minimum, 
#'            columns Pos_markers and Neg_markers
#' @return message indicating that there is a marker combo conflict between 
#'         cell types with the same identifier
duplicateClassConflicts <- function(df){
    pos <- unique(unlist(strsplit(df %>% pull(Pos_markers), ",")))
    neg <- unique(unlist(strsplit(df %>% pull(Neg_markers), ",")))
    if(length(intersect(pos,neg)) > 0){
        return(paste0("Duplicated Subtypes with conflicting marker combinations: ",
                      paste(stDups, collapse=", ")))
    }
    return(NULL)
}

#' Perform various checks to ensure cell types spreadsheet has been filled out correctly.
#' 
#' Criteria to be checked include:
#'    * Values in Subtype and Tag columns must be unique 
#'    * No value may occur in more than one column
#'    * All markers listed in columns Pos_markers and Neg_markers must be in markers XLSX file
#'    * All Pos_markers+Neg_markers combinations must be unique to avoid ambiguous cell type definitions
validateCellTypes <- function(cellTypesXLSX, markerDescFile){

    errs <- wrns <- list()
    
    cts   <- getCellTypes(cellTypesXLSX)
    mDesc <- read.xlsx(markerDescFile, 1, check.names=F) %>% as_tibble()
    errs <- checkRequiredCols(names(cts), 'cell_types')
    errs <- c(errs, checkRequiredCols(names(mDesc), 'markers'))

    ## 
    ## unique subtypes and tags
    ##
    for(col in c("Subtype")){ #,"Tag")){
        #colDups <- duplicatesWithinCol(cts %>% filter(!is.na(!!as.name(col)), Pos_required == 'all') %>% pull(col), col, rtrn='values')
        colDups <- cts %>% 
                   filter(!is.na(!!as.name(col)), Pos_required == 'all') %>% 
                   pull(col) %>%
                   duplicatesWithinCol(col, rtrn = 'values')
        ## if there are subtype duplicates, check that their marker combinations do not conflict
        if(!is.null(colDups)){
            colErrs <- c()
            for(dup in colDups){
                cnfl <- duplicateClassConflicts(cts %>% 
                                                filter_at(all_of(col), any_vars(. == dup)) %>% 
                                                filter(Pos_required == 'all'))
                colErrs <- c(colErrs, cnfl) 
            } 
            if(length(colErrs) > 0){ 
                errs <- c(errs, colErrs) 
            } else {
                wrns <- c(wrns, paste0("Duplicate but apparently valid ",col,"(s) found: ",
                                       paste(colDups, collapse=", ")))
            }
        }
    }

    ##
    ## no duplicates between Category, Cell_type, Subtype
    ##
    errs[[length(errs) + 1]] <- duplicatesAcrossCols(cts %>% select(Category, Cell_type))
    errs[[length(errs) + 1]] <- duplicatesAcrossCols(cts %>% select(Category, Subtype))
    errs[[length(errs) + 1]] <- duplicatesAcrossCols(cts %>% select(Cell_type, Subtype))

    ##
    ## all markers are in MarkerDesc
    ##
    allM <- unique(unlist(strsplit(c(cts$Pos_markers, cts$Neg_markers),",")))
    allM <- allM[!is.na(allM)]
    if(!all(allM %in% mDesc$Marker_name)){
        e <- paste0("Unrecognized markers (markers not in Markers file): ",
                    paste(allM[!allM %in% mDesc$Marker_name], collapse=", "))
        errs[[length(errs) + 1]] <- e
    }

    ##
    ## unique marker combinations (one set of classifications per combination)
    ##
    fullCombos <- cts %>%
                  mutate(Negs = sapply(Neg_markers, function(x){ 
                                  paste(sort(paste0(unlist(strsplit(x,",")),"-")), 
                                         collapse=",") 
                                }), 
                         Pos = sapply(Pos_markers, function(x){ 
                                  paste(sort(unlist(strsplit(x,","))),
                                         collapse=",") 
                                })
                        )

    if(any(duplicated(fullCombos))){
        e <- paste0("Duplicate marker combinations: \n",
                    paste0(fullCombos %>% 
                           unite("Marker combo", Pos, Negs, sep=",") %>% 
                           pull(`Marker combo`) %>% 
                           duplicated, collapse="\n"))
        errs[[length(errs) + 1]] <- e
    }
 
    if(length(wrns) > 0){
        for(w in unlist(wrns)){ log_warn(w) }
    }
    if(length(errs) > 0){
        for(e in unlist(errs)){ log_error(e) }
        return(FALSE)
    }
    return(TRUE)    
}

#' Check for invalid and duplicate cell types and states
#'
#'
validateCellStates <- function(csXLSX, ctXLSX, markerXLSX){

    errs <- list()

    states  <- read.xlsx(csXLSX, 1, check.names = F) %>% as_tibble
    types   <- read.xlsx(ctXLSX, 1, check.names = F) %>% as_tibble
    markers <- read.xlsx(markerXLSX, 1, check.names = F) %>% as_tibble

    validCols <- types %>%
                 select(Category, Cell_type, Subtype) %>%
                 unlist %>%
                 unique
    if(!all(names(states)[-1] %in% validCols)){
        e <- paste0("Invalid cell type(s): ", 
                    paste(setdiff(names(states)[-1], validCols), 
                          collapse = ", "))
        errs[[length(errs) + 1]] <- e
    }    
    
    sts <- unlist(strsplit(states$state, ","))
    if(!all(sts %in% markers$Marker_name)){
        e <- paste0("Invalid state(s): ", 
                    paste(setdiff(sts, markers$Marker_name),
                          collapse = ", "))
        errs[[length(errs) + 1]] <- e                  
    }    

    if(any(duplicated(names(states)))){  ## error should be thrown by as_tibble above
        dups <- names(states)[duplicated(names(states))]
        e <- paste0("Cell type(s) duplicated: ", paste(dups, collapse = ", "))
        errs[[length(errs) + 1]] <- e    
    }

    if(any(duplicated(states$state))){
        dups <- states$state[duplicated(states$state)]
        e <- paste0("Cell state(s) duplicated: ", paste(dups, collapse = ", "))
        errs[[length(errs) + 1]] <- e
    }

    if(length(errs) > 0){
        for(e in unlist(errs)){ log_error(e) }
        return(FALSE)
    }
    return(TRUE)
}

#' Check for duplicate Patient_ID
#' 
#' @param patientAnnotationXLSX    XLSX file containing columns Patient_ID and Patient_response
#'                                 at least
#' @return list of errors, if any
validatePatientAnnotation <- function(patientAnnotationXLSX){
    errors <- list()
    pa <- read.xlsx(patientAnnotationXLSX, 1, check.names=F) %>% as_tibble() %>% filter(!is.na(Patient_ID))
    errors <- checkRequiredCols(names(pa), 'patients')
    errors[[length(errors) + 1]] <- duplicatesWithinCol(pa %>% pull(Patient_ID), "Patient_ID")
    return(list(errs=errors))
}


#' Check for errors in sample annotation
#' 
#' @param sampleAnnotationXLSX    XLSX file containing all sample annotation
#' @param patientAnnotationXLSX   XLSX file containing all patient annotation
#' @return  list of errors and warnings 
validateSampleAnnotation <- function(sampleAnnotationXLSX){ #, patientAnnotationXLSX){
    errors <- warnings <- list() 

    sa <- read.xlsx(sampleAnnotationXLSX, 1, check.names=F) %>% as_tibble() %>%
          filter(!is.na(CellDive_ID))
    errors <- checkRequiredCols(names(sa), 'samples') 

    ##
    ## unique CellDive_ID and Patient_ID + Sample_number (Sample_ID) combinations
    ##
    errors[[length(errors) + 1]] <- duplicatesWithinCol(sa %>% pull(CellDive_ID), "CellDive_ID")
    errors[[length(errors) + 1]] <- duplicatesWithinCol(sa %>% pull(Patient_ID), "Patient_ID")

    return(list(errs=errors, wrns=warnings))
}

#' Check for errors in FOV annotation
#' 
#' @param fovAnnotationXLSX    XLSX file containing all FOV-level annotation
#' @param sampleAnnotationXLSX XLSX file containing all Sample-level annotation
#' @return  list of errors and warnings
validateFOVannotation <- function(fovAnnotationXLSX, sampleAnnotationXLSX){

    errors <- warnings <- messages <- list()

    sa <- as_tibble(read.xlsx(sampleAnnotationXLSX, 1, check.names=F)) %>%
          filter(!is.na(CellDive_ID))
    fa <- as_tibble(read.xlsx(fovAnnotationXLSX, 1, check.names=F)) %>%
          filter(!is.na(FOV_number), !is.na(CellDive_ID)) %>%
          mutate(CDID_FOV_combo = paste0(CellDive_ID,"_", FOV_number)) %>%
          filter(!grepl("_WS$", CellDive_ID))
    errors <- checkRequiredCols(names(fa), 'fovs')

    ##
    ## 1:1 CellDiveID to FOV_number (unique fovs within each sample)
    ##
    errors[[length(errors) + 1]] <- duplicatesWithinCol(fa %>% pull(CDID_FOV_combo), "CellDive_ID + FOV combo")

    ##
    ## warn of any differences between sample annotation and fov annotation sample lists
    ##
    if(!identical(sort(unique(fa$CellDive_ID)), sort(unique(sa$CellDive_ID)))){
        inSAnotFA <- setdiff(sa$CellDive_ID, fa$CellDive_ID) %>% unique
        inFAnotSA <- setdiff(fa$CellDive_ID, sa$CellDive_ID) %>% unique
        if(length(inSAnotFA) > 0 & !all(grepl("_WS", inSAnotFA))){
            ## this is a warning, not an error because some samples are whole slides (no FOVs)
            warnings[[length(warnings) + 1]] <- 
                     paste0("CellDive_ID(s) found in *SampleAnnotations.xlsx but not *FOVannotations.xlsx: ",
                            paste(inSAnotFA, collapse=", "))
        }
        if(length(inFAnotSA) > 0){
            errors[[length(errors) + 1]] <- 
                     paste0("CellDive_ID(s) found in *FOVannotations.xlsx but not *SampleAnnotations.xlsx: ",
                            paste(inFAnotSA, collapse=", "))
        }
    }
 
    return(list(errs=errors, wrns=warnings))
}

#' Check for errors in all study annotation files
#' 
#' @param fovAnnotationXLSX 
#' @param sampleAnnotationXLSX
#' @param FALSE if any errors found, TRUE if not
validateStudyAnnotation <- function(fovAnnotationXLSX, sampleAnnotationXLSX){
    log_debug("Validating sample annotation...\n")
    resS <- validateSampleAnnotation(sampleAnnotationXLSX)
    log_debug("Validating FOV annotation...")
    resF <- validateFOVannotation(fovAnnotationXLSX, sampleAnnotationXLSX)
   
    errors   <- c(unlist(resS$errs), unlist(resF$errs))
    warnings <- c(unlist(resS$wrns), unlist(resF$wrns))
    messages <- c(unlist(resS$msgs), unlist(resF$msgs))

    cat("\n\n")
    if(length(warnings) > 0){
        for(w in unlist(warnings)){ log_warn(w) }
    }
    if(length(messages) > 0){
        for(m in unlist(messages)){ log_info(m) }
    }
    if(length(errors) > 0){
        for(e in unlist(errors)){ log_error(e) }
        return(FALSE)
    }
    return(TRUE)
}

validateStudyConfig <- function(sCfg){

    errors <- warnings <- messages <- c()

    allVals <- c(validateMetaFileList, validatePipelineFiles)

    for(x in seq(allVals)){
        validation <- allVals[x][[1]]
        res <- validation(sCfg)
        errors <- c(errors, unlist(res$errs))
        warnings <- c(warnings, unlist(res$wrns))
        messages <- c(messages, unlist(res$msgs))
    }

    cat("\n")
    if(length(warnings) > 0){
        for(w in warnings){ log_warn(w) }
    }
    if(length(messages) > 0){
        for(m in messages){ log_info(m) }
    }
    if(length(errors) > 0){
        for(e in errors){ log_error(e) }
        return(FALSE)
    }
    return(TRUE)
}

validateExclusions <- function(dat, studyAnn){

    errors <- warnings <- messages <- list()

    ## full FOV exclusions (FOVs should NOT appear in data AT ALL)
    exFOVs <- studyAnn %>% 
              select(FOV, dplyr::matches("FOV_exclusion")) %>%
              filter_all(any_vars(str_detect(., pattern = "X"))) %>%
              pull(FOV)
    errors <- warnings <- messages <- list()
    if(any(exFOVs %in% dat$FOV)){
        errors[[length(errors) + 1]] <- paste0("These FOVs are still in the data: ", paste(exFOVs[exFOVs %in% dat$FOV], collapse=", "))
    }

    ## marker exclusions
    mEx <- studyAnn %>% 
           select(FOV, Marker_exclusion) %>%
           filter(!is.na(Marker_exclusion))
    for(r in 1:nrow(mEx)){
        ex <- mEx[r,]
        for(m in unlist(strsplit(ex$Marker_exclusion))){
            mDat <- dat %>% 
                    filter(FOV == ex$FOV, grepl(getClassifierPattern(m, delim=","), FullPosMarkerStr))
            if(!is.null(mDat) && nrow(mDat) > 0){
                errors[[length(errors) + 1]] <- paste0("Data still contains marker ",m," in FOV ",ex$FOV)
            }
        }
    }

    return(list(errs=errors, wrns=warnings))
             
}       

samplesInDir <- function(datDir, ext){
    gsub(ext,"",dir(datDir))
}

getStudySampleList <- function(samplesXLSX, imageDir, csvDir){
    meta <- read.xlsx(samplesXLSX, 1, check.names=F) %>% as_tibble() %>% pull(CellDive_ID)
    imgs <- samplesInDir(imageDir, "_scanplan_AllFOVs.tif")
    csvs <- samplesInDir(csvDir, "_v1.csv.gz")

    errs <- wrns <- c()

    if(!all(meta %in% imgs)){
        for(m in setdiff(meta,imgs)){ errs <- c(errs, paste0("Sample in meta data but NOT in image dir: ",m)) }    
    }
    if(!all(meta %in% csvs)){
        for(m in setdiff(meta,csvs)){ errs <- c(errs, paste0("Sample in meta data but NOT in CSV dir: ",m)) }   
    }
    if(!all(imgs %in% meta)){
        for(i in setdiff(imgs,meta)){ wrns <- c(wrns, paste0("Sample in image dir but NOT in meta data: ",i)) }
    }
    if(!all(csvs %in% meta)){
        for(c in setdiff(csvs,meta)){ wrns <- c(wrns, paste0("Sample in CSV dir but NOT in meta data: ",c)) }
    }

    if(length(wrns) > 0){
        for(w in wrns){ log_warn(w) }
    }
    if(length(errs) > 0){
        for(e in errs){ log_error(e) }
    }
    return(meta)
}

validateQuestionDesign <- function(question, validMeta, markers = NULL){

    errs    <- c()

    nonAnnoCompVars <- list("Cell Region" = c(NA, "All", "HRS_nbhd")) 

    if(!is.null(markers)){
        names(markers) <- paste0("HRS,", markers, "_env")
        nonAnnoCompVars <- c(nonAnnoCompVars,
                             lapply(markers, function(x){ 
                                 c("pos.env","neg.env","mixed.env") 
                             }))
    }

    ## what is to be compared?
    compVar  <- unique(question$ComparisonVariable)
    if(is.na(compVar)){
        errs <- c(errs, paste0("Question has no ComparisonVariable."))
        return(errs)
    }
    if(!compVar %in% names(question)){
        errs <- c(errs, paste0("Question has invalid ComparisonVariable. '",compVar ,"' is not a column in questions file."))
        return(errs)
    }
    grpNames <- question[[compVar]]
    if(any(is.na(grpNames))){
        errs <- c(errs, paste0("Question is missing one or more comparison group"))
        return(errs)
    }
    grpVals  <- grpNames %>% strsplit(",") %>% unlist %>% trimws

    ## check for invalid cell regions
    if(compVar %in% names(nonAnnoCompVars)){
        if(!all(tolower(grpNames) %in% nonAnnoCompVars[[compVar]])){
            err  <- paste0("The following group names are not valid for comparison variable ",
                           compVar, ": ",
                           paste(grpVals[ ! grpVals %in% nonAnnoCompVars], collapse=", "))
            #errs <- c(errs, err)
        }
    } else {
        ## check for invalid annotation values
        invalidGrps <- grpVals[ !grpVals %in% c(NA, validMeta[[compVar]]) ]

        if(length(invalidGrps) > 0){
            err  <- paste0("The following group names are not valid for comparison variable ",
                           compVar, ": ",
                           paste(invalidGrps, collapse=", "))
            #errs <- c(errs, err)
        }
    }
 
    ## check that groups are actually different
    if(grpNames[1] == grpNames[2]){
        grpStr <- paste(paste0(c("Group 1: [", "Group 2: ["), 
                        paste0(grpNames, "]")), collapse = ", ")
        errs <- paste0("ComparisonVariable seems to be the same in both groups: ", 
                        grpStr, ". Is ComparisonVariable wrong?")
    }
    return(errs)

}

validateQuestions <- function(questionsXLSX, validMeta, markers = NULL){

    errs <- list()

    qs <- read.xlsx(questionsXLSX, 1, check.names = F) 
    qs <- qs %>% as_tibble() %>% filter(!is.na(QuestionNumber))
    errs <- checkRequiredCols(names(qs), 'questions')
    
    ## make sure all column labels match those in meta data
    filtCols  <- setdiff(names(qs), req_cols$questions)
    validCols <- c(names(validMeta), "Cell Region")
    if(!is.null(markers)){ validCols <- c(validCols, paste0("HRS,", markers, "_env")) }

    if(!all(filtCols %in% validCols)){
        invalid <- setdiff(filtCols, validCols)
        msg <- paste0("Filter column(s) in questions file not found in meta data: ", 
                      paste(invalid, collapse = ", "))
        #errs[[length(errs) + 1]] <- msg 
    }

    ## make sure all values in filter columns match those in meta data
    cols2check <- intersect(filtCols, names(validMeta))
    for(fc in cols2check){
        obsVals <- trimws(unique(unlist(strsplit(as.character(qs[[fc]]), ","))))
        expVals <- c(unique(validMeta[[fc]]), NA)
        if(!all(obsVals %in% expVals)){
            errs[[length(errs) + 1]] <- paste0("Invalid value(s) in column ",fc," of questions file: ",
                                               paste(obsVals[!obsVals %in% expVals], collapse=", "))
        }
    }

    ## check each question
    for(qu in unique(qs$QuestionNumber)){
        log_debug("Validating question ", qu)
        quest <- qs %>% filter(QuestionNumber == qu)
        qErrors <- validateQuestionDesign(quest, validMeta, markers = markers)
        if(!is.null(qErrors)){
            errs[[length(errs) + 1]] <- paste0(qu, ": ",qErrors)
        }
    }
 
    if(length(errs) > 0){
        for(e in unlist(errs)){ log_error(e) }
        return(FALSE)
    }
    return(TRUE)

}

#' Check for errors in all input data including questions, study annotation and
#' cell type definitions
#' 
#' @param study_config    study configuration in list format, with at minimum 
#'                        'meta_dir' defined
#' @return nothing; quit with exit status 1 if any errors found
validateInput <- function(study_config){

    # TODO: refactor this to be less redundant

    log_info("Validating meta files...\n")
    mFiles <- getCurrentMetaFiles(metaDir = study_config$meta_dir)
    saValid <- validateStudyAnnotation(getMetaFile("FOVs", mFiles),
                                       getMetaFile("Samples", mFiles))
    qValid <- FALSE
    if(saValid){
        log_info("Validating questions file...\n")
        meta <- loadStudyAnnotations(metaFiles = mFiles)$flat
        markers <- read.xlsx(getMetaFile("Markers", mFiles), 1, check.names = F) %>%
                   as_tibble() %>%
                   getAllMarkers()
        qValid <- validateQuestions(getMetaFile("Questions", mFiles), meta, markers = markers)
    }

    log_info("Validating markers and cell types...\n")
    #ctValid <- validateCellTypes(getMetaFile("CellTypes", mFiles), getMetaFile("Markers", mFiles))

    log_info("Validating cell states...\n")
    csValid <- validateCellStates(getMetaFile("CellStates", mFiles),
                                  getMetaFile("CellTypes", mFiles),
                                  getMetaFile("Markers", mFiles))

    if(!all(saValid, qValid, csValid)){ #ctValid, csValid)){
        log_fatal("Validation FAILED.")
        cat("\n")
        q(save="no", status=1)
    }

    log_info("All input is VALID!")
}

