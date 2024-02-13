#' Get numbers of cells belonging to a list of populations of interest
#' 
#' For each population of interest, count the number of cells determined to be in that population
#' 
#' @param dat           tibble of annotated cell data
#' @param populations   vector of comma-delimited character strings describing the cell populations to be 
#'                      counted;
#'                      each element of a single delimited string may be a cell Classifier or a positive or
#'                      negative marker; a positive marker is indicated simply by the marker name and a negative; IMPORTANT: this function assumes all populations can, in theory exist in every FOV (i.e., populations do not include functional markers that are missing in data)
#'                      marker is indicated by the marker name with '-' appended (e.g., "PD1-")
#' @param markers       vector of all marker names, used to distinguish a marker from a classifier in a
#'                      population string
#' @param calcUnit      data column name containing values for which a single count should be made
#'                      (default: FOV_ID)
#' @param numThreads    max number of threads to use when available
#' @param outFile       path to file that contains pre-counted population data OR the file to which data
#'                      should be saved
#' @param forceRecalc   logical, default = FALSE; count populations from data even if outFile
#'                      exists from a previous run
#'
#' @return table with columns 'Population', [calcUnit], and Count
getPopulationCounts <- function(dat, populations, markers, calcUnit = "FOV_ID", 
                                 numThreads = 1, outFile = NULL, forceRecalc = FALSE){

    if(fileDone(outFile) && !forceRecalc){
        log_info(paste0("Loading pre-computed counts from file ",outFile))
        return(readRDS(outFile))
    }
    log_info("No pre-computed counts file found and/or forced rerun turned ON. Counting now...")

    populations <- unique(populations[!is.na(populations)])
    popCounts <- lapply(1:length(populations), function(x){
                     popdat <- filterDataForCondition(dat, populations[x], markers)
                     if(is.null(popdat) || nrow(popdat) == 0){ 
                         dat %>% 
                         select_at(calcUnit) %>% 
                         unique %>%
                         mutate(Population = populations[x], Count = ifelse(is.null(popdat), NA_real_, 0)) 
                     } else {
                         popdat %>%
                         group_by_at(calcUnit) %>%
                         dplyr::summarize(Count = n()) %>%
                         mutate(Population = populations[x])
                     }
                 }) %>%
                 bind_rows() #%>%
                 #spread(calcUnit, Count, fill = 0)

    #allCounts <- tibble(Population = unique(populations)) %>% 
    #             left_join(popCounts, by = "Population") %>%
    allCounts <- popCounts %>%
                 spread(calcUnit, Count, fill = 0) %>%
                 gather(2:ncol(.), key = !!as.name(calcUnit), value = "Count")

    if(!is.null(outFile)){
        saveRDS(allCounts, file = outFile)
    }

    allCounts
}

#' Get fractions of cell populations
#'
#' Given a tibble of cell population counts and a set of population fractions to calculate,
#' generate a table of said fractions 
#'
#' @param counts      tibble of counts for all cell populations included in list of fractions
#' @param conditions  tibble describing population fractions to calculate, including at minimum
#'                    columns Subpopulation and Population, where Population is a character string
#'                    describing the overall cell population (denominator) and Subpopulation is a
#'                    character string describing the cell state (numerator)
#' @param calcUnit    data column name containing values for which a single fraction
#'                    should be calculated (default: FOV_ID)
#' @param outFile     optional; path to RDA file where results should be saved
#' @param forceRecalc logical, default = FALSE; calculate population fractions from data even if
#'                    outFile exists from a previous run
#'
#' @return tibble of cell population fractions
getPopulationFractions <- function(counts, conditions, calcUnit = "FOV_ID", 
                                    outFile = NULL, forceRecalc = FALSE){

    if(fileDone(outFile) && !forceRecalc){
        log_info(paste0("Loading pre-computed fractions from file ",outFile))
        return(readRDS(outFile))
    }

    log_info("No pre-computed fractions file found and/or forced rerun turned ON. Calculating fractions now...")

    subCounts <- counts %>% select(Subpopulation = Population, SubCount = Count, everything())
    popCounts  <- counts %>% select(PopCount = Count, everything())

    fracs <- conditions %>%
             left_join(subCounts, by = intersect(names(.), names(subCounts))) %>%
             left_join(popCounts, by = intersect(names(.), names(popCounts))) %>%
             mutate(Fraction = SubCount/PopCount)

    if(!is.null(outFile)){
        saveRDS(fracs, file = outFile)
    }

    fracs
}

