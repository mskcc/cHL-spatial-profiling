#' Generate one or more fraction condition labels by pasting pairwise the elements
#' in two vectors with a specified delimiter in between.
#'
#' Given two vectors, one for numerator values and one for denominator values,
#' paste together the values at each corresponding position, with any specified
#' delimiter in between the two.
#'
#' @param numerator    vector of character strings representing numerator cell 
#'                     states (subpopulations)
#' @param denominator  vector of character strings representing denominator cell
#'                     states (populations)
#' @param delim        character to use as delimiter between numerator strings and 
#'                     denominator strings
#' 
#' @return  vector of character strings (complete fraction labels)
getFractionConditionLabel <- function(numerator, denominator, delim = "|"){
    if(!length(numerator) == length(denominator)){
        log_warn("length differ between input vectors. results may be incorrect.")
    }
    sapply(1:length(numerator), function(x){
        #cnd <- gsub(paste0(denominator[x], ","), "", numerator[x]) 
        cnd <- numerator[x]
        pop <- denominator[x]
        paste(cnd, delim, pop)
    }) %>% unlist()
}


#' 
getDensityConditionLabel <- function(value){
    sapply(1:length(value), function(x){
        ifelse(grepl(",", value[x]),
            paste0(unlist(strsplit(value[x], ","))[-1], collapse = ","),
            "All")
    }) %>% unlist()
}

getCountConditionLabel <- function(value){
    value
}



assignClasses <- function(dat, facets, cellTypes, facetOrder = NULL){
    if(is.null(facetOrder)){
        facetOrder <- lapply(facets, function(x){ unique(cellTypes[[x]]) })
        names(facetOrder) <- facets
    }
    grps <- dat %>%
            select_at(unique(c(facets, "Subpopulation", "Population"))) %>%
            mutate(ct = extractCT(Subpopulation))
    for(classType in facets[facets %in% names(cellTypes)]){
        grps <- grps %>%
                mutate(!!as.name(classType) := getClassifier(ct, cellTypes, classType))
        grps[[classType]][is.na(grps[[classType]])] <- ""
    }
    grps <- grps %>% select(-ct) %>% unique()

    for(classType in facets){
        grps[[paste0(classType, "_order")]] <- customOrder(grps[[classType]], facetOrder[[classType]])
        grps[[classType]] <- factor(grps[[classType]], levels = facetOrder[[classType]])
    }
    dat %>%
    select_at(names(dat)[!names(dat) %in% facets]) %>%
    left_join(grps, by=intersect(names(.), names(grps))) %>%
    arrange_at(paste0(facets, "_order"))
}




getCondOrderByPopulation <- function(dat, popOrder = NULL){
    if(!is.null(popOrder)){
        vals <- extractCT(dat$Population)
        dat$popOrder <- customOrder(vals, popOrder)
    } else {
        dat$popOrder <- as.numeric(dat$`Cell State ID`)
    }
    dat
}



setConditionOrder <- function(conds, ids, cellTypes, statsFile = NULL, sheet = NULL,
                              orderBy = "Cell State ID", calcType = "fractions",
                              facetY = NULL, facetOrder = NULL, popOrder = NULL, idOrder = NULL){

    ids <- as.numeric(ids)

    ## eliminate any irrelevant columns
    cnds <- conds %>%
            mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
            filter(`Cell State ID` %in% ids) %>%
            select(all_of(c("Cell State ID", "Subpopulation", "Population")), 
                   any_of(c("Cell_type", "Subtype", facetY))) %>%
            mutate(Subpopulation = formatMarkers(Subpopulation),
                   Population = formatMarkers(Population)) %>%
            left_join(cellTypes %>% select(Category, dplyr::matches("type")) %>% unique,
                      by = c("Cell_type", "Subtype"))

    ## data already contains Cell_type column associated with NUMERATOR condition, 
    ## but we need to change that in this case to match the DENOMINATOR condition
    if(!is.null(facetY)){
        cnds <- assignClasses(cnds, facetY, cellTypes, facetOrder = facetOrder) %>% unique()
    } else {
        cnds$facet <- NA
    }

    if(calcType == "fractions"){
        cnds <- cnds %>%
                mutate(Cond = getFractionConditionLabel(Subpopulation, Population))
    } else if(calcType == "densities"){
        cnds <- cnds %>%
                mutate(Cond = getDensityConditionLabel(Subpopulation))
    }

    if(!is.null(idOrder)){
        cnds <- cnds %>%
                mutate(ID_order = customOrder(as.character(`Cell State ID`), as.character(idOrder))) %>%
                arrange(desc(ID_order)) %>%
                mutate(labelY = row_number()) %>%
                group_by_at(facetY) %>%
                mutate(y = row_number()) %>%
                ungroup() %>%
                select_at(c("Cell State ID", "Subpopulation", "Population", "Cond", facetY, "y",
                         "Cell_type_figures", "Cell_type_figures_subscript",
                         "Subtype_figures", "Subtype_figures_subscript", "labelY"))
        cnds$y <- factor(cnds$y, levels = unique(cnds$y))
        return(cnds)
    }

    ## if rows are to be ordered based on denominator population, we need a new column
    ## on which we can sort
    if(!is.null(popOrder)){
        cnds <- cnds %>% getCondOrderByPopulation(popOrder = popOrder)
    }
  
    if(!orderBy %in% names(cnds)){
        if(any(is.null(c(statsFile, sheet)))){
            stop(paste0("Can not order by ", orderBy, ". Did you forget to provide statsFile or sheet?"))
        }
        dat <- read.xlsx(statsFile, sheet, check.names = F) %>% as_tibble()
        cnds <- cnds %>%
                left_join(dat %>%
                          mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
                          select_at(c("Cell State ID", orderBy)),
                          by = intersect(names(.), c("Cell State ID", orderBy)))
    }

orderCols <- orderBy
if(!is.null(facetY)){
    orderCols <- c(paste0(facetY, "_order"), orderCols)
}
    cnds <- cnds %>%
            arrange_at(vars(orderCols), list(~ desc(.))) #%>%
if(!is.null(facetY)){
    cnds <- cnds %>%
            group_by_at(facetY)
}
cnds <- cnds %>%
            mutate(y = row_number()) %>%
            ungroup() %>%
            mutate(labelY = row_number()) 

    cnds$labelY <- factor(cnds$labelY, levels = unique(cnds$labelY))
    cnds$y <- factor(cnds$y, levels = unique(cnds$y))
    cnds

}



formatConditionsForPlotting <- function(conds, includeIDs, calcType, cellTypes, orderBy = NULL,
                                        statsFile = NULL, sheet = NULL, idOrder = NULL,
                                        facetY = NULL, facetOrder = NULL, popOrder = NULL){

    cnds <- setConditionOrder(conds, includeIDs, cellTypes, orderBy = orderBy,
                              statsFile = statsFile, sheet = sheet, idOrder = idOrder,
                              facetY = facetY, facetOrder = facetOrder) %>%
            gather(c("Cond", "Cell State ID"), key = "Column", value = "Value") %>%
            mutate(hjust = ifelse(Column == "Cond", 0.5, 1))

    if(calcType == "densities"){
        cnds$Value = gsub(" \\|.*", "", cnds$Value)
    }

    cnds$Column <- factor(cnds$Column, levels = c("Cond", "Cell State ID"))
    cnds
}



formatForOverallvsIndiv <- function(overallFile, indivFiles, idsToPlot, max.fdr = 0.05,
                                   dataCol = "Fraction", effectCol = "Odds Ratio",
                                   statUnit = "Patient"){

    sheet <- "all_fractions_and_densities"
    signifCol <- paste(dataCol, " adjusted p.value")
    pvalCol   <- paste(dataCol, " p.value")
    overall <- read.xlsx(overallFile, sheet, check.names = F) %>% as_tibble() %>%
               mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
               filter(`Cell State ID` %in% idsToPlot)

    countCols <- names(overall)[grepl("total cell count", names(overall))]
    medCols   <- names(overall)[grepl(paste("median",tolower(dataCol)), names(overall))]

    overall <- overall %>%
               select(`Cell State ID`, `Cell State`, Population,
                      OverallSubpopulationCount = !!as.name(countCols[1]),
                      OverallPopulationCount = !!as.name(countCols[2]),
                      OverallGroup1Median = !!as.name(medCols[1]),
                      OverallGroup2Median = !!as.name(medCols[2]),
                      OverallEffect = !!as.name(effectCol),
                      OverallPval = !!as.name(pvalCol),
                      OverallFDR = !!as.name(signifCol))


    match_status <- function(indivEffect, allEffect, indivFDR, max.fdr = 0.05){
        lapply(1:length(indivEffect), function(x){
            if(is.na(indivEffect[x])){ return("stats unavailable") }
            if(is.na(allEffect[x])){ return(NA) }

            if(all(c(indivEffect[x], allEffect[x]) > 1) || all(c(indivEffect[x], allEffect[x]) < 1)){
                ifelse(indivFDR[x] < max.fdr,
                       "LOR same direction, FDR < 0.05",
                       "LOR same direction, FDR \u2265 0.05")
            } else {
                "LOR opposite direction"
            }

         }) %>% unlist()
    }
    get_value <- function(indivEffect, allEffect){
        lapply(1:length(indivEffect), function(x){

            if(is.na(indivEffect[x])){ return(NA) }
            if(is.na(allEffect[x])){ return(NA) }

             if(all(c(indivEffect[x], allEffect[x]) > 1) || all(c(indivEffect[x], allEffect[x]) < 1)){
                 indivEffect[x]
             } else {
                 -1
             }

        }) %>% unlist()
    }

    dat <- tibble()
    for(ind in indivFiles){
        #print(ind)
        id <- gsub(".xlsx", "", basename(ind))
        prfx <- "Indiv"
        indivEffect <- "IndivEffect"
        indivSignif <- "IndivFDR"
        indivPval   <- "IndivPval"
        overallEffect <- "OverallEffect"

        stats <- read.xlsx(ind, sheet, check.names = F) %>%
                 as_tibble() %>%
                 mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
                 filter(`Cell State ID` %in% idsToPlot) %>%
                 rename(!!as.name(indivEffect) := effectCol,
                        !!as.name(indivPval) := pvalCol,
                        !!as.name(indivSignif) := signifCol) %>%
                 left_join(overall, by = intersect(names(overall), names(.))) %>%
                 mutate(!!as.name(statUnit) := id,
                        Status = match_status(!!as.name(indivEffect),
                                              !!as.name(overallEffect),
                                              !!as.name(indivSignif),
                                              max.fdr = max.fdr),
                        Value = get_value(!!as.name(indivEffect),
                                          !!as.name(overallEffect)),
                        signif = ifelse(!!as.name(indivSignif) < max.fdr, "*", NA)
                 ) %>%
                 select_at(c("Cell State ID", "Cell State", "Population",
                             names(overall)[grepl("Overall", names(overall))],
                             indivEffect, indivPval, indivSignif,
                             statUnit, "Status", "Value", "signif"))

        dat <- dat %>% bind_rows(stats)
    }

    dat$Status <- factor(dat$Status,
                         levels = c("LOR same direction, FDR < 0.05",
                                    "LOR same direction, FDR \u2265 0.05",
                                    "LOR opposite direction",
                                    "stats unavailable"))
    dat

}

formatStatsForPlottingEffects <- function(statsFile, sheet, calcType, dataCol, allConds, ids,
                                          facetY = NULL, facetOrder = NULL, cellTypes = NULL,
                                          popOrder = NULL, orderBy = NULL, idOrder = NULL){

    cnds <- setConditionOrder(allConds, ids, cellTypes, orderBy = orderBy,
                              statsFile = statsFile, sheet = sheet, idOrder = idOrder,
                              facetY = facetY, facetOrder = facetOrder)

    effect <- ifelse(dataCol == "Fraction", "Odds Ratio", "Fold Change")
    stats <- read.xlsx(statsFile, sheet, check.names = F) %>%
             as_tibble() 

    fdrCol  <- ifelse("adjusted p.value" %in% names(stats), "adjusted p.value", paste0(dataCol, " adjusted p.value"))
    ci.high <- ifelse("CI.high" %in% names(stats), "CI.high", paste0(dataCol, " CI.high"))
    ci.low  <- ifelse("CI.low" %in% names(stats), "CI.low", paste0(dataCol, " CI.low"))
    stats <- stats %>%
             select(`Cell State ID`, `Cell State`, Population, dplyr::matches(dataCol),
                    !!as.name(effect), !!as.name(fdrCol), 
                    !!as.name(ci.low), !!as.name(ci.high)) %>%
             mutate(Condition = `Cell State`, `Cell State ID` = as.numeric(`Cell State ID`)) %>%
             filter(`Cell State ID` %in% ids)

    if(is.null(stats)){ stop("Couldn't find appropriate sheet in stats file") }

    names(stats) <- gsub(paste0("^",dataCol," "),"",names(stats))

    if(nrow(stats) < 1){
        log_warn("    WARNING: no statistics available for this question.")
        return(NULL)
    }

    ## format y labels
    stats <- stats %>%
             unique() %>%
             ungroup() %>%
             select(-Condition, -Population) %>%   ## do this because they may have changed for facetting purposes
             mutate(signif = signifStars(`adjusted p.value`))

    stats %>% left_join(cnds, by = intersect(names(.), names(cnds)))
}



formatFOVdataForPlottingDetail <- function(question, sGrps, condsToPlot, allConds, analyses,
                                           calcType, dataCol, markers, metricsDir, orderBy = NULL,
                                           statsFile = NULL, sheet = NULL,
                                           nbhdCounts = NULL, calcUnit = "FOV_ID",
                                           facetY = NULL, facetOrder = NULL, cellTypes = NULL,
                                           popOrder = NULL, idOrder = NULL){

    log_debug("Setting condition order")
    condsToPlot <- as.numeric(condsToPlot)
    cnds <- setConditionOrder(allConds, condsToPlot, cellTypes, orderBy = orderBy,
                              statsFile = statsFile, sheet = sheet, idOrder = idOrder,
                              facetY = facetY, facetOrder = facetOrder)

    log_debug("Getting comparison data")
    dat <- sGrps %>%
           getComparisonData(question, metricsDir, filePattern = calcType, 
                             calcUnit = calcUnit) 
    if(is.null(dat)){
        log_warn("   WARNING: One or more groups in this question contains no data. Skipping.\n")
        return(NULL)
    }

    if(calcType == "densities"){
        dat <- dat %>%
               full_join(analyses$densities %>% select(`Cell State ID`, Population), by = "Population")
    }

    log_debug("Calculating Q1,Q3,Median & IQR")

    maxx <- ifelse(calcType == "fractions", 1, 10^4)
    minx <- ifelse(calcType == "densities", 10^-1, 0)

    selCols <- c("Cell State ID", facetY, "Subpopulation", "Population",
                 question$groupVar, calcUnit, dataCol)

    fovDat <- dat %>% 
              ungroup() %>%
              mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
              filter(`Cell State ID` %in% condsToPlot) %>%
              select_at(selCols[selCols %in% names(.)]) %>%
              group_by_at(c("Cell State ID", question$groupVar)) %>%
              mutate(Median = median(!!as.name(dataCol), na.rm=T),
                     Q1 = summary(!!as.name(dataCol), na.rm=T)[2],
                     Q3 = summary(!!as.name(dataCol), na.rm=T)[4],
                     IQR = Q3 - Q1,
                     min = max(Q1 - (1.5*IQR), minx, na.rm=T),
                     max = min(Q3 + (1.5*IQR), maxx, na.rm=T))

    if(!"Subpopulation" %in% names(fovDat)){
        fovDat <- fovDat %>% mutate(Subpopulation = Population)
    }

    ## adjust group labels
    grpLbls <- fovDat %>%
               ungroup() %>%
               select_at(c(calcUnit, question$groupVar)) %>%
               unique() %>%
               group_by_at(question$groupVar) %>%
               summarize(Count = n()) %>%
               mutate(GroupLabel = paste0(!!as.name(question$groupVar), " (FOVs=", Count, ")")) %>%
               select(-Count)

    grpLbls$GroupLabel <- factor(grpLbls$GroupLabel, levels = unique(grpLbls$GroupLabel))

    fovDat <- fovDat %>%
              left_join(grpLbls, by = c(question$grpVar)) %>%
              ungroup() %>%
              select(-any_of(c("Subpopulation", "Population", "facetY"))) %>%   ## do this because they may have changed for facetting purposes
              left_join(cnds) %>% #, by = intersect(names(.), names(cnds))) %>%
              filter(!is.na(!!as.name(dataCol)))

    fovDat

}



formatMedians <- function(dat, c2p, groups, calcCol, layout){
    meds <- dat %>%
            filter(`Cell State ID` %in% c2p) %>%
            select(`Cell State ID`, dplyr::matches("Overall Median")) %>%
            rename_at(vars(names(.)[grepl("Overall Median", names(.))]),
                      list(~ gsub(" Overall Median.*", "", .))) %>%
            gather(2:ncol(.), key = "Group", value = "raw median") 

    signif <- dat %>%
              filter(`Cell State ID` %in% c2p) %>%
              select(`Cell State ID`, dplyr::matches("Overall FDR")) %>%
              rename_at(vars(names(.)[grepl("Overall FDR", names(.))]),
                        list(~ gsub(" vs.*", "", .))) %>%
              gather(2:ncol(.), key = "Group", value = "FDR") %>%
              mutate(signif = ifelse(FDR < 0.05, "*", ""))

    meds <- meds%>%
            left_join(signif, by = c("Cell State ID", "Group")) %>%
            left_join(layout, by = "Cell State ID") %>%
            group_by(`Cell State ID`) %>%
            mutate(maxMed = max(`raw median`, na.rm=T),
                   median = `raw median`/maxMed,
                   y2 = as.numeric(y),
                   bottom = y2 - 0.5,
                   top = (y2 - 0.5) + median*0.7)

    medIdxs <- grep("Overall Median", names(meds))
    names(meds)[medIdxs] <- gsub(" Overall Median.*", "", names(meds)[medIdxs])

    meds$x     <- customOrder(meds$Group, groups)
    meds$Group <- factor(meds$Group, levels = groups)
    meds$y     <- factor(as.numeric(meds$y), levels = seq(0.5, max(as.numeric(meds$y), na.rm=T) + 0.5, 0.5))
    meds$signif <- factor(meds$signif, levels = c("*", ""))
    meds
}

formatIndivDatForPlotting <- function(dat, idMap, layout, xVar = "Sample_ID", xOrder = NULL){
    iDat <- dat %>%
            left_join(idMap, by = intersect(names(.), names(idMap))) %>%
            left_join(layout, by = intersect(names(.), names(layout))) %>%
            arrange(desc(labelY)) %>%
            mutate(status = ifelse(Value == -1, "opposite",
                                   ifelse(Value > 1, "up",
                                          ifelse(Value < 1, "down",
                                                 ifelse(Value == 1, "no change", NA)))))

    xVals <- intersect(unique(iDat[[xVar]]), xOrder)
    if(!all(xVals %in% xOrder)){
       log_warn("Some x variables included in xOrder are missing.")
    } else if(!all(xOrder %in% iDat[[xVar]])){
       log_warn("Some x variables included in data are NOT included in xOrder.")
    }
    xOrder <- xOrder[xOrder %in% xVals]

    iDat$x <- customOrder(iDat[[xVar]], xOrder) 
    iDat$x <- factor(iDat$x, levels = seq(xOrder)) 

    iDat
}


