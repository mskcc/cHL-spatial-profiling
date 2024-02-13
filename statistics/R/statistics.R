
na_population <- function(populations, na_func_markers){                        
    pattern = paste0("(^|,)(", paste(na_func_markers, collapse = ("|")), ")(,|\\-|$)")
                                                                                
    sapply(populations, function(x){                                            
        grepl(pattern, x)                                                       
    })                                                                          
}       

#' Adjust p value to 1 when median of both groups in a two-group comparison
#' is 0
#' 
#' For any condition where means/medians of both groups is 0, set pvals to 1
#' 
#' @param dat     data table including columns for 'p.value' and 'adjusted p.value'
#' @param column  vector of 2 column names, one for median of group 1 and one for median
#'                of second group
fixPvals <- function(dat, columns){
    idxs <- which(dat %>% ungroup() %>% select_at(columns) %>% rowSums() == 0)
    dat$`p.value`[idxs] <- 1
    dat$`adjusted p.value`[idxs] <- 1
    dat
}

#' Convert value to log odds
#' 
#' Convert value to log odds
#' 
#' @param x  numeric value to convert
#' @return log odds of x
logOdds <- function(x){
    log(x/(1-x))
}

#' Convert log odds to a fraction
#' 
#' Convert log odds to a fraction
#' 
#' @param lo   numeric value; log odds to be converted
#' @return fraction
logOddsToFraction <- function(lo){
    odds <- exp(lo)
    odds/(1+odds)
}

#' Convert log odds ratio to fraction
#' 
#' Convert log odds ratio to fraction
#'
#' @param x  numeric; value between 0 and 1
#' @return log odds ratio
fractionToLogOddsRatio <- function(x){
    log((1/x) -1)
}

expandFunctionalStates <- function(cellStates){

    states <- cellStates
    combos <- states %>% filter(grepl(",", state)) %>% pull(state)
    for(combo in combos){
        applyTo  <- states %>% filter(state == combo) %>% select(-state)
        toAdd    <- tibble(state = rev(unlist(expandMarkerCombos(combo, emptyCombo = TRUE))))
        expanded <- cbind(toAdd, applyTo)
        states <- states %>%
                  bind_rows(expanded) %>%
                  unique() ## need this because the original combo is now in tibble twice
    }

    states
}

addUniversalMarkers <- function(cellStates, markers = NULL){
    if(is.null(markers)){ return(cellStates) }
    markers <- unlist(strsplit(markers, ","))
    states <- cellStates
    for(marker in markers){
        states[[marker]] <- marker
    }
    states
}


configureConditions <- function(types, states, configYAML = NULL, force = FALSE){

    if(fileDone(configYAML) && force == FALSE){
        return(read_yaml(configYAML))
    }
    condCfg <- list()
    condCfg$annotation_order <- c("Cell State ID", "Cell_type", "Subtype", "Proliferation", "Function")
    condCfg$spatial_annotation_order <- c("Cell State ID", "Center Population A", "Neighborhood Population A",
                                          "Center Population B", "Neighborhood Population B")
    condCfg$arrange_annotation <- list(Category = c(NA, "All", "DAPI", unique(types$Category)),
                                       Cell_type = c(NA, "All", unique(types$Cell_type)),
                                       Subtype = c(NA, "All", unique(types$Subtype[types$Subtype != "All"])),
                                       Proliferation = c(NA,"","KI67"),
                                       Function = c(NA,"",unique(states$state)))
    if(!is.null(configYAML)){
        write_yaml(condCfg, configYAML)
    }
    condCfg
}


getAllSampleFilters <- function(questions){
    lapply(stDat$allQuestions, function(x){ 
        allF <- unique(c(names(x$`Group 1`), names(x$`Group 2`))) 
    }) %>% 
    unlist %>% 
    unique
}

getAllDataAnnotation <- function(dataFiles, sampAnn, allFilters, calcUnit = "FOV_ID", threads = 1,
                                   filter_stats_exclusions = FALSE){

    if(filter_stats_exclusions & any(!is.na(sampAnn$Exclude_from_stats))){
        excl <- sampAnn %>% 
                filter(Exclude_from_stats == "X") %>%
                select(CellDive_ID, FOV_number)
        for(ex in unique(excl$CellDive_ID)){
            fovs <- excl %>% filter(CellDive_ID == ex) %>% pull(FOV_number)
            log_warn(paste0("Excluding FOVs [", paste(fovs, collapse = ","), "] from sample [", ex, "] from all stats!")) 
        }
        sampAnn <- sampAnn %>%
                   filter(is.na(Exclude_from_stats))
    }

    cl <- makeCluster(threads, type = "FORK", outfile = "")
    clusterExport(cl, c("sampAnn", "allFilters", "calcUnit"), envir = environment())

    ann <- parLapply(cl, dataFiles, function(df){
             suppressMessages( 
               readRDS(df) %>%
               select_at(c("CellDive_ID", calcUnit, intersect(names(.), allFilters))) %>%
               unique %>%
               left_join(sampAnn %>% select_at(c(calcUnit, intersect(names(.), allFilters)))) 
             )
           }) %>%
           bind_rows %>%
           unique

    stopCluster(cl)

    missing <- setdiff(allFilters, names(ann))
    if(length(missing) > 0){
        tmp <- lapply(missing, function(m){
                   log_debug(paste("No data for annotation [", m, "]."))
               })
    }

    ann
}


getSampleGroups <- function(question, sampAnn, calcUnit = "FOV_ID"){

    var <- question$groupVar

    sGrps = list()
    for(grp in c("Group 1", "Group 2")){
        filt <- question[[grp]]
        gSmps <- sampAnn %>% 
                 select_at(unique(c("CellDive_ID", calcUnit))) %>% 
                 distinct
        if(length(setdiff(names(filt), names(sampAnn))) > 0){
            log_warn("STILL NEED TO FILTER ON [", paste0(setdiff(names(filt), names(sampAnn)), collapse=", "), "]")
        }
        for(x in intersect(names(filt), names(sampAnn))){
            gSmps <- sampAnn %>% 
                     filter_at(vars(all_of(x)), all_vars(. %in% filt[[x]])) %>% 
                     select_at(unique(c("CellDive_ID", calcUnit))) %>%
                     filter_at(vars(calcUnit), all_vars(. %in% gSmps[[calcUnit]]))
        }
        sGrps[[paste(filt[[var]], collapse = "_or_")]] <- gSmps
    }

    sGrps

}

getGroupMicroenvData <- function(dataDir, groupInfo, nbhd, envir, calcType = "fractions", calcUnit = "FOV_ID"){

    fp <- paste0(calcType, ".*___(", paste0(nbhd, "|", gsub(",", "__", nbhd)), ").rda")
    files <- dir(dataDir, full.names = T, pattern = fp)

    cdids <- unique(groupInfo$CellDive_ID)                                      

    dat <- tibble()
    for(cdid in cdids){
        flNm <- files[grepl(paste0("^", cdid, "_"), basename(files))]
        if(!fileDone(flNm)){ log_error("No microenvironment file for ", cdid, ".") }
        gDat <- readRDS(flNm) %>%
                mutate(CellDive_ID = cdid) %>%
                separate(paste0(calcUnit, "___", nbhd), c(calcUnit, nbhd), sep="___", remove = F) %>%
                filter_at(vars(all_of(calcUnit)), all_vars(. %in% groupInfo[[calcUnit]])) %>%
                filter_at(all_of(nbhd), all_vars(. == envir))
        if(nrow(gDat) == 0){
            log_debug("No [", envir, "] found in sample ", cdid, ".")
            next
        }
        dat <- dat %>% bind_rows(gDat)
    }

    dat

}

getGroupSampleData <- function(dataDir, groupInfo, calcUnit = "FOV_ID", filePattern = ".rda"){

    files <- dir(dataDir, full.names = T, pattern = filePattern)

    cdids <- unique(groupInfo$CellDive_ID)
    
    dat <- tibble()
    for(cdid in cdids){
        flNm <- files[grepl(paste0("^",cdid,"_"), basename(files))]
        if(!fileDone(flNm)){ 
            log_warn("No file for ",cdid,".") 
            next
        }
        cdid_dat <- readRDS(flNm) %>%
                    mutate(CellDive_ID = cdid) %>%
                    filter_at(vars(all_of(calcUnit)), all_vars(. %in% groupInfo[[calcUnit]]))
        dat <- dat %>%
               bind_rows(cdid_dat)
    }

    dat

}


groupNameToNumber <- function(question, grpName){
    vals <- unique(unlist(strsplit(grpName, "_or_")))
    quVec <- unlist(question)
    unique(gsub("\\..*", "", names(quVec[quVec %in% vals])))
}


getComparisonData <- function(groups, question, dataDir, filePattern = ".rda",
                               calcUnit = "FOV_ID"){

    groupVar <- question$groupVar
    lapply(names(groups), function(x){
        grp <- groupNameToNumber(question, x)
        nbhdGrp <- "Neighborhood" %in% names(question[[grp]])
        grpDat <- NULL
        if(nbhdGrp){
            nbhd <- question[[grp]]$Neighborhood
            env <- question[[grp]]$Environment
            grpDat <- getGroupMicroenvData(dataDir, groups[[x]], nbhd, env, 
                                           calcUnit = calcUnit) %>% 
                      mutate(!!as.name(groupVar) := x) 
        } else {
            grpDat <- getGroupSampleData(dataDir, groups[[x]], calcUnit = calcUnit, 
                                         filePattern = filePattern) %>%
                      mutate(!!as.name(groupVar) := x)
        }
        log_debug(paste0(nrow(grpDat), " cells in group [", x, "]"))
        grpDat
    }) %>%
    bind_rows %>%
    mutate(!!as.name(groupVar) := factor(!!as.name(groupVar), levels = rev(names(groups))))

}


#' Run statistical test on one or more columns of data in a table
#'
#' Given a table containing a column or columns of groups (e.g., 'Lesion Response') and 
#' one or more columns of data, run specified statistical test on one or more
#' of the data columns
#' 
#' @param fun        statistical function to run (default=t.test)
#' @param df         data table
#' @param vars       names of data column(s) on which to run the statistical test
#' @param groupVar   name of column containing the groups each value in data column belongs
#'                   to; column must contain exactly two unique values and will be converted
#'                   to a factor whose levels are sorted alphanumerically
#'
#' @return a list of stats results where each element contains results for one item in {vars}
#' @export
multi.tests <- function(fun = t.test, df, vars, groupVar, ...) {
    if(!is.factor(df[[groupVar]])){
        df[[groupVar]] <- factor(df[[groupVar]], levels=unique(sort(df[[groupVar]])))
        log_debug(paste0("Factor levels and order not set for ", groupVar,". Setting order to : ", 
                  paste0(levels(df[[groupVar]]), collapse=", ")))
    }
    groupVar <- paste0("`",groupVar,"`")

    lst <- list()
    fails <<- 0
    for(x in 1:length(vars)){
        var <- vars[x]
        newName <- paste0("`",var,"`")

        formula <- as.formula(paste(newName, "~", groupVar))
        lst[[var]] <- tryCatch({
                         fun(data = df %>% filter(!is.na(!!as.name(var))), formula, conf.int=TRUE)
                      }, error = function(e){
                         fails <<- fails + 1
                         log_debug(paste0("[",fails,"] Could not calculate stats for ",var))
                         NULL
                      })
    }
    return(lst)
}


#' Format numbers for XLSX output
#' 
#' Format numbers for XLSX output
#' 
#' @param df   data frame to be formatted
#'
#' @return data frame with all numeric values formatted
formatFinalNumbers <- function(df){
    types <- sapply(df, class)

    ints <- unique(c(names(df)[grepl("count",names(df),ignore.case=TRUE)],
                     names(types)[which(types == "integer")]))
    nums <- unique(c(names(df)[grepl("median|density$|area|estimate|low|high|p\\.value|ratio|mean",
                               names(df),
                               ignore.case=TRUE)],
                     names(types)[which(types == "numeric")]))

    nums <- nums[!nums %in% ints]
    for(n in nums){
        df[[n]] <- as.numeric(df[[n]])
        df[[n]] <- as.numeric(formatC(df[[n]], digits=3))
    }

    return(df)
}


#' Write stats results for a single question and a single cell region to XLSX file
#' 
#' Write stats results for a single question and a single cell region to XLSX file
#' 
#' @param tblList                list of stats results where each item is a table of results 
#'                               for one calculation type (e.g., densities or fractions), and/or 
#'                               sample information; each item will be written to its own tab
#' @param fileName               XLSX file to which results will be saved
#' @param filterName             name of filter used (simply for naming tabs of filtered data)
#'
#' @return nothing
writeStatsQuestionXLSX <- function(tblList, fileName){
    countStyle <- openxlsx::createStyle(numFmt = 'COMMA')
    wb = openxlsx::createWorkbook()
    for(i in seq(tblList)){
        tblName <- names(tblList)[i]
        tbl <- tblList[[tblName]]
        if(nrow(tbl) > 0){
            tbl <- formatFinalNumbers(tbl)
        }
        openxlsx::addWorksheet(wb, tblName)
        openxlsx::writeData(wb, i, tbl)
        openxlsx::freezePane(wb, tblName, firstRow = TRUE)
        countColIdxs <- grep("count", names(tbl), ignore.case=T)
        if(length(countColIdxs) > 0){
            for(cc in countColIdxs){
                openxlsx::addStyle(wb, i, style=countStyle, cols=cc, rows=1:nrow(tbl)+1)
            }
        }
    }
    openxlsx::saveWorkbook(wb, fileName, overwrite=TRUE)
}


pseudoScale <- function(dat, groupVar, valueCol){
    tmp <- dat
    ## assume fractions at first
    tmp$OM <- 1 - tmp[[valueCol]]
    tmp <- tmp %>%
           ungroup() %>%
           group_by_at(groupVar) %>%
           filter_at(vars(all_of(valueCol)), any_vars(. > 0 & . < 1))
    if(nrow(tmp) > 0){
        ps <- tmp %>%
              summarize_at(vars(all_of(valueCol),"OM"), min) %>%
              select_at(c(valueCol,"OM")) %>%
              min()/2

        scaled <- dat[[valueCol]]
        scaled[scaled == 0] <- ps
        scaled[scaled == 1] <- 1 - ps
    } else {  ## nothing in between 0 and 1, assume densities
        scaled <- dat[[valueCol]] + 0.001
    }
    scaled
}


getStatsLO <- function(loDat, groupVar, varCols, func = wilcox.test){
    stats <- multi.tests(fun = func, df = loDat, vars = varCols, groupVar = groupVar) %>%
             wilcoxResToTable(., log2exp = TRUE) %>%
             dplyr::rename(`Odds Ratio` = estimate) %>%
             mutate(`log(Odds Ratio)` = log(as.numeric(`Odds Ratio`)),
                    `Inverse Odds Ratio` = 1/as.numeric(`Odds Ratio`))#,
    stats$`adjusted p.value` <- NA
    tum <- grep('HRS', stats$CondTitle)
    stats$`adjusted p.value`[tum] <- p.adjust(stats$`p.value`[tum], method = "bonferroni")

    imm <- which(!grepl('HRS', stats$CondTitle))
    stats$`adjusted p.value`[imm] <- p.adjust(stats$`p.value`[imm], method = "bonferroni")     

    stats
}

#' Convert fractions to log odds
#' 
#' Given a table of fractions and group assignments, first pseudoscale
#' fractions to handle values of 0 and 1, then convert pseudoscaled
#' values to log odds
#' 
#' @param fracDat  table of fractions and sample group assignments
#' @param fracCol  name of column containing fractions
#' @param groupVar name of column containing sample group assignments
#'
#' @return fracDat table with additional columns, psFraction and LO
fractionsToLO <- function(fracDat, fracCol, groupVar){
    dat <- fracDat
    dat$psFraction <- pseudoScale(dat, groupVar, fracCol)
    dat %>% mutate(LO = logOdds(psFraction))
}


#' Get degrees of freedom for each cell state in table of fractions or densities
#'
#' Given a table of cell state IDs, calculation units (e.g., samples or FOVs), and
#' a data value for each (e.g., fractions or densities), count the number of
#' calculation units that have a non-NA value and subtract one to get the degrees
#' of freedom.
#'
#' @param dat      table of cell state IDs, calculation units and data values
#' @param datCol   name of column containing data value
#' @param calcUnit name of column containing calculation unit; default = 'FOV_ID'
#'
#' @return two column table with Cell State ID and DF
getCellStateDegreesOfFreedom <- function(dat, datCol, calcUnit = "FOV_ID"){
    dat %>%
    filter_at(vars(all_of(datCol)), all_vars(!is.na(.))) %>%
    select_at(c("Cell State ID", calcUnit)) %>%
    unique %>%
    group_by(`Cell State ID`) %>%
    summarize(DF = n() - 1)
}

#' Run wilcox test and format log odds statistics 
#' 
#' Run wilcox test and format log odds statistics
#'
#' @param fracDat      pre-filtered table of fraction data with, at minimum, columns: Condition, Population,
#'                     calculation unit ("FOV_ID" by default), Fraction and groupVariable (e.g., Patient_response) 
#' @param groupVar     string representing name of variable to be compared; will be either a column name from sampAnnot
#'                     or arbitrary string when comparing non-clinical variable (e.g., inside/outside, pos/neg microenv)
#' @param calcUnit     data column name containing values for which a single area value
#'                     should be calculated (default: Sample)
#' @param max.fdr      maximum adjusted p.value to consider a condition 'Passed'
#' @param min.or       minimum odds ratio to consider a condition 'Passed'
#' @param  full fraction statistics table
reportLogOdds <- function(fracDat, groupVar, calcUnit = "FOV_ID", fracCol = "Fraction"){

    dat <- fracDat %>% mutate(CondTitle = fractionCond(Subpopulation, Population))

    df     <- dat %>% getCellStateDegreesOfFreedom(fracCol, calcUnit = calcUnit)
    counts <- dat %>% getTotalCountsForLOReport()
    meds   <- dat %>% getGroupMedianFractions(groupVar, fracCol)

    ## transform fractions to log odds for stats
    loTbl <- dat %>%
             fractionsToLO(fracCol, groupVar) %>%
             select_at(c("CondTitle", calcUnit, groupVar, "LO")) %>%
             spread(CondTitle, LO)

    stats  <- loTbl %>% getStatsLO(groupVar, unique(dat$CondTitle))

    ## compile report
    counts %>%
    left_join(df, by = "Cell State ID") %>%
    left_join(meds, by = intersect(names(.), names(meds))) %>%
    left_join(stats, by = intersect(names(.), names(stats))) %>%
    fixPvals(., paste(unique(dat[[groupVar]]), "median fraction")) %>%
    select(`Cell State ID`, `Cell State` = Subpopulation, Population,
           DF, `Cell State total cell count`, `Population total cell count`,
           contains("median fraction"),
           `log(Odds Ratio)`, `Odds Ratio`, `Inverse Odds Ratio`,
           everything()) %>%
    select(-CondTitle) %>%
    arrange(as.numeric(`Cell State ID`)) %>%
    ungroup()

}


#' Convert list of results from multiple wilcox tests to a table
#' 
#' Results from multi.tests() are returned in list form, each element containing
#' raw test results for one condition. When test is wilcox.test, convert that
#' list to a table with columns CondTitle, estimate, CI.low, CI.high and p.value.
#' 
#' @param statsRes   results of multiple wilcox.test runs in list form
#' @param log2exp    logical; if tests were run on log values, set this value to TRUE to
#'                   convert estimate and confidence interval values back to normal scale
#' @return table of all wilcox test results, one row per test 
wilcoxResToTable <- function(statsRes, log2exp=FALSE){
    res <- tibble()
    for(cond in names(statsRes)){
        sr <- statsRes[[cond]]
        row <- tibble(CondTitle = cond,
                      estimate = sr$estimate,
                      CI.low = sr$`conf.int`[1],
                      CI.high = sr$`conf.int`[2],
                      p.value = sr$`p.value`)
        if(log2exp){
            row <- row %>%
                   mutate(estimate = exp(estimate),
                          CI.low = exp(CI.low),
                          CI.high = exp(CI.high))
        }
        res <- res %>% bind_rows(row)
    }
    res
}


fractionCond <- function(conds, pops){
    lapply(1:length(conds), function(x){
        gsub("/| ","_",paste(conds[x], pops[x], sep="___"))
    }) %>%
    unlist()
}

getTotalCountsForLOReport <- function(allCountDat){
    allCountDat %>%
    group_by(`Cell State ID`, Subpopulation, Population, CondTitle) %>%
    summarize(`Cell State total cell count` = sum(SubCount, na.rm=T),
              `Population total cell count` = sum(PopCount, na.rm=T)) %>%
    unique()
}

getGroupMedianFractions <- function(allCountDat, groupVar, fracCol){
    groupNames <- unique(as.vector(allCountDat[[groupVar]]))

    allCountDat %>%
    group_by_at(c("CondTitle", groupVar)) %>%
    summarize(Median = median(!!as.name(fracCol), na.rm=T)) %>%
    spread_(groupVar, "Median", fill = 0) %>%
    dplyr::rename_at(all_of(groupNames), list(~ paste(., "median fraction")))
}


filterLOreport <- function(report, min.subpop.count = 300, min.pop.count = 1000,
                           min.median.frac.diff = 0.1, max.fdr = 0.05){

    report %>%
    bioFilter(calc = "fractions",
              min.subpop.count = min.subpop.count,
              min.pop.count    = min.pop.count,
              min.median.diff  = min.median.frac.diff,
              max.fdr          = max.fdr,
              subpop.col       = "Cell State total cell count",
              pop.col          = "Population total cell count",
              fdr.col          = "adjusted p.value",
              median.cols      = names(report)[grepl("median", names(report))])

}

getStatsFC <- function(logDat, groupVar, varCols, func = wilcox.test){
    stats <- multi.tests(fun = func, df = logDat, vars = varCols, groupVar = groupVar) %>%
             wilcoxResToTable(., log2exp = TRUE) %>%
             dplyr::rename(`Fold Change` = estimate) %>%
             mutate(`log(Fold Change)` = log(as.numeric(`Fold Change`)),
                    `Inverse Fold Change` = 1/as.numeric(`Fold Change`)) #,

    ## adjust pvals for HRS and immune cells separately since there are so 
    ## few HRS cells compared to immune
    stats$`adjusted p.value` <- NA
    tum <- grep('HRS', stats$CondTitle)
    stats$`adjusted p.value`[tum] <- p.adjust(stats$`p.value`[tum], method = "bonferroni")

    imm <- which(!grepl('HRS', stats$CondTitle))
    stats$`adjusted p.value`[imm] <- p.adjust(stats$`p.value`[imm], method = "bonferroni")

    stats

}


getTotalCountsForDensityFCReport <- function(dat){
    dat %>%
    select(`Cell State ID`, CondTitle, Count, TotalArea, Density) %>%
    group_by(`Cell State ID`, CondTitle) %>%
    summarize(`Cell State total cell count` = sum(Count, na.rm = T),
              `Total Area` = sum(TotalArea, na.rm = T),
              `Overall Median Density` = median(Density, na.rm = T))
}


#' Get median density of each sample group for every cell state
#'
#' Given a table of cell state density values for each calculation 
#' unit (generally FOV), and a group assignment for each unit, 
#' calculate the median density value in each group for every cell state.
#'
#' @param dat       table containing columns 'Cell State ID', 'CondTitle', and 
#'                  a column containing sample group assignments
#' @param groupVar  name of column that contains sample group assignments
#' 
#' @return table in 'dat' with two columns added, each containing median
#'         density values for one of the two sample groups
getGroupMedianDensities <- function(dat, groupVar){
    groupNames <- unique(as.vector(dat[[groupVar]]))

    dat %>%
    group_by_at(c("Cell State ID", "CondTitle", groupVar)) %>%
    summarize(Median = median(Density, na.rm=T)) %>%
    spread(groupVar, Median, fill=0) %>%
    dplyr::rename_at(all_of(groupNames), list(~ paste(., "median density")))
}


#' Remove rows of duplicate conditions from conditions index
#'
#' Conditions index contains some duplicates in the case of densities
#' and neighborhood averages as a result of joining with multiple
#' fractions conditions; Here we find and remove those duplicates
#'
#' @param dat         tibble that possibly contains duplicate conditions
#' @param uniqueCols  vector of column names required to distinguish unique conditions
#' 
#' @return tibble with duplicate conditions removed
removeDuplicateRows <- function(dat, uniqueCols){
    tmp <- dat %>% select_at(uniqueCols)
    if(any(duplicated(tmp))){
        log_debug(paste0("Removing ",length(which(duplicated(tmp))),
                        " rows with duplicate conditions before analysis."))
        dupIDs <- sort(as.numeric(unique(dat$`Cell State ID`[which(duplicated(tmp))])))
        log_debug(paste("   IDs removed: ", paste(dupIDs, collapse=",")))
        dat <- dat[-which(duplicated(tmp)),]
    }
    dat
}

getDensityConditions <- function(condIndex){
    condIndex %>% 
    removeDuplicateRows("Subpopulation") %>% 
    select(`Cell State ID`, Population = Subpopulation)
}

#' Run wilcox test and format log fold change statistics 
#' 
#' Run wilcox test and format log fold change statistics
#'
#' @param denDat       table of data with following columns at minimum: 
#'                     Population, Samp1, Samp2, ..., SampN
#'                        OR
#'                     Population, Sample, Count, Area, Density, [Group] 
#' @param allSamples   a vector containing the union of group1 samples and group2 samples ONLY
#' @param groupVar     string representing name of variable to be compared; will be either a column name from sampAnnot
#'                     or arbitrary string when comparing non-clinical variable (e.g., inside/outside, pos/neg microenv)
#' @param calcUnit     data column name containing values for which a single area value
#'                      should be calculated (default: Sample)
#' @param max.fdr      maximum adjusted p.value for 'passing' results
#' @param min.fc       minimum absolute fold change for 'passing' results
#' 
#' @return tibble of density stats including mean densities per group, fold changes, confidence intervals,
#'         p values and FDRs for all conditions
reportDensityStats <- function(denDat, groupVar, calcUnit = "FOV_ID"){

    dat <- denDat %>%
           rename(CondTitle = Population) %>%
           removeDuplicateRows(c("CondTitle", calcUnit, groupVar))

    dat$psDensity <- pseudoScale(dat, groupVar, "Density")
    dat <- dat %>% mutate(logDen = log(psDensity))

    logDens <- dat %>%
               select_at(c(groupVar, "CondTitle", calcUnit, "logDen")) %>%
               spread(CondTitle, logDen)

    df     <- dat %>% getCellStateDegreesOfFreedom("Density", calcUnit = calcUnit)
    counts <- dat %>% getTotalCountsForDensityFCReport()
    stats  <- logDens %>% getStatsFC(groupVar, unique(dat$CondTitle))
    meds   <- dat %>% getGroupMedianDensities(groupVar)

    counts %>%
    left_join(df, by = "Cell State ID") %>%
    left_join(meds, by = intersect(names(.), names(meds))) %>%
    left_join(stats, by = intersect(names(.), names(stats))) %>%
    fixPvals(., paste(unique(dat[[groupVar]]), "median density")) %>%
    select(`Cell State ID`,
           `Cell State` = CondTitle,
           DF,
           `Cell State total cell count`,
           `Total Area`,
           `Overall Median Density`,
           contains("median density"),
           `log(Fold Change)`,
           `Fold Change`,
           `Inverse Fold Change`,
           everything()) %>%
    ungroup() %>%
    arrange(as.numeric(`Cell State ID`))

}


filterFCreport <- function(report, min.subpop.count = 300, max.fdr = 0.05,
                           min.median.diff = 0.1){

    report %>%
    bioFilter(calc = "densities",
              min.subpop.count = min.subpop.count,
              min.pop.count    = 0,
              min.median.diff  = min.median.diff,
              max.fdr          = max.fdr,
              subpop.col       = "Cell State total cell count",
              pop.col          = "Population total cell count",
              fdr.col          = "adjusted p.value",
              median.cols      = names(report)[grepl("median", names(report))])

}




#' Filter statistics results for biologically significant results
#' 
#' Remove from statistics results any conditions that do not pass
#' criteria including minimum population and subpopulation counts,
#' minimum difference in group medians and FDR.
#'
#' @param dat               tibble of statistics to be filtered
#' @param calc              calulation type ['fractions'|'densities'] (default: fractions)
#' @param min.subpop.count  minimum subpopulation cell count (default: 300)
#' @param min.pop.count     minimum population cell count (default: 1000)
#' @param min.median.diff   minimum difference between group median fraction or 
#'                          density (default = 0.1 or 10%)
#' @param pop.col           column header of population cell count 
#'                          (default = "Population total cell count")
#' @param subpop.col        column header of subpopulation cell count
#'                          (default = "Cell State total cell count")
#' @param fdr.col           column header of FDR column (default = "adjusted p.value")
#' @param median.cols       vector of column headers of median values
#' 
#' @return tibble of conditions that pass all criteria
bioFilter <- function(dat, calc = "fractions",
                      min.subpop.count = 300, min.pop.count = 1000,
                      min.median.diff = 0.1, max.fdr = 0.05,
                      subpop.col = 'Cell State total cell count',
                      pop.col = 'Population total cell count',
                      fdr.col = 'adjusted p.val', median.cols = NULL){

    if(calc == "fractions"){
        return(dat %>%
               filter(!!as.name(subpop.col) >= min.subpop.count,
                      !!as.name(pop.col) >= min.pop.count,
                      !!as.name(fdr.col) < max.fdr,
                      abs(!!as.name(median.cols[2]) - !!as.name(median.cols[1])) >= min.median.diff))
    } else if(calc == "densities"){
        return(dat %>%
               filter(!!as.name(subpop.col) >= min.subpop.count,
                      !!as.name(fdr.col) < max.fdr,
                      abs((!!as.name(median.cols[2]) - !!as.name(median.cols[1]))/
                           ((!!as.name(median.cols[2]) + !!as.name(median.cols[1]))/2)) >=
                              min.median.diff))
    }

    msg <- paste0("Unrecognized calc: [", calc, "]")
    log_error(msg)
    stop(msg)
}

filterListToTable <- function(filters){
    filtTbl <- tibble()
    for(f in names(filters)){
        filtTbl <- filtTbl %>%
                   bind_rows(filters[[f]] %>% as_tibble() %>% mutate(Calculation = f))
    }
    filtTbl %>% 
    gather(1:3, key='filter', value='value') %>% 
    spread(Calculation, value, fill=0) %>%
    select(filter, everything()) %>%
    mutate(densities = ifelse(filter == "minimum_median_difference",
                              paste0(densities*100, "%"), densities))
}

neighborhoodInQuestion <- function(question){

    nbhd <- NULL

    nbhd <- unique(unlist(question)[grep("\\.Neighborhood$", names(unlist(question)))])
    if(length(nbhd) > 0){
        nbhd <- gsub("_env", "", nbhd)
        log_debug("Neighborhood(s): ", paste(nbhd, collapse = " & "))
    }

    nbhd

}

compareSampleGroups <- function(question, analyses, annotation, resultsFilters, metricsDir, 
                                calcUnit = "FOV_ID", filterID = NULL){

    filtName <- ifelse(is.null(filterID), "filters", filterID)

    qRes <- list()
    grps <- getSampleGroups(question, annotation, calcUnit = calcUnit)
    nbhd <- neighborhoodInQuestion(question)
    statsCalcUnit <- calcUnit
    if(!is.null(nbhd) & length(nbhd) > 0){
        statsCalcUnit <- paste(calcUnit, nbhd, sep = "___")
    }

    tmp <- lapply(seq(grps), function(x){ 
        log_debug(paste0(nrow(grps[[x]]), " ", calcUnit, 
                         " in group [", names(grps)[x], "]"))
    })

    if("fractions" %in% names(analyses)){
        log_debug("Calculating log odds ratios...")
        fltr <- resultsFilters$fractions

        lo <- grps %>%
              getComparisonData(question, metricsDir, filePattern = "fractions", 
                                calcUnit = calcUnit) %>%
              reportLogOdds(question$groupVar, calcUnit = statsCalcUnit)

        lof  <- lo %>%
                filterLOreport(min.subpop.count     = fltr$minimum_subpopulation_count,
                               min.pop.count        = fltr$minimum_population_count,
                               min.median.frac.diff = fltr$minimum_median_difference,
                               max.fdr              = fltr$fdr_cutoff) %>%
                arrange(desc(abs(`log(Odds Ratio)`)), `adjusted p.value`)

        qRes$all_fractions <- lo %>% 
                              mutate(`Fraction Passed` = ifelse(`Cell State ID` %in% lof$`Cell State ID`, "X", NA)) %>%
                              select(`Fraction Passed`, everything())
        qRes[[paste0("fractions_", filtName)]] <- lof
    }

    if("densities" %in% names(allAnalyses)){
        log_debug("Calculating density fold changes...")
        fltr <- resultsFilters$densities
        fc <- grps %>%
              getComparisonData(question, metricsDir, filePattern = "densities") %>%
              right_join(analyses$densities %>% select(`Cell State ID`, Population), by = "Population") %>%
              reportDensityStats(question$groupVar, calcUnit = statsCalcUnit)

        fcf <- fc %>%
               filterFCreport(min.subpop.count = fltr$minimum_subpopulation_count,
                              min.median.diff  = fltr$minimum_median_difference,
                              max.fdr          = fltr$fdr_cutoff) %>%
               arrange(desc(abs(`log(Fold Change)`)), `adjusted p.value`)

        qRes$all_densities <- fc  %>% 
                              mutate(`Density Passed` = ifelse(`Cell State ID` %in% fcf$`Cell State ID`, "X", NA)) %>%
                              select(`Density Passed`, everything())

        qRes[[paste0("densities_", filtName)]] <- fcf
    }

    if(all(c("fractions", "densities") %in% names(allAnalyses))){
        log_debug("Generating report...")
        overlaps <- c("DF", "CI.high", "CI.low", "p.value", "adjusted p.value")
        lo <- qRes$all_fractions %>% rename_at(vars(any_of(overlaps)), list(~ paste("Fraction", .)))
        fc <- qRes$all_densities %>% rename_at(vars(any_of(overlaps)), list(~ paste("Density", .)))
        fin <- full_join(lo, fc) %>%
               select(`Fraction Passed`, `Density Passed`,
                      `Cell State ID`, `Cell State`, Population, `Cell State total cell count`,
                      `Population total cell count`, `Fraction DF`, dplyr::matches("median fraction"),
                      `log(Odds Ratio)`, `Odds Ratio`, `Inverse Odds Ratio`,
                      `Fraction CI.low`, `Fraction CI.high`, `Fraction p.value`, `Fraction adjusted p.value`,
                      `Density DF`, `Total Area`, `Overall Median Density`, dplyr::matches("median density"),
                      `log(Fold Change)`, `Fold Change`, `Inverse Fold Change`,
                      `Density CI.low`, `Density CI.high`, `Density p.value`, `Density adjusted p.value`)
        qRes$all_fractions_and_densities <- fin
        qRes$all_fractions <- NULL
        qRes$all_densities <- NULL
    }

    log_debug("  reporting filters used")
    qRes[[filtName]] <- filterListToTable(resultsFilters)

    log_debug("  summarizing samples included in comparison")
    qRes$samples <- allAnnot %>% filter(`CellDive_ID` %in% unlist(grps))
    if(question$groupVar %in% names(qRes$samples)){
        qRes$sample_summary <- qRes$samples %>% 
                               group_by_at(question$groupVar) %>% 
                               summarize(!!as.name(paste(calcUnit, "(n=)")) := n())
    }

    tabs <- c("all_fractions_and_densities",
              "all_fractions",
              paste0("fractions_", filtName),
              "all_densities",
              paste0("densities_", filtName),
              filtName,
              "samples",
              "sample_summary")
    tabs <- tabs[tabs %in% names(qRes)]
    
    qRes[tabs]

}
