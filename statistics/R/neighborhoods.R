loadNeighbors <- function(nFile, maxRadius = 30, debug = FALSE){
    if(debug){
        nFileSub = gsub(".rda", "___SubSample_100k.rda", nFile)
        if(!file.exists(nFileSub)){
            nn        <- readRDS(nFile)$neighbors %>% filter(C.UUID != N.UUID, Dij.micron <= maxRadius)
            randCells <- nn %>% distinct(C.UUID) %>% sample_n(10000) %>% pull
            n1        <- filter(nn, C.UUID %in% randCells)
            saveRDS(n1, nFileSub, compress=T)
        }
        return(readRDS(nFileSub))
    } 
    readRDS(nFile)$neighbors %>% filter(C.UUID != N.UUID, Dij.micron <= maxRadius)
}


getMicroEnv <- function(nDat, cellAnn, nCellClass, nCellMarkerCombo = NULL, allMarkers = NULL, 
                         posCutoff = 0.95, minCellsInMicroEnv = 1){

    ## keep list of all cells of nCells class, regardless of marker positivity (e.g., HRS)
    nbrCells <- cellAnn %>% filterDataForCondition(nCellClass, markers = allMarkers) %>% pull(UUID)
    mrkrPos  <- nbrCells

    ## if marker combo given, get list of cells that match marker combo (e.g., B2M,MHCI,PD1-_
    if(!gtools::invalid(nCellMarkerCombo)){
        if(is.null(allMarkers)){
            msg <- "Must provide vector of all study markers."
            log_error(msg)
            stop(msg)
        }
        mrkrPos <- cellAnn %>% 
                   filterDataForCondition(paste0(nCellClass, ",", nCellMarkerCombo), markers = allMarkers) %>% 
                   pull(UUID) 
    }

    ## here, C.UUIDs refer to "center" cells and N.UUIDs are neighbor (e.g., HRS,B2M,MHCI,PD1-) cells
    ## for each center - neighbor cell pair, if neighbor is in list of cells that match nCellClass+nCellMarkerCombo,
    ## assign "Pos"; if neighbor is of class nCellClass but does not match marker combo, assign "Neg"; if
    ## the cell is not of nCellClass, assign "non_[class]_neighbor". 
    ## for each center cell, count how many Pos and Neg neighbors and then the fraction of Pos
    dat <- nDat %>%
           mutate(meCell = ifelse(N.UUID %in% mrkrPos, "Pos", 
                              ifelse(N.UUID %in% nbrCells, "Neg", paste0("non_", nCellClass, "_neighbors")))) %>%
           group_by(C.UUID, meCell) %>%
           summarize(count = n()) %>%
           spread(meCell, count, fill = 0)

    if(!'Pos' %in% names(dat)){
        dat <- dat %>% mutate(Pos = 0)
    }
    if(!'Neg' %in% names(dat)){
        dat <- dat %>% mutate(Neg = 0)
    }

    meCol <- paste0(nCellClass, "_microEnv")
    if(!gtools::invalid(nCellMarkerCombo)){
        meCol <- paste0(nCellClass, ",", nCellMarkerCombo, "_microEnv")
    } 
                    
    dat %>%
    mutate(Total = Pos + Neg, 
           fPos = Pos/Total,
           !!as.name(meCol) :=
                 case_when(
                     Total >= minCellsInMicroEnv & fPos > posCutoff ~ "Pos.Env",
                     Total >= minCellsInMicroEnv & fPos < (1 - posCutoff) ~ "Neg.Env",
                     Total >= minCellsInMicroEnv ~ "Mixed.Env",
                     Total < minCellsInMicroEnv ~ "noEnv"
                 )
           ) %>%
    ungroup

}

get_cell_clusters <- function(n_dat, uuids, radius = 30){

    cells <- n_dat %>% 
             filter(C.UUID %in% uuids, 
                    N.UUID %in% uuids,
                    Dij.micron <= radius)

    aggs <- tibble(UUID = cells$C.UUID[1],
                   SPOT = cells$SPOT[1],
                   AggregateID = 1)

    ids <- unique(c(cells$C.UUID, cells$N.UUID))
    while(length(ids) > 0) {
        agg <- cells %>% 
               filter(C.UUID == ids[1] | N.UUID == ids[1]) %>%
               gather(C.UUID, N.UUID, key = 'tmp', value = 'UUID') %>%
               select(-Dij.micron, -tmp) %>%
               mutate(AggregateID = max(aggs$AggregateID) + 1) %>%
               unique

        add_to <- aggs %>%
                  filter(UUID %in% agg$UUID) %>%
                  pull(AggregateID) %>%
                  unique

        if(length(add_to) > 0){
            aggs <- aggs %>%
                    mutate(AggregateID = ifelse(AggregateID %in% add_to, 
                                                 min(add_to), AggregateID))
            agg <- agg %>% mutate(AggregateID = min(add_to))
        }
        aggs <- bind_rows(aggs, agg) %>%
                unique

        ids <- setdiff(ids, aggs$UUID)
    }

    aggs %>%
    group_by(AggregateID) %>%
    mutate(AggregateSize = n())
}

get_tumor_cell_aggregates <- function(ann_dat, n_dat, radius = 30, 
                                       min_agg_size = 20){

    uuids <- ann_dat %>%
             filter(Cell_type == "HRS") %>%
             pull(UUID)

    if(length(uuids) == 0){ 
        return(tibble())
    }

    aggs <- get_cell_clusters(n_dat, uuids, radius = radius)

    ann_dat %>%
    filter(Category == "HRS") %>%
    select(UUID, SPOT = FOV_number, dplyr::matches("_ID"), x0, y0) %>% 
    left_join(aggs, by = c("UUID", "SPOT")) %>%
    mutate(HRS_spatial_class = ifelse(!is.na(AggregateID) &  
                                  AggregateSize >= min_agg_size, 
                                   'clustered', 'isolated'))
}
