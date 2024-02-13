suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript annotate_cells.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --annotated_cells_file  path to RDA file (output) that will contain a 
                                    single table where rows are cells and columns 
                                    are all data for that single cell
            --cell_dive_id          CellDive_ID of sample being analyzed
            --neighbor_cell_class   Category, Cell_type or Subtype of neighbor cells 
                                    of interest 
            --neighbor_file         path to RDA file containing all cell-to-cell distances
            --meta_dir              path to meta data files
            --microenvironment_dir  output directory; RDA file with microenvironment 
                                    assignments will be saved here


          [OPTIONAL]
            --manifest                    YAML file containing one or more parameter; NOTE: 
                                          arguments on command line override manifest 
                                          arguments!!! default = NULL        
            --min_microenv_cells          minimum number of total cells in a cell's 
                                          environment to be assigned positive, negative, or 
                                          mixed environment; default = 1
                                          interest
            --neighborhood_radius         radius in microns defining a cell neighborhood; 
                                          default = 30
            --number_threads              number of threads to use for parallel processes; 
                                          default = 4
            --positive_cutoff             minimum fraction of marker positive tumor cells in 
                                          a 'center' cell's environment required to be considered 
                                          a positive environment; the inverse of this value will 
                                          be considered the maximum value of negative environments; 
                                          all values in between are considered a mixed 
                                          environment; default = 0.95
            --testing                     run code on subset of data for testing purposes
                                  
        \n"
    )
}

minReq   <- list("annotated_cells_file",
                 "meta_dir",
                 "microenvironment_dir",
                 "neighbor_cell_class",
                 "neighbor_file",
                 "number_threads")
defaults <- list(number_threads = 4, debug = TRUE, neighborhood_radius = 30, 
                 neighbor_cell_marker_combo = "aggregated",
                 min_microenv_cells = 1, positive_cutoff = 0.95, testing = FALSE) 
allUsed  <- unique(c(unlist(minReq), names(defaults), "manifest"))

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- read_yaml("input/config/study_config.yaml")
    args$cell_dive_id <- "HL_8"
    args$neighbor_cell_class = "HRS"
    args$neighbor_cell_marker_combo = "aggregated"
    args$pos_cutoff = 0.75
    args$testing <- FALSE
    args$neighbor_file = "preprocessing/neighbors/HL_8___HaloObj_v10.9__211112_Exclusions_ReassignAll_NeighborTbl_50.rda" 
    args$annotated_cells_file = "preprocessing/annotated/HL_8_annotated_cells.rda" 
    args <- processCMD(args, defaults, minReq, usage)
}

setLogThreshold(debug = args$debug)
logParams(args, allUsed)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## prep data
log_info("Loading cell annotation from file: ", args$annotated_cells_file)      
annCells  <- readRDS(args$annotated_cells_file)                                 
hrs_uuids <- annCells %>% filter(Cell_type == "HRS") %>% pull(UUID)             

## get ALL HRS neighbor cells, including non-clustered, which are not in aggregate table
log_info("Loading neighbors from file: ", args$neighbor_file)
neighbors <- loadNeighbors(args$neighbor_file, 
                           maxRadius = args$neighborhood_radius, 
                           debug = args$testing) %>%
             select(C.UUID, 
                    HRS.UUID = N.UUID, 
                    FOV_number = SPOT, 
                    everything()) %>%
             filter(HRS.UUID %in% hrs_uuids) %>% 
             mutate(CellDive_ID = args$cell_dive_id)

## read file with HRS cell cluster assignments & aggregate calls
aggFile <- "hrsCellAggregationTable_ceef61e2_.csv.gz"                           
log_info("Loading HRS cell aggregates from file: ", aggFile)
aggs <- read.csv(aggFile) %>% 
        as_tibble %>% 
        select(CellDive_ID = CellDiveID, 
               FOV_number = SPOT, 
               HRS.UUID = UUID, 
               everything()) %>%
        filter(CellDive_ID == args$cell_dive_id) ## this is not really necessary

## for each immune cell, label all HRS cells in its neighborhood as aggregated or not  
lbls <- c("isolated", "clustered, non-aggregated", "aggregated")
names(lbls) <- lbls

log_debug("assigning aggregation labels to each HRS cell in immune cell neighborhoods")
hrs_nbhd <- neighbors %>%   
            left_join(aggs, by = c("HRS.UUID", "CellDive_ID", "FOV_number")) %>%
            mutate(HRSneighbor = 
                      case_when(is.na(AggregateGe20) ~ lbls['isolated'],
                                !AggregateGe20 ~ lbls['clustered, non-aggregated'], 
                                AggregateGe20 ~ lbls['aggregated'])) 

## summarize fractions of HRS neighbor types in neighborhood of each immune cell
log_debug("summarizing counts and fractions of aggregated HRS cells for each immune cell")
hrs_counts <- hrs_nbhd %>% 
              group_by(CellDive_ID, FOV_number, C.UUID) %>%
              mutate(Total = n()) %>% ## total HRS cells in each cell's neighborhood
              group_by(CellDive_ID, FOV_number, C.UUID, Total, HRSneighbor) %>%         
              summarize(Count = n()) %>%   
              spread(HRSneighbor, Count, fill = 0) 

missing <- setdiff(lbls, names(hrs_counts))
if(length(missing) > 0){
    for(col in missing){
        hrs_counts[[col]] <- 0
    }
}


log_debug("assigning final aggregated/non-aggregated labels to each immune cell") 
## assign aggregated/non-aggregated labels depending on cutoffs
min_frac = args$positive_cutoff
min_hrs_count = args$min_microenv_cells

hrs_microenv <- 
    hrs_counts %>%
    mutate(FractionAgg = !!as.name(lbls['aggregated']) / Total,
           `HRS,aggregated_microEnv` = 
             case_when(Total >= min_hrs_count & 
                         FractionAgg >= min_frac ~ "HRS aggregated neighborhood",
                       Total >= min_hrs_count & 
                         FractionAgg <= (1 - min_frac) ~ "HRS non-aggregated neighborhood",
                       Total >= min_hrs_count ~ "HRS mixed neighborhood",
                       Total < min_hrs_count ~ "no HRS neighbors"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#                       
## save                                                                         
of <- file.path(args$microenvironment_dir,                                      
                gsub("_annotated_cells",                                        
                     paste0("___HRS_aggregated___microenvironments__frAgg_", 
                             min_frac, "__minCells_", min_hrs_count),
                     basename(args$annotated_cells_file)))                      
saveRDS(hrs_microenv, of, compress = T)                                                  
log_info("Done!")



##### QC
## quick test to see if table makes sense                                   
if(interactive()){

agg_coords <- aggs %>% left_join(annCells %>% select(HRS.UUID = UUID, x0, y0)) %>% filter(FOV_number == 15)

    test <- annCells %>% 
            filter(FOV_number == 16) %>% 
            left_join(hrs_microenv %>% select(UUID = C.UUID, `HRS,aggregated_microEnv`))

    clrs = c('magenta', 'lightblue', 'orange')
    names(clrs) = c("HRS aggregated neighborhood", "HRS non-aggregated neighborhood", "HRS mixed neighborhood")
    pointSize = 0.5

    ## quick test to see if table makes sense
    ggplot(test, aes(x = x0, y = y0)) +
    #geom_point(size = 0.25, color = "lightgray") + 
    geom_point(test %>% filter(Cell_type == "HRS"), 
                 mapping = aes(x = x0, y = y0), size = pointSize, color = "red") +
    geom_point(test %>% filter(Cell_type != "HRS", !is.na(`HRS,aggregated_microEnv`)),
                 mapping = aes(x = x0, y = y0, color = `HRS,aggregated_microEnv`), size = pointSize) +
geom_point(data = agg_coords, mapping = aes(x = x0, y = y0, color = AggregateGe20, shape = AggregateGe20), size = 1) +
    #scale_x_continuous(limits = c(1000, 1550)) +
    #scale_y_continuous(limits = c(-1000, -750)) +
    theme_minimal() +
    scale_color_manual(values = clrs, name = "Center Immune Cell") 
}




