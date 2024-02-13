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
            --neighbor_cell_marker_combo  comma-delimited marker combination that together 
                                          with neighbor_cell_class define neighborhood of 
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
defaults <- list(number_threads = 4, debug = TRUE, neighborhood_radius = 30, neighbor_cell_marker_combo = NULL, 
                  min_microenv_cells = 1, positive_cutoff = 0.95, testing = FALSE)
allUsed  <- unique(c(unlist(minReq), names(defaults), "manifest"))

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- read_yaml("input/config/study_config.yaml")
    args$neighbor_cell_class = "HRS"
    #args$neighbor_cell_marker_combo = "B2M,MHCI,PDL1-"
    #args$neighbor_file = "preprocessing/neighbors/HL_27___HaloObj_v10.9__211112_Exclusions__Normalized_NeighborTbl_50.rda"
    #args$annotated_cells_file = "preprocessing/annotated/HL_27_annotated_cells.rda"
    args$testing <- FALSE
    args$neighbor_file = "preprocessing/neighbors/HL_8___HaloObj_v10.9__211112_Exclusions_ReassignAll_NeighborTbl_50.rda"  
    args$annotated_cells_file = "preprocessing/annotated/HL_8_annotated_cells.rda"
    args <- processCMD(args, defaults, minReq, usage)
}

setLogThreshold(debug = args$debug)
logParams(args, allUsed)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## prep data
allMarkers <- getAllMarkers(markerFile = getMetaFile("Markers", metaDir = args$meta_dir))
me <- args$neighbor_cell_class
if(!is.na(args$neighbor_cell_marker_combo)){
    me <- paste0(args$neighbor_cell_class, ",", args$neighbor_cell_marker_combo)
}

log_info("Loading cell annotation from file: ", args$annotated_cells_file)
annCells  <- readRDS(args$annotated_cells_file) %>% 
             select(UUID, Category, Cell_type, Subtype, PositiveMarkers)

log_info("Loading neighbors from file: ", args$neighbor_file)
neighbors <- loadNeighbors(args$neighbor_file, 
                           maxRadius = args$neighborhood_radius, 
                           debug = args$testing) %>%
             mutate(SPLIT = substr(C.UUID, 1, 2)) 

parts      <- unique(neighbors$SPLIT)
threads    <- min(args$number_threads, length(parts))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## make assignments
log_info("Assigning [", me, "] microenvironment statuses") 
cl <- makeCluster(threads, type="FORK", outfile="")
clusterExport(cl, c("neighbors", "annCells", "args"), envir = environment())

env <- parLapply(cl, parts, function(x){
           ndat <- neighbors %>% filter(SPLIT == x)
           getMicroEnv(ndat,
                       annCells,
                       args$neighbor_cell_class,
                       nCellMarkerCombo = args$neighbor_cell_marker_combo,
                       allMarkers = allMarkers, 
                       posCutoff = args$positive_cutoff,
                       minCellsInMicroEnv = args$min_microenv_cells)
       }) %>%
       bind_rows

stopCluster(cl)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## save
of <- file.path(args$microenvironment_dir, 
                gsub("_annotated_cells", 
                     paste0("___", me, "___microenvironments__fPos_", 
                     args$positive_cutoff, "__minCells_", 
                     args$min_microenv_cells), 
                     basename(args$annotated_cells_file)))
saveRDS(env, of, compress = T)


