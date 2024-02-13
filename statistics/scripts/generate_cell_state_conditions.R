suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_debug(paste0("loading source files from: ",sdir))

log_threshold(DEBUG)

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript generate_cell_state_conditions.R

          [REQUIRED (may be defined on command line OR in manifest file)]
            --meta_dir                    directory containing all meta files, including cell types
                                          and cell states

          [OPTIONAL]
            --manifest                    YAML file containing custom config; default = NULL
            --statistics_conditions_file  XLSX output file of all conditions to be analyzed, including
                                          fractions and densities; default = 'input/config/cell_state_conditions.xlsx'
            --statistics_conditions_index XLSX output file to which indexed conditions should be saved; 
                                          default = 'results/statistics/cell_state_conditions_index.xlsx'
            --annotation_config_file      YAML file to which annotation config should be saved; 
                                          default = 'input/config/annotation_config.yaml'
            --universal_markers           comma-delimited string (no spaces) of markers to be combined 
                                          with every other combination created; default = 'KI67'
            --index_conditions            logical; when TRUE, all cell states will be converted to 
                                          cell fraction and density conditions and indexed according
                                          to annotation config; default = FALSE
        \n"
    )
}

minReq <- list("meta_dir")

defaults <- list(statistics_conditions_file = "input/config/cell_state_conditions.xlsx", 
                 statistics_conditions_index = "results/statistics/cell_state_conditions_index.xlsx",
                 annotation_config_file = "input/config/annotation_config.yaml",
                 universal_markers = "KI67",
                 index_conditions = FALSE)

args <- processCMD(commandArgs(asValue = TRUE), defaults, minReq, usage)

ctFile <- getMetaFile("CellTypes", metaDir = args$meta_dir)
types <- getCellTypes(ctFile) %>%
         select(Category, Cell_type, Subtype, Cell_type_figures, Cell_type_figures_subscript, 
                Subtype_figures, Subtype_figures_subscript) %>%
         unique() 

csFile <- getMetaFile("CellState", metaDir = args$meta_dir)
states <- read.xlsx(csFile, 1, check.names = F) %>% 
          as_tibble() %>%
          select(!dplyr::matches("NA")) %>%
          filter(!state %in% args$universal_markers) %>%
          expandFunctionalStates() 

states <- states %>%
          gather(2:ncol(.), key = "Subtype", value = "Apply") %>%
          filter(Apply == "X") %>%
          select(Subtype, state, Apply) %>%
          left_join(types, by = "Subtype") %>%
          mutate(DAPI = "DAPI") %>% 
          select(DAPI, Category, Cell_type, Subtype, state, 
                 Cell_type_figures, Cell_type_figures_subscript, Subtype_figures, Subtype_figures_subscript) %>%
          addUniversalMarkers(markers = args$universal_markers) %>%
          mutate(Cell_type_figures = ifelse(is.na(Cell_type_figures), Cell_type, Cell_type_figures),
                 Subtype_figures = ifelse(is.na(Subtype_figures), Subtype, Subtype_figures))

treeOrder <- list("DAPI", "Category", 
                  c("Cell_type"), #"Cell_type_figures", "Cell_type_figures_subscript"), 
                  c("Subtype"), #"Subtype_figures", "Subtype_figures_subscript"), 
                  "state")

fracTypes  <- list(list(nm = "Category",  dnm = "DAPI"),
                   list(nm = "Cell_type", dnm = "DAPI"),
                   list(nm = "Subtype",   dnm = "DAPI"),
                   list(nm = "Cell_type", dnm = "Category"),
                   list(nm = "Subtype",   dnm = "DAPI"),
                   list(nm = "Subtype",       dnm = "Category"),
                   list(nm = "Subtype",       dnm = "Cell_type"),
                   list(nm = c("Category","state"),  dnm = "Category"),
                   list(nm = c("Cell_type","state"), dnm = "Cell_type"),
                   list(nm = c("Subtype","state"),   dnm = "Subtype"))

for(mrkr in args$universal_markers){
    mrkrTypes <- list(list(nm = c("Category", mrkr),         dnm = "Category"),
                      list(nm = c("Cell_type", mrkr),        dnm = "Cell_type"),
                      list(nm = c("Subtype",mrkr),           dnm = "Subtype"),
                      list(nm = c("Category","state",mrkr),  dnm = c("Category","state")),
                      list(nm = c("Cell_type","state",mrkr), dnm = c("Cell_type","state")),
                      list(nm = c("Cell_type","state",mrkr), dnm = c("Cell_type")),
                      list(nm = c("Subtype","state",mrkr),   dnm = c("Subtype","state")),
                      list(nm = c("Subtype","state",mrkr),   dnm = c("Subtype")))
    fracTypes <- c(fracTypes, mrkrTypes)
}

fractions <- tibble()

for(ft in fracTypes){
    lvl <- lapply(treeOrder[1:4], function(x) any(x %in% ft$nm)) %>% unlist() %>% which() %>% max()
    parentCols <- unique(c(unlist(treeOrder[1:(lvl-1)])))
    cols <- unique(c(unlist(treeOrder[1:lvl]), unlist(ft, use.names = F)))

    conds <- states %>% 
             select(all_of(cols)) %>%
             unite("Subpopulation", ft$nm, sep=",", remove = F, na.rm = T) %>%
             unite("Population", ft$dnm, sep=",", remove = F, na.rm = T) %>%
             unique() %>%
             mutate(FracType = paste0(paste0(ft$nm, collapse = ","), "/", paste0(ft$dnm, collapse = ",")), 
                    level = lvl) %>%
             filter(Subpopulation != Population) 

    if("Cell_type" %in% cols){
        conds <- conds %>%
                 left_join(types %>% select(dplyr::matches("Cell_type")) %>% unique)
    }
    if("Subtype" %in% cols){
        conds <- conds %>%
                 left_join(types %>% select(dplyr::matches("Subtype")))
    }

    if(nrow(conds %>% filter(is.na(Category))) != 0){
        log_error("INVALID CONDITIONS CREATED!")
    }

    fractions <- fractions %>% bind_rows(conds) %>% unique
    
}

## remove duplicates
fractions <- fractions %>%
             group_by(Subpopulation, Population) %>% 
             mutate(Count = n(), maxLevel = max(level)) %>% 
             filter(level == maxLevel) %>%
             select(FracType, Subpopulation, Population, Category, Cell_type, Subtype,
                    Function = state, Proliferation = KI67,
                    Cell_type_figures, Cell_type_figures_subscript,
                    Subtype_figures, Subtype_figures_subscript) 

densities <- fractions %>% 
             select(Subpopulation, Category, Cell_type, Subtype, Function, Proliferation, 
                    Cell_type_figures, Cell_type_figures_subscript,
                    Subtype_figures, Subtype_figures_subscript) %>%
             unique() %>%
             bind_rows(tibble(Category = "DAPI", Subpopulation = "DAPI"))

## save file
log_debug(paste0("Saving conditions to file: ", args$statistics_conditions_file))
analyses <- list(fractions = fractions, densities = densities)
write.xlsx(analyses, file = args$statistics_conditions_file, check.names = F)

## index conditions
if(args$index_conditions){
    if(fileDone(args$statistics_conditions_index)){
        unlink(args$statistics_conditions_index)
    }
    annotConfig <- configureConditions(types, states, 
                                       configYAML = args$annotation_config_file, 
                                       force = TRUE)
    condIdx <- getConditionsIndex(args$statistics_conditions_index, 
                                  args$statistics_conditions_file, 
                                  annotConfig$arrange_annotation) 
}
