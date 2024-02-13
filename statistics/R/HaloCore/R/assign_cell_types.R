#' Identify cell types in marker data talbe
#'
#' For each cell in marker data table, identify its cell type based on
#' the positive/negative requirements in cell types definition table
#'
#' @param m_dat      tibble; marker data table that is part of the halo object 
#'                   (obj$marker.data) joined with Sample and FOV by UUID
#' @param cell_types tibble; cell type definitions read directly from cell 
#'                   types meta file
#' @param markers    vector of all identity markers used in study; IMPORTANT: 
#'                   do NOT included non-identity markers here
#' @param flexible   logical; when TRUE, definitions will be considered
#'                   'flexible' in that if a marker is a required negative
#'                   but is missing in the data, ignore that negative
#'                   requirement; note: if a required positive marker is
#'                   missing, requirement can NOT be ignored and cell type
#'                   will not be assigned to any cells
#'
#' @return m_dat tibble with assigned cell type identities. Keyed with UUID
#' @export
assign_cell_types_by_fov <- function(m_dat, cell_types, markers, 
                                      flexible = FALSE){
    label_cols <- setdiff(names(cell_types), markers)

    if(!all(c("Sample","FOV") %in% colnames(m_dat))){
        log_fatal("Sample and/or FOV not in m_dat")
        rlang::abort("FATAL ERROR: assign_cell_types_by_fov - missing columns")
    }

    dat <- m_dat %>% 
           filter(!is.na(Positive_Classification)) %>%
           select(UUID, Sample, FOV, Marker, Positive_Classification) %>%
           spread(Marker, Positive_Classification)

    fovs <- dat %>% 
            select(Sample, FOV) %>%
            unique

    cells <- lapply(1:nrow(fovs), function(fov_row){
               smp <- fovs[fov_row,]$Sample
               fov <- fovs[fov_row,]$FOV 

               log_info("Sample ", smp, ", FOV ", fov)

               fov_dat <- dat %>%
                           filter(Sample == smp, FOV == fov) 
    
               lapply(1:nrow(cell_types), function(ct_row){
                   cell_labels <- cell_types[ct_row, label_cols]

                   log_info(paste(cell_labels %>% 
                                   select(Category, 
                                          Cell_type, 
                                          Subtype), 
                                   collapse = "/"))

                   cell_def <- cell_types[ct_row, markers] %>% 
                               select_if(!is.na(.))        
                   ct_dat <- fov_dat %>%
                             filter_by_cell_type_definition(cell_def, 
                                                            flexible = flexible)

                   if(!is.null(ct_dat) & nrow(ct_dat) > 0){
                      ct_dat <- ct_dat %>%
                                select(UUID, Sample, FOV) %>%
                                bind_cols(cell_labels) %>%
                                unique  ## in case a cell statisfies more
                                        ## than one of the expanded cell type 
                                        ## definitions
                      log_debug("Found [ ", nrow(ct_dat), " ] ",
                                 paste0(unlist(cell_labels %>% 
                                               select(Category,
                                                      Cell_type,
                                                      Subtype)), 
                                        collapse = "/"), 
                                  " cells defined by [ ROW ", 
                                  ct_row, " ] of cell type file.")
                      return(ct_dat)
                   }
                   return(NULL)
               }) %>%
               bind_rows
             }) %>%
             bind_rows

    ## cells with no positive ID markers are super neg; ignore missing markers
    mrkrs <- intersect(markers, names(dat))
    sn <- dat %>%
          filter_at(vars(all_of(mrkrs)), all_vars(. == 0)) 
    if(nrow(sn) > 0){
        sn <- sn %>%
              select(UUID, Sample, FOV) %>%                                        
              unique %>%                                                           
              mutate(Category = "superNeg",                                         
                     Cell_type = "superNeg",                                        
                     Subtype = "superNeg")   
    }
    unk <- dat %>%
           select(UUID, Sample, FOV) %>%
           unique %>%
           mutate(Category = "UNKNOWN", 
                  Cell_type = "UNKNOWN", 
                  Subtype = "UNKNOWN")
    if(nrow(cells) > 0){
        unk <- unk %>%
               filter(!UUID %in% cells$UUID, !UUID %in% sn$UUID)
    }
    log_debug(paste0("Assigning [ ", nrow(unk),
                     " ] unidentified cells to ", 
                     "Category, Cell_type & Subtype 'UNKNOWN'"))

    cells <- bind_rows(cells, sn, unk)

    log_debug("Checking for identity conflicts...")
    errsFound <- cell_identity_conflicts_found(cells)
    if(errsFound){
        log_error("Found conflicting cell types. Can not identify cells.")
        stop()
    }
    log_debug("All good!")

    return(cells)
}


#' Identify cell types in marker data talbe
#'
#' For each cell in marker data table, identify its cell type based on
#' the positive/negative requirements in cell types definition table
#'
#' @param m_dat       tibble; marker data table that is part of the halo object 
#'                   (obj$marker.data)
#' @param cell_types  tibble; cell type definitions read directly from cell 
#'                   types meta file
#' @param markers    vector of all identity markers used in study; IMPORTANT: 
#'                   do NOT included non-identity markers here
#'
#' @return m_dat tibble with assigned cell type identities. Keyed with UUID
#' @export
assign_cell_types <- function(m_dat, cell_types, markers){ 
    .Deprecated("assign_cell_types_by_fov", 
                msg = paste("'assign_cell_types' is deprecated and",
                            "may not work with the latest cell type",
                            "definitions. Use",
                            "'assign_cell_types_by_fov' instead."))
    
    label_cols <- setdiff(names(cell_types), markers)

    dat <- m_dat %>% 
           select(UUID, Marker, Positive_Classification) %>%
           spread(Marker, Positive_Classification)

    c_types <- cell_types
    na_cell_types <- get_na_cell_types(m_dat, cell_types, markers)
    if(!is.null(na_cell_types) && nrow(na_cell_types) > 0){
        c_types <- c_types %>% 
                   select_at(label_cols) %>%
                   bind_rows(na_cell_types) %>%
                   group_by_at(label_cols) %>%
                   summarize(Count = n()) %>%
                   filter(Count == 1) %>%
                   left_join(cell_types, 
                             by = intersect(names(.), names(cell_types))) %>%
                   ungroup
        if(nrow(c_types) == 0){
            log_error("No cell types left after removing NA cell types!")
            return(NULL)
        }
    }

    ids <- lapply(1:nrow(c_types), function(i){
                  cell_labels <- c_types[i, label_cols]
                  cell_def <- c_types[i, markers] %>% select_if(!is.na(.))
                  if(!all(names(cell_def) %in% names(dat))){
                      missing <- setdiff(names(cell_def), names(dat))
                      log_warn(paste0("[", 
                               paste(c(cell_labels$Cell_type, 
                                       cell_labels$Subtype),
                                     collapse = "/"),
                               "] missing markers: ", 
                               paste(setdiff(names(cell_def), names(dat)), 
                                     collapse = ",")))
                      cell_def <- remove_neg_requirements(cell_def, missing)
                  }
                  ct_dat <- dat %>%
                            filter_by_cell_type_definition(cell_def) %>%
                            select(UUID) %>%
                            bind_cols(cell_labels) %>%
                            unique  ## important in case a cell statisfies more
                                    ## than one of the expanded cell type 
                                    ## definitions

                  log_debug("  Found [ ", nrow(ct_dat), " ] ",
                             paste0(unlist(cell_labels %>% 
                                           select(Category,Cell_type,Subtype)), 
                                    collapse = "/"), " cells")

                  ct_dat

            }) %>%
            bind_rows

    ## label all other cells unknown
    unk_cells <- dat %>% filter(!UUID %in% ids$UUID) %>% pull(UUID) %>% unique
    unk <- tibble(UUID = unk_cells) %>%
           mutate(Category = "UNKNOWN", Cell_type = "UNKNOWN")

    log_debug(paste0("Assigning [ ", nrow(unk), 
                     " ] unidentified cells to Category & Cell_type 'UNKNOWN'"))
    ids <- bind_rows(ids, unk)

    log_debug("Checking for identity conflicts...")
    errsFound <- cell_identity_conflicts_found(ids)
    if(errsFound){
        log_error("Cell types invalid. Returning original data without cell ", 
                  "type annotation.")
        return(m_dat)
    }
    log_debug("All good!")

    return(ids)
}

###############################################################################
# Internal code
###############################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(logger)
    library(openxlsx)
    library(purrr)
})

remove_neg_requirements <- function(cell_def, missing_markers){

    if(length(missing_markers) == 0){ return(cell_def) }

    log_debug("  POS:[", 
              paste0(names(cell_def)[cell_def == 1], collapse = ","), "]; NEG:[", 
              paste0(names(cell_def)[cell_def == 0], collapse = ","), "]")

    ok_to_rmv <- cell_def %>%
                 select_at(missing_markers) %>%
                 unique %>%
                 select_if(as.numeric(.) == 0) %>%
                 names

    if(!all(missing_markers %in% ok_to_rmv)){
        mm <- setdiff(missing_markers, ok_to_rmv)
        log_warn("Missing required positive marker(s) [",  
                 paste0(mm, collapse = ","), 
                 "]. Can not remove marker requirment.")
        return(NULL)
    }

    log_warn("Removing negative requirement of ", 
              paste0(ok_to_rmv, collapse = ","))

    cell_def %>% select(-any_of(ok_to_rmv))
}


## m_dat here is the marker data that has already been SPREAD by marker
get_na_cell_types_by_fov <- function(m_dat, cell_types, markers){
    label_cols <- setdiff(names(cell_types), markers)

    fovs <- m_dat %>%
             select(Sample, FOV) %>%
             unique

    ### NOTE: right now the warnings printed may be redundant
    ### due to expanded definitions in the case where a definition
    ### in the meta data contains a 'one or more of these markers'
    ### portion 
    lapply(1:nrow(fovs), function(x){
        sample_fov <- fovs[x,]
        dat <- m_dat %>%
               filter(Sample == sample_fov$Sample,
                      FOV == sample_fov$FOV) %>%
               select_if(~!all(is.na(.)))

        lapply(1:nrow(cell_types), function(i){
            cell_labels <- cell_types[i, label_cols]
            cell_def <- cell_types[i, markers] %>%
                        select_if(!is.na(.) & . != 0)

            if(!all(names(cell_def) %in% names(dat))){
                log_warn(paste0("[",
                         paste(c(cell_labels$Cell_type, cell_labels$Subtype),
                                 collapse = "/"),
                         "] missing marker(s) [",
                         paste(setdiff(names(cell_def), names(dat)),
                               collapse = ",")), 
                         "] in sample [",
                         sample_fov$Sample, "], FOV [",
                         sample_fov$FOV, "]."
                               )
                return(bind_cols(sample_fov, cell_types[i,]))
            }
            return(tibble())
        })
    }) %>%
    bind_rows()
}



get_na_cell_types <- function(m_dat, cell_types, markers){
    ### should we deprecate this function as well or leave it 
    ### and assume that m_dat contains data for only one Sample,FOV?

    label_cols <- setdiff(names(cell_types), markers)

    dat <- m_dat %>%
           select(UUID, Sample, FOV, Marker, Positive_Classification) %>%
           spread(Marker, Positive_Classification)

    lapply(1:nrow(cell_types), function(i){
        cell_labels <- cell_types[i, label_cols]
        cell_def <- cell_types[i, markers] %>% 
                    select_if(!is.na(.) & . != 0)
        if(!all(names(cell_def) %in% names(dat))){
            log_warn(paste0("[",
                     paste(c(cell_labels$Cell_type, cell_labels$Subtype),
                             collapse = "/"),
                     "] missing markers: ",
                     paste(setdiff(names(cell_def), names(dat)),
                           collapse = ",")))
            return(cell_labels)
        }
        return(tibble())
    }) %>%
    bind_rows()
}

#' Expand cell type definition if it has a set of markers of which
#' one or more must be positive 
#'
#' For cell type definitions that include a set of conditionally positive 
#' markers, create the multiple definitions that actually work. 
#' 
#' @param cell_def  tibble; a one-row tibble containing only the marker columns 
#'                 one cell type definition. a cell type is defined by positive 
#'                 markers indicated with a 1, negative markers indicated with 
#'                 a 0. To indicate groups of markers, of which one or more must 
#'                 be positive assign a group identifier, e.g., 'A' and append 
#'                 the minimum number of positive markers within the group that 
#'                 will satisfy the definition. For example, if Marker1, Marker2, 
#'                 and Marker3 are all marked with 'A_2', any cell with at least 
#'                 two of those three markers can be assigned to that cell type. 
#'                 Multiple groups can be defined within the same definition by 
#'                 changing the group identifier (e.g., 'B_1'). If, in the 
#'                 definition requiring at least two markers in Group A to be 
#'                 positive other markers marked with 'B_1', cells must ALSO 
#'                 have at least one positive marker in group B.
#' @return tibble with all possible marker combinations defining a cell type
expand_cell_definition <- function(cell_def){

    to_exp <- cell_def %>% select_if(grepl("_", .))
    if(ncol(to_exp) == 0){ return(cell_def) }

    exp_grps <- unique(unlist(to_exp))

    lapply(exp_grps, function(i){
        grp_exp <- names(to_exp)[to_exp == i]
        min_pos <- as.integer(gsub(".*_", "", i))
        combos  <- combn(grp_exp, min_pos, simplify = FALSE)
        lapply(combos, function(combo){
            x <- as.list(rep(1, min_pos))
            names(x) <- combo
           
            cell_def %>% 
            select(-names(to_exp)) %>%
            bind_cols(as_tibble(x)) %>%
            mutate(TMP = 'X') # in case there are no markers separate from grouped
                              # markers, need a tmp column to join multiple expanded 
                              # groups
    
        }) %>%
        bind_rows
    }) %>%
    reduce(full_join) %>% ## finally, get all combos of group combos
    select(-TMP)
}


#' Filter cell marker data for cells fitting a single cell type definition
#' 
#' Given a table of cell marker data and a single cell type definition 
#' consisting of one row from the cell types meta table, extract all cells 
#' matching the definition.
#'
#' @param dat      tibble; marker data table from halo object
#' @param cell_def  tibble; a one-row tibble containing only the marker columns  
#'                 one cell type definition. a cell type is defined by positive 
#'                 markers indicated with a 1, negative markers indicated with  
#'                 a 0, and/or  a set of markers, of which one or more must be  
#'                 positive, each  indicated with a '+'        
#' @param flexible logical; when TRUE, negative requirement of a missing marker
#'                 will be ignored; default = FALSE 
#' 
#' @return filtered data table with only cells satisfying given cell type 
#'         definition
filter_by_cell_type_definition <- function(dat, cell_def, flexible = FALSE){

    defs <- cell_def %>% expand_cell_definition 
    cell_dat <- dat %>% 
                select_if(~!all(is.na(.)))

    lapply(1:nrow(defs), function(x){
        def <- defs[x,] %>% 
               select_if(!is.na(.)) %>% 
               mutate_if(is.character, as.numeric)
        missing <- setdiff(names(def), names(cell_dat))

        if(length(missing) > 0 && !flexible){
            log_warn(paste0("[", paste(missing, collapse = ","), 
                     "] missing from data. Can not assign cell type."))
            return(tibble())
        }

        def <- remove_neg_requirements(def, missing)
        if(!is.null(def)){
            def %>%
            left_join(cell_dat, by = names(def)) %>% 
            filter(!is.na(UUID))  ## UUID is NA if there are zero cells fitting def
        }
    }) %>% 
    bind_rows 
}


#' Validate cell type assignments
#'
#' Check for cells assigned to multiple cell types. If any conflicts are found, 
#' log them as errors and return FALSE. Otherwise, return TRUE.
#'
#' @param dat  tibble; data table with UUIDs and Category, Cell_type, Subtype 
#'             assignments
#'
#' @return logical
cell_identity_conflicts_found <- function(dat){
    dat <- dat %>% 
           select(UUID, Category, Cell_type, Subtype) %>% 
           unique 

    dups <- dat %>%
            group_by(UUID) %>%
            mutate(Count = n()) %>%
            filter(Count > 1) 

    if(nrow(dups) > 0){
        conf <- dups %>%
                mutate(Conf = paste(Subtype, collapse = ", ")) %>%
                pull(Conf) %>%
                unique
        log_error(paste0(nrow(dups), " cells were assigned multiple cell ",
                         "types. Conflicting definitions between:"))
        log_error(paste(names(conf), collapse="\t"))
        lapply(1:length(conf), function(x){ 
            log_error(conf[x])
        })
        return(TRUE)
    }
    FALSE
}

get_all_expanded_cell_types <- function(cell_types_xlsx, markers_xlsx){

    markers <- read.xlsx(markers_xlsx, 1, check.names = F) %>% 
               as_tibble %>% 
               filter(Identity == "Y") %>%
               pull(Marker_name)
    cell_types <- read.xlsx(cell_types_xlsx, 1, check.names = F) %>% 
                  as_tibble %>% 
                  mutate_all(as.character)
    classes <- setdiff(names(cell_types), markers)
    if(!"superNeg" %in% cell_types$Cell_type){
        super_neg <- c("superNeg", "superNeg", rep("0", length(markers))) %>% 
                    as.list
        names(super_neg) <- c("Category", "Cell_type", markers)
        super_neg <- super_neg %>% 
                    as_tibble 
        cell_types <- cell_types %>% 
                      bind_rows(super_neg)
    }

    lapply(1:nrow(cell_types), function(i){
        cell_types[i, classes] %>% 
        bind_cols(cell_types[i, markers] %>% 
                  select_if(!is.na(.)) %>%
                  expand_cell_definition %>%
                  mutate_if(is.character, as.numeric)) 
    }) %>%
    bind_rows

}

                   
