
validate_cell_type_definitions <- function(cell_types_file, markers_file){
    exp_cts <- get_all_expanded_cell_types(cell_types_file, markers_file)
    markers <- exp_cts %>% select_if(is.numeric) %>% names
    pos_combos <- lapply(markers, function(x){ 
                      tmp <- list(0:1); names(tmp) <- x; tmp 
                  }) %>% 
                  unlist(recursive = F) %>%
                  expand.grid %>%
                  as_tibble 

    m_dat <- pos_combos %>% 
             mutate(Sample = "TEST", FOV = 0, UUID = row_number()) %>%
             gather(markers, key = "Marker", value = "Positive_Classification")

    tryCatch({
        test <- assign_cell_types_by_fov(m_dat, exp_cts, markers) 
    }, error = function(e){
        log_error("Cell type definition(s) invalid.")
    })
}
