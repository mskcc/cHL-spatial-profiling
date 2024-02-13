microenv_file_pattern <- function(cdid = NULL, neighbor = NULL, state = NULL, 
                                  frac_pos = NULL, min_cells = NULL){
                                                                                
    pat <- ""                                                                   
    if(!is.null(cdid)){                                                         
        pat <- paste0(pat, "^", cdid, "___.*")                                  
    }                                                                           
    pat <- paste0(pat, neighbor, "(,|_)")
    if(!is.null(state)){
        pat <- paste0(pat, state, "___microenvironments__.*")                    
    }
    if(!is.null(frac_pos)){                                                     
        pat <- paste0(pat, "(fPos|frAgg)_", frac_pos, "__.*")                           
    }                                                                           
    if(!is.null(min_cells)){                                                    
        pat <- paste0(pat, "minCells_", min_cells)                              
    }                                                                           
    pat <- paste0(pat, "\\.rda")                                                
    pat                                                                         
                                                                                
}                                                                               
                                        
                                        
annotated_cell_file_pattern <- function(cdid){                                  
    paste0("^", cdid, "_")                                                      
}                                                                               

                                                                                
population_metrics_file_pattern <- function(cdid = NULL, metric = NULL, calc_unit = NULL){
    pat <- ""                                                                   
    if(!is.null(cdid)){                                                         
        pat <- paste0(pat, "^", cdid, "_.*")                                    
    }                                                                           
    if(!is.null(metric)){                                                       
        pat <- paste0(pat, "population_", metric, "_per__.*")                   
    }                                                                           
    if(!is.null(calc_unit)){                                                    
        pat <- paste0(pat, calc_unit)                                           
    }                                                                           
    pat <- paste0(pat, "\\.rda")                                                
    pat                                                                         
}                                       
