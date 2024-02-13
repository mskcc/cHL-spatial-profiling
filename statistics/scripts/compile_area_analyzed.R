suppressMessages(library(funr))                                                 
sdir <- dirname(get_script_path())                                              
source(file.path(sdir, "source_all.R"))                                         
log_info(paste0("Loaded source files from: ",sdir))                             
                                                                                
#####################################                                           
#####        SOURCE CODE        #####                                           
#####################################                                           
                                                                                
usage <- function(){                                                            
                                                                                
    cat("\nUsage:  Rscript calculate_area.R                                     
                                                                                
          [REQUIRED (may be defined on command line OR in manifest file)]       
            --fov_area_dir     output directory for FOV area files (one per sample)
            --meta_dir         path to meta files in XLSX format, required IF meta_data_file is NULL           
            --out_file         XLSX output file 
                                                                                
          [OPTIONAL]                                                            
            --manifest         YAML file containing one or more parameter;      
                               NOTE: arguments on command line override manifest arguments!!!
        \n"                                                                     
    )                                                                           
}                                                                               
                                                                                
### process command line                                                        
minReq   <- list("fov_area_dir",                                                  
                 c("meta_dir","meta_files"),
                 "out_file")                     
defaults <- list()
args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)          
                                                                                
logParams(args, names(args))                                                    
                                                                                
                                                                                
### get list of samples by cell dive ID and patient ID                                              
ids <- loadStudyAnnotations(metaFiles = getCurrentMetaFiles(metaDir = args$meta_dir))$flat %>%
       select(Patient_ID, CellDive_ID) %>%                                      
       distinct                                                                 

fovAreas <- lapply(ids$CellDive_ID, function(cdid){
                af <- file.path(args$fov_area_dir, paste0(cdid, "_total_area_per_FOV_ID.rda"))
                readRDS(af) 
            }) %>%
            bind_rows

cdidAreas <- ids %>%
             left_join(fovAreas %>%
                       group_by(CellDive_ID, Area_unit) %>%
                       summarize(`Total Area Analyzed` = sum(Area)),
                       by = "CellDive_ID") %>%
             select(Patient_ID, CellDive_ID, `Total Area Analyzed`, unit = Area_unit)

tbls <- list(area_by_sample = cdidAreas,
             area_by_fov = fovAreas)
write.xlsx(tbls, args$out_file) 
