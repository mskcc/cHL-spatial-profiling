suppressMessages(library(funr))
suppressMessages(library(ComplexHeatmap))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))

library(ggforce)
library(showtext)

### TEMPORARY
htmpAnnClrs <- c('Positive' = '#6290C8', 'Negative' = '#f0f0f0')

usage <- function(){

    cat("\nUsage:  Rscript plot_drift.R 

          [REQUIRED (may be defined on command line OR in manifest file)]
            --cell_data_dir         path to annotated cell files 
            --meta_dir              path to meta files in XLSX format
            --study_figure_dir      output directory, where figure should be saved

          [OPTIONAL]
            --exclude_samples     comma-delimited string of cell dive IDs of samples
                                  to exclude from heatmap
            --annotation_columns  comma-delimited string of meta data columns with
                                  which to annotate the heatmap
            --manifest            YAML file containing one or more parameter; NOTE: arguments on command
                                  line override manifest arguments!!! default = NULL
            --meta_files          comma-delimited list of meta files; default = NULL
            --plot_config_file    YAML file containing plot config, specifically color
                                  assignments for annotations, if any
            --height              PDF height in inches
            --width               PDF width in inches
        \n"
    )
}

minReq   <- list("cell_data_dir",
                 c("meta_dir","meta_files"),
                 "number_threads",
                 "study_figure_dir")
defaults <- list(number_threads = 4, debug = TRUE, 
                 height = 6.5, width = 8,
                 annotation_columns = NULL,
                 plot_config_file = NULL,
                 exclude_samples = NULL)
allUsed  <- unique(c(unlist(minReq), names(defaults), "manifest"))

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(args, allUsed)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### load meta data
stDat <- loadStudyData(args, annotatedCells = F)

### formatting 
yOrder <- c("Leukocyte_Other", "Plasma",
            stDat$cellTypes %>%
            filter(Category == "Immune_All", 
                  !Subtype %in% c("Plasma", "Leukocyte_Other")) %>%
            pull(Subtype) %>% unique %>% rev)
yLabels <- make_figure_labels(stDat$cellTypes, "Subtype")

htmpAnnCols <- NULL                                                  
#htmpAnnClrs <- NULL
if(!is.null(args$annotation_columns) & !is.na(args$annotation_columns)){                                
#    if(is.null(args$plot_config_file)){
#        log_error("Please provide --plot_config_file when requesting heatmap annotation")
#        stop()
#    } 
    htmpAnnCols <- unlist(strsplit(args$annotation_columns, ","))                
#    htmpAnnClrs <- read_yaml(args$plot_config_file)$clinical[htmpAnnCols]
}                       
                                                                               

cl <- makeCluster(args$number_threads, type="FORK", outfile="")
clusterExport(cl, c("args", "stDat"), envir=environment())

dat <- parLapply(cl, unique(stDat$sampAnn$CellDive_ID), function(cdid){

           if(!is.null(args$exclude_samples) & !is.na(args$exclude_samples)){
               ex_smp <- unlist(strsplit(args$exclude_samples, ","))
               if(cdid %in% ex_smp){ 
                   log_warn(paste0("Excluding sample ", cdid))
                   return(NULL) 
               }               
           }

           af <- file.path(args$cell_data_dir, 
                           paste0(cdid, "_annotated_cells.rda"))

           if(!file.exists(af)){
               log_warn("Annotated cells file not found for sample ", cdid)
               return(NULL)
           }

           loadAnnotatedCells(stDat$sampAnn, annFile = af) %>%
           filter(Category == "Immune_All") %>%
           group_by_at(c("Patient_ID", "FOV_ID", "Subtype", 
                         all_of(htmpAnnCols))) %>%
           summarize(Count = n()) %>%
           ungroup %>%  
           mutate(Total = sum(Count), Fraction = Count/Total) %>%
           select(Patient_ID, FOV_ID, Subtype, all_of(htmpAnnCols), 
                  Count, Total, Fraction)

       }) %>%
       bind_rows %>%
       mutate(Patient_ID = factor(Patient_ID, levels = sort(unique(Patient_ID))),
              FOV_ID = factor(FOV_ID, levels = unique(FOV_ID)),
              Subtype = factor(Subtype, levels = yOrder))

stopCluster(cl)


fn <- file.path(args$study_figure_dir, "fov_variation_heatmap.pdf")
log_info("Saving plot to file ", fn)

#print(htmpAnnCols)
#print(htmpAnnClrs)
#print(yLabels)

grDevices::cairo_pdf(fn, height = args$height, width = args$width)
plot_fov_variation_annotated_heatmap(dat, 
                                     cell_category = "Subtype", 
                                     avg_within = "Patient_ID", 
                                     cell_type_labels = yLabels,
                                     data_col = "Fraction",
                                     annotation_columns = htmpAnnCols, 
                                     annotation_colors = htmpAnnClrs,
                                     data_file = gsub(".pdf", ".xlsx", fn))
dev.off()


