suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))
log_threshold(DEBUG)

library(showtext)
font <- "Roboto"
font_add_google(font, font)
showtext_auto()

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript annotate_cells.R 

          [REQUIRED (may be defined on command line OR in manifest file)]
            --meta_dir            path to meta files in XLSX format
            --study_figure_dir    output directory

          [OPTIONAL]
            --plot_config_file    YAML file containing color assignments for values in
                                  study overview meta file 
            --xVar                column in meta data to be plotted on x-axis
        \n"
    )
}

## set up required args & defaults 
minReq   <- list(c("meta_dir","meta_files"))
defaults <- list(number_threads = 1, study_figure_dir = getwd(), xVar = "Patient_ID")
args     <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

###################
mFiles <- getCurrentMetaFiles(metaDir = args$meta_dir)

dat <- loadStudyAnnotations(metaFiles = mFiles)$flat %>%
       select(dplyr::matches("_ID")) %>%
       mutate(!!as.name(args$xVar) := factor(!!as.name(args$xVar), levels = unique(!!as.name(args$xVar))))


log_debug("Plotting number of FOVs in each sample...")
fovs <- dat %>% select(FOV_ID, args$xVar) %>% unique
fov_plot  <- plot_fov_count_per_sample(fovs,
                                       sampCol = args$xVar) +
             xlab(gsub("_", " ", args$xVar)) +
             theme(plot.margin = margin(t = 0, b = 0.5, r = 0.5, l = 0.5, "in"))

fn <- file.path(args$study_figure_dir, paste0("fov_counts_per_", args$xVar, ".pdf"))
ggsave(fn,                                                                       
       device = cairo_pdf, bg = "white",                                        
       height = 4, width = 8, units = "in")                  

cps <- fovs %>%
       group_by_at(args$xVar) %>%
       summarize(`FOV count` = n())
write.xlsx(cps, gsub(".pdf", ".xlsx", fn))
