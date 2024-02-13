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
            --height              PDF height in inches
            --width               PDF width in inches
        \n"
    )
}

## set up required args & defaults 
minReq   <- list(c("meta_dir","meta_files"))
defaults <- list(number_threads = 1, study_figure_dir = getwd(), xVar = "Patient_ID", height = 5, width = 11)
args     <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)


###################
mFiles <- getCurrentMetaFiles(metaDir = args$meta_dir)
figCfg <- read.xlsx(getMetaFile("Samples", metaFiles = mFiles), 2, check.names = F) %>% as_tibble %>% 
          select(Column, Section, Order, Sub_order, figure_label)
ordr   <- figCfg %>%
          select(-Column) %>% unique %>%
          arrange(Order, desc(Sub_order))

dat <- loadStudyAnnotations(metaFiles = mFiles)$flat %>%
       select(dplyr::matches("_ID"), figCfg$Column) %>%
       gather(figCfg$Column, key = "Column", value = "Val") %>%
       mutate(Val = ifelse(tolower(Val) == "n/a", NA, Val)) %>%
       full_join(figCfg, by = "Column") %>%
       filter(!is.na(!!as.name(args$xVar))) %>%
       mutate(!!as.name(args$xVar) := factor(!!as.name(args$xVar), levels = unique(!!as.name(args$xVar))),
              Section = factor(Section, levels = unique(ordr$Section)),
              figure_label = factor(figure_label, levels = unique(ordr$figure_label)),
              Val = factor(Val))

allClrs <- NULL
if(!is.null(args$plot_config_file)){
    pcfg <- read_yaml(args$plot_config_file)
    if('clinical' %in% names(pcfg)){
        allClrs <- unlist(pcfg$clinical)
    }
}

log_debug("Making heatmap of clinical data...")
meta_plot <- plot_study_overview(dat, x = args$xVar, clrs = allClrs) +
             theme(#axis.text.x = element_blank(),  
                   plot.margin = margin(t = 0.5, b = 0, r = 0.5, l = 0.5, "in"),
                   legend.key.size = unit(0.1, 'in'),
                   legend.text = element_text(size = 10))

pdf(file.path(args$study_figure_dir, "study_overview.pdf"), height = args$height, width = args$width)
grid.draw(meta_plot)
dev.off()


