suppressMessages(library(funr))
sdir <- dirname(get_script_path())
source(file.path(sdir, "source_all.R"))
log_info(paste0("Loaded source files from: ",sdir))

library(ggforce)
library(showtext)


usage <- function(){

    cat("\nUsage:  Rscript plot_drift.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --drift_summary_dir     path to drift summary files
            --meta_dir              path to meta files in XLSX format
            --study_figure_dir      output directory, where figure should be saved

          [OPTIONAL]
            --manifest            YAML file containing one or more parameter; NOTE: arguments on command 
                                  line override manifest arguments!!! default = NULL        
            --meta_files          comma-delimited list of meta files; default = NULL
        \n"
    )
}

minReq   <- list("drift_summary_dir", 
                 c("meta_dir","meta_files"),
                 "number_threads",
                 "study_figure_dir")
defaults <- list(number_threads = 4, debug = TRUE)
allUsed  <- unique(c(unlist(minReq), names(defaults), "manifest"))

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

setLogThreshold(debug = args$debug)
logParams(args, allUsed)

font <- "Roboto"
font_add_google(font, font)
showtext_auto()


mFiles <- getCurrentMetaFiles(metaDir = args$meta_dir, metaFiles = args$meta_files)
print(mFiles)

ids <- loadStudyAnnotations(metaFiles = mFiles)$IDs %>%
       select(CellDive_ID, Patient_ID) %>%
       distinct

smryFiles <- dir(args$drift_summary_dir, full.names=TRUE, pattern="_summary.txt.stats")

dat <- lapply(ids$CellDive_ID, function(cdid){
           smry <- smryFiles[basename(smryFiles) == paste0(cdid, "_drift_summary.txt.stats.gz")]
           
           read.csv(smry, header=T, sep="\t") %>%
           as_tibble %>%
           mutate(CellDive_ID = cdid)
       }) %>%
       bind_rows %>%
       left_join(ids)

pTheme <- theme(text = element_text(family = font, size = 16, color = "black"),
                axis.text.x = element_text(size = 14, color = "black"),
                axis.text.y = element_text(size = 14, color = "black"),
                legend.text = element_text(size = 12, margin = unit(c(0, 0, 0, 0.15), "in"), hjust = 0),
                legend.key.height = unit(0.15, "in"),
                legend.title = element_blank(),
                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in"))

bpTheme <- theme_bw() +
           pTheme +
           theme(text = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title.y = element_text(margin = unit(c(0, 0.2, 0, 0), "in")),
            axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
            legend.position = "top",
            legend.justification = "left",
            legend.title = element_blank(),
            legend.text = element_text(margin = unit(c(0, 0.5, 0, 0), "in")),
            panel.grid.major.x = element_blank(),
            panel.border = element_blank(),
            axis.line.x = element_line(color = "black", size = 0.5),
            axis.line.y = element_line(color = "black", size = 0.5))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
qual_col_pals = qual_col_pals[c("Dark2", "Set1", "Set2", "Set3", "Paired", "Accent", "Pastel1", "Pastel2"),]

### cell loss
loss <- plot_accumulated_cell_loss(dat, ptheme = pTheme, col_pals = qual_col_pals)
pct_loss <- plot_accumulated_cell_loss_percent(dat, ptheme = pTheme, col_pals = qual_col_pals)

### cells remaining
zoom <- c(round_any(1 - max(dat$loss_pct), 0.1, floor), round_any(1 - min(dat$loss_pct), 0.1, ceiling))
if(zoom[2] - zoom[1] > 0.5){ zoom <- NULL } ## only zoom if it's worth it
remaining <- plot_percent_cells_remaining(dat, ptheme = pTheme, col_pals = qual_col_pals, zoom_y = zoom) 

fn <- file.path(args$study_figure_dir, "drift_by_cycle.pdf")
pdf(fn, height = 8.5 * 0.75, width = 11)
print(loss)
print(pct_loss)
grid.draw(remaining[[1]]) ## plot
grid.newpage()
grid.draw(remaining[[2]]) ## legend
dev.off()

write.xlsx(dat, gsub(".pdf", ".xlsx", fn), check.names = F)


######### DRIFT BY SAMPLE (NO CYCLE DATA)

dbsCounts <- plot_total_drift(dat, groups = "Patient_ID", bptheme = bpTheme, sort_by = 'cells remaining')
dbsPcts   <- plot_total_drift(dat, groups = "Patient_ID", bptheme = bpTheme, sort_by = 'cells remaining', pct = TRUE)
dbsPctsRmn <- plot_total_drift(dat, groups = "Patient_ID", bptheme = bpTheme, sort_by = 'pct cells remaining', pct = TRUE)

fn <- file.path(args$study_figure_dir, "drift_by_sample.pdf")
pdf(fn, height = 8.5 * .6, width = 11*.75)
print(dbsCounts)
print(dbsPcts)
print(dbsPctsRmn)
dev.off()

write.xlsx(dat, gsub(".pdf", ".xlsx", fn), check.names = F)
