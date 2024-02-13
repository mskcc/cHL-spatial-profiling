source("/home/byrne/halo/dev/hodgkins_dev/source_all.R")                        

suppressMessages(library(funr))                                                 
sdir <- dirname(get_script_path())                                              
source(file.path(sdir, "source_all.R"))                                         
log_debug(paste0("loading source files from: ",sdir))                           
                                                                                
#####################################                                           
#####       GET USER INPUT      #####                                           
#####################################                                           
                                                                                
usage <- function(){                                                            
    cat("\nUsage:  Rscript calculate_metrics.R                                  
                                                                                
          [REQUIRED (may be defined on command line OR in manifest file)]      
            --cell_types_file         CellTypes meta file 
            --conditions_file         XLSX file with a column containing the cell state IDs
                                      of the conditions to be included in the plot; the name
                                      of the column must match exactly the value passed to 
                                      '--question'.  
            --question                QuestionNumber of question to plot
            --statistics_tables_dir   path to XLSX stats tables
            --statistics_figures_dir  output directory  

          [OPTIONAL]
            --width       width of output file in inches; default = 7
            --height      height of output file in inches; default = 7.5 
            --rel_wdiths  comma-delimited string of two numbers representing
                          the relative widths of the plot and the labels; 
                          default = '1,1'
        \n"                                                                     

    )                                                                           
}                                                                               
                                                                                
## names of equired args   
minReq <- c("cell_types_file", "conditions_file", "question", 
            "statistics_tables_dir", "statistics_figures_dir")

defaults <- list(height = 7.5, width = 7, rel_widths = "1,1") # store rel_widths as it would be on the command line

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- list()
    args$cell_types_file  <- "input/meta/HodgkinLymphoma_CellTypes__V0.xlsx"
    args$question <- "12_UT_EBV-_agg"
    args$conditions_file <- "Figure_5_states.xlsx"
    args$statistics_tables_dir <- "results/statistics/spatial/tables"
    args$statistics_figures_dir <- "results/statistics/spatial/figures/odds_ratios"

    args <- processCMD(args, defaults, minReq, usage)          
}

mkdir(args$statistics_figures_dir)

quest <- args$question                                              
conds <- rev(read.xlsx(args$conditions_file, 1)[,quest])                         
conds <- conds[!is.na(conds)]
condOrder <- setNames(1:length(conds), conds)                                   

cellTypes <- getCellTypes(args$cell_types_file)
                                                                                
formatStatsResultsForOddsRatioPlots <- function(fl){                              
                                                                                 
    dt <- read.xlsx(fl, 1, check.names = F) %>%                                       
          as_tibble 

    grps <- names(dt)[grep(" median fraction", names(dt))]
    grps <- gsub(" median fraction", "", grps)
    comp <- paste0(grps, collapse = " vs ")

    dt <- dt %>%
          select(`Cell State ID`, `Cell State`, Population, 
                 dplyr::matches("Odds Ratio"), dplyr::matches("adj"),
                  `CI.low`, `CI.high`) %>%
          filter(`Cell State ID` %in% as.character(conds)) %>%                    
          mutate(Comparison = comp, 
                 y = condOrder[`Cell State ID`],                                  
                 signif = factor(case_when(`adjusted p.value` < 0.001 ~ "***",
                                           `adjusted p.value` < 0.01 ~ "**",       
                                           `adjusted p.value` < 0.05 ~ "*"), 
                                 levels = c("***", "*")))
    dt
}   

resFile <- file.path(args$statistics_tables_dir, paste0(quest, ".xlsx"))

dat <- formatStatsResultsForOddsRatioPlots(resFile)      
comps <- unique(dat$Comparison)                                                 


## ODDS RATIO PLOT 
or_plot <- dat %>%
           mutate(labelY = y) %>%
           plotEffectSize(effectCol = "Odds Ratio", yVar = "y",
                          stripeBG = FALSE, separateFacets = FALSE, 
                          facetY = NULL, build_plot = FALSE, fontsize = 12) +
           theme(plot.margin = unit(c(1, 0, 0.5, 1), "in"))

or_plot <- ggplot_gtable(ggplot_build(or_plot))

### LABELS
lblDat <- dat %>%                                                              
          ungroup() %>%                                                         
          mutate(y = y - 0.5) %>%                                               
          select(`Cell State ID`, `Cell State`, Population, y) %>%              
          unite("Cell State Label", `Cell State`, Population, sep=" | ") %>%          
          unique()                                                              
                                                                                
lbls <- ggplot(lblDat, aes(x = 0, y = y)) +                                     
        scale_x_continuous(limits = c(-0.05, 2), expand = c(0,0)) +             
        scale_y_continuous(limits = c(0.5, nrow(lblDat) + 0.5), expand = c(0,0.1)) #, breaks = seq(length(ids))) 

for(rw in 1:nrow(lblDat)){                                                      
    row <- lblDat[rw,]                                                          
    numSp <- (4 - nchar(row$`Cell State ID`)) * 2                               
    leadSp <- paste(rep(" ", numSp), collapse = "")                             
    id <- paste0(leadSp, row$`Cell State ID`, " ")                              
    label <- conditionLabel(row$`Cell State Label`[1], cellTypes, returnExpression = FALSE)
    lbls <- lbls +                                                              
            geom_text(dat = row,                                                
                      aes(x = 0, y = y, hjust = 0),                             
                      vjust = 0.0,                                              
                      label = as.expression(bquote(.(id) ~ .(label))),          
                      size = 4.25)                                              
}                                                                               
                                                                                
lbls <- lbls +
        labs(caption = paste0(quest, "\n", comps[1])) +                                                     theme_minimal() +   
        theme(panel.grid.major = element_blank(),     
              panel.grid.minor = element_blank(),                               
              axis.text = element_blank(),                                      
              axis.title = element_blank(),                                     
              plot.margin = unit(c(1, 0.5, 0.5, 0), "in"),
              plot.caption = element_text(size =10 , hjust = 1, vjust = -3))  
 
gglbls <- ggplot_gtable(ggplot_build(lbls))
 
outPDF     <- file.path(args$statistics_figures_dir, paste0(quest, "__OR.pdf"))
outXLSX    <- gsub(".pdf", ".xlsx", outPDF)                                          
rel_widths <- as.numeric(unlist(strsplit(args$rel_widths, ",")))                                                        
print(rel_widths)

                                                                                
plotList <- matchPanelHeights(or_plot, gglbls)                                
pg <- plot_grid(plotList[[1]], plotList[[2]], align = 'h', axis = 'bt', rel_widths = rel_widths, rel_heights = c(1,1), nrow = 1)

write.xlsx(dat, outXLSX)
ggsave(pg, height = args$height, width = args$width, units = "in",                        
       device = cairo_pdf, filename = outPDF)                                   
               
    
