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
            --cell_data_dir       path to anntoated cells RDA files             
            --meta_dir            path to meta files in XLSX format             
            --study_figure_dir    output directory                              
            --plot_config_file    YAML file containing color assignments for values in
                                  study overview meta file 
                                                                                
          [OPTIONAL]                                                            
            --x_var                column in meta data to be plotted on x-axis  
        \n"                                                                     
    )                                                                           
}                                                                               
                                                                                
## set up required args & defaults                                              
minReq   <- list(c("meta_dir","meta_files"))                                    
defaults <- list(number_threads = 1, study_figure_dir = getwd(), x_var = "Patient_ID")
used     <- c(unlist(minReq), names(defaults), "cell_data_dir", "meta_dir", "study_figure_dir", "plot_config_file")
args     <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)      
                                                                                
logParams(args, used)                                                           
                                                                                
stDat <- loadStudyData(args, annotatedCells = FALSE)                            
threads   <- min(args$number_threads, detectCores() - 2)                        
                                                                        

library(showtext)
font <- "Roboto"
font_add_google(font, font)
showtext_auto()


#####################################                                           
#####       PREPARE  DATA       #####                                           
#####################################                                           
 
##### load all cell counts ####
cl <- makeCluster(threads, type="FORK", outfile="")
clusterExport(cl, c("args", "stDat"), envir=environment())

dat <- parLapply(cl, unique(stDat$sampAnn$CellDive_ID), function(cdid){
             
             sampdat <- tryCatch({
                 log_info("Loading data for ", cdid, ",")
                 ann_cells <- loadAnnotatedCells(stDat$sampAnn,
                                                 annDir = args$cell_data_dir,
                                                 cellDiveID = cdid)                                
                 if(is.null(ann_cells)){ return(NULL) }
                 lblcols <- grep("Category|Subtype|Cell_type", names(ann_cells))
                 ann_cells %>% 
                 group_by_at(c(args$x_var, names(.)[lblcols])) %>%
                 summarize(Count = n()) 
             }, error = function(){
                 log_error(paste0("Error summarizing cells for ", cdid))
                 return(NULL)
             })
        }) %>%
        bind_rows() %>%
        mutate(!!as.name(args$x_var) := factor(!!as.name(args$x_var), 
                                        levels = unique(!!as.name(args$x_var))))
stopCluster(cl)


#### TODO: MAKE THIS DYNAMIC BY ALLOWING USER TO SPECIFY OTHER CLINICICAL CHARS
####         TO FACET BY
plot_dat <- dat %>% left_join(stDat$sampAnn %>% select(Patient_ID, EBV_final))

####  prep data for plotting ####
clrs <- NULL
if(!is.null(args$plot_config_file)){
    clrs <- read_yaml(args$plot_config_file)
}

## fill in labels specified in cell_types table 
labels <- stDat$cellTypes %>%                                             
          select(dplyr::matches("Category|type")) %>%                
          unique
other <- c("superNeg", "UNKNOWN")
plot_dat <- dat %>% 
            left_join(labels,
                      by = intersect(names(.), names(labels))) %>%
            mutate(Category = ifelse(is.na(Category) | Category %in% other,
                                     "Other", Category),
                   Cell_type = ifelse(is.na(Cell_type) | Cell_type %in% other, 
                                  "Other", Cell_type),
                   Subtype = ifelse(is.na(Subtype) | Subtype %in% other, 
                                      Cell_type, Subtype))

#### TODO: MAKE THIS DYNAMIC BY ALLOWING USER TO SPECIFY OTHER CLINICICAL CHARS 
####         TO FACET BY                                                        
plot_dat <- plot_dat %>%
            left_join(stDat$sampAnn %>% select(Patient_ID, EBV_final)) %>%
            unique
plot_dat$EBV_final <- factor(plot_dat$EBV_final, levels = c("Positive", "Negative"))

facet_col <- "EBV_final"

## create list of count summaries by each level of cell identities
## Category, Cell_type, Subtyp
plist <- list()
pdata <- list()
for(cell_subset in c('All', 'Known', 'Immune_All')){
    lvls <- c("Category", "Cell_type", "Subtype")
    for(lvl in lvls){
        pdat <- plot_dat
        
        if(cell_subset == 'Known'){
            pdat <- pdat %>%
                    filter_at(vars(all_of(lvls)), all_vars(. != 'Other'))
        } else if(cell_subset != "All"){                                               
            pdat <- pdat %>%                                                    
                    filter_at(vars(any_of(lvls)), any_vars(. == cell_subset))   
        }                                                                       

        ct_order <- c("Other", rev(unique(stDat$cellTypes[[lvl]])))                  
        ct_order <- ct_order[ct_order %in% pdat[[lvl]]]
        pdat <- pdat %>%
                mutate(!!as.name(lvl) := 
                       factor(!!as.name(lvl), levels = ct_order)) %>%
                group_by_at(c(args$x_var, facet_col, lvl)) %>%
                summarize(Count = sum(Count, na.rm = TRUE)) %>%
                group_by_at(args$x_var) %>%
                mutate(Total = sum(Count, na.rm = TRUE), PCT = Count/Total)

        lbls <- levels(pdat[[lvl]])
        names(lbls) <- lbls
        if(lvl %in% c("Cell_type", "Subtype")){
            lbls <- make_figure_labels(stDat$cellTypes, lvl)
            lbls['Other'] <- 'Other'
        }

        pcolors <- clrs[ct_order]                             
        for(pct in c(FALSE, TRUE)){
print(paste(cell_subset, lvl, pct, sep=";"))
            if(pct & length(ct_order) == 1){ next } 
            p <- plot_cell_types(pdat,                                              
                                 pcolors,                                           
                                 xvar = args$x_var,                                 
                                 xvar_facet = facet_col,
                                 yvar = ifelse(pct, "PCT", "Count"),                         
                                 color_by = lvl,                                    
                                 clr_lbls = lbls,                   
                                 pct = pct) +                                     
                  ylab(ifelse(pct, 
                              paste("Percent", cell_subset, "Cells"), 
                              "Cell Count"))       

            plist[[length(plist) + 1]] <- p
###
### SCALED AXES DON'T GO WITH STACKED CHARTS
###
            ## add a plot with y axis scaled
            #if(!pct & lvl == "Category"){
            #    p_scaled <- p + scale_y_sqrt(labels = 
            #    plist[[length(plist) + 1]] <- p_scaled 
            #}
                                         
            if(pct & lvl == "Category"){
                ## add another plot with percentage labels just for the Category
                ## plot
                lbl_colors <- lapply(clrs[levels(pdat$Category)], function(x){
                                  ifelse(is_dark(x), 'white', 'black')
                              }) %>% unlist
                pdat <- pdat %>%
                        select(-Count, -Total) %>%
                        spread(Category, PCT, fill = 0) %>%
                        gather((length(facet_col) + 2):ncol(.), key = "Category", value = "PCT") %>%  
                        mutate(label = paste0(sprintf("%0.1f", PCT * 100), "%"), 
                               hjust = ifelse(Category == 
                                               ct_order[length(ct_order)] & 
                                                PCT < 0.1, -0.1, 1.1),
                               color = lbl_colors[Category])
                pdat$Category <- factor(pdat$Category, levels = ct_order) 
    
                if(!is.null(facet_col)){                                        
                    p <- p +                                                    
                         facet_grid(facet_col,                                  
                                    space = "free_y",                           
                                    scales = "free_y",                          
                                    switch = "both") +                          
                        theme(strip.text.y = element_text())                    
                }                                                               

                p <- p + 
                     geom_text(pdat, 
                               mapping = aes(label = label, 
                                             hjust = hjust,
                                             color = color), 
                               position = 'fill',
                               size = 3) +
                     coord_flip() +
                     scale_color_manual(values = c(white = 'white', 
                                                   black = 'black'),
                                        guide = NULL) +
                     theme(axis.text.y = element_text(size = 10)) 

                plist[[length(plist) + 1]] <- p           
            } 
        }                                                                           
        pdata[[paste0(cell_subset, "__", lvl, "__by__", args$x_var)]] <- pdat
    }
}

fn <- file.path(args$study_figure_dir, 
                paste0("cell_type__counts_by__", args$x_var, ".pdf"))
pdf(fn, height = 8 * 0.75, width = 11 * 0.75)
for(p in plist){
    print(p)
}
dev.off()

write.xlsx(pdata, gsub(".pdf", ".xlsx", fn), check.names = F)

