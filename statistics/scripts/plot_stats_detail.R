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
            --detail_figure_config_file      YAML file of figure configuration including
                                             facet order & biological filters
            --fov_area_dir                   path to RDA files, each containing table of FOVs and 
                                             total FOV area for all FOVs in a single sample 
            --fov_metrics_dir                root directory of all precalculated metrics
            --meta_dir                       path to meta files in XLSX format
            --plot_config_file               YAML file containing all color assignments 
            --statistics_conditions_file     XLSX file listing all cell states/conditions to 
                                             compare between two sample groups
            --statistics_conditions_index    XLSX file with pre-indexed cell states/conditions
            --statistics_config_file         YAML file containing stats config including results filters
            --statistics_detail_dir          directory where plot data and figures will be saved
            --statistics_tables_dir          output directory where XLSX files of results should 
                                             be written

          [OPTIONAL]
            --detail_figure_conditions_file  XLSX file containing a column for each question, 
                                             each column containing any number of Cell State IDs
                                             that should be included in the figure
            --neighbor_dir                   directory containing RDA files of all pairwise distances 
                                             between cells, at least those <= 30 microns. generally 
                                             each file will contain a table of all 'center' cells of
                                             a certain type, but this is not a requirement as all files
                                             will be loaded together; required if question data is restricted
                                             to tumor or macrophage neighborhoods
            --question                       question number of comparison to be plotted
            --plot_calculations               calculation to be plotted (densities or fractions)
            --manifest                       YAML file containing one or more parameter; NOTE: 
                                             arguments on command line override manifest arguments!!!         
            --microenvironment_metrics_dir   directory containing RDA files of microenvironment counts and fractions 
                                             not required if no questions involve microenvironment data 
        \n"
    )
}

## names of required args
minReq <- c("detail_figure_config_file", "fov_area_dir", "fov_metrics_dir", 
            "meta_dir", "plot_config_file",
            "statistics_conditions_file", "statistics_conditions_index",
            "statistics_config_file", "statistics_detail_dir", "statistics_tables_dir") 

used <- c(minReq, "detail_figure_conditions_file", "neighbor_dir", 
          "question", "plot_calculations", "manifest",
          "microenvironment_metrics_dir")

defaults <- list(detail_figure_conditions_file = NULL)

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

###################################
##   CONFIGURE & INITIALIZE DATA ##
###################################
cfg <- resolveConfig(read_yaml(args$plot_config_file), 
                     read_yaml(args$detail_figure_config_file), 
                     read_yaml(args$statistics_config_file), 
                     args)

if(!is.null(cfg$plot_calculations)){ 
    cfg$calculations <- cfg$calculations[cfg$plot_calculations]
}

logParams(cfg, used)

calcUnit   <- "FOV_ID"
statsFiles <- getFiles(path = cfg$statistics_tables_dir) 
rel_widths <- list("fractions" = c(3.6, 1.25, 1.25),
                   "densities" = c(2.5, 1.25, 1.25))

## we don't need to load any macro neighborhood counts because
## here we are only plotting general fractions/densites (no spatial
## conditions)
stDat <- loadStudyData(cfg, 
                       annotatedCells = F,
                       analyses = T, 
                       conditions = T, 
                       questions = T) 

###################################
##              PLOT             ##
###################################
for(calc in names(cfg$calculations)){

    analyses <- stDat$analysisList[names(stDat$analysisList) == calc]

    conds2plot <- list()
    if(!is.null(cfg$detail_figure_conditions_file) && !is.na(cfg$detail_figure_conditions_file)){ 
        if(!file.exists(cfg$detail_figure_conditions_file)){
            log_error(paste0("Could not find file: ", cfg$detail_figure_conditions_file))
            q(status = 1, save = "no")
        }
        conds2plot <- read.xlsx(cfg$detail_figure_conditions_file, 1, check.names=F)
        conds2plot <- conds2plot[,which(colnames(conds2plot) != "NA")] %>% as_tibble()
    }

    for(q in names(stDat$allQuestions)){

        outXLSX <- file.path(cfg$statistics_detail_dir, paste(q, calc, "detail.xlsx", sep = "_"))
        outPDF  <- file.path(cfg$statistics_detail_dir, paste(q, calc, "detail.pdf", sep = "_"))

        metDir <- cfg$fov_metrics_dir
        sheetName <- "all_fractions_and_densities"

        if(any(grepl("Neighborhood", names(unlist(stDat$allQuestions[[q]]))))){
            metDir <- cfg$microenvironment_metrics_dir
            sheetName <- "all_fractions"

            if(calc == "densities"){
                ## create empty files to allow pipeline to continue/complete
                file.create(outXLSX)
                file.create(outPDF)
                next ## questions with data restricted to tumor neighborhoods only 
                     ## have stats on fraction values
                ### TODO: FIGURE OUT A BETTER WAY TO HANDLE THIS
            }
        } 

        log_info(paste0(q, ": Formatting data for plotting ", toupper(calc), "..."))
        ids <- NULL
        orderBy <- cfg$calculations[[calc]]$effect_col
        if(q %in% names(conds2plot)){
            ids <- conds2plot[[q]][!is.na(conds2plot[[q]])] %>% as.numeric %>% rev
            orderBy <- NULL
            log_debug(paste0("Including ", length(ids), " conditions specified in file ",
                             cfg$detail_figure_conditions_file))
        }
        plotDat <- getStatsDetailData(statsFiles[basename(statsFiles) == paste0(q, ".xlsx")], 
                                      stDat$sampAnn, 
                                      stDat$allQuestions[[q]], 
                                      stDat$conds, 
                                      stDat$cellTypes, 
                                      analyses,
                                      stDat$markers, 
                                      metDir,
                                      sheetName      = sheetName, 
                                      calcUnit       = "FOV_ID",
                                      calculation    = calc, 
                                      calcColumn     = cfg$calculations[[calc]]$calc_column, 
                                      cellStateIDs   = ids,
                                      facets         = c("Cell_type","Subtype"), 
                                      statsFilters   = cfg$results_filters[[cfg$use_filter]][[calc]], 
                                      idOrder        = ids,
                                      orderBy        = orderBy, 
                                      nbhdCounts     = stDat$nbhdCounts) 

        if(is.null(plotDat) || length(plotDat) == 0 || 
              any(unlist(lapply(plotDat, function(x) is.null(x) || nrow(x) == 0)))){ 
            log_warn(paste0("Data insufficient to generate figure."))
            ## create empty files to allow pipeline to continue/complete
            file.create(outXLSX)
            file.create(outPDF)
            next
        }

### TEMP FIX
ct_lvls <- levels(plotDat$labels$Cell_type)
st_lvls <- levels(plotDat$labels$Subtype)

plotDat$fovVals$Cell_type <- factor(plotDat$fovVals$Cell_type, levels = ct_lvls)
plotDat$stats$Cell_type <- factor(plotDat$stats$Cell_type, levels = ct_lvls) 
#plotDat$labels$Cell_type <- factor(plotDat$stats$Cell_type, levels = unique(plotDat$labels$Cell_type))

plotDat$fovVals$Subtype <- factor(plotDat$fovVals$Subtype, levels = st_lvls)
plotDat$stats$Subtype <- factor(plotDat$stats$Subtype, levels = st_lvls) 
#plotDat$labels$Subtype <- factor(plotDat$stats$Subtype, levels = unique(plotDat$labels$Subtype))


        xlsxout <- file.path(cfg$statistics_detail_dir, paste(q, calc, "detail.xlsx", sep="_"))
        write.xlsx(plotDat, xlsxout, check.names = F)

        clrs <- c("#009ed8", "#ff8a30") # blue, orange
        names(clrs) <- levels(plotDat$fovVals$GroupLabel)

        pdfHeight <- (0.10 * nrow(plotDat$labels)) + 2.5
        pdfWidth  <- ifelse(is.null(cfg$pdfWidth), 13, as.numeric(cfg$pdfWidth))

        log_info("  Generating figure...")
        tryCatch({

#source("/home/byrne/halo/dev/hodgkins_dev/source_all.R")
            p <- plotQuestionResults(plotDat$fovVals, 
                                     plotDat$stats, 
                                     plotDat$labels, 
                                     clrs, 
                                     "GroupLabel",
                                     calcType    = calc,
                                     cellTypes   = stDat$cellTypes,
                                     xVar        = cfg$calculations[[calc]]$calc_column,
                                     yVar        = "y",
                                     yNudge      = 0.1,
                                     fontsize    = 12,
                                     facetY = c("Cell_type", "Subtype"),
                                     separateFacets = FALSE,
                                     stripeBG = FALSE,
                                     spacerColor = "#e0e0e0",
                                     rel_widths  = rel_widths[[calc]], 
                                     stripWidth = unit(4,"cm"))

            log_info(paste0("  Saving figure to file: ", outPDF))
            ggsave(p, height = pdfHeight, width = pdfWidth, units = "in",
                   device = cairo_pdf, filename = outPDF)

          }, error = function(e){
                print(e)
        })
    }
}
