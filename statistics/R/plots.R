library(ComplexHeatmap)

is_dark <- function(colr) { 
              (sum( col2rgb(colr) * c(299, 587,114))/1000 < 123) 
           }

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

count_bar_theme <- theme_minimal() +
                   theme(strip.background = element_blank(),
                       strip.text.y = element_blank(),
                       #panel.background = element_rect(fill = "white", color = "darkgray"),
                       panel.grid.major.x = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line.x = element_line(size = 0.5, color = "black"),
                       axis.line.y = element_line(size = 0.5, color = "black"),
                       axis.text.y = element_text(size = 11, color = "black"),
                       axis.text.x = element_text(size = 11, color = "black", 
                                                  angle = 90, vjust=0.5),
                       legend.title = element_blank(),
                       legend.text = element_text(size = 10, color = "black", 
                                                  hjust = 0))


make_figure_labels <- function(cell_types, label_type){
    
    main <- paste0(label_type, "_figures")
    subs <- paste0(label_type, "_figures_subscript")

    dat <- cell_types %>% 
           select(all_of(c(label_type, main, subs))) %>% 
           unique %>%
           mutate(!!as.name(main) := ifelse(is.na(!!as.name(main)), 
                                             !!as.name(label_type), 
                                               !!as.name(main)),
                  !!as.name(subs) := ifelse(is.na(!!as.name(subs)), "", 
                                              !!as.name(subs)))

    lbls <- lapply(1:nrow(dat), function(x){
               ss <- dat[[subs]][x]
               mn <- gsub("\\+", "", gsub(ss, "", dat[[main]][x]))
               bquote(.(mn)[.(ss)]) 
            }) %>% 
            as.expression()
    names(lbls) <- dat[[label_type]]
    lbls
}


plot_fov_variation_heatmap <- function(counts, cell_category = "Subtype", 
                                       avg_within = "CellDive_ID", 
                                       cell_type_labels = NULL, 
                                       data_col = "Fraction", max_val = 3){

    pdat <- counts %>%
            group_by_at(all_of(cell_category)) %>%
            mutate(log_val = log(!!as.name(data_col) + 0.001), 
                   stdev = sd(log_val)) %>%
            group_by_at(c(avg_within, all_of(cell_category))) %>%
            mutate(mn = mean(log_val),
                   val = (log_val - mn)/stdev) %>%
            mutate(val = ifelse(abs(val) > max_val, 
                                 ifelse(val < 0, max_val*-1, max_val), 
                                         val))


   ggplot(pdat, aes_string(x = "FOV_ID", y = cell_category, fill = "val")) +
   geom_tile() +
   xlab("FOV") +
   ylab("") +
   scale_y_discrete(labels = y_labels, position = "right", expand = c(0,0)) +
   scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", 
                         midpoint = 0, limits = c(-max_val,max_val)) +
   theme_minimal() +
   theme(panel.border = element_rect(color = "black", fill = NA),
         axis.text.x = element_blank(),
         axis.title.x = element_text(size = 14, color = "black"),
         axis.text.y = element_text(size = 14, color = "black"),
         legend.title = element_blank(),
         legend.position = "bottom",
         legend.direction = "horizontal")
}

plot_fov_variation_annotated_heatmap <- function(counts, cell_category = "Subtype",
                                                 avg_within = "CellDive_ID", 
                                                 data_col = "Fraction",
                                                 max_val = 3,                   
                                                 annotation_columns = NULL, 
                                                 annotation_colors = NULL,
                                                 cell_type_labels = NULL,       
                                                 data_file = NULL){

    pdat <- counts %>%
            group_by_at(all_of(cell_category)) %>%
            mutate(log_val = log(!!as.name(data_col) + 0.001),
                   stdev = sd(log_val)) %>%
            group_by_at(c(avg_within, 
                         all_of(cell_category), 
                         all_of(annotation_columns))) %>%
            mutate(mn = mean(log_val),
                   val = (log_val - mn)/stdev) %>%
            mutate(val = ifelse(abs(val) > max_val, 
                            ifelse(val < 0, max_val*-1, max_val), 
                               val)) %>%
            ungroup() %>%
            select(all_of(cell_category), FOV_ID, val) %>%
            spread(FOV_ID, val, fill = 0)

    if(!is.null(data_file)){
        write.xlsx(pdat, data_file)
    }

    rn <- pdat[[cell_category]]
    pdat <- pdat %>% select(-all_of(cell_category)) %>% as.matrix
    rownames(pdat) <- rn

    ann <- counts %>% 
           ungroup %>% 
           select(all_of(c("FOV_ID", annotation_columns, avg_within))) %>%
           unique %>%
           select(-FOV_ID) %>%
           as.data.frame

    avg_within_clrs = gg_color_hue(length(unique(counts[[avg_within]]))) %>%
                       setNames(levels(counts[[avg_within]]))

    ann_clrs <- lapply(setdiff(annotation_columns, avg_within), function(x){
                    annotation_colors[unique(ann[[x]])]
                }) %>%
                setNames(setdiff(annotation_columns, avg_within))
    clrs <- ann_clrs
    clrs[[avg_within]] <- avg_within_clrs                 

    ha <- HeatmapAnnotation(df = ann,  
                            #annotation_label = avg_within,
                            annotation_legend_param = 
                                 list(direction = "vertical", 
                                      ncol = 4),
                            col = clrs)

    htmp <- Heatmap(pdat, 
                    cluster_rows = F, 
                    row_order = rev(levels(counts[[cell_category]])), 
                    row_labels = cell_type_labels[rownames(pdat)],
                    show_column_names = F,
                    top_annotation = ha,
                    column_title = paste0("FOV (n=", 
                                           length(levels(counts$FOV_ID)), 
                                          ")"),
                    column_title_side = "bottom",
                    heatmap_legend_param = list(title = "z-score"))

    ComplexHeatmap::draw(htmp, annotation_legend_side = "bottom")
}

plot_cell_types <- function(pdat, fill_clrs, xvar = "Patient_ID", yvar = "Count", 
                            xvar_facet = NULL,
                             color_by = "Cell_type", clr_lbls = NULL, 
                             ptheme = count_bar_theme, pct = FALSE){

    ymax <- pdat %>% 
            group_by_at(xvar) %>%
            summarize(tot = sum(!!as.name(yvar))) %>% 
            pull(tot) %>% 
            max

    yscale <- scale_y_continuous(limits = c(0, ymax), expand = c(0,0.5), 
                                 labels = scales::number_format(big.mark = ",")) 
    if(pct){
        yscale = scale_y_continuous(limits = c(0,ymax), expand = c(0,0), 
                                    labels = scales::percent)
    }

    p <- ggplot(pdat, 
                aes_string(x = xvar, 
                           y = paste0("`", yvar, "`"), 
                           fill = color_by)) +
         geom_bar(stat = "identity", position = position_stack()) + # "stack") +
         scale_fill_manual(values = fill_clrs, 
                           labels = clr_lbls, 
                           na.value = "grey80",
                           guide = guide_legend(reverse = TRUE) ) +
         yscale +
         ptheme +
         theme(axis.title.x = element_blank())

    if(!is.null(xvar_facet)){
        p <- p +
             facet_grid(as.formula(paste0(". ~ `",xvar_facet,"`")),
                        space = "free_x", 
                        scales = "free_x", 
                        switch = "both") + 
             theme(strip.text.x = element_text())
    }

    p
}


plot_fov_count_per_sample <- function(fovs, sampCol = "Patient_ID", 
                                        ptheme = count_bar_theme){

    fovs <- fovs %>% group_by_at(sampCol) %>% mutate(y = row_number()) 
   
    ggplot(fovs, aes(x = !!as.name(sampCol), y = y)) +
    geom_point(size = 1) +
    scale_y_continuous(limits = c(0, max(fovs$y) + 5), expand = c(0,0)) +
    ylab("FOV count") +
    ptheme 

}


plot_study_overview <- function(sampledat, clrs, x = "Patient_ID", 
                                 ptheme = NULL){

    if(is.null(ptheme)){
       ptheme <- theme_minimal() +
                 theme(strip.background = element_blank(),
                       strip.text.y = element_blank(),
                       panel.background = element_rect(fill = "white", 
                                                       color = "darkgray"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.y = element_text(size = 11, color = "black"),
                       axis.text.x = element_text(size = 11, color = "black", 
                                                   angle = 90, vjust=0.5),
                       legend.title = element_blank(),
                       legend.text = element_text(size = 10, color = "black"))
    }

    #if(is.null(clrs)){
    #    yes <- '#023E8A'
    #    no  <- '#f0f0f0'
    #    pos <- yes
    #    neg <- '#f0f0f0'

    #    clrs <- c('F' = '#5e548e', 'M' = '#ffd97d',     ## sex
    #              'NS' = '#8f2d56', 'MC' = '#19647e', 'LR' = '#d8973c', 'NLPHL' = '#455972', ## tumor subtype
    #              'I' = '#A3B18A', 'II' = '#588157', 'III' = '#3A5A40', 'IV' = '#344E41', ## clinical stage
    #              'LN' = '#ad2831',                   ## site
    #              'Positive' = pos, 'Negative' = neg, #'Negative' = '#1f487e', ## EBV status
    #              'Untreated' = no, 'Treated' = yes,  ## Treatment
    #              'Yes' = yes, 'No' = no,             ## Imm before, other before, imm after, other after        
    #              '1' = yes, '0' = no,                ## we have data or not
    #              'T-cell' = '#2073ac', 'B-cell' = '#33a547', 'indeterminate' = 'grey15', ## main cell population
    #              'pos' = pos, 'neg' = neg, 'cytoplasmic' = '#e8871e',  ## IHC
    #              'unknown' = 'grey75', 'Other' = 'gray50')
    #}

    ggplot(sampledat, 
           aes(x = !!as.name(x), y = figure_label, fill = Val)) +
    geom_tile(aes(width = 0.99, height = 0.99), 
              color = "gray85") +
    geom_point(sampledat %>% filter(is.na(Val)),
               mapping = aes(x = !!as.name(x), 
                             y = figure_label, 
                             shape = '4', 
                             color = 'darkgray'), 
               fill = NA) +
    facet_grid(Section ~ ., scales = "free", space = 'free') +
    scale_shape_manual(values = c('4' = 4), 
                       breaks = c('4'), 
                       name = "NA", 
                       labels = "NA") +
    scale_color_manual(values = c('darkgray' = 'darkgray'), 
                       breaks = c('darkgray'), 
                       name = "x", 
                       guide = "none", 
                       na.value = 'white') +
    scale_fill_manual(values = clrs, 
                      name = "value", 
                      breaks = names(clrs), 
                      na.value = 'white') +
    xlab("") +
    ptheme

}


plot_accumulated_cell_loss <- function(drift_summary, groups = "Patient_ID", 
                                        justify_labels = TRUE, 
                                        ptheme = NULL, clrs = NULL, 
                                        col_pals = NULL){

    ds <- drift_summary %>% 
          group_by_at(groups) %>%
          mutate(MaxLoss = max(accumulated_cell_loss, na.rm=T),
                 Label = !!as.name(groups)) %>%
          arrange(MaxLoss) %>%
          ungroup()    

    if(justify_labels){
        ds <- ds %>%
              mutate(longestID = max(nchar(!!as.name(groups))),
                     str_num = paste0("(", 
                                     formatC(accumulated_cell_loss, 
                                             big.mark=","), 
                                     ")"),
                     longestLabel = max(nchar(str_num) + longestID),
                     spacesNeeded = longestLabel - (nchar(!!as.name(groups)) + 
                                     nchar(str_num)) + 1)
        ds$spaces <- sapply(1:nrow(ds), function(x){ 
                        paste(rep(" ", ds[x,] %>% pull(spacesNeeded)), 
                              collapse = "") 
                     })
        ds$Label <- paste0(ds[[groups]], ds$spaces, ds$str_num)
    }
    ds$Label <- factor(ds$Label, levels = unique(ds$Label))
    ds[[groups]] <- factor(ds[[groups]], levels = unique(ds[[groups]]))

    lbls <- ds %>% 
            filter(MaxLoss == accumulated_cell_loss) %>% 
            pull(Label) %>% 
            unique %>% 
            as.character
    names(lbls) <- levels(ds[[groups]])

    if(is.null(clrs)){
        if(is.null(col_pals)){
            col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
            col_pals = col_pals[c("Dark2", "Set1", "Set2", "Set3", "Paired", 
                                  "Accent", "Pastel1", "Pastel2"),]
        }
        clrs = unlist(mapply(brewer.pal, col_pals$maxcolors, 
                              rownames(col_pals)))[1:length(lbls)]
        names(clrs) <- names(lbls)
    }

    ggplot(ds, 
           aes(x = stage, 
               y = accumulated_cell_loss/1000, 
               color = !!as.name(groups))) +
    geom_line() +
    ylab("Accumulated Cell Loss (k)") +
    xlab("Stage") +
    scale_color_manual(values = clrs, 
                       labels = lbls, 
                       breaks = rev(levels(ds[[groups]]))) +
    theme_minimal() +
    ptheme +
    guides(color = guide_legend(ncol = 2,
                                label.position = "left",
                                label.hjust = 1))

}

plot_accumulated_cell_loss_percent <- function(drift_summary, 
                                                groups = "Patient_ID", 
                                                clrs = NULL, 
                                                justify_labels = TRUE, 
                                                ptheme = NULL, col_pals = NULL){

    ds <- drift_summary %>%
          group_by_at(groups) %>%
          mutate(TotalPct = max(loss_pct, na.rm=T),
                 Label = !!as.name(groups)) %>%
          arrange(TotalPct) %>%
          ungroup()

    if(justify_labels){
        ds <- ds %>%
              mutate(longestID = max(nchar(!!as.name(groups))),
                     MaxPct = round(TotalPct * 100, 2),
                     str_num = paste0(" (", 
                                     trimws(format(round(MaxPct, 2), nsmall=2)), 
                                            "%)"),
                     longestLabel = max(nchar(str_num) + longestID),
                     spacesNeeded = longestLabel - (nchar(!!as.name(groups)) + 
                                     nchar(str_num)))
        ds$spaces <- sapply(1:nrow(ds), function(x){ 
                       paste(rep(" ", ds[x,] %>% pull(spacesNeeded)), 
                             collapse = "") })
        ds$Label <- paste0(ds[[groups]], ds$spaces, ds$str_num)
    }
    ds$Label <- factor(ds$Label, levels = unique(ds$Label))
    ds[[groups]] <- factor(ds[[groups]], levels = unique(ds[[groups]]))

    lbls <- ds %>% 
            filter(TotalPct == loss_pct) %>% 
            pull(Label) %>% 
            unique %>% 
            as.character()
    names(lbls) <- levels(ds[[groups]])

    if(is.null(clrs)){
        if(is.null(col_pals)){
            col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
            col_pals = col_pals[c("Dark2", "Set1", "Set2", "Set3", "Paired", 
                                  "Accent", "Pastel1", "Pastel2"),]
        }
        clrs = unlist(mapply(brewer.pal, col_pals$maxcolors, 
                              rownames(col_pals)))[1:length(lbls)]
        names(clrs) <- names(lbls)
    }

    ggplot(ds, aes(x = stage, y = loss_pct, color = !!as.name(groups))) +
    geom_line() +
    scale_y_continuous(limits = c(0,1), label=scales::percent) +
    ylab("Accumulated Cell Loss") +
    xlab("Stage") +
    scale_color_manual(values = clrs, 
                       labels = lbls, 
                       breaks = rev(levels(ds[[groups]]))) +
    theme_minimal() +
    ptheme +
    guides(color = guide_legend(ncol = 2,
                                label.position = "left",
                                label.hjust = 1))

}


plot_percent_cells_remaining <- function(drift_summary, groups = "Patient_ID", 
                                         clrs = NULL, justify_labels = TRUE, 
                                         ptheme = NULL, col_pals = NULL, 
                                         zoom_y = NULL){

    ds <- drift_summary %>%
          group_by_at(groups) %>%
          mutate(TotalPct = max(loss_pct, na.rm=T),
                 Label = !!as.name(groups),
                 RmnPct = 1 - TotalPct,
                 MinPct = round(RmnPct * 100, 2)) %>%
          ungroup() %>%
          arrange(desc(RmnPct))

    if(justify_labels){
        ds <- ds %>%
              mutate(longestID = max(nchar(!!as.name(groups))),
                     str_num = paste0(" (", 
                                      trimws(format(round(MinPct, 2), nsmall=2)), 
                                      "%)"),
                     longestLabel = max(nchar(str_num) + longestID),
                     spacesNeeded = longestLabel - (nchar(!!as.name(groups)) + 
                                     nchar(str_num))) %>%
          ungroup() %>%
          arrange(desc(RmnPct))

          ds$spaces <- sapply(1:nrow(ds), function(x){ 
                         paste(rep(" ", ds[x,] %>% pull(spacesNeeded)), 
                               collapse = "") })
          ds$Label <- paste0(ds[[groups]], ds$spaces, ds$str_num)
    }
    ds$Label <- factor(ds$Label, levels = unique(ds$Label))
    ds[[groups]] <- factor(ds[[groups]], levels = unique(ds[[groups]]))

    lbls <- ds %>% 
            filter(RmnPct == 1 - loss_pct) %>%                                                          
            pull(Label) %>% 
            unique %>% 
            as.character()
    names(lbls) <- levels(ds[[groups]])

    if(is.null(clrs)){
        if(is.null(col_pals)){
            col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
            col_pals = col_pals[c("Dark2", "Set1", "Set2", "Set3", "Paired", 
                                   "Accent", "Pastel1", "Pastel2"),]
        }
        clrs = unlist(mapply(brewer.pal, col_pals$maxcolors, 
                       rownames(col_pals)))[1:length(lbls)]
        names(clrs) <- names(lbls)
    }

    pzoom <- ggplot(ds,
                    aes(x = stage, y = 1 - loss_pct, color = !!as.name(groups))) +
             geom_line() +
             scale_y_continuous(limits = zoom_y, 
                                expand = c(0,0), 
                                label=scales::percent) +
             scale_x_continuous(expand = c(0, 0)) +
             ylab("Cells Remaining") +
             xlab("Stage") +
             scale_color_manual(values = clrs, 
                                labels = lbls, 
                                breaks = levels(ds[[groups]])) +
             theme_bw() +
             ptheme +
             theme(legend.text = element_text(size = 10, 
                                              margin = unit(c(0,0,0,0.5), "in"), 
                                              hjust = 0)) +
             guides(color = guide_legend(ncol = 2,
                                         label.position = "left",
                                         label.hjust = 1)) 
    
    pfull <- pzoom + 
             scale_y_continuous(limits = c(0,1), 
                                expand = c(0,0), 
                                label = scales::percent)
    legend <- getLegend(pfull) 

    if(is.null(zoom_y)){
        return(list(pfull))
    }
 
    tmarg <- bmarg <- 1
    pzoom <- pzoom + 
             theme(plot.margin = margin(t = tmarg, r = 1, b = bmarg, l = 0, "in"),
                   legend.position = "none",
                   axis.title.y = element_blank())

    pfull <- pfull +
             geom_hline(yintercept = min(zoom_y), size = 0.5) +
             geom_hline(yintercept = max(zoom_y), size = 0.5) +
             geom_ribbon(aes(ymin = min(zoom_y), ymax = max(zoom_y)), 
                         size = 0.5, fill = "#e0e0e0", alpha = 0.5, color = NA) +
             geom_line() +
             theme(legend.position = "none",
                   plot.margin = margin(t = tmarg, r = 0, b = bmarg, l = 1, "in"))

    pz <- ggplot_build(pzoom)
    ylbls <- pz$layout$panel_params[[1]]$y$get_labels()
    lblY  <- 1/(length(ylbls) - 1) * c(0, seq(ylbls[-1]))

    magdat <- tibble(x = c(0, 0, 1, 1), y = c(zoom_y, 1, 0), order = c(1,2,2,1))  
    axisDat <- tibble(y = 1/length(ylbls) * seq(ylbls), label = ylbls) 
    mag <- ggplot(magdat, aes(x = x, y = y)) +
           geom_line(aes(group = order), size = 0.5) + 
           geom_polygon(fill = "#e0e0e0", alpha = 0.5, color = NA) +
           scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
           annotate("text", 
                    x = 0.85, y = lblY, 
                    label = ylbls, 
                    size = 14 * 0.352) +
           coord_cartesian(clip = "off") +
           theme_void() +
           theme(plot.margin = 
                   margin(r = -0.15, l = -0.1, t = tmarg, b = bmarg, "in"))
    pzoom <- pzoom + theme(axis.text.y = element_blank())

    pg <- plot_grid(pfull, mag, pzoom, 
                    rel_widths = c(1, 0.5, 1), 
                    axis = 'bt', 
                    align = 'h', 
                    nrow = 1)
    #pg <- plot_grid(pg, legend, rel_heights = c(1, 1), align = 'v', ncol = 1)

    list(pg, legend)
}


plot_total_drift <- function(drift_summary, groups = "Patient_ID", 
                              bptheme = NULL, pct = FALSE, 
                               sort_by = 'remaining'){

    clrs <- c('cells lost' = 'black', 
              'pct cells lost' = 'black',
              'cells remaining' = 'gray', 
              'pct cells remaining' = 'gray')

    ds <- drift_summary %>%
          group_by_at(groups) %>%
          mutate(`cells lost` = max(accumulated_cell_loss), 
                 `cells remaining` = total_count - `cells lost`,
                 `pct cells lost` = `cells lost`/total_count, 
                 `pct cells remaining` = `cells remaining`/total_count) %>%
          select(all_of(groups), dplyr::matches("cells"), total_count) %>%
          unique %>%
          arrange(desc(!!as.name(sort_by))) 
    ds[[groups]] <- factor(ds[[groups]], levels = unique(ds[[groups]]))

    if(pct){
        dat <- ds %>% 
               select(all_of(groups), dplyr::matches('pct')) %>%
               gather(dplyr::matches('pct'), key = 'cat', value = 'pct') %>%
               mutate(cat = factor(cat, 
                                   levels = c('pct cells remaining', 
                                              'pct cells lost')))

        bp <- ggplot(dat, aes(x = !!as.name(groups), y = pct, fill = cat)) +
              geom_bar(stat = "identity", 
                       position = "stack", 
                       color = "black", 
                       size = 0.25) +
              scale_y_continuous(limits = c(0, 1), 
                                 expand = c(0,0.01), 
                                 labels = scales::percent) +
              scale_fill_manual(values = clrs[levels(dat$cat)], 
                                labels = gsub("pct ", "", levels(dat$cat))) +
              xlab("") +
              ylab("Cell Percentage") +
              bptheme
        return(bp)
    }

    dat <- ds %>%
           select(all_of(groups), !dplyr::matches('pct')) %>%
           gather(dplyr::matches('cells'), key = 'cat', value = 'count') %>%
           mutate(cat = factor(cat, levels = c('cells remaining', 'cells lost')))

    ggplot(dat, 
           aes(x = !!as.name(groups), y = count/1000000, fill = cat)) +
    geom_bar(stat = "identity", 
             position = "stack", 
             color = "black", 
             size = 0.25) +
    scale_y_continuous(limits = c(0, ceiling(max(dat$count)/1000000)), 
                       expand = c(0,0.01)) +
    scale_fill_manual(values = clrs[levels(dat$cat)], 
                      labels = levels(dat$cat)) +
    xlab("") +
    ylab("Cell Count (millions)") +
    bptheme

}



blank <- function(dat, blockClr = "#e8e8e8", blockH = 0.98, blockW = 0.98){
    ggplot(dat,
           aes(x = `Patient ID`, y = Var),
           color = blockClr) +
    geom_tile(width = blockW, height = blockH, fill = blockClr)
}

study_sample_heatmap <- function(dat, scale = "discrete", 
                                facet_cols = "Patient ID",
                                clrs = NULL, palette = "YlGnBu",
                                gradientL = NULL, gradientH = NULL,
                                blockClr = "#e8e8e8", naClr = "#e8e8e8",
                                blockH = 0.98, blockW = 0.98, clrH = 0.92, 
                                clrW = 0.92, plot_theme = NULL, legendNcol = 1){

    if(is.null(plot_theme)){
        plot_theme <- theme(axis.text.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     legend.position = "right",
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank())
    }

    p <- blank(dat, blockClr = blockClr, blockH = blockH, blockW = blockW) +
         geom_tile(aes(x = `Patient ID`, y = Var, fill = Val),
                   color = blockClr,
                   width = clrW,
                   height = clrH) +
         scale_y_discrete(expand = c(0,0)) +
         theme_minimal() +
         plot_theme 

    if(!is.null(facet_cols)){
        p <- p + facet_grid(as.formula(paste0(". ~ `",facet_cols,"`")), 
                            scales = "free", space = "free")
    }

    colorFunc <- scale_fill_manual(values = clrs, 
                                   na.value = naClr, 
                                   drop = FALSE)
    if(is.null(clrs)){
        colorFunc <- scale_fill_brewer(palette = palette, drop = FALSE)
    }

    if(scale == "discrete"){
        if(is.null(clrs)){
            p <- p + scale_fill_brewer(palette = palette, 
                                       drop = FALSE, 
                                       na.value = naClr)
        }  else {
            p <- p + scale_fill_manual(values = clrs, 
                                       na.value = naClr,  
                                       drop = FALSE)
        }
        if(legendNcol > 0){
            p <- p + guides(fill = guide_legend(ncol = legendNcol))
        }
    } else if(scale == "continuous"){
       p <- p + scale_fill_gradient(low = gradientL, 
                                    high = gradientH, 
                                    guide = "colourbar", 
                                    na.value = naClr)
    } else {
        stop(paste0("Unrecognized scale: ",scale))
    }
    p
}

addNoTumor <- function(plotDat, plot, noTumorVal = 'No tumor', 
                       lineColor = "#e8e8e8"){
    dat <- plotDat %>%
           group_by(`Patient ID`, Sample_ID) %>%
           mutate(top = row_number() + 0.5, bot = row_number() -0.5) %>%
           group_by(`Patient ID`, Var) %>%
           mutate(left = row_number() - 0.5, right = row_number() + 0.5)
    plot +
    geom_tile(data = dat %>% filter(Val == noTumorVal), width = 0.85, height = 0.85, fill = "white") +
    geom_segment(data = dat %>% filter(Val == noTumorVal),
                 aes(x = left + 0.1, xend = right - 0.05, y = top - 0.05, yend = bot + 0.05),
                 size = 0.25, color = "gray") +
    geom_segment(data = dat %>% filter(Val == noTumorVal),
                 aes(x = left + 0.05, xend = right - 0.05, y = bot + 0.05, yend = top - 0.05),
                 size = 0.25, color = "gray")
}


pseudo_color_fov <- function(fovDat, clrs, colorBy = "Cell_type", pointSize = 0.05){
    
    ggplot(fovDat, aes_string(x = "x0", y = "y0", color = colorBy)) +
    geom_point(size = pointSize) +
    scale_color_manual(values = clrs) + 
    theme_minimal() + 
    facet_wrap(FOV_number ~ .)


}
