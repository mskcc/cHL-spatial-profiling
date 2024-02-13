source("/home/byrne/halo/dev/hodgkins_dev/source_all.R")


#dat <- read.xlsx("/juno/res/bic/shared/Multiomyx/SocciN/Results/2021-10-18/markerComboCounts_Hodgkins_PhaseI_v10.7__210714_b_CTDv____14858b32_20211018_150717_.xlsx", 1, check.names = F) %>%
#       as_tibble()

#dat <- read.xlsx("markerComboCounts_Hodgkins_PhaseI_v10.7__210714_b_CTDv____3887616d_20211025_125306_.xlsx", 1, check.names = F) %>%
#       as_tibble() 

ids <- read.xlsx("input/meta/HodgkinLymphoma_Samples__V0.xlsx", 1, check.names = F) %>%
       as_tibble() %>%                                                          
       select(CellDive_ID, Patient_ID, EBV_final)

ctOrder <- read.xlsx("input/meta/HodgkinLymphoma_CellTypes__V0.xlsx", 1, check.names = F) %>%
           as_tibble() %>%                                                      
           select(Cell_type, Cell_type_figures, Cell_type_figures_subscript, Subtype, Subtype_figures, Subtype_figures_subscript)

dat <- read.xlsx("results/counts/marker_combination_cell_type_counts.xlsx", 4, check.names = F) %>% 
       as_tibble 

pdat <- dat %>% 
       left_join(ctOrder %>% select(Cell_type, Subtype), by = "Subtype") %>%
       select(-Total, -DF) %>%
       select(Subtype, Cell_type, everything()) %>%
       gather(3:ncol(.), key = CellDive_ID, value = `Cell Count`) %>% 
       left_join(ids) %>% 
       select(-CellDive_ID) %>%
       left_join(ctOrder) %>%                                                  
       mutate(Patient_ID = factor(Patient_ID, levels = sort(unique(Patient_ID))),
              Subtype = factor(Subtype, levels = rev(c(ctOrder$Subtype, "superNeg", "UNKNOWN"))),
              EBV_final = factor(EBV_final, levels = c("Positive", "Negative")))
         

clrs <- read_yaml("input/config/global_plot_colors.yaml")
vals <- unique(c(as.character(pdat$Cell_type), as.character(pdat$Subtype)))
clrs <- clrs[vals] %>% unlist
names(clrs)[names(clrs) == "Unknown"] <- "UNKNOWN"


lbls <- lapply(1:nrow(ctOrder %>% filter(!is.na(Subtype))), function(x){
             sbs <- ctOrder$Subtype_figures_subscript[x]
             lbl <- ifelse(!is.na(sbs), gsub(sbs, "", ctOrder$Subtype_figures[x]), ctOrder$Subtype_figures[x])
             if(is.na(sbs)){ return(lbl) }
             bquote(.(lbl)[.(sbs)])
        })
names(lbls) <- ctOrder %>% filter(!is.na(Subtype)) %>% pull(Subtype) 
lbls$superNeg <- "superNeg"
lbls$UNKNOWN <- "Unknown"


ctTheme <-      theme_minimal() +
     theme(text = element_text(size = 12, color = "black"),
           panel.grid.major.x = element_blank(),
           axis.line.x = element_line(size = 0.5, color = "black"),
           axis.line.y = element_line(size = 0.5, color = "black"),
           axis.text.x = element_text(color = "black", angle = 90),
           axis.title.x = element_blank(),
           legend.title = element_blank(), 
           legend.key.height = unit(0.17, "in")) 


maxY <- pdat %>% group_by(Patient_ID) %>% summarize(tot = sum(`Cell Count`)) %>% pull(tot) %>% max()/1000000
p <- ggplot(pdat, aes(x = Patient_ID, y = `Cell Count`/1000000, fill = Subtype)) +
     geom_bar(position = "stack", stat = "identity") +
     scale_y_continuous(limits = c(0, round_any(maxY, 0.5, f = ceiling)), expand = c(0,0),
                        breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3.0), labels = c("0M", "0.5M", "1.0M", "1.5M", "2.0M", "2.5M", "3.0M")) +
     scale_fill_manual(values = clrs, breaks = c(ctOrder$Subtype,  "superNeg", "UNKNOWN"), labels = lbls) +
     facet_grid(~EBV_final, space = "free_x", scales = "free_x", switch = "both") +
     ylab("Cell Count") +
     ctTheme +
     guides(fill= guide_legend(ncol = 1))


immDat <- pdat %>% filter(!Subtype %in% c("superNeg", "UNKNOWN", "HRS"))
maxY <- immDat %>% group_by(Patient_ID) %>% summarize(tot = sum(`Cell Count`)) %>% pull(tot) %>% max()/100000

p2 <- ggplot(immDat, aes(x = Patient_ID, y = `Cell Count`/100000, fill = Subtype)) +
     geom_bar(position = "stack", stat = "identity") +
     scale_y_continuous(limits = c(0, round_any(maxY, 0.5, f = ceiling)), expand = c(0,0),
                        breaks = c(0, 2.5, 5, 7.5, 10, 12.5), labels = c("0K", "250K", "500K", "750K", "1.00M", "1.25M")) + 
     scale_fill_manual(values = clrs, breaks = c(ctOrder$Subtype,  "superNeg", "UNKNOWN"), labels = lbls) +
     ylab("Cell Count") +
     ctTheme +
     facet_grid(~EBV_final, space = "free_x", scales = "free_x", switch = "both") +
     guides(fill= guide_legend(ncol = 1))


catDat <- pdat %>% 
          mutate(Category = factor(ifelse(Cell_type %in% c("UNKNOWN", "superNeg"), "superNeg+Unknown", 
                                     ifelse(Cell_type == "HRS", "HRS", "Immune")), 
                                   levels = rev(c('HRS', 'Immune', 'superNeg+Unknown')))) %>%
          group_by(Patient_ID, Category, EBV_final) %>%
          summarize(`Cell Count` = sum(`Cell Count`, na.rm = T))

maxY <- catDat %>% group_by(Patient_ID) %>% summarize(tot = sum(`Cell Count`)) %>% pull(tot) %>% max()
catClrs <- c(clrs['HRS'], c('Immune' = '#1261A0', 'superNeg+Unknown' = 'lightgray'))

catP <- ggplot(catDat %>% mutate(`Cell Count` = `Cell Count` + 1), aes(x = Patient_ID, y = `Cell Count`, fill = Category)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_fill_manual(values = rev(catClrs), breaks = names(catClrs)) +
        scale_y_continuous(expand = c(0, 0), breaks = c(250000, 1e+06, 2e+06, 3e+06), labels = c("250K", "1M", "2M", "3M")) +
        ylab("Cell Count") +
        ctTheme +
        facet_grid(~EBV_final, space = "free_x", scales = "free_x", switch = "both") +
        theme(#legend.position = "top",
              #legend.direction = "horizontal",
              plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "in")) +
        coord_trans(y = "sqrt") 

catPctDat <- catDat %>% group_by(Patient_ID, EBV_final) %>%
              mutate(`Cell Percentage` = `Cell Count`/sum(`Cell Count`, na.rm=T),
                     labelY = ifelse(Category == "HRS", `Cell Percentage` + 0.02, ifelse(Category == 'Immune', `Cell Percentage` - 0.025, 0.97)),
                     label = paste0(format(round(`Cell Percentage` * 100, 1), nsmall = 1), "%"),
                     labelColor = ifelse(Category == "superNeg+Unknown", "black", "white"))
                         

catPct <- ggplot(catPctDat, aes(x = Patient_ID, y = `Cell Percentage`, fill = Category)) +
          geom_bar(position = "stack", stat = "identity") +
          scale_y_continuous(limits = c(0, 1.01), expand = c(0,0), label = scales::percent) +
          scale_fill_manual(values = catClrs, breaks = names(catClrs), labels = lbls) +
          ctTheme  +
         facet_grid(~EBV_final, space = "free_x", scales = "free_x", switch = "both") +
          theme(#legend.position = "top",
                #legend.direction = "horizontal",
                plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "in")) +
          coord_trans(y = "sqrt")


catPctAlt <- catPct +
          xlab("") +
          coord_flip() +
          theme(legend.position = "top",
                legend.direction = "horizontal",
                plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "in")) +
          facet_grid(EBV_final ~ ., space = "free_y", scales = "free_y", switch = "both") +
          geom_text(aes(x = Patient_ID, y = labelY, label = label, color = labelColor), size = 3) +
          scale_color_manual(values = c('white' = 'white', 'black' = 'black'), guide = FALSE)


#pctDat <- dat %>% 
#        rename(Cell_type = CellType, Subtype = SubType) %>%
#        select(Cell_type, Subtype, dplyr::matches("\\.N\\.PCT$")) %>% 
#        gather(3:ncol(.), key = CellDive_ID, value = `Cell Percentage`) %>%
#        mutate(CellDive_ID = gsub("\\.N.*", "", CellDive_ID)) %>%
#        left_join(ids) %>% select(-CellDive_ID) %>%
#        left_join(ctOrder) %>%
#        mutate(Patient_ID = factor(Patient_ID, levels = sort(unique(Patient_ID))),
#               Subtype = factor(Subtype, levels = rev(c(ctOrder$Subtype, "superNeg", "UNKNOWN"))))

pctDat <- pdat %>%
          group_by(Patient_ID, EBV_final) %>%
          mutate(total = sum(`Cell Count`, na.rm = T), `Cell Percentage` = `Cell Count`/total) 

pct <- ggplot(pctDat, aes(x = Patient_ID, y = `Cell Percentage`, fill = Subtype)) +
     geom_bar(position = "stack", stat = "identity") +
     scale_y_continuous(limits = c(0, 1.01), expand = c(0,0), label = scales::percent) +
     scale_fill_manual(values = clrs, breaks = c(ctOrder$Subtype,  "superNeg", "UNKNOWN"), labels = lbls) +
     facet_grid(~EBV_final, space = "free_x", scales = "free_x", switch = "both") +
     ctTheme  +
     guides(fill= guide_legend(ncol = 1))


immPct <- immDat %>% 
          filter(!Subtype %in% c("superNeg", "UNKNOWN")) %>%
          group_by(Patient_ID) %>%
          mutate(`Cell Percentage` = `Cell Count`/sum(`Cell Count`, na.rm=T))          
 
pct2 <- ggplot(immPct, aes(x = Patient_ID, y = `Cell Percentage`, fill = Subtype)) +
     geom_bar(position = "stack", stat = "identity") +
     scale_y_continuous(limits = c(0, 1.01), expand = c(0,0), label = scales::percent) +
     scale_fill_manual(values = clrs, breaks = c(ctOrder$Subtype,  "superNeg", "UNKNOWN"), labels = lbls) +
     facet_grid(~EBV_final, space = "free_x", scales = "free_x", switch = "both") +
     ctTheme  +
     guides(fill= guide_legend(ncol = 1))






pdf("cell_type_counts.pdf", height = 4, width = 8)
print(p)
print(pct)
print(p2)
print(pct2)
print(catP)
print(catPct)
dev.off()

pdf("ctCategory_counts_alt.pdf", height = 8, width = 8)
print(catPctAlt)
dev.off()

