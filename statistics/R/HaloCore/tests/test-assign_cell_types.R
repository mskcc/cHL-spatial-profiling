context("assign cell types")

source("../R/assign_cell_types.R")

meta       <- read.xlsx('data/FOVs.xlsx', 1, check.names = F) %>% as_tibble
cell_types <- read.xlsx("data/good_cell_types.xlsx", 1, check.names = F) %>% 
              as_tibble
markers    <- read.xlsx('data/markers.xlsx', 1, check.names = F) %>% 
              as_tibble %>%
              filter(Identity == 'Y') %>%
              pull(Marker_name)
m_dat      <- readRDS("data/test_marker_data_rand.rda")


ct_strict <- assign_cell_types_by_fov(m_dat, cell_types, markers)
ct_flex   <- assign_cell_types_by_fov(m_dat, cell_types, markers, flexible = TRUE)


test_that("missing markers prevent cell type identification on FOV level only (not entire samples)", { 
    ## get the cell types that can be assigned when removal
    ## of negative requirements is NOT allowed
    olig2_na <- cell_types %>% filter(is.na(OLIG2)) 
    obs <- ct_strict %>% filter(Subtype %in% olig2_na$Subtype)
    expect_true(all(m_dat$FOV %in% obs$FOV))

    ## the only subtype assign-able in FOVs 4 and 7 is Tumor (and NA for UNKNOWN)
    # olig2_req <- cell_types %>% filter(!is.na(OLIG2))
    # Tumor only
    obs <- ct_strict %>% filter(FOV %in% c(4,7)) %>% pull(Subtype) %>% unique
    expect_equal(obs, c("Tumor","UNKNOWN"))

    ## all other subtypes (except T_Reg8, T_Reg4, T_RegDN & Neuron which all happen to have count zero, 
    ## confirmed manually) should be assigned in all other FOVs
    obs <- ct_strict %>% filter(!FOV %in% c(4,7)) %>% pull(Subtype) %>% unique
    exp <- cell_types %>% 
           filter(!Subtype %in% c("T_Reg4", "T_Reg8", "T_RegDN", "Neuron")) %>% 
           pull(Subtype) %>% 
           unique
    expect_equal(sort(obs), c(sort(exp), "UNKNOWN")) 
})


test_that("cell type definitions that do NOT include a missing marker are not affected", {
    ## counts should not differ between strict/flex for definitions
    ## that do not include the missing marker
    olig2_na <- cell_types %>% filter(is.na(OLIG2))
    counts_strict <- ct_strict %>% 
                     filter(Subtype %in% olig2_na$Subtype) %>% 
                     group_by(Subtype) %>% 
                     summarize(Count = n())
    counts_flex <- ct_flex %>% 
                   filter(Subtype %in% olig2_na$Subtype) %>% 
                   group_by(Subtype) %>% 
                   summarize(Count = n())
    expect_equal(counts_strict, counts_flex)
})

test_that("cell type definitions are NEVER assigned when missing marker is a required positive", {
    ## NPC is the only cell def with OLIG2 pos required
    obs <- ct_flex %>% filter(Cell_type == "NPC")
    expect_false(any(c(4,7) %in% obs$FOV))
    expect_true(all(c(5,8,9) %in% obs$FOV)) ## the other FOVs have count zero in this subset (confirmed manually) 
})


test_that("cell type definitions ARE modified when missing marker is only a required negative and 'flexible' is TRUE", {
    olig2_neg <- cell_types %>% filter(OLIG2 == 0)
    obs <- ct_flex %>% filter(FOV %in% c(4,7))
    expect_gt(length(unique(obs$Subtype)), 2)
    
    unk_flex <- ct_flex %>% filter(Cell_type == "UNKNOWN")
    expect_equal(nrow(unk_flex), 4885)
})

test_that("cell type definitions are NEVER modified when 'flexible' is FALSE", {
    olig2_neg <- cell_types %>% filter(OLIG2 == 0)

    ## olig2 neg cell types still should not be in ct_strict FOVs 4 & 7
    neg_strict <- ct_strict %>% filter(Subtype %in% olig2_neg$Subtype, FOV %in% c(4,7))
    expect_equal(nrow(neg_strict), 0)

    unk_strict <- ct_strict %>% filter(Cell_type == "UNKNOWN")
    expect_equal(nrow(unk_strict), 5124)   ## this should be more than when flexible is TRUE
})

test_that("cell type definitions are expanded completely", {

    cell_def <- tibble(m1 = 0, m2 = 1, m3 = 1)
    obs <- expand_cell_definition(cell_def)
    expect_identical(obs, cell_def)

    cell_def <- tibble(m1 = 0, m2 = 0, m3 = 'A_1', m4 = 'A_1')
    exp <- tibble(m1 = c(0,0), 
                  m2 = c(0,0),
                  m3 = c(1,NA),
                  m4 = c(NA,1))

    obs <- expand_cell_definition(cell_def)
    expect_identical(obs, exp)

    cell_def <- tibble(m1 = 0, m2 = 'B_2', m3 = 'B_2', m4 = 'B_2', m5 = 1)
    exp <- tibble(m1 = c(0,0,0),
                  m2 = c(1,1,NA),
                  m3 = c(1,NA,1),
                  m4 = c(NA,1,1),
                  m5 = c(1,1,1))
    obs <- expand_cell_definition(cell_def) %>% select(m1,m2,m3,m4,m5)
    expect_equal(obs, exp)  

    cell_def <- tibble(m1 = 'A_1', m2 = 'A_1',  m3 = 'B_2', m4 = 'B_2', m5 = 'B_2')
    exp <- tibble(m1 = c(1,1,1,NA,NA,NA),
                  m2 = c(NA,NA,NA,1,1,1),
                  m3 = c(1,1,NA,1,1,NA),
                  m4 = c(1,NA,1,1,NA,1),
                  m5 = c(NA,1,1,NA,1,1))
    obs <- expand_cell_definition(cell_def) %>% select(m1,m2,m3,m4,m5)
    expect_equal(obs,exp)

})

test_that("conflicting cell types are caught", {
    cell_types <- read.xlsx('data/conflicting_cell_type_definitions.xlsx', 1, check.names = F) %>% as_tibble
    expect_error(assign_cell_types_by_fov(m_dat, cell_types, markers)) 
})

test_that("data is filtered properly by cell type definitions", {

    mdat <- m_dat %>% 
            select(Sample, FOV, UUID, Marker, Positive_Classification) %>%
            spread(Marker, Positive_Classification)

    ## olig2 is the only pos
    cell_def <- cell_types %>% filter(Cell_type == "NPC") %>% select(markers)
    
    ## when a marker missing in an FOV, empty tibble returned
    cells <- mdat %>% 
             filter(FOV == 4) %>%
             filter_by_cell_type_definition(cell_def, flexible = F)
    expect_equal(nrow(cells), 0)

    ## no markers missing
    cells <- mdat %>%
             filter(FOV == 0) %>%
             filter_by_cell_type_definition(cell_def, flexible = F)
    expect_true(all(cells$OLIG2 == 1))

    neg <- cell_def %>% select_if(~all(. == 0)) %>% names
    dvnl <- lapply(neg, function(x){
        expect_true(all(cells[[x]] == 0))
    })

    ## now allow removal of neg req for missing marker(s) 
    cell_def <- cell_types %>% filter(Subtype == "Other_macro") %>% select(markers)
    neg <- cell_def %>% select_if(~all(. == 0)) %>% names

    cells <- mdat %>% 
             filter(FOV == 4) %>%
             filter_by_cell_type_definition(cell_def, flexible = T)
    expect_equal(nrow(cells), 4)
    expect_true(all(cells$CD14 == 1))
    expect_false("OLIG2" %in% names(cells))

    neg <- cell_def %>% select_if(~all(. == 0)) %>% names
    dvnl <- lapply(setdiff(neg, "OLIG2"), function(x){
        expect_true(all(cells[[x]] == 0))
    })    
})

test_that("error is thrown when marker data is missing columns", {

    mdat <- m_dat %>% 
            select(UUID, Marker, Positive_Classification) %>%
            spread(Marker, Positive_Classification)

    expect_error(assign_cell_types_by_fov(mdat, cell_types, markers))   

})

