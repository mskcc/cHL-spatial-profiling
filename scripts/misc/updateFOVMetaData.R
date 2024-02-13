# “rules” for updating the Hodgkin meta data:
#
# - Use:
#       HodgkinLymphoma_FOVs__V0.xlsx
#
# - For FOVs that contain <0.1% HRS cell (over DAPI), place an ‘X’ in the
#   ‘FOV_exclusion’ column AND write “post-halo” in the ‘FOV_exclusion_stage’
#   column

get_category_PCT_by_fov <- function(rdaFile) {

    obj=readRDS(rdaFile)

    allMarkers <- obj$marker.data %>% distinct(Marker) %>% pull
    m_dat=obj$marker.data %>%
            left_join(obj$geom.data %>% select(UUID,Sample,SPOT,Exclude) %>% mutate(FOV=SPOT)) %>%
            filter(!Exclude) %>%
            select(-Exclude)

    if(nrow(m_dat)==0) {
        cat("\n   No non-excluded cells\n   Exiting\n\n")
        return(NULL)
    }

    ct=assign_cell_types_by_fov(
                        m_dat,
                        cellTypeTbl,
                        intersect(cellTypeMarkers,allMarkers),
                        flexible=T
                    )

    ct %>%
        group_by(Sample,FOV) %>%
        count(Category) %>%
        mutate(PCT=n/sum(n)) %>%
        select(-n) %>%
        spread(Category,PCT,fill=0) %>%
        ungroup

}

##############################################################################
require(tidyverse)
require(readxl)
require(openxlsx)
##############################################################################
source("HaloX/tools.R")
source("HaloX/HaloCore/R/assign_cell_types.R")
##############################################################################
if(exists("SOURCED") && SOURCED) halt("INCLUDE")
##############################################################################


if(interactive()) {
    PARALLEL=TRUE
} else {
    PARALLEL=TRUE
}

if(PARALLEL) {
     library(parallel)
     nCores=24
     cat("Running in PARALLEL nCores=",nCores,"\n")
} else {
     cat("Running in SERIAL\n")
}

rdaFiles=fs::dir_ls(regex="_Reassign_01_.rda$",recur=T)

projParams=read_yaml("config/study.yaml")
cellTypeTbl = read_metadata_as_tibble(projParams,"CELLTYPES")
cellTypeMarkers = read_metadata_as_tibble(projParams,"MARKERS") %>%
    filter(Identity=="Y") %>%
    pull(Marker_name)

if(PARALLEL) {

    aa=mclapply(rdaFiles,get_category_PCT_by_fov,mc.cores = nCores)

} else {

    aa=lapply(rdaFiles,get_category_PCT_by_fov)

}

tbl1=bind_rows(aa) %>% mutate_all(~replace(., is.na(.), 0)) %>% arrange(HRS)
tbl2=tbl1 %>% filter(HRS<0.001)

FOVMETADATA="data/meta/raw/221024/HodgkinLymphoma_FOVs__V0.xlsx"

fovs=read_xlsx(FOVMETADATA)

extbl=tbl2 %>%
    select(CellDive_ID=Sample,FOV_number=FOV,PCT.HRS=HRS) %>%
    mutate(FOV_number=as.character(FOV_number))

fovsNew=fovs %>%
    left_join(extbl) %>%
    mutate(FOV_exclusion=ifelse(is.na(PCT.HRS),FOV_exclusion,"X")) %>%
    mutate(FOV_exclusion_stage=ifelse(is.na(PCT.HRS),FOV_exclusion_stage,"post-halo")) %>%
    mutate(Comment.RA01=case_when(
                        is.na(PCT.HRS) ~ "",
                        T ~ paste0("PCT.HRS:",prettyNum(round(100*PCT.HRS,2)),"%"))) %>%
    select(-PCT.HRS)

countFovs.orig=fovs %>% group_by(CellDive_ID) %>% count(CellDive_ID,name="NoExc")
countFovs.PreHRS=fovs %>% filter(is.na(FOV_exclusion)) %>% group_by(CellDive_ID) %>% count(CellDive_ID,name="PreHRSExc")
countFovs.PostHRS=fovsNew %>% filter(is.na(FOV_exclusion)) %>% group_by(CellDive_ID) %>% count(CellDive_ID,name="PostHRS")
fovCounts=full_join(countFovs.orig,countFovs.preHRS) %>% full_join(countFovs.PostHRS) %>% mutate_all(~replace(., is.na(.), 0)) %>% arrange(PostHRS)


write.xlsx(fovsNew,"HodgkinLymphoma_FOVs__V1.xlsx")
write.xlsx(list(FOVCounts=fovCounts,Exclude=tbl2,Full=tbl1),"rptFOVExclusionsUpdateHRSCells.xlsx")

SOURCED=T
