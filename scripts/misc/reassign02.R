####################################################################################
args=commandArgs(trailing=T)

if(len(grep("=",args,invert=T))<1) {
    cat("\n")
    cat("   usage: reassign02.R STUDY.yaml RDAFILE\n")
    cat("\n")
    cat("      - STUDY.yaml              = Study config file\n")
    cat("      - RDAFILE                 = HALO OBJ RDAFILE\n")
    cat("\n")
    quit()
}
names(args)=args

suppressPackageStartupMessages(require(stringr))

ii=grep("=",args)
if(len(ii)>0) {
    parseArgs=str_match(args[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){opt.args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

if(len(ii)>0) {
    args=args[-ii]
}

####################################################################################
source("HaloX/parseAI.R")
source("HaloX/tools.R")
source("HaloX/HaloCore/utils.R")
source("HaloX/HaloCore/R/assign_cell_types.R")
source("HaloX/VERSION.R")
####################################################################################

require(tidyverse, quietly = T, warn.conflicts=F)

####################################################################################
if(exists("SOURCED") && SOURCED) halt(".INCLUDE")
####################################################################################

projParams=read_yaml(args[1])
rdaFile=args[2]
obj=readRDS(rdaFile)

#
# Get cell types
#

cellTypeTbl = read_metadata_as_tibble(projParams,"CELLTYPES")

cellTypeMarkers = read_metadata_as_tibble(projParams,"MARKERS") %>%
    filter(Identity=="Y") %>%
    pull(Marker_name)

allMarkers <- obj$marker.data %>% distinct(Marker) %>% pull

m_dat=obj$marker.data %>%
        left_join(obj$geom.data %>% select(UUID,Sample,SPOT,Exclude) %>% mutate(FOV=SPOT)) %>%
        filter(!Exclude) %>%
        select(-Exclude)

if(nrow(m_dat)==0) {
    cat("\n   No non-excluded cells\n   Exiting\n\n")
    quit(status=1)
}

ct=assign_cell_types_by_fov(
                    m_dat,
                    cellTypeTbl,
                    intersect(cellTypeMarkers,allMarkers),
                    flexible=T
                )

#
#
#

cells_unknown=obj$geom.data %>%
    filter(!Exclude) %>%
    left_join(ct,by=c("UUID","Sample")) %>%
    filter(Category=="UNKNOWN") %>%
    pull(UUID)

cells_dblPos=obj$marker.data %>%
    filter(Marker %in% c("CD3","CD20") & Positive_Classification==1) %>%
    group_by(UUID) %>% filter(n()==2) %>%
    pull(UUID)

rr=getRanksPerMarkerByFOV(obj)
rankSelect=rr %>%
    filter(UUID %in% cells_unknown) %>%
    filter(UUID %in% cells_dblPos) %>%
    filter(Marker %in% c("CD3","CD20")) %>%
    arrange(desc(Pr)) %>%
    distinct(UUID,.keep_all=T)

cells_CD20=rankSelect %>% filter(Marker=="CD20") %>% pull(UUID)
cells_CD3=rankSelect %>% filter(Marker=="CD3") %>% pull(UUID)

marker_new=obj$marker.data %>%
    mutate(Positive_Classification.orig=Positive_Classification) %>%
    mutate(Positive_Classification=case_when(
        UUID %in% cells_CD3 & Marker=="CD20" ~ 0,
        UUID %in% cells_CD20 & Marker %in% c("CD3","CD4","CD8","FOXP3") ~ 0,
        T ~ Positive_Classification.orig)
)


#
#
#

m_dat_new=marker_new %>%
    left_join(obj$geom.data %>% select(UUID,Sample,SPOT,Exclude) %>% mutate(FOV=SPOT)) %>%
    filter(!Exclude) %>%
    select(-Exclude)


ct_new=assign_cell_types_by_fov(
                    m_dat_new,
                    cellTypeTbl,
                    intersect(cellTypeMarkers,allMarkers),
                    flexible=T
                )

stats=full_join(ct,ct_new,by="UUID") %>%
    select(UUID,Sample=Sample.x,FOV=FOV.x,matches("Cell_type")) %>%
    group_by(Sample,FOV) %>%
    count(Cell_type.x, Cell_type.y) %>%
    unite(Change,Cell_type.x,Cell_type.y,sep=":") %>%
    group_by(Sample,FOV) %>%
    mutate(Total=sum(n)) %>%
    filter(grepl("UNKNOWN",Change)) %>%
    ungroup

before=full_join(ct,ct_new,by="UUID") %>%
    select(UUID,Sample=Sample.x,FOV=FOV.x,matches("Cell_type")) %>%
    group_by(Sample,FOV) %>%
    count(Cell_type.x) %>%
    rename(Cell_type=Cell_type.x,Before=n)
after=full_join(ct,ct_new,by="UUID") %>%
    select(UUID,Sample=Sample.x,FOV=FOV.x,matches("Cell_type")) %>%
    group_by(Sample,FOV) %>%
    count(Cell_type.y) %>%
    rename(Cell_type=Cell_type.y,After=n)

stats2=full_join(before,after) %>%
    ungroup %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate(`R.A/B`=(After+1)/(Before+1)) %>%
    filter(Cell_type %in% c("T_All","B_All")) %>%
    gather(Metric,Value,Before,After,`R.A/B`) %>%
    unite(CM,Cell_type,Metric,sep=".") %>%
    spread(CM,Value)


write_csv(stats,cc("stats_Reassign02_A",obj$sample.data$CellDive_ID,".csv"))
write_csv(stats2,cc("stats_Reassign02_B",obj$sample.data$CellDive_ID,".csv"))

obj$marker.data=marker_new

obj$metadata$reassign=obj$metadata$reassign[["Set02"]]="CD3|CD20"

oFile=gsub("Reassign_01_.rda","Reassign_01,02_.rda",basename(rdaFile))
saveRDS(obj,oFile,compress=T)

SOURCED=T
