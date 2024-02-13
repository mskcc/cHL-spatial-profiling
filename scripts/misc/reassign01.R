####################################################################################
args=commandArgs(trailing=T)

if(len(grep("=",args,invert=T))<1) {
    cat("\n")
    cat("   usage: [DATAROOT] STUDY.yaml RDAFILE\n")
    cat("\n")
    cat("      - RDAFILE                 = HALO OBJ RDAFILE\n")
    cat("\n   Optionals:\n")
    cat("      - DATAROOT [default=\"\"]   = Root of data file folder\n")
    cat("\n")
    quit()
}
names(args)=args

suppressPackageStartupMessages(require(stringr))

opt.args=list(DATAROOT="")

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
obj=parse_AnalysisInput(obj,opt.args$DATAROOT)

thetas=obj$thetas %>%
    select(Marker,matches("threshold_weak")) %>%
    filter(Marker %in% c("CD30","MUM1")) %>%
    gather(threshold,theta,-Marker) %>%
    filter(
        Marker=="CD30" & threshold=="dye_cyto_positive_threshold_weak" |
        Marker=="MUM1" & threshold=="dye_nuclei_positive_threshold_weak"
    ) %>%
    select(-threshold)

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

tbl=obj$geom.data %>%
    filter(!Exclude) %>%
    left_join(ct,by=c("UUID","Sample")) %>%
    filter(Category=="UNKNOWN") %>%
    select(UUID,Category)

tbl2=obj$marker.data %>%
    filter(Marker %in% c("CD30","MUM1")) %>%
    select(UUID,Marker,TIntensity,Positive_Classification) %>%
    left_join(thetas) %>%
    inner_join(tbl) %>%
    mutate(New.Pos=ifelse(TIntensity>theta,1,0)) %>%
    select(UUID,Marker,New.Pos)

marker_new=obj$marker.data %>%
    left_join(tbl2) %>%
    mutate(Positive_Classification.orig=Positive_Classification) %>%
    mutate(Positive_Classification=ifelse(is.na(New.Pos),Positive_Classification.orig,New.Pos)) %>%
    select(-New.Pos)


m_dat=marker_new %>%
        left_join(obj$geom.data %>% select(UUID,Sample,SPOT,Exclude) %>% mutate(FOV=SPOT)) %>%
        filter(!Exclude) %>%
        select(-Exclude)

ct_new=assign_cell_types_by_fov(
                    m_dat,
                    cellTypeTbl,
                    intersect(cellTypeMarkers,allMarkers),
                    flexible=T
                )

stats=full_join(ct,ct_new,by="UUID") %>%
    select(UUID,Sample=Sample.x,FOV=FOV.x,matches("Category")) %>%
    group_by(Sample,FOV) %>%
    count(Category.x, Category.y) %>%
    unite(CategoryChange,Category.x,Category.y) %>%
    group_by(Sample,FOV) %>%
    mutate(Total=sum(n)) %>%
    filter(grepl("HRS|UNKNOWN",CategoryChange)) %>%
    ungroup

write_csv(stats,cc("stats_Reassign01",obj$sample.data$CellDive_ID,".csv"))

obj$marker.data=marker_new

obj$metadata$reassign=list(Set01="CD30,MUM1")

oFile=gsub(".rda","Reassign_01_.rda",basename(rdaFile))
saveRDS(obj,oFile,compress=T)

SOURCED=T
