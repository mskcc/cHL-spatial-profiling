####################################################################################
args=commandArgs(trailing=T)

if(len(grep("=",args,invert=T))<1) {
    cat("\n")
    cat("   usage: reassign03.R STUDY.yaml RDAFILE\n")
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

obj=parse_AnalysisInput(obj,"")

thetas=obj$thetas %>%
    select(Marker,matches("threshold_weak")) %>%
    filter(Marker %in% c("CD4","CD8")) %>%
    gather(threshold,theta,-Marker) %>%
    filter(
        Marker=="CD4" & threshold=="dye_cyto_positive_threshold_weak" |
        Marker=="CD8" & threshold=="dye_cyto_positive_threshold_weak"
    ) %>%
    select(-threshold)


#
# Get cell types
#

getCellTypes<-function(markerData) {

    cellTypeTbl = read_metadata_as_tibble(projParams,"CELLTYPES")

    cellTypeMarkers = read_metadata_as_tibble(projParams,"MARKERS") %>%
        filter(Identity=="Y") %>%
        pull(Marker_name)

    allMarkers <- obj$marker.data %>% distinct(Marker) %>% pull

    m_dat=markerData %>%
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

}

compPCTSubtype<-function(ctbl) {

    ctcounts=ctbl %>%
        group_by(Sample,FOV) %>%
        count(Cell_type) %>%
        rename(CT.n=n)
    stcounts=ctbl %>%
        group_by(Sample,FOV) %>%
        count(Cell_type,Subtype) %>%
        mutate(TotalAll=sum(n)) %>%
        rename(ST.n=n)

    left_join(stcounts,ctcounts) %>%
        filter(Cell_type=="T_All") %>%
        mutate(PCT.ST.TotalA=ST.n/TotalAll,PCT.ST.T=ST.n/CT.n) %>%
        mutate(PCT=PCT.ST.T) %>%
        ungroup

}

ct=getCellTypes(obj$marker.data)
ct_orig=ct

#
# Step 3-A: Resolution of T cell/null phenotype population
#
# MUST BE DONE AFTER Step 2 and after type reassignment
#
# 1. Take T cell/null phenotype population (double neg cells CD4/CD8)
#

cells_Tnull=obj$geom.data %>%
    filter(!Exclude) %>%
    left_join(ct,by=c("UUID","Sample")) %>%
    filter(Subtype=="T_Null") %>% pull(UUID)

# 2. If PCT[CD4-,CD8-]>5%

pct_DN=compPCTSubtype(ct)

fovsGt5pctDN=pct_DN %>%
    filter(Subtype=="T_Null" & PCT>0.05) %>%
    pull(FOV)

if(len(fovsGt5pctDN)>0) {

    #
    # 2. then drop CD4,CD8 completeness threshold for cells that meet cytoplasmic intensity threshold
    #

    tblRA=obj$marker.data %>%
        filter(Marker %in% c("CD4","CD8")) %>%
        filter(UUID %in% cells_Tnull) %>%
        left_join(ct %>% select(UUID,FOV)) %>%
        filter(FOV %in% fovsGt5pctDN) %>%
        select(UUID,Marker,TIntensity,Positive_Classification) %>%
        left_join(thetas) %>%
        mutate(New.Pos=ifelse(TIntensity>theta,1,0)) %>%
        select(UUID,Marker,New.Pos)

    marker_new=obj$marker.data %>%
        left_join(tblRA) %>%
        mutate(Positive_Classification.orig=Positive_Classification) %>%
        mutate(Positive_Classification=ifelse(is.na(New.Pos),Positive_Classification.orig,New.Pos)) %>%
        select(-New.Pos)

    ct_new=getCellTypes(marker_new)
    obj$marker.data=marker_new

    cells_Tnull_new=obj$geom.data %>%
        filter(!Exclude) %>%
        left_join(ct_new,by=c("UUID","Sample")) %>%
        filter(Subtype=="T_Null") %>% pull(UUID)

    pct_DN_new=compPCTSubtype(ct_new)
    fovsGt5pctDN=pct_DN_new %>%
        filter(Subtype=="T_Null" & PCT>0.05) %>%
        pull(FOV)

    if(len(fovsGt5pctDN)>0) {
        #
        # If PCT[CD4-,CD8-]>5% then take remainder of cells (cells that both fail the CD4 and CD8 intensity thresholds)
        #   Rank reassignment for CD4/CD8 (5% T cell/null phenotype cells allowed)
        #       If CD4 “wins”: CD4 positivity switches to positive
        #           Cell becomes a CD4+ T cell
        #       If CD8 “wins”: CD8 positivity switches to positive
        #           Cell becomes a CD8+ T cell
        #
        rr=getRanksPerMarkerByFOV(obj)
        tblRA2=rr %>%
            filter(Marker %in% c("CD4","CD8")) %>%
            filter(UUID %in% cells_Tnull_new) %>%
            left_join(ct %>% select(UUID,FOV)) %>%
            filter(FOV %in% fovsGt5pctDN) %>%
            select(-R) %>%
            spread(Marker,Pr) %>%
            mutate(lOR=log((CD4/(1-CD4))/(CD8/(1-CD8)))) %>%
            arrange(desc(abs(lOR)))

        numCells=pct_DN_new %>%
            filter(Subtype=="T_Null" & PCT>0.05) %>%
            mutate(N.PCT.5=round(ST.n-CT.n*0.05)) %>%
            select(FOV,N.PCT.5)

        cellsToFlip=tblRA2 %>%
            arrange(desc(abs(lOR))) %>%
            left_join(numCells) %>%
            group_by(FOV) %>%
            mutate(R=rank(-abs(lOR))) %>%
            filter(R<=N.PCT.5) %>%
            mutate(Marker=ifelse(lOR>0,"CD4","CD8")) %>%
            select(UUID,lOR,R,Marker)

        marker_new=obj$marker.data %>%
            left_join(cellsToFlip) %>%
            mutate(Positive_Classification=ifelse(!is.na(lOR),1,Positive_Classification)) %>%
            select(-FOV,-lOR,-R)

        obj$marker.data=marker_new

    }



}

rpt0=pct_DN %>% select(Sample,FOV,Subtype,PCT.ST.T) %>% spread(Subtype,PCT.ST.T,fill=0) %>% arrange(desc(T_Null))

#
# Step 3-B: Resolution of CD4+/CD8+ T cells
#
# 0. MUST BE DONE AFTER Step 2, Step 3-A and after type reassignment
# 1. Take CD4+/CD8+ T cell population
# 2. Rank assignment for CD4/CD8 (5% CD4+/CD8+ T cells allowed)

ct=getCellTypes(obj$marker.data)
pct_ST=compPCTSubtype(ct)

rpt1=pct_ST %>% select(Sample,FOV,Subtype,PCT.ST.T) %>% spread(Subtype,PCT.ST.T,fill=0) %>% arrange(desc(T_Null))

fovsGt5pctDP=pct_ST %>%
    filter(Subtype=="T_CD4/CD8" & PCT>0.05) %>%
    pull(FOV)

if(len(fovsGt5pctDP)>0) {
    #
    # 1. Take CD4+/CD8+ T cell population
    #
    cells_Tdp=obj$geom.data %>%
        filter(!Exclude) %>%
        left_join(ct,by=c("UUID","Sample")) %>%
        filter(Subtype=="T_CD4/CD8") %>% pull(UUID)

    #
    #   i. If CD4 “wins”: CD8 positivity switches to negative
    #       - Cell becomes a CD4+ T cell or CD4+ regulatory T cell
    #   ii. If CD8 “wins”: CD4 positivity switches to negative
    #       - Cell becomes a CD8+ T cell or CD8+ regulatory T cell
    #

    rr=getRanksPerMarkerByFOV(obj)
    tblDP=rr %>%
        filter(Marker %in% c("CD4","CD8")) %>%
        filter(UUID %in% cells_Tdp) %>%
        left_join(ct %>% select(UUID,FOV)) %>%
        filter(FOV %in% fovsGt5pctDP) %>%
        select(-R) %>%
        spread(Marker,Pr) %>%
        mutate(lOR=log((CD4/(1-CD4))/(CD8/(1-CD8)))) %>%
        arrange(desc(abs(lOR)))

    numCells=pct_ST %>%
        filter(Subtype=="T_CD4/CD8" & PCT>0.05) %>%
        mutate(N.PCT.5=round(ST.n-CT.n*0.05)) %>%
        select(FOV,N.PCT.5)

    cellsToFlip=tblDP %>%
        arrange(desc(abs(lOR))) %>%
        left_join(numCells) %>%
        group_by(FOV) %>%
        mutate(R=rank(-abs(lOR))) %>%
        filter(R<=N.PCT.5) %>%
        mutate(Marker=ifelse(lOR<0,"CD4","CD8")) %>%
        select(UUID,lOR,R,Marker)

    marker_new=obj$marker.data %>%
        left_join(cellsToFlip) %>%
        mutate(Positive_Classification=ifelse(!is.na(lOR),0,Positive_Classification)) %>%
        select(-FOV,-lOR,-R)

    obj$marker.data=marker_new


}

if(exists("PCT_DN")) rm("pct_DN")
if(exists("PCT_DN_new")) rm("pct_DN_new")
if(exists("PCT_ST")) rm("pct_ST")

ct_final=getCellTypes(obj$marker.data)
pct_Final=compPCTSubtype(ct_final)
rpt2=pct_Final %>% select(Sample,FOV,Subtype,PCT.ST.T) %>% spread(Subtype,PCT.ST.T,fill=0) %>% arrange(desc(T_Null))

stats=full_join(ct_orig,ct_final,by="UUID") %>%
    filter(Cell_type.x=="T_All") %>%
    select(UUID,Sample=Sample.x,FOV=FOV.x,matches("Subtype")) %>%
    group_by(Sample,FOV) %>%
    count(Subtype.x, Subtype.y) %>%
    unite(Change,Subtype.x,Subtype.y,sep=":") %>%
    group_by(Sample,FOV) %>%
    mutate(Total=sum(n)) %>%
    ungroup %>%
    filter(grepl("T_Null:|T_CD4/CD8:",Change))

before=full_join(ct_orig,ct_final,by="UUID") %>%
    select(UUID,Sample=Sample.x,FOV=FOV.x,matches("Subtype")) %>%
    group_by(Sample,FOV) %>%
    count(Subtype.x) %>%
    rename(Subtype=Subtype.x,Before=n)

after=full_join(ct_orig,ct_final,by="UUID") %>%
    select(UUID,Sample=Sample.x,FOV=FOV.x,matches("Subtype")) %>%
    group_by(Sample,FOV) %>%
    count(Subtype.y) %>%
    rename(Subtype=Subtype.y,After=n)

stats2=full_join(before,after) %>%
    ungroup %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate(`R.A/B`=(After+1)/(Before+1)) %>%
    filter(Subtype %in% c("T_CD4","T_CD8")) %>%
    gather(Metric,Value,Before,After,`R.A/B`) %>%
    unite(CM,Subtype,Metric,sep=".") %>%
    spread(CM,Value)

write_csv(stats,cc("stats_Reassign03_A",obj$sample.data$CellDive_ID,".csv"))
write_csv(stats2,cc("stats_Reassign03_B",obj$sample.data$CellDive_ID,".csv"))

obj$metadata$reassign=obj$metadata$reassign[["Set03"]]="CD4/CD8 Dn,Dp"

oFile=gsub("Reassign_01,02_.rda","ReassignAll.rda",basename(rdaFile))
saveRDS(obj,oFile,compress=T)

SOURCED=T
