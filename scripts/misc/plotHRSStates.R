require(tidyverse)
require(yaml)

gColors=yaml::read_yaml("Halo-Hodgkins/assets/graphics/global_plot_colors.yaml")
ptColors=gColors$Patient_ID_FINAL %>% unlist

ptMap=readxl::read_xlsx("heatMapSamplePatientMap.xlsx")
reMap=ptMap$Samples
names(reMap)=ptMap$Patient_ID_FINAL

sampColors=ptColors
names(sampColors)=reMap[names(sampColors)]

config=read_yaml("config/study.yaml")
samplesMetaFile=fs::dir_ls("data/meta",recur=T,regex=paste0(config$SAMPLES,".*Sheet1.yaml"))
clin=read_yaml(samplesMetaFile) %>% bind_rows

VERSION="v2.2"
#cbjRDAFolder="/juno/work/bic/byrne/celldive/Hodgkins/2022-12-08/processed/metrics/fovs"
#cbjRDAFolder="/juno/work/bic/byrne/celldive/Hodgkins/current/processed/metrics/fovs/WITH_ADDITIONAL_HRS_STATES"
cbjRDAFolder="/juno/work/bic/byrne/celldive/Hodgkins/current/processed/metrics/fovs/WITH_EVEN_MORE_HRS_STATES"

dc=fs::dir_ls(cbjRDAFolder,regex="population_counts_per_FOV_ID") %>% map(readRDS) %>% bind_rows

totalHRS=dc %>% filter(Population=="HRS") %>% pull(Count) %>% sum
totalHRSByFOV=dc %>% filter(Population=="HRS") %>% group_by(FOV_ID) %>% summarize(FTotal=sum(Count))

idMap=dc %>%
    distinct(FOV_ID) %>%
    mutate(PP=str_match(FOV_ID,"(.*)_(\\d\\d)$")) %>%
    mutate(PP=as.data.frame(PP)) %>%
    unnest_wider(PP) %>%
    rename(Sample=V2,FOV=V3) %>%
    select(-V1)

totalHRSBySample=dc %>% left_join(idMap) %>% filter(Population=="HRS") %>% group_by(Sample) %>% summarize(STotal=sum(Count))

pc=dc %>%
    filter(grepl("^HRS,",Population)) %>%
    group_by(Population) %>%
    summarize(PCount=sum(Count,na.rm=T)) %>%
    arrange(desc(PCount)) %>%
    mutate(PCT=PCount/totalHRS)
oo=pull(pc,Population)
pc=pc %>% mutate(Population=factor(Population,levels=oo))

pb0=pc %>%
    ggplot(aes(Population,PCT)) +
    theme_light() +
    geom_bar(stat="identity") +
    scale_y_continuous(labels = scales::percent_format(scale = 100),limits=c(0,1)) +
    ggtitle("Population Fraction of Total HRS Cells") +
    xlab("") + ylab("")

pb1=pb0 + scale_x_discrete(guide=guide_axis(angle=45))
pb2=pb0 + coord_flip()

write_csv(pc,cc("hrsStatesHeatmap06.barplot",VERSION,".csv"))

xx=dc %>%
    filter(grepl("^HRS,",Population)) %>%
    left_join(idMap) %>%
    group_by(Sample,Population) %>%
    summarize(N=sum(Count,na.rm=T)) %>%
    left_join(totalHRSBySample) %>%
    mutate(PCT=(N+1)/(STotal+2)) %>%
    mutate(lOdds=log10(PCT/(1-PCT)))

# ml=xx %>%
#     select(-N,-STotal,-PCT) %>%
#     spread(Sample,lOdds) %>%
#     column_to_rownames("Population")

mm=xx %>% select(Sample,Population,PCT) %>% left_join(ptMap,by=c(Sample="Samples")) %>% ungroup %>% select(Patient_ID_FINAL,Population,PCT) %>% spread(Patient_ID_FINAL,PCT) %>% column_to_rownames("Population") %>% as.matrix
require(ComplexHeatmap)

#mm1=mm[apply(mm,1,\(x){max(x)})>.1,]
mm1=mm
ro=levels(factor(rownames(mm1),levels=intersect(oo,rownames(mm1))))

annotePt=clin %>% filter(Patient_ID %in% colnames(mm1)) %>% select(Patient_ID,EBV_final) %>% arrange(EBV_final,Patient_ID)

mm1=mm1[,annotePt$Patient_ID]

ebvColors=gColors$EBV_final
names(ebvColors)=sort(unique(annotePt$EBV_final))

ha0=HeatmapAnnotation(
    Patient=annotePt$Patient_ID,
    EBVFinal=annotePt$EBV_final,
    col=list(Patient=ptColors,EBVFinal=ebvColors)
)


pg=Heatmap(mm1,col = circlize::colorRamp2(c(0, 1), c("white", "Darkblue")),row_order=ro,row_names_gp = gpar(fontsize = 9),top_annotation=ha0, na_col="grey90", cluster_columns=F)

mm1 %>%
    data.frame(check.names=F) %>%
    rownames_to_column("Population") %>%
    write_csv(cc("hrsStatesHeatmap06.PatientHeatmap",VERSION,".csv"))

xf=dc %>%
    filter(grepl("^HRS,",Population)) %>%
    group_by(FOV_ID,Population) %>%
    summarize(N=sum(Count,na.rm=T)) %>%
    left_join(totalHRSByFOV) %>%
    mutate(PCT=(N+1)/(FTotal+2)) %>%
    mutate(Samples=gsub("_\\d\\d$","",FOV_ID)) %>% 
    left_join(ptMap)


mf=xf %>% select(FOV_ID,Population,PCT) %>% spread(FOV_ID,PCT) %>% column_to_rownames("Population") %>% as.matrix

annote=tibble(FOV_ID=colnames(mf)) %>% mutate(Samples=gsub("_\\d\\d$","",FOV_ID)) %>% left_join(ptMap)
ca=annote %>% distinct(Patient_ID_FINAL,COLOR)

ha=HeatmapAnnotation(
    Patient=annote$Patient_ID_FINAL,
    col=list(Patient=ptColors))

pf0=Heatmap(mf,col = circlize::colorRamp2(c(0, 1), c("white", "Darkblue")),row_order=ro,row_names_gp = gpar(fontsize = 9), show_column_names = FALSE, top_annotation=ha, cluster_columns=F, na_col="grey90")
#pf1=Heatmap(mf,col = circlize::colorRamp2(c(0, 1), c("white", "Darkblue")),row_order=ro,row_names_gp = gpar(fontsize = 9), show_column_names = FALSE, top_annotation=ha, na_col="grey90")

xf %>% select(FOV_ID,Population,PCT) %>% spread(FOV_ID,PCT) %>%
    write_csv(cc("hrsStatesHeatmap06.fovHeatmap",VERSION,".csv"))

S=xf %>% 
    ungroup %>% 
    group_by(Population) %>% 
    mutate(PCT=(N)/(sum(N))) %>% 
    summarize(S=-sum(ifelse(PCT>0,PCT*log2(PCT),0))) %>% 
    arrange(desc(S))

ns=xf %>% distinct(FOV_ID) %>% nrow

S=S %>% mutate(Snorm=S/log2(ns))

write_csv(S,cc("hrsStatesHeatmap06.entropy",VERSION,".csv"))

S=left_join(S,pc)

ps=ggplot(S,aes(reorder(Population,-PCount),Snorm)) + theme_light() + geom_bar(stat="identity") + scale_x_discrete(guide=guide_axis(angle=45)) + xlab("") + ylab("Entropy (normalized)") + ggtitle("Entropy per HRS-state")

vpct=xf %>% 
    group_by(Patient_ID_FINAL,Population) %>% 
    filter(n()>1) %>% 
    summarize(V=var(PCT)) %>% 
    spread(Patient_ID_FINAL,V) %>% 
    data.frame %>% 
    column_to_rownames("Population") %>%
    as.matrix

maxV=quantile(vpct,.99,na.rm=T)
ha1=HeatmapAnnotation(Patient=colnames(vpct),col=list(Patient=ptColors))

ph=Heatmap(vpct,row_order=ro,row_names_gp = gpar(fontsize = 9), cluster_columns=F,col = circlize::colorRamp2(c(0, maxV), c("white", "Darkgreen")), top_annotation=ha1, na_col="grey90")


pct_pop_in_pt=xf %>% select(Population,Patient_ID_FINAL,N) %>% group_by(Patient_ID_FINAL,Population) %>% summarize(Count=sum(N)) %>% mutate(PCT=Count/sum(Count))
pct_pt_in_pop=xf %>% select(Population,Patient_ID_FINAL,N) %>% group_by(Population,Patient_ID_FINAL) %>% summarize(Count=sum(N)) %>% mutate(PCT=Count/sum(Count))


maxCov=xf %>%
    group_by(Patient_ID_FINAL,Population) %>%
    filter(n()>1) %>%
    summarize(COV=sd(PCT)/mean(PCT)) %>%
    group_by(Population) %>%
    summarize(MaxVar=max(COV)) %>%
    arrange(MaxVar)

covPerPt=xf %>%
     group_by(Patient_ID_FINAL,Population) %>%
     filter(n()>1) %>%
     summarize(COV=sd(PCT)/mean(PCT)) %>%
     group_by(Population)

openxlsx::write.xlsx(list(MaxCov=maxCov,CovPerPt=covPerPt),cc("intraPatientCoV_PCTHRSStates",VERSION,".xlsx"))

pfMv=Heatmap(mf,col = circlize::colorRamp2(c(0, 1), c("white", "Darkblue")),row_order=ro,row_names_gp = gpar(fontsize = 9), show_column_names = FALSE, top_annotation=ha, cluster_columns=F, na_col="grey90")

pdf(file=cc("hrsStatesHeatmapCoVOrder",VERSION,".pdf"),width=18.5,height=11)
print(pfMv)
dev.off()


pdf(file=cc("hrsStatesHeatmap06",VERSION,".pdf"),width=18.5,height=11)
print(pb1)
print(pb2)
print(pg)
print(pf0)
#print(pf1) #do not plot with col-clusters
print(ph)
print(ps)
dev.off()
