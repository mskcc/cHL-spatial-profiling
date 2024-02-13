require(tidyverse)

fovAreasRDAFolder="/juno/work/bic/byrne/celldive/Hodgkins/2022-12-08/processed/metrics/fovs"

fovAreas=fs::dir_ls(fovAreasRDAFolder,regex="area") %>% map(readRDS) %>% bind_rows

atlas=readRDS("cellAtlas_HodgkinsV10s_v10.9__211112_b_CTD_V0___0576fa92_20221031_231605.rda")

aa=atlas %>% filter(!Exclude)

counts=aa %>%
    group_by(Sample,SPOT) %>%
    count(Cell_type) %>%
    mutate(Cell_type=paste0("ct:",Cell_type)) %>%
    spread(Cell_type,n,fill=0) %>%
    mutate(Total=rowSums(across(matches("^ct")))) %>%
    rename_all(~gsub("^ct:","",.)) %>%
    ungroup %>%
    rename(CellDive_ID=Sample,FOV_number=SPOT)

counts=left_join(counts,fovAreas) %>% filter(!is.na(Area)) %>% mutate(N=HRS,f=HRS/Total,rho=HRS/Area)

agg=readxl::read_xlsx("stats_Aggregation_HRS_N2_AggSize__20_CLIN.xlsx") %>% rename(CellDive_ID=Sample,FOV_number=SPOT)

pg=list()

dd=left_join(counts,agg) %>% arrange(desc(HRSStatus))

minFracClustered=min(dd$f[dd$HRSStatus=="Clustered"])

pg[[len(pg)+1]]=dd %>% ggplot(aes(N,f,color=HRSStatus)) + theme_light() + geom_point(alpha=.75,size=3) + scale_x_log10() + scale_y_log10() + scale_color_manual(values=c("darkred","grey60")) + xlab("Number HRS Cells per FOV") + ylab("Fraction HRS (HRS/DAPI)") + ggtitle("Fraction HRS vs Number HRS per FOV")
pg[[len(pg)+1]]=dd %>% ggplot(aes(N,rho,color=HRSStatus)) + theme_light() + geom_point(alpha=.75,size=3) + scale_x_log10() + scale_y_log10() + scale_color_manual(values=c("darkred","grey60")) + xlab("Number HRS Cells per FOV") + ylab("Density HRS Cells") + ggtitle("Density HRS vs Number HRS per FOV")
pg[[len(pg)+1]]=dd %>% ggplot(aes(f,rho,color=HRSStatus)) + theme_light() + geom_point(alpha=.75,size=3) + scale_x_log10() + scale_y_log10() + scale_color_manual(values=c("darkred","grey60")) + xlab("Fraction HRS (HRS/DAPI)") + ylab("Density HRS Cells") + ggtitle("Density HRS vs Fraction HRS per FOV") + geom_vline(xintercept=minFracClustered,alpha=.4) + annotate(geom="text",x=minFracClustered,y=10,label=paste0(" ",round(minFracClustered,4)),hjust=0,vjust=0,size=5)


pg[[len(pg)+1]]=dd %>% ggplot(aes(HRSStatus,rho,color=HRSStatus)) + theme_light() + geom_violin() + geom_jitter() + scale_color_manual(values=c("darkred","grey60")) + ylab("Density HRS Cells") + ggtitle("Density HRS Cells vs HRS Status")


pg[[len(pg)+1]]=filter(dd,f>minFracClustered) %>% ggplot(aes(HRSStatus,rho,color=HRSStatus)) + theme_light() + geom_violin() + geom_jitter() + scale_color_manual(values=c("darkred","grey60")) + ylab("Density HRS Cells") + ggtitle("Density HRS Cells vs HRS Status // FOV Filter frac > min(frac-Clustered)")

pdf(file="phaseDiagramA.pdf",width=11,height=8.5)
print(pg)
dev.off()
