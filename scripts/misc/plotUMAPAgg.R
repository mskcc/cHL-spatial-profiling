suppressPackageStartupMessages({
    require(stringr)
    require(tidyverse)
    require(yaml)
    require(fs)
    require(ggpubr)
})

###############################################################################
source("HUtils/toolsPlotUMAP.R")
###############################################################################

cArgs=commandArgs(trailing=T)

#
# This code will parse command line args in the form of
#    KEY=VAL
# and sets
#    args[[KEY]]=VAL
#

# Set defaults first

args=list(DEBUG=FALSE)

config=read_yaml("config/study.yaml")

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

args$DEBUG=as.logical(args$DEBUG)

pArgs=grep("=",cArgs,invert=T,value=T)
atlasFile=pArgs[1]
umapObjFile=pArgs[2]

atlas=readRDS(atlasFile)
uuObj=readRDS(umapObjFile)

geom=readRDS(dir_ls(regex="cellGeom.*\\.rda"))

hrsAgg=read_csv("./hrsCellAggregationTable_ceef61e2_.csv.gz")

umap=data.frame(uuObj$mm) %>%
    rownames_to_column("UUID") %>%
    tibble %>%
    bind_cols(as_tibble(uuObj$uu$embedding)) %>%
    left_join(atlas) %>%
    left_join(hrsAgg)

umap.o=umap

cat("\n\n    Ncells =",nrow(umap),"\n\n")

samplesMetaFile=fs::dir_ls("data/meta",recur=T,regex=paste0(config$SAMPLES,".*Sheet1.yaml"))
clin=read_yaml(samplesMetaFile) %>% bind_rows

umap=umap %>% left_join(clin,by=c(Sample="CellDive_ID"))
umap=umap %>% left_join(geom,by="UUID")

colorPal=read_yaml("Halo-Hodgkins/assets/graphics/global_plot_colors.yaml")

pt.size=sqrt(1e4/(2*nrow(umap)))
nCells.title=ggtitle(paste(" N.cells =",formatC(nrow(umap.o),big.mark=",")))

#apal=RColorBrewer::brewer.pal(7,"OrRd")[c(3,7)]
apal=c("#8faad8", "#990000")
cat("apal =",apal,"\n")

umap=umap %>% mutate(Aggregate=case_when(AggregateGe20 ~ "Yes", !AggregateGe20 ~ "No", T ~ "No"))

pa1=umap %>% ggplot(aes(V1,V2,color=Aggregate)) + theme_umap() + geom_point(size=pt.size) + scale_color_manual(values=apal) + guides(color = guide_legend(override.aes = list(size=5)))

source("HaloX/plotTools.R")

uuid.md5=substr(digest::digest(sort(rownames(uuObj$mm))),1,8)

oTag=cc(
    uuObj$params$cellSet,
    uuObj$params$args$SAMPLES,
    "_",
    uuObj$params$args$paramTag,
    paste0("UUID-",uuid.md5)
    )

oFileFactory<-function(base) {
    ODIR=file.path("plots/UMAPS",oTag)
    dir.create(ODIR,recursive=T,showWarnings=F)
    file.path(ODIR,sub(",",paste0("__",oTag,"__"),base))
}

png(oFileFactory("umapPlotAggV2,%03d.png"),14,11,pointsize=12,res=400)
print(pa1)
dev.off()

