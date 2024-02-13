get_cluster_table<-function(gi) {
    graph=igraph::graph_from_edgelist(select(gi,C.UUID,N.UUID)%>%as.matrix)
    components=igraph::components(graph)
    tibble(Sample=sampleId,SPOT=gi$SPOT[1],Cell_type="HRS",ClusterSize=components$csize)
}

if(interactive() && exists("SOURCED") && SOURCED) {
  stop(".INCLUDE")
}

require(tidyverse)

args=commandArgs(trailing=T)
graphTbls=readRDS(args[1])
sampleId=gsub("neighborGraph_","",basename(args[1])) %>% gsub("_HRS_.rda","",.)

clusterStats=map(graphTbls,get_cluster_table) %>% bind_rows %>% arrange(desc(ClusterSize))
#clusterStats=map(graphTbls,filter,DelaunayNeighbor) %>% map(get_cluster_table) %>% bind_rows %>% arrange(desc(ClusterSize))

ODIR=file.path("db","HRSAggregates")
dir.create(ODIR,recursive=T,showWarnings=F)
write_csv(clusterStats,file.path(ODIR,cc("stats","HRSAggregatges",sampleId,".csv")))


SOURCED=T
