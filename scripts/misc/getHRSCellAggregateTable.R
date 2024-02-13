#
# Dump a table with the aggregation status of every cell
# Input is the neighborGraph_ RDA files
# Loop over all samples
#

get_cluster_table_by_fov<-function(gi) {

    graph=igraph::graph_from_edgelist(select(gi,C.UUID,N.UUID)%>%as.matrix)
    components=igraph::components(graph)

    cSize=tibble(Cluster=seq(components$csize),ConnCompSize=components$csize)

    components$membership %>%
        data.frame %>%
        rownames_to_column("UUID") %>%
        tibble %>%
        rename(Cluster=2) %>%
        left_join(cSize,by="Cluster") %>%
        mutate(AggregateGe20=ConnCompSize>=20) %>%
        mutate(SPOT=gi$SPOT[1])

}

get_cluster_table_by_sample<-function(nfile) {
    graphTbl=readRDS(nfile)
    cellDiveID=gsub(".*neighborGraph_","",nfile) %>% gsub("_HRS_.rda$","",.)

    map(graphTbl,get_cluster_table_by_fov) %>%
        bind_rows %>%
        mutate(CellDiveID=cellDiveID)

}

require(tidyverse)

if(interactive() && exists("SOURCED") && SOURCED) {
  stop(".INCLUDE")
}

args=commandArgs(trailing=T)
tbl=map(args,get_cluster_table_by_sample) %>% bind_rows

write_csv(tbl,cc("hrsCellAggregationTable",substr(digest::digest(args),1,8),".csv.gz"))

SOURCED=T
