suppressMessages(library(funr))                                                 
srcRoot <- get_script_path()                                              

#srcRoot <- "/home/byrne/halo/dev/hodgkins_dev"

#suppressMessages(library(assertthat))
suppressMessages(library(contoureR))
suppressMessages(library(cowplot))
suppressMessages(library(digest))
#suppressMessages(library(dplyr))
suppressMessages(library(effsize))
suppressMessages(library(egg))
suppressMessages(library(logger))
suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(kimisc))
suppressMessages(library(lemon))
suppressMessages(library(magrittr))
suppressMessages(library(parallel))
suppressMessages(library(pheatmap))
suppressMessages(library(plotrix))
suppressMessages(library(plyr))  ## load BEFORE tidyverse
suppressMessages(library(randtoolbox))
suppressMessages(library(raster))
suppressMessages(library(RColorBrewer))
suppressMessages(library(R.utils))
suppressMessages(library(reshape))
suppressMessages(library(rgeos))
suppressMessages(library(rJava))
suppressMessages(library(rlang)) ## load BEFORE tidyverse
suppressMessages(library(scales))
suppressMessages(library(SearchTrees))
suppressMessages(library(sp))
suppressMessages(library(stringr))
suppressMessages(library(tidyverse))
#suppressMessages(library(tidyr))
#suppressMessages(library(testthat))
suppressMessages(library(tools))
suppressMessages(library(xlsx))
suppressMessages(library(openxlsx))
suppressMessages(library(xlsxjars))
suppressMessages(library(XML))
suppressMessages(library(yaml))


sourceDir <- file.path(srcRoot, "R") 
for(f in dir(sourceDir, full.names = T, recursive = T, pattern = "*.R$")){      
    if(!grepl("test", f)){                                                      
        source(f)                                                               
    }                                                                           
}       
