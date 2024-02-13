require(yaml)
require(tidyverse)

paletteFile="global_plot_colors.yaml"

pal=read_yaml(paletteFile)

xx=tibble(X=1,Y=len(pal)-seq(len(pal)),Color=names(pal),Hex=unlist(pal))

xScale=0.01

pg2=ggplot(xx,aes(X,Y,fill=Color,label=Hex)) +
    scale_fill_manual(values=unlist(pal)) +
    xlim(1+xScale*c(-1,1)) +
    geom_point(size=8,shape=21) +
    geom_label(hjust=0,nudge_x=xScale/6,fill="#FEFEFE",family="mono") +
    ggtitle(paletteFile) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank()
        )

#print(pg2)

pdf(file="testPalette.pdf",width=11,height=8.5)
print(pg2)
dev.off()

