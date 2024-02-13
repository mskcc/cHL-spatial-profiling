require(tidyverse)
require(openxlsx)

dd=fs::dir_ls("out/reassign/01",regex="stats.*Reassign01_.*.csv") %>% map(read_csv) %>% bind_rows

rpt=dd %>%
    spread(CategoryChange,n,fill=0) %>%
    rename(HRS.Orig=HRS_HRS) %>%
    mutate(HRS.New=HRS.Orig+UNKNOWN_HRS) %>%
    mutate(PCT.HRS.Orig=HRS.Orig/Total,PCT.HRS.New=HRS.New/Total) %>%
    mutate(`R.New/Orig`=(HRS.New+1)/(HRS.Orig+1)) %>%
    select(-matches("UNKNOWN")) %>%
    arrange(desc(`R.New/Orig`))

tbl2=tibble(LEVEL="FOV",
    Num.Pct_0=sum(rpt$PCT.HRS.New>0),Pct.Pct_0=mean(rpt$PCT.HRS.New>0),
    Num.Pct_0.1=sum(rpt$PCT.HRS.New>0.001),Pct.Pct_0.1=mean(rpt$PCT.HRS.New>0.001),
    Num.Pct_0.2=sum(rpt$PCT.HRS.New>0.002),Pct.Pct_0.2=mean(rpt$PCT.HRS.New>0.002),
    Num.Pct_0.5=sum(rpt$PCT.HRS.New>0.005),Pct.Pct_0.5=mean(rpt$PCT.HRS.New>0.005),
    Num.Pct_1.0=sum(rpt$PCT.HRS.New>0.01),Pct.Pct_1.0=mean(rpt$PCT.HRS.New>0.010),
)

rpt1=rpt %>%
    group_by(Sample) %>%
    summarize(Total=sum(Total),HRS.Orig=sum(HRS.Orig),HRS.New=sum(HRS.New)) %>%
    mutate(PCT.HRS.Orig=HRS.Orig/Total,PCT.HRS.New=HRS.New/Total) %>%
    mutate(`R.New/Orig`=(HRS.New+1)/(HRS.Orig+1)) %>%
    arrange(desc(`R.New/Orig`))

tbl2s=tibble(LEVEL="SAMPLE",
    Num.Pct_0=sum(rpt1$PCT.HRS.New>0),Pct.Pct_0=mean(rpt1$PCT.HRS.New>0),
    Num.Pct_0.1=sum(rpt1$PCT.HRS.New>0.001),Pct.Pct_0.1=mean(rpt1$PCT.HRS.New>0.001),
    Num.Pct_0.2=sum(rpt1$PCT.HRS.New>0.002),Pct.Pct_0.2=mean(rpt1$PCT.HRS.New>0.002),
    Num.Pct_0.5=sum(rpt1$PCT.HRS.New>0.005),Pct.Pct_0.5=mean(rpt1$PCT.HRS.New>0.005),
    Num.Pct_1.0=sum(rpt1$PCT.HRS.New>0.01),Pct.Pct_1.0=mean(rpt1$PCT.HRS.New>0.010),
)

write.xlsx(list(Summary=bind_rows(tbl2,tbl2s),FOV=rpt,Sample=rpt1),"rpt_Reassign01_C.xlsx")
