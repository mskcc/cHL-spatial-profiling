require(tidyverse)
require(openxlsx)

dA=fs::dir_ls("out/reassign/02",regex="stats.*Reassign02_A.*.csv") %>% map(read_csv) %>% bind_rows
dB=fs::dir_ls("out/reassign/02",regex="stats.*Reassign02_B.*.csv") %>% map(read_csv) %>% bind_rows

rpt=dB %>%
    left_join(dA %>% select(Sample,FOV,Total) %>% distinct) %>%
    mutate(PCT.B_All.After=B_All.After/Total,PCT.T_All.After=T_All.After/Total)

write.xlsx(list(FOV=rpt),"rpt_Reassign02_A.xlsx")
