require(tidyverse)
require(openxlsx)

dA=fs::dir_ls("out/reassign/03",regex="stats.*Reassign03_A.*.csv") %>% map(read_csv) %>% bind_rows
dB=fs::dir_ls("out/reassign/03",regex="stats.*Reassign03_B.*.csv") %>% map(read_csv) %>% bind_rows

rpt=dB %>%
    left_join(dA %>% select(Sample,FOV,Total) %>% distinct)

write.xlsx(list(FOV=rpt),"rpt_Reassign03_A.xlsx")
