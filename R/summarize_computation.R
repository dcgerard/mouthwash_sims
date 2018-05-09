## Summarize computation time results
suppressMessages(library(tidyverse))
library(stringr)
library(xtable)
compdat <- as_data_frame(readRDS(file = "./Output/computation/comp_time.Rds"))

str_replace(names(compdat), "_time", "") %>%
  str_replace("ols", "OLS") %>%
  str_replace("ruv2", "RUV2") %>%
  str_replace("ruv3", "RUV3") %>%
  str_replace("catenc", "CATEnc") %>%
  str_replace("sva", "SVA") %>%
  str_replace("caterr", "CATErr") %>%
  str_replace("mouthwash", "MOUTHWASH") %>%
  str_replace("backwash", "BACKWASH") %>%
  str_replace("Nsamp", "n") %>%
  str_replace("seed", "Seed") ->
  names(compdat)

compdat %>%
  select(-contains("BACKWASH")) %>%
  gather(key = "Method", value = "Time", OLS:MOUTHWASH) %>%
  group_by(n, Method) %>%
  summarize(Median = median(Time), Lower = quantile(Time, 0.025), Upper = quantile(Time, 0.975)) %>%
  ungroup() %>%
  mutate(Time = paste0(round(Median, digits = 2), " (",
                       round(Lower, digits = 2), ", ", round(Upper, digits = 2), ")")) ->
  sumdat

sumdat %>%
  arrange(Median) %>%
  select(Method, Time) ->
  tabdat


print(file = "./Output/computation/comp_tab.txt",
      xtable(tabdat), include.rownames = FALSE)
