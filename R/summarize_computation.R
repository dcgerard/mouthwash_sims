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
  gather(key = "Method", value = "Time", OLS:BACKWASH) %>%
  group_by(n, Method) %>%
  summarize(Median = median(Time), Lower = quantile(Time, 0.025), Upper = quantile(Time, 0.975)) %>%
  ungroup() %>%
  mutate(Time = paste0(round(Median, digits = 2), " (",
                       round(Lower, digits = 2), ", ", round(Upper, digits = 2), ")")) ->
  sumdat

ggplot(data = sumdat, mapping = aes(x = n, y = Median)) +
  facet_wrap(~ Method, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_segment(mapping = aes(x = n, xend = n, y = Lower, yend = Upper), lty = 2)

sumdat %>% filter(n == 20) %>%
  arrange(Median) %>%
  select(Method, Time) ->
  tabdat

tabchar <- print(xtable(tabdat), include.rownames = FALSE, file = "./Output/computation/comp_tab.txt")
