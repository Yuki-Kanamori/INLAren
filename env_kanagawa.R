require(tidyverse)
require(abind)
require(openxlsx)


# メモ ------------------------------------------------------------
# temp = env %>% select(temp, tag)
# sal = env %>% select(salinity, tag)
# do = env %>% select(DO_mg, tag)
# ph = env %>% select(pH, tag)

# st27, 134, 136 ------------------------------------------------
dir = "/Users/Yuki/Dropbox/eDNA_INLA/env_kanagawa/"
st = c("st27", "st134", "st136")
t2 = NULL
t3 = NULL

for(i in 1:length(st)){
  path = paste0(dir, st[i])
  setwd(path)
  files = list.files()
  
  for(j in 1:length(files)){
    t = read.table(files[j], header = T) %>% mutate(file = paste0(files[j]))
    min = min(t$DepSM)
    max = max(t$DepSM)
    s = t %>% filter(DepSM == min) %>% mutate(SorB = "S")
    b = t %>% filter(DepSM == max) %>% mutate(SorB = "B")
    sb = rbind(s, b)
    t2 = abind(t2, sb, along = 1)
  }
  #t3 = abind(t3, t2, along = 1)
}
envk1 = data.frame(t2) %>% mutate(year = 2018, month = as.numeric(str_sub(file, 5, 6)), day = NA, Point = as.numeric(str_sub(file, 8, 10)))
(colnames(envk1))
envk1 = envk1 %>% select("T090C", "Ph", "Sal00", "Sbeox0Mg.L","SorB", "year", "month", "day", "Point")
colnames(envk1) = c("temp", "pH", "salinity", "DO_mg","SorB", "year", "month", "day", "Point")

# st129 ---------------------------------------------------------
st2 = "st129"
path2 = paste0(dir, st2)
setwd(path2)
(files2 = list.files())
t2 = NULL

for(i in 1:length(files2)){
  t = read.xlsx(files2[i], sheet = "129", startRow = 44, colNames = T) %>% dplyr::rename(dep = "深度.[m]")
  min = min(t[, 3])
  max = max(t[, 3])
  s = t %>% filter(dep == min) %>% mutate(SorB = "S")
  b = t %>% filter(dep == max) %>% mutate(SorB = "B")
  sb = rbind(s, b)
  t2 = abind(t2, sb, along = 1)
}
envk2 = data.frame(t2) %>% mutate(date = as.Date(日付情報)) %>% mutate(year = 2018, month = as.numeric(str_sub(date, 6, 7)), day = as.numeric(str_sub(date, 9, 10)), Point = "129")
colnames(envk2) #pHがない！！
envk2 = envk2[, c(4,5,14,17,19:22)]
colnames(envk2)
colnames(envk2) = c("temp", "salinity", "DO_mg","SorB", "year", "month", "day", "Point")
envk2 = envk2 %>% mutate(pH = NA)

env = rbind(envk1, envk2)

setwd("/Users/Yuki/Dropbox/eDNA_INLA")
write.csv(env, "env2018kanagawa.csv", fileEncoding = "CP932")
