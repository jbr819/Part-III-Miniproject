## ribbon plots of lesioning 
library(tidyverse)
library(ggplot2)
setwd("~/Course_Materials/Part-III-Miniproject/Stochastic_temporal/data")
getwd()
lesion_keratinocyte_data <- (read.csv("3000_kerat_lesion_pop_data_for_R_analysis"))
lesion_stemcell_data <- (read.csv("3000_stem_lesion_pop_data_for_R_analysis"))

colnames(lesion_keratinocyte_data) <- paste0("run", 1:ncol(lesion_keratinocyte_data))
colnames(lesion_stemcell_data) <- paste0("run",  1:ncol(lesion_stemcell_data))

tidy_stem_data <- lesion_stemcell_data %>% mutate(week = (1-200):(nrow(lesion_stemcell_data)-200)) %>% pivot_longer(cols=-week, names_to="run")
tidy_kerat_data <- lesion_keratinocyte_data %>% mutate(week = (1-200):(nrow(lesion_keratinocyte_data)-200)) %>% pivot_longer(cols=-week, names_to="run")

#tidy_stem_data <- lesion_stemcell_data %>% mutate(week = (1-200):(nrow(lesion_stemcell_data)-200)) %>% pivot_longer(cols=-week, names_to="run")
#tidy_kerat_data <- lesion_keratinocyte_data %>% mutate(week = (1-200):(nrow(lesion_keratinocyte_data)-200)) %>% pivot_longer(cols=-week, names_to="run")

tidy_stem_data$value <- as.numeric(tidy_stem_data$value)
tidy_kerat_data$value <- as.numeric(tidy_kerat_data$value)

tidy_stem_data$run <- as.factor(tidy_stem_data$run)
tidy_kerat_data$run <- as.factor(tidy_kerat_data$run)

write_csv(tidy_kerat_data, "kerat_data.csv")

ribbon_kerat_summary <- tidy_kerat_data %>% group_by(week) %>% 
  summarise(mean = mean(value), 
            usd = mean(value) + sd(value),
            lsd = mean(value) - sd(value),
            usd2 = mean(value) + 2*sd(value),
            lsd2 = mean(value) - 2*sd(value),
            se = mean(value)/sqrt(length(value)),
            uci = mean(value) + 1.96*sd(value)/sqrt(length(value)),
            lci = mean(value) - 1.96*sd(value)/sqrt(length(value)),
            sd_by_mean = sd(value)/mean(value))

ribbon_stem_summary <- tidy_stem_data %>% group_by(week) %>%
  summarise(mean = mean(value), 
            usd = mean(value) + sd(value),
            usd2 = mean(value) + 2*sd(value),
            lsd = mean(value) - sd(value),
            lsd2 = mean(value) - 2*sd(value),
            se = mean(value)/sqrt(length(value)),
            uci = mean(value) + 1.96*sd(value)/sqrt(length(value)),
            lci = mean(value) - 1.96*sd(value)/sqrt(length(value)),
            sd_by_mean = sd(value)/mean(value))


plot_kerat <- ggplot(ribbon_kerat_summary, aes(x=week)) + 
  geom_ribbon(data=ribbon_kerat_summary, aes(ymin=lsd, ymax=usd), fill="blue", alpha=0.2) +
  geom_line(aes(y=mean), color="blue4") +
  geom_line(aes(y=usd), color="blue4", linetype=2) +
  geom_line(aes(y=lsd), color="blue4", linetype=2) +
  geom_line(aes(y=(usd2)), color="blue4", linetype=3)+
  geom_line(aes(y=(lsd2)), color="blue4", linetype=3)+
  geom_ribbon(data = ribbon_stem_summary %>% filter(week >= 50 & week <= 60), 
              aes(ymin = 0, ymax =3900), 
              fill = 'antiquewhite4', 
              alpha = 0.30) +
  ylim(0, 3900) +
  xlim(0, 200) +
  labs(x = "Time (weeks)", y = "Population of Keratinocytes", title = "") +
  theme_bw()

plot_kerat

plot_stem_variance_to_mean <- ggplot(ribbon_kerat_summary, aes(x=week)) +
  geom_line(aes(y=sd_by_mean), color="blue") +
  labs(x = "Time (weeks)", y = "Ratio of SD to Mean", title = "Ratio of Standard Variation to Mean Population of Keratinocyte", parse=T) +
  xlim(0,200) +
  theme_bw() +
  geom_ribbon(data = ribbon_kerat_summary %>% filter(week >= 50 & week <= 60), 
              aes(ymin = 0, ymax =0.3), 
              fill = 'antiquewhite4', 
              alpha = 0.30) 
plot_stem_variance_to_mean

plot_stem  <- ggplot(ribbon_stem_summary, aes(x=week)) + 
  geom_ribbon(aes(ymin=lsd, ymax=usd), fill="orange", alpha=0.2) +
  geom_line(aes(y=mean), color="orange3") +
  geom_line(aes(y=usd), color="orange3", linetype=2) +
  geom_line(aes(y=lsd), color="orange3", linetype=2) +
  geom_line(aes(y=(usd2)), color="orange4", linetype=3)+
  geom_line(aes(y=(lsd2)), color="orange4", linetype=3)+
  labs(x = "Week", y = "Population of Stem Cells", title = "") +
  geom_ribbon(data = ribbon_stem_summary %>% filter(week >= 50 & week <= 60), 
              aes(ymin = 0, ymax =60), 
              fill = 'antiquewhite4', 
              alpha = 0.30) +
  xlim(0, 200) +
  theme_bw()
  
  plot_stem

  plot_stem_variance_to_mean <- ggplot(ribbon_stem_summary, aes(x=week)) +
    geom_line(aes(y=sd_by_mean), color="orange") +
    labs(x = "Time (weeks)", y = "Ratio of Variance to Mean", title = "Ratio of Standard Variation to Mean Population of Stem Cells", parse=T) +
    xlim(0,200) +
    theme_bw() +
    geom_ribbon(data = ribbon_stem_summary %>% filter(week >= 50 & week <= 60), 
                aes(ymin = 0, ymax =2.5), 
                fill = 'antiquewhite4', 
                alpha = 0.30) 
  
  plot_stem_variance_to_mean
  
plot_kerat + plot_stem 





