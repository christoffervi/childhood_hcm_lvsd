##### age-specific incidences of composite outcome according to LVSD status
#This plot gives the age-specific incidence of the composite outcome, according to whether LVSD has been identified (white) in the subject or not (dark purple).
######

# This code filters the data frame `df` to obtain three new data frames,
# `no_lvef_50`, `pre_lvef_50`, and `post_lvef_50`, which correspond to
# different groups of patients based on their LVSD status and the timing
# of when LVSD was diagnosed.

# `no_lvef_50` contains patients who did not have LVSD at any point during
# the study period.
no_lvef_50 <- df %>% 
  filter(event_lvef50==0) %>% 
  mutate(lvef50 = 0,
         fu_comp_event = t2_d_htx_vad-echo_age0)

# `pre_lvef_50` contains patients who had LVSD diagnosed at some point
# during the study period, but not at the time of the final follow-up
# visit (`t2_d_htx_vad`), and therefore have a period of follow-up time
# without LVSD. The `t2_d_htx_vad-t2_lvef50>0` condition ensures that
# only patients who had LVSD diagnosed before the final follow-up visit
# are included in this group.
pre_lvef_50 <- df %>% 
  filter(event_lvef50==1 & t2_d_htx_vad-t2_lvef50>0) %>%
  mutate(lvef50 = 0, 
         fu_comp_event = t2_lvef50-echo_age0,
         d_htx_vad_c =0)

# `post_lvef_50` contains patients who had LVSD diagnosed at some point
# during the study period and had LVSD at the final follow-up visit.
post_lvef_50 <- df %>% 
  filter(event_lvef50==1) %>%
  mutate(lvef50 = 1, 
         fu_comp_event = t2_d_htx_vad-t2_lvef50)

# The three data frames above are combined using `bind_rows` and the
# resulting data frame is assigned to `lvef50`. The `age_comp_event`
# column is created to represent the age at which the composite outcome
# occurs, which is the sum of the patient's age at the time of the final
# follow-up visit (`echo_age0`) and the duration of follow-up time
# (`fu_comp_event`).
lvef50 <- bind_rows(no_lvef_50, pre_lvef_50, post_lvef_50) %>% 
  mutate(age_comp_event = echo_age0+fu_comp_event)

# `lvef50_df` contains only the rows of `lvef50` where the age at which
# the composite outcome occurs is greater than or equal to the patient's
# age at the time of the final follow-up visit. This ensures that the
# composite outcome occurs after the final follow-up visit.
lvef50_df <- lvef50 %>% filter(age_comp_event-echo_age0>=0)


#Create survival object using echo_age0, age_comp_event, and d_htx_vad_c
fit1 <- Surv(lvef50_df$echo_age0, lvef50_df$age_comp_event+.001, event = lvef50_df$d_htx_vad_c)~event_lvef50

#Split data by age group using the cut points and episode name
spl <- as_tibble(survSplit(fit1, data = lvef50_df, cut = c(12,20,30,40), episode = "agegroup")) %>%
  #Calculate time for each age group and event type
  mutate(time = tstop-tstart) %>% 
  group_by(agegroup, event_lvef50) %>% 
  summarise(nyear = sum(time), ncas = sum(event)) %>% 
  ungroup() 

#Create table of observed events by age group and LVSD status
obs_tar <- spl %>% select(agegroup,event_lvef50, ncas) %>% pivot_wider(names_from = agegroup, values_from = ncas) %>% bind_rows(spl %>% select(-ncas) %>% pivot_wider(names_from = agegroup, values_from = nyear))

#Create matrix of expected events using standard population rates
std <- matrix(data = c(887+115,1773+202,1518+338,625+307,390+259), nrow = 1, dimnames = list(c(""), c(1,2,3,4,5)))

#Use epiR package to calculate age-standardized rates and relative risk
dir <- epiR::epi.directadj(obs_tar[1:2,2:6] %>% as.matrix(rownames.force=T), obs_tar[3:4,2:6] %>% as.matrix(rownames.force=T), std = std %>% as.matrix())

#Round relative risks and confidence intervals to one decimal place
label_paste <- dir$adj.strata %>% tibble() %>% mutate(across(.cols = 4:6,~round(.x*100, 1)))

#Create a data frame with observed and expected events, calculate the complement of the age-standardized rate (for plotting purposes), and bind to the original data
dat.df <- spl %>% select(ncas,nyear) %>% 
  as.matrix() %>% 
  epiR::epi.conf(ctype = "inc.rate", method = "byar") %>%  
  bind_cols(spl) %>% 
  mutate(sur = 1-est)

#### In this final part we make the plot based on the data above
dat.df %>% 
  mutate(agegroup = factor(agegroup, labels = c("<12","12-19", "20-29","30-39", ">40"))) %>% 
  ggplot(aes(x= agegroup, y= est, group = event_lvef50, ymin = lower, ymax= upper, fill = event_lvef50))+
  geom_errorbar(position = position_dodge(width = .7), width = .1)+
  geom_line(position = position_dodge(width = .7))+
  geom_point(position = position_dodge(width = .7), shape = 21, size = 4, show.legend = F)+
  scale_fill_scico(palette = "davos")+
  scale_y_continuous(breaks = seq(0,20,.02))+
annotate("text", x= c(.9,1.25), y = c(0.0085895339,0.069574897),
         label = c("No LVSD","LVSD"), hjust = 0, vjust= 0)+
  theme_pubclean()+
  labs(y ="Incidence Rate of Death, Cardiac Txp or LVAD per Person-year",
       x = "Age groups")+
  theme(axis.title.x = element_text(family = "Helvetica", size = 10, hjust = 1),
        axis.title.y = element_text(family = "Helvetica", size = 10, hjust = 1),
        legend.position = "none",
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(family = "Helvetica", size = 10),
        axis.text.y = element_text(family = "Helvetica", size = 10))
#ggsave("new_Figure S5.pdf", dpi = 1200, units = "cm", height = 12, width = 18)
#ggsave("new_Figure S5.tiff",compression= "lzw", dpi = 1200, units = "cm", height = 12, width = 18)

