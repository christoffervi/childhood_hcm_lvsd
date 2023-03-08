library(tidyverse);library(lubridate);library(janitor);library(readr)
library(ggforce); library(scico); library(ggtext);library(ggthemes);library(gtsummary)
library(survival); library(survminer); library(gt);library(patchwork);library(rms)
library(gghalves)

#Time to LVSD from first SHaRe visit
df1<- df %>% filter(t2_d_htx_vad>echo_age0 & event_lvef50==1) %>% 
  filter(is_proband==1) %>% 
  select(pid, t2_srt, t2_d_htx_vad,event_srt, d_htx_vad_c, echo_age0,
         t2_lvef35, event_lvef35,t2_lvef50,
         lvsd_timing,sex,genetics, primary_diagnosis_age) %>% 
  mutate(fu_comp_event = t2_d_htx_vad,
         fu_lvef35_event = t2_lvef35,
         fu_srt = if_else(t2_srt-t2_lvef50<0,0,t2_srt),
         genetics = fct_relevel(genetics, "Non-Sarcomeric HCM"),
         fu_lvsd = (t2_lvef50-primary_diagnosis_age
         )/5
  ) %>% 
  filter(fu_comp_event>t2_lvef50)


df1_new<- 
  survival::tmerge(df1, df1, id= pid, srt = tdc(t2_srt, event_srt,0), 
                   lvef35 = tdc(t2_lvef35, event_lvef35,0),
                   endpt = event(t2_d_htx_vad, d_htx_vad_c),
                   
                   tstart = t2_lvef50, tstop = fu_comp_event)

cox1<-
  df1_new %>% tibble() %>% 
  mutate(sex = fct_relevel(sex, "Male"),
         genetics = case_when(genetics== "Not tested"~"Not Genotyped",
                              genetics== "Genonegative"~"Non-Sarcomeric HCM",
                              genetics== "Non-Sarcomeric HCM"~"Non-Sarcomeric HCM",
                              genetics== "VUS"~"VUS",
                              T~"Sarcomeric HCM"),
         genetics = fct_relevel(genetics, "Non-Sarcomeric HCM")) %>% 
  survival::coxph(
    survival::Surv(tstart, tstop, endpt)~
      fu_lvsd+sex+lvef35+genetics+srt, data = .) %>% 
  broom::tidy(exponentiate=T, conf.int=T)

cox1
cox1_df<-
  cox1 %>% 
  select(term, .est = estimate, .low = conf.low, .high = conf.high, p= p.value) %>% 
  mutate(term = case_when(str_detect(term, "lvsd")~"Per 5-year increase",
                          str_detect(term, "Not Genotyped")~"Not Genotyped",
                          str_detect(term, "Thick")~"Thick filament P/LP",
                          str_detect(term, "VUS")~"VUS",
                          str_detect(term, "Non-Sarcomeric HCM")~"Non-Sarcomeric HCM",
                          str_detect(term, "Sarcomeric HCM")~"Sarcomeric HCM",
                          
                          str_detect(term, "Female")~"Female",
                          str_detect(term, "35")~"Yes",
                          str_detect(term, "Multiple P/LP")~"Multiple P/LP",
                          T~term
  )) %>% rbind(
    tibble(term = c("Male", "Non-Sarcomeric HCM", "No","No srt",
                    "**Time to LVSD from diagnosis**", "**Sex**", "**Genetics**", "**LVEF < 35%**","**Septal reduction therapy**"),
           .est = c(1,1,1,1,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_),
           .low= c(NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_),
           .high = c(NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_),
           p = c(NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_))) %>%
  mutate(term = fct_relevel(term,"**Time to LVSD from diagnosis**", "Per 5-year increase", "**Sex**",  "Male","Female",
                            "**Genetics**", "Non-Sarcomeric HCM", "Sarcomeric HCM",  "VUS","Not Genotyped", 
                            "**LVEF < 35%**",
                            "No", "Yes", "**Septal reduction therapy**", "No srt", "srt")) %>% 
  arrange(term)

# the plot-part of the figure
plot_fig<-
  cox1_df %>% 
  mutate(term = fct_rev(term)) %>% 
  ggplot()+
  geom_hline(aes(yintercept = 1), color = "#C87970", linetype = 3)+
  geom_errorbar(aes(x = term, ymin= .low, ymax = .high), width = .2)+
  geom_point(aes(x= term, y = .est), size = 3, shape = 21, 
             fill = "#60A5DF")+ 
  #fill ="#C87970")+
  #geom_segment(aes(x = .7, xend=.7, y= .2, yend=4))+
  labs(y= "Hazard ratio",
       x = "")+
  scale_y_log10(breaks = c(0.2,.5,1,2,4,8))+
  theme_void()+theme(axis.text.x = element_markdown(size = 10),
                     #axis.line.x = element_line(),
                     axis.ticks.length.x = unit(.3, "cm"),
                     axis.ticks.x = element_line(),
                     axis.line.x = element_blank()
  )+ coord_flip()+
  geom_segment(aes(x = 0.4, xend = .4, y=0.2, yend = 8))
plot_fig

# p-values for the plot
plot_p<-
  cox1_df %>% 
  mutate(term = fct_rev(term),
         p = case_when(str_detect(term, "LVSD|Genetics|LVEF|Sex|Septal")~NA_character_,
                       is.na(p)~"ref",
                       p<0.001~"<0.001",
                       T~as.character(round(p,3)))) %>% 
  ggplot()+
  geom_text(aes(x = term, y=1, label = p), size = 3)+
  theme_void()+
  labs(title = "**p-value**")+
  coord_flip()+ theme(plot.title = element_markdown(hjust = .5, size = 10))
plot_hr<-
  cox1_df %>% 
  mutate(term = fct_rev(term),
         p = case_when(str_detect(term, "LVSD|Genetics|LVEF|Sex|Septal")~NA_character_,
                       is.na(p)~"ref",
                       T~paste(round(.est,2), "(",round(.low,2),"-",round(.high,2),")" ))) %>% 
  ggplot()+
  geom_text(aes(x = term, y=1, label = p), size = 3)+
  theme_void()+
  labs(title = "**Hazard ratio <br> (95% CI)**")+
  coord_flip()+ theme(plot.title = element_markdown(hjust = .5, size = 10))

# Variable names for the plot
plot_terms<-
  cox1_df %>% 
  mutate(term = fct_rev(term),
         term_bold = case_when(str_detect(term, "Time to")~"Time to LVSD from diagnosis",
                               str_detect(term, "Sex")~"Sex",
                               str_detect(term, "Genet")~"Genetics",
                               str_detect(term, "LVEF")~"LVEF < 35%",
                               str_detect(term, "Septal")~"Septal reduction therapy",
                               T~NA_character_),
         term_plain = case_when(str_detect(term, "Sex|Genetics|LVEF|Septal|Time")~NA_character_,
                                term=="No srt"~"No",
                                term=="srt"~"Yes",
                                T~as.character(term))) %>% 
  ggplot()+
  geom_text(aes(x = term, y=.5, label = term_plain), size = 3, hjust = 0)+
  geom_text(aes(x = term, y=0, label = term_bold), size = 3, hjust = 0, fontface = "bold")+
  #geom_richtext(aes(x = term, y=.5, label = term), fill = NA, label.color = NA, 
  #            size = 3, hjust = 0)+
  theme_void()+
  labs(title = "**Variable**")+
  coord_flip(ylim = c(0,3))+ theme(plot.title = element_markdown(hjust = .47, size = 10),
  )

# This is for non-timevarying
plot_events<-
  df %>% filter(event_lvef50==1) %>%
filter(is_proband==1) %>% 
    select(t2_lvef50, sex, genetics, lvef35_at_onset, event_srt ,d_htx_vad_c) %>%
  mutate(genetics = case_when(genetics== "Not tested"~"Not Genotyped",
                              genetics== "Genonegative"~"Non-Sarcomeric HCM",
                              genetics== "Non-Sarcomeric HCM"~"Non-Sarcomeric HCM",
                              genetics== "VUS"~"VUS",
                              T~"Sarcomeric HCM"),
         genetics = fct_relevel(genetics, "Non-Sarcomeric HCM"),
         event_srt = as.character(event_srt),
         t2_lvef50 = "Per 5-year increase") %>% 
  pivot_longer(1:5) %>% group_by(value) %>%
  summarise(combined = sum(d_htx_vad_c), n=n())%>%
  rbind(
    tibble(value = c("**Time to LVSD from diagnosis**", "**Sex**", "**Genetics**", "**LVEF < 35%**", "**Septal reduction therapy**"),
           combined = c(NA_real_,NA_real_,NA_real_,NA_real_,NA_real_),
           n =c(NA_real_,NA_real_,NA_real_,NA_real_,NA_real_))) %>%
  mutate(term = fct_relevel(value,"**Time to LVSD from diagnosis**", "Per 5-year increase", "**Sex**",  "Male","Female",
                            "**Genetics**", "Non-Sarcomeric HCM", "Sarcomeric HCM",  "VUS","Not Genotyped",
                            "**LVEF < 35%**",
                            "no", "yes", "**Septal reduction therapy**", "0", "1")) %>% 
  select(term, combined, n) %>%
  arrange(term) %>%
  mutate(term = fct_rev(term),
         prop = round(combined/n *100,0),
         plot_label = case_when(str_detect(term, "Time to|Sex|Genet|LVEF|Sept")~NA_character_,
                                T~paste(combined, "/", n, " (",prop, "%)"))) %>%
  ggplot()+
  geom_text(aes(x = term, y=.5, label = plot_label), size = 3)+
  theme_void()+
  labs(title = "**Death, Txp, LVAD <br> n (%)**")+
  coord_flip()+ theme(plot.title = element_markdown(hjust = .5, size = 10),
  )
plot_fig
plot_terms+plot_events+plot_hr+plot_fig+plot_p+plot_layout(widths = c(4,3,3,6,1))+
  plot_annotation(caption = "**Sensitivity analysis**: Only probands have been included in the analysis",
                  theme = theme(plot.caption = element_markdown()))
ggsave("proband_only_figure 4.png",units = "cm", width = 20, height = 14, dpi = 1200)
ggsave("proband_only_figure 4.pdf",units = "cm", width = 20, height = 14, dpi = 1200)










##################
#########################################################
# FIGURE X CAUSE-specific hazard
##########################################################

#all patients without srt
split_df1 <- df_dcc2 %>% filter(event_srt==0) %>% 
  mutate(fu_lvsd1 =0,
         fu_lvsd2 = fu_lvsd_cmp)

#patients with SRT after first share visit in the period prior to srt
split_df2 <- df_dcc2 %>% 
  filter(event_srt==1 & t2_srt-echo_age0>0) %>% 
  mutate(fu_lvsd = case_when(t2_srt>=t2_lvef50~fu_lvsd_cmp,
                             T~t2_srt-echo_age0),
         event_srt = 0,
         event_lvef50 = case_when(t2_srt>=t2_lvef50~event_lvef50,
                                  T~0),
         fu_lvsd1 =0,
         fu_lvsd2 = fu_lvsd_cmp)

#patients with SRT after first share visit in the period after to srt
split_df3 <- df_dcc2 %>% 
  filter(event_srt==1 & t2_srt-echo_age0>0) %>% 
  mutate(fu_lvsd = case_when(t2_srt>=t2_lvef50~NA_real_,
                             T~t2_lvef50-t2_srt),
         event_srt = 1,
         event_lvef50 = case_when(t2_srt>=t2_lvef50~0,
                                  T~event_lvef50),
         fu_lvsd1 =t2_srt-echo_age0,
         fu_lvsd2 = fu_lvsd_cmp-fu_lvsd1)

#patients with SRt at baseline visit
split_df4<- df_dcc2 %>% filter(event_srt==1 & t2_srt<=echo_age0) %>% 
  mutate(fu_lvsd1 =0,
         fu_lvsd2 = fu_lvsd_cmp)

df_dcc3 <- rbind(split_df1,split_df2, split_df3, split_df4)
df_dcc3 <- df_dcc3 %>% mutate(sex = fct_relevel(sex, "Female"))




df_dcc3_a<- 
  df_dcc3 %>% 
  filter(is_proband==1) %>% 
  mutate(
    genetics = case_when(genetics== "Not tested"~"Not Genotyped",
                         genetics== "Genonegative"~"Non-Sarcomeric HCM",
                         genetics== "Non-Sarcomeric HCM"~"Non-Sarcomeric HCM",
                         genetics== "VUS"~"VUS",
                         T~"Sarcomeric HCM"),
                     genetics = fct_relevel(genetics, "Non-Sarcomeric HCM", "Sarcomeric HCM","VUS")
  )


cox2<- 
  coxph(Surv(fu_lvsd1, fu_lvsd2, cmpr_event==1)~
          hcm_pub+sex+genetics+lvef5+event_srt, data = df_dcc3_a)
cox2_df<-
  cox2 %>% 
  broom::tidy(conf.int = T, exponentiate = T) %>% 
  select(term, .est = estimate, .low = conf.low, .high = conf.high, p= p.value) %>% 
  mutate(term = case_when(str_detect(term, "hcm")~"<12",
                          str_detect(term, "Thin")~"Thin filament P/LP",
                          str_detect(term, "Thick")~"Thick filament P/LP",
                          str_detect(term, "VUS")~"VUS",
                          str_detect(term, "Non-Sarcomeric HCM")~"Non-Sarcomeric HCM",
                          str_detect(term, "Not Genotyped")~"Not Genotyped",
                          str_detect(term, "Male")~"Male",
                          str_detect(term, "lvef")~"LVEF (per 5% decrease)",
                          str_detect(term, "geneticsSarcomeric")~"Sarcomeric HCM",
                          str_detect(term, "event_srt")~"Yes",
                          T~term
  )) %>% rbind(
    tibble(term = c("12-18", "Female", "Non-Sarcomeric HCM", "No",
                    "**Age at HCM diagnosis**", "**Sex**", "**Genetics**", "**Septal reduction theraphy**"),
           .est = c(1,1,1,1,NA_real_,NA_real_,NA_real_,NA_real_),
           .low= c(NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_),
           .high = c(NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_),
           p = c(NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_,NA_real_))) %>% 
  mutate(term = fct_relevel(term,"**Age at HCM diagnosis**", "12-18", "<12", "**Sex**", "Female", "Male",
                            "**Genetics**", "Non-Sarcomeric HCM", "Sarcomeric HCM", "VUS","Not Genotyped",  
                            
                            "**Septal reduction theraphy**", "No", "Yes", 
                            "LVEF (per 5% decrease)")) %>% 
  arrange(term)

plot_fig<-
  cox2_df %>% 
  mutate(term = fct_rev(term)) %>% 
  ggplot()+
  geom_hline(aes(yintercept = 1), color = "#C87970", linetype = 3)+
  geom_errorbar(aes(x = term, ymin= .low, ymax = .high), width = .2)+
  geom_point(aes(x= term, y = .est), size = 3, shape = 21, 
             fill = "#60A5DF")+ 
  #fill ="#C87970")+
  #geom_segment(aes(x = .7, xend=.7, y= .2, yend=4))+
  labs(y= "Hazard ratio",
       x = "")+
  scale_y_log10(breaks = c(0.2,.5,1,2,4))+
  theme_void()+theme(axis.text.x = element_markdown(size = 10),
                     axis.line.x = element_blank(),
                     axis.ticks.length.x = unit(.3, "cm"),
                     axis.ticks.x = element_line()
                     #     axis.text.y = element_markdown()
  )+coord_flip()+
  geom_segment(aes(x= .4, xend=.4, y=.2, yend=4))
  #lemon::coord_capped_flip(bottom = "both", ylim = c(.5,6))

plot_p<-
  cox2_df %>% 
  mutate(term = fct_rev(term),
         p = case_when(str_detect(term, "at HCM|Genetics|Sex|Septal")~NA_character_,
                       is.na(p)~"ref",
                       p<0.001~"<0.001",
                       T~as.character(round(p,3)))) %>% 
  ggplot()+
  geom_text(aes(x = term, y=1, label = p), size = 3)+
  theme_void()+
  labs(title = "**p-value**")+
  coord_flip()+ theme(plot.title = element_markdown(hjust = .5, size = 10))

plot_hr<-
  cox2_df %>% 
  mutate(term = fct_rev(term),
         across(.cols = c(.est, .low, .high),~if_else(.x>3, round(.x,1), round(.x,2))),
         
         p = case_when(str_detect(term, "at HCM|Genetics|Sex|Septal")~NA_character_,
                       is.na(p)~"ref",
                       T~paste(sep ="",round(.est,2), " (",round(.low,2),"-",round(.high,2),")" ))) %>% 
  ggplot()+
  geom_text(aes(x = term, y=1, label = p), size = 3)+
  theme_void()+
  labs(title = "**Hazard ratio <br> (95% CI)**")+
  coord_flip()+ theme(plot.title = element_markdown(hjust = .5, size = 10))

plot_terms<-
  cox2_df %>% 
  mutate(term = fct_rev(term),
         term_bold = case_when(str_detect(term, "Age")~"Age at HCM diagnosis",
                               str_detect(term, "Sex")~"Sex",
                               str_detect(term, "Genet")~"Genetics",
                               str_detect(term, "LVEF")~"LVEF (per 5% decrease)",
                               str_detect(term, "Septal")~"Septal reduction therapy",
                               T~NA_character_),
         term_plain = case_when(str_detect(term, "Age|Sex|Genet|LVEF|Septal")~NA_character_,
                                T~as.character(term))) %>% 
  ggplot()+
  geom_text(aes(x = term, y=.5, label = term_plain), size = 3, hjust = 0)+
  geom_text(aes(x = term, y=0, label = term_bold), size = 3, hjust = 0, fontface = "bold")+
  #geom_richtext(aes(x = term, y=.5, label = term), fill = NA, label.color = NA, 
  #            size = 3, hjust = 0)+
  theme_void()+
  labs(title = "**Variable**")+
  coord_flip(ylim = c(0,3))+ theme(plot.title = element_markdown(hjust = .47, size = 10),
  )

plot_events<-
  df_dcc2 %>% 
  filter(is_proband==1) %>% 
  select(hcm_pub, sex, genetics, primary_diagnosis,event_srt, event_lvef50,t2_lvef50,t2_srt) %>%
  mutate(primary_diagnosis= "LVEF (per 5% decrease)",
         event_srt = case_when(event_srt==1 & event_lvef50==1 & t2_lvef50<t2_srt~0,
                               T~event_srt),
         event_srt = factor(event_srt)) %>% 
  mutate(genetics = case_when(genetics== "Not tested"~"Not Genotyped",
                              genetics== "Genonegative"~"Non-Sarcomeric HCM",
                              genetics== "Non-Sarcomeric HCM"~"Non-Sarcomeric HCM",
                              genetics== "VUS"~"VUS",
                              T~"Sarcomeric HCM"),
         genetics = fct_relevel(genetics, "Non-Sarcomeric HCM", "Sarcomeric HCM","VUS")
  ) %>% 
  pivot_longer(1:5) %>% group_by(value) %>% 
  summarise(combined = sum(event_lvef50), n=n())%>% 
  rbind(
    tibble(value = c("**Age at HCM diagnosis**", "**Sex**", "**Genetics**", "**Septal reduction therapy**"),
           combined = c(NA_real_,NA_real_,NA_real_,NA_real_),
           n =c(NA_real_,NA_real_,NA_real_,NA_real_))) %>% 
  mutate(term = fct_relevel(value,"**Age at HCM diagnosis**", "12-18", "<12", "**Sex**", "Female", "Male",
                            "**Genetics**", "Non-Sarcomeric HCM", "Sarcomeric HCM", "VUS", "Not Genotyped",  
                            
                            "**Septal reduction therapy**", "0", "1", 
                            "LVEF (per 5% decrease)")) %>% 
  select(term, combined, n) %>%  
  arrange(term) %>% 
  mutate(term = fct_rev(term),
         prop = round(combined/n *100,0),
         plot_label = case_when(str_detect(term, "Age|Sex|Genet|Septal")~NA_character_,
                                T~paste(sep ="",combined, "/", n, " (",prop, "%)"))) %>% 
  ggplot()+
  geom_text(aes(x = term, y=.5, label = plot_label), size = 3)+
  theme_void()+
  labs(title = "**LVSD <br> n (%)**")+
  coord_flip()+ theme(plot.title = element_markdown(hjust = .5, size = 10),
  )

plot_terms+plot_events+plot_hr+plot_fig+plot_p+plot_layout(widths = c(2.7,2,2,5,1))+
  plot_annotation(caption = "**Sensitivity analysis**: Only probands have been included in the analysis",
                  theme = theme(plot.caption = element_markdown()))

ggsave("proband_analysis 2.pdf",units = "cm", width = 20, height = 14, dpi = 1200)
ggsave("proband_analysis 2.png", units = "cm", width = 20, height = 14, dpi = 1200)
#ggsave("new_figure 2_.tiff", units = "cm", width = 18, height = 12, dpi = 1800, compression = "lzw")






plot_terms+plot_events+plot_hr+plot_fig+plot_p+plot_layout(widths = c(4,3,3,6,1))+
  plot_annotation(caption = "**Sensitivity analysis**: Only probands have been included in the analysis",
                  theme = theme(plot.caption = element_markdown()))
