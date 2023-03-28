# Load required libraries
library(tidyverse);library(lubridate);library(janitor);library(readr)
library(ggforce); library(scico); library(ggtext);library(ggthemes);library(gtsummary)
library(survival); library(survminer); library(gt);library(patchwork);library(rms)
library(gghalves)

# Time to LVSD from first SHaRe visit
# This section filters and selects columns from a data frame called 'df' and assigns them to a new data frame called 'df1'
df1<- df %>% filter(t2_d_htx_vad>echo_age0 & event_lvef50==1) %>% 
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

# Merge data to create time-dependent covariates
# This section uses the 'tmerge' function from the survival package to merge the 'df1' data frame with itself and create time-dependent covariates
df1_new<- 
  survival::tmerge(df1, df1, id= pid, srt = tdc(t2_srt, event_srt,0), 
                   lvef35 = tdc(t2_lvef35, event_lvef35,0),
                   endpt = event(t2_d_htx_vad, d_htx_vad_c),
                   
                   tstart = t2_lvef50, tstop = fu_comp_event)

# Run Cox regression model
# This section fits a Cox regression model using the 'coxph' function from the survival package and calculates hazard ratios and confidence intervals using the 'broom::tidy' function
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

# Prepare data for table and plot
# This section selects columns from the 'cox1' data frame, renames them, and assigns them to a
# a new data frame named "cox1_df" by selecting specific columns from "cox1" data frame and then modifying the 
# "term" column using case_when and str_detect functions
cox1_df <- cox1 %>% 
  select(term, .est = estimate, .low = conf.low, .high = conf.high, p= p.value) %>% 
  mutate(term = case_when(
    str_detect(term, "lvsd") ~ "Per 5-year increase",
    str_detect(term, "Not Genotyped") ~ "Not Genotyped",
    str_detect(term, "Thick") ~ "Thick filament P/LP",
    str_detect(term, "VUS") ~ "VUS",
    str_detect(term, "Non-Sarcomeric HCM") ~ "Non-Sarcomeric HCM",
    str_detect(term, "Sarcomeric HCM") ~ "Sarcomeric HCM",
    str_detect(term, "Female") ~ "Female",
    str_detect(term, "35") ~ "Yes",
    str_detect(term, "Multiple P/LP") ~ "Multiple P/LP",
    T ~ term # if none of the above conditions are met, keep the original "term" value
  )) %>%
  # Append a new row to the data frame containing additional variable names and default values
  rbind(
    tibble(
      term = c("Male", "Non-Sarcomeric HCM", "No", "No srt",
               "**Time to LVSD from diagnosis**", "**Sex**", "**Genetics**", "**LVEF < 35%**", "**Septal reduction therapy**"),
      .est = c(1, 1, 1, 1, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
      .low = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
      .high = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
      p = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
    )
  ) %>%
  # Reorder the rows of the data frame by releveling the "term" column in the desired order
  mutate(term = fct_relevel(term, "**Time to LVSD from diagnosis**", "Per 5-year increase", "**Sex**", "Male", "Female",
                            "**Genetics**", "Non-Sarcomeric HCM", "Sarcomeric HCM", "VUS", "Not Genotyped",
                            "**LVEF < 35%**", "No", "Yes", "**Septal reduction therapy**", "No srt", "srt")) %>%
  arrange(term)


#######
# In the following section I am creating the plot in 5 individual pieces which are later pasted together:
# A plot section (plot_fig) 
# A section containing the names of the variable (plot_terms)
# A section containing the p-values (plot_p)
# A section containing the number of events and overall counts in each subgroup (plot_events)
# A section containgin the hazard ratios (plot_hr)

# plot section of the figure
plot_fig<- 
  # create the plot using ggplot
  cox1_df %>% 
  # reverse the order of the y-axis
  mutate(term = fct_rev(term)) %>% 
  ggplot()+
  # add horizontal dashed line at y=1
  geom_hline(aes(yintercept = 1), color = "#C87970", linetype = 3)+
  # add error bars for confidence intervals
  geom_errorbar(aes(x = term, ymin= .low, ymax = .high), width = .2)+
  # add points for estimates
  geom_point(aes(x= term, y = .est), size = 3, shape = 21, 
             fill = "#60A5DF")+ 
  # label the y-axis
  labs(y= "Hazard ratio",
       x = "")+
  # use a logarithmic scale for the y-axis
  scale_y_log10(breaks = c(0.2,.5,1,2,4,8))+
  # remove all non-data elements of the theme
  theme_void()+
  # customize the appearance of the x-axis labels
  theme(axis.text.x = element_markdown(size = 10),
        axis.ticks.length.x = unit(.3, "cm"),
        axis.ticks.x = element_line(),
        axis.line.x = element_blank()
  )+ 
  # flip the x and y axes
  coord_flip()+
  # add a vertical line at x=0.4
  geom_segment(aes(x = 0.4, xend = .4, y=0.2, yend = 8))

# display the plot
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

#Using the patchwork package, they are now pasted together
plot_terms+plot_events+plot_hr+plot_fig+plot_p+plot_layout(widths = c(4,3,3,6,1))
# and then saved as tiff and pdf-files
ggsave("new_figure 4.tiff", compression = "lzw",units = "cm", width = 20, height = 14, dpi = 1200)
ggsave("new_figure 4.pdf",units = "cm", width = 20, height = 14, dpi = 1200)

