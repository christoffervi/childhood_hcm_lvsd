# Create a new dataframe with additional variables
df_dcc2 <- df %>% 
  mutate(
    lvef = lvef_baseline,
    lvef_cat = case_when(
      lvef >= 80 ~ "≥80%",
      lvef >= 60 ~ "60-80",
      lvef > 55 ~ "55-60",
      lvef >= 50 ~ "50-55",
      lvef >= 30 ~ "30-50",
      T ~ NA_character_
    ),
    genetics = case_when(
      str_detect(genes_plp, ";") ~ "Multiple P/LP",
      str_detect(genes_plp, "TNN|TPM|ACT") ~ "Thin filament P/LP",
      str_detect(genes_plp, "MY") ~ "Thick filament P/LP",
      sarc_status == "SARC(U)" ~ "VUS",
      sarc_status == "SARC(-)" ~ "Non-Sarcomeric HCM",
      T ~ "Not tested"
    ),
    genetics = fct_infreq(genetics),
    genetics = fct_relevel(genetics, "Non-Sarcomeric HCM"),
    sex = fct_relevel(sex, "Female"),
    hcm_pub = case_when(
      primary_diagnosis_age >= 12 ~ "12-18",
      is.na(primary_diagnosis_age) ~ NA_character_,
      T ~ "<12"
    ),
    hcm_pub = fct_relevel(hcm_pub, "12-18"),
    lvef5 = -lvef / 5,
    lvef10 = lvef / 10,
    hcm_age5 = primary_diagnosis_age / 5
  ) %>% 
  filter(fu_lvsd_cmp > 0)

# Note: This code creates a new dataframe `df_dcc2` based on an existing dataframe `df`
# The code adds new variables to `df_dcc2` using `mutate()`. The new variables are:
#   - `lvef`: the LVEF at baseline
#   - `lvef_cat`: categorical variable based on LVEF values
#   - `genetics`: categorizes genetic mutations as "Multiple P/LP", "Thin filament P/LP", 
#                 "Thick filament P/LP", "VUS", "Non-Sarcomeric HCM", or "Not tested"
#   - `sex`: reorders sex so "Female" is first
#   - `hcm_pub`: categorical variable based on age at diagnosis
#   - `lvef5`: LVEF divided by 5
#   - `lvef10`: LVEF divided by 10
#   - `hcm_age5`: age at diagnosis divided by 5
# The code then filters out rows where `fu_lvsd_cmp` is not greater than 0


############
# Create a new dataframe with patients split based on the event_srt variable

# all patients without srt
split_df1 <- df_dcc2 %>% filter(event_srt==0) %>% 
  mutate(fu_lvsd1 =0,
         fu_lvsd2 = fu_lvsd_cmp)

# patients with SRT after first share visit in the period prior to srt
split_df2 <- df_dcc2 %>% 
  filter(event_srt==1 & t2_srt-echo_age0>0) %>% 
  mutate(fu_lvsd = case_when(t2_srt>=t2_lvef50~fu_lvsd_cmp,
                             T~t2_srt-echo_age0),
         event_srt = 0,
         event_lvef50 = case_when(t2_srt>=t2_lvef50~event_lvef50,
                                  T~0),
         fu_lvsd1 =0,
         fu_lvsd2 = fu_lvsd_cmp)

# patients with SRT after first share visit in the period after srt
split_df3 <- df_dcc2 %>% 
  filter(event_srt==1 & t2_srt-echo_age0>0) %>% 
  mutate(fu_lvsd = case_when(t2_srt>=t2_lvef50~NA_real_,
                             T~t2_lvef50-t2_srt),
         event_srt = 1,
         event_lvef50 = case_when(t2_srt>=t2_lvef50~0,
                                  T~event_lvef50),
         fu_lvsd1 =t2_srt-echo_age0,
         fu_lvsd2 = fu_lvsd_cmp-fu_lvsd1)

# patients with SRt at baseline visit
split_df4<- df_dcc2 %>% filter(event_srt==1 & t2_srt<=echo_age0) %>% 
  mutate(fu_lvsd1 =0,
         fu_lvsd2 = fu_lvsd_cmp)

# Combine the split dataframes into a single dataframe
df_dcc3 <- rbind(split_df1,split_df2, split_df3, split_df4)

# Rearrange the factor levels of the sex variable
df_dcc3 <- df_dcc3 %>% mutate(sex = fct_relevel(sex, "Female"))

#The code splits the df_dcc2 dataframe into four different dataframes (split_df1, split_df2, split_df3, and split_df4) 
# based on the value of the event_srt variable. 
#The df_dcc3 dataframe is then created by concatenating these four dataframes together.
#Each of the four split dataframes is created by filtering df_dcc2 based on different conditions. 
# The split_df1 dataframe includes all patients without SRT, 
# the split_df2 dataframe includes patients with SRT after their first SHARE visit but before the first LVEF < 50% event, 
# the split_df3 dataframe includes patients with SRT after their first SHARE visit and after the first LVEF < 50% event, 
# and the split_df4 dataframe includes patients with SRT at their baseline visit.
#The df_dcc3 dataframe has the same columns as the df_dcc2 dataframe, but the fu_lvsd1 and fu_lvsd2 columns have 
# been added to indicate the time periods before and after the Split


#####
# Fit Cox regression model with time-to-event outcomes and explanatory variables
cox2 <- coxph(Surv(fu_lvsd1, fu_lvsd2, cmpr_event == 1) ~ hcm_pub + sex + genetics + lvef5 + event_srt, data = df_dcc3)

# Prepare summary table of results from Cox regression model
cox2_df <- cox2 %>%
  broom::tidy(conf.int = T, exponentiate = T) %>% 
  # Select relevant columns and relabel them
  select(term, .est = estimate, .low = conf.low, .high = conf.high, p = p.value) %>%
  # Relabel terms based on their original names
  mutate(term = case_when(
    str_detect(term, "hcm") ~ "<12",
    str_detect(term, "Thin") ~ "Thin filament P/LP",
    str_detect(term, "Thick") ~ "Thick filament P/LP",
    str_detect(term, "VUS") ~ "VUS",
    str_detect(term, "Non-Sarcomeric HCM") ~ "Non-Sarcomeric HCM",
    str_detect(term, "Not tested") ~ "Not tested",
    str_detect(term, "Male") ~ "Male",
    str_detect(term, "lvef") ~ "LVEF (per 5% decrease)",
    str_detect(term, "Multiple P/LP") ~ "Multiple P/LP",
    str_detect(term, "event_srt") ~ "Yes",
    T ~ term
  )
  ) %>%
  # Append table with additional rows for explanatory variables that aren't in original model
  rbind(
    tibble(
      term = c(
        "12-18", "Female", "Non-Sarcomeric HCM", "No",
        "**Age at HCM diagnosis**", "**Sex**", "**Genetics**", "**Septal reduction theraphy**"
      ),
      .est = c(1, 1, 1, 1, NA_real_, NA_real_, NA_real_, NA_real_),
      .low = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
      .high = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
      p = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
    )
  ) %>% 
  # Reorder rows by explanatory variable
  mutate(term = fct_relevel(
    term, "**Age at HCM diagnosis**", "12-18", "<12", "**Sex**", "Female", "Male",
    "**Genetics**", "Non-Sarcomeric HCM", "Not tested", "VUS", "Thick filament P/LP", "Thin filament P/LP",
    "Multiple P/LP",
    "**Septal reduction theraphy**", "No", "Yes", "LVEF (per 5% decrease)"
  )) %>%
  arrange(term) # Sort rows by explanatory variable

######
#This code creates four plots: plot_fig, plot_p, plot_terms and plot_hr.
#####
#plot_fig is a forest plot displaying the hazard ratios of various variables from the cox regression analysis. 
#The variables are arranged in decreasing order of hazard ratio, and each variable is represented by a point. 
#The x-axis displays the variables, and the y-axis displays the hazard ratios on a logarithmic scale. 
#The width of the error bars represents the confidence intervals, and the color of the points represents thevariables.

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
  scale_y_log10(breaks = c(0.2,.5,1,2,4,8,16))+
  theme_void()+theme(axis.text.x = element_markdown(size = 10),
                     axis.line.x = element_line(),
                     axis.ticks.length.x = unit(.3, "cm"),
                     axis.ticks.x = element_line()
                     #     axis.text.y = element_markdown()
  )+
  lemon::coord_capped_flip(bottom = "both")

# plot_p is a plot displaying the p-values of the variables from the cox regression analysis. 
# The p-values are displayed as text on the plot, and the variables are arranged in decreasing order of 
# hazard ratio. The x-axis displays the variables, and the y-axis is empty.
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


# plot_hr is a plot displaying the hazard ratios and their 95% confidence intervals of the variables 
# from the cox regression analysis. The hazard ratios are displayed as text on the plot, 
# and the variables are arranged in decreasing order of hazard ratio. 
# The x-axis displays the variables, and the y-axis is empty.
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


# plot_terms is a plot displaying the variable names in a tabular form. 
# The bold text represents the variable category, and the plain text represents the variable name. 
# The x-axis displays the variable names, and the y-axis is empty. 
# The variable names are arranged in decreasing order of hazard ratio.
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

# plot_events creates a plot of the number and percentage of patients who developed LVSD for 
# each level of the categorical variables included in the model. 
# The pivot_longer() function is used to reshape the data from wide to long format, and summarise() is 
# used to calculate the total number of patients who developed LVSD and the total number of patients for 
# each level of the categorical variables. 
# The ggplot() function is used to plot the results as text on the y-axis, 
# with the categorical variables on the x-axis.

plot_events<-
  df_dcc2 %>% 
  select(hcm_pub, sex, genetics, primary_diagnosis,event_srt, event_lvef50,t2_lvef50,t2_srt) %>%
  mutate(primary_diagnosis= "LVEF (per 5% decrease)",
         event_srt = case_when(event_srt==1 & event_lvef50==1 & t2_lvef50<t2_srt~0,
                               T~event_srt),
         event_srt = factor(event_srt)) %>% 
  pivot_longer(1:5) %>% group_by(value) %>% 
  summarise(combined = sum(event_lvef50), n=n())%>% 
  rbind(
    tibble(value = c("**Age at HCM diagnosis**", "**Sex**", "**Genetics**", "**Septal reduction therapy**"),
           combined = c(NA_real_,NA_real_,NA_real_,NA_real_),
           n =c(NA_real_,NA_real_,NA_real_,NA_real_))) %>% 
  mutate(term = fct_relevel(value,"**Age at HCM diagnosis**", "12-18", "<12", "**Sex**", "Female", "Male",
                            "**Genetics**", "Non-Sarcomeric HCM", "Not tested", "VUS", "Thick filament P/LP", "Thin filament P/LP", 
                            "Multiple P/LP",
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

# Finally, the plot_layout() function from the patchwork package is used to combine all of the plots into a 
# single grid.
plot_terms+plot_events+plot_hr+plot_fig+plot_p+plot_layout(widths = c(2.7,2,2,5,1))