library(tidyverse); library(gt); library(gtsummary); library(flextable)

# Creating a massive table containing 6 groups of patients in the study
df_dcc %>%  
  select(sex, primary_diagnosis_age, first_encounter_age, is_proband, nyha, race, 
         genetics,genes_table1,
         lvef_baseline, #event_lvoto50, 
         lvot_obs,event_lvef50,
         fu_lvsd_table1, lad, echo_max_lvt) %>% 
  mutate(genetics2 = case_when(genetics== "Not tested"~NA_character_,
                              genetics== "Genonegative"~"Non-Sarcomeric HCM",
                              genetics== "Non-Sarcomeric HCM"~"Non-Sarcomeric HCM",
                              genetics== "VUS"~"VUS",
                              T~"Sarcomeric HCM"),
         genetics2 = fct_infreq(genetics2),
         fu_lvsd_table1 = case_when(fu_lvsd_table1 == "prevalent lvsd"& primary_diagnosis_age>18~"Prevalent LVSD, adult-onset",
                             fu_lvsd_table1 == "prevalent lvsd"& primary_diagnosis_age<=18~"Prevalent LVSD, pediatric-onset",
                             event_lvef50 ==1 & primary_diagnosis_age>18~"Incident LVSD, adult-onset",
                             event_lvef50 ==1 & primary_diagnosis_age<=18~"Incident LVSD, pediatric-onset",
                             primary_diagnosis_age>18 ~"No LVSD, adult-onset",
                             primary_diagnosis_age<= 18~"No LVSD, pediatric-onset",
                             T~NA_character_)) %>% 
  select(!event_lvef50) %>% 
  gtsummary::tbl_summary(by = fu_lvsd_table1,
                         value = list(),
                         statistic = list(primary_diagnosis_age~"{median} ({p25} to {p75})",
                                          #'Total follow-up, yrs'~"{median} ({p25} to {p75})",
                                          echo_max_lvt~"{median} ({p25} to {p75})",
                                          #  'Max LV wall thickness, z-score'~"{median} ({p25} to {p75})",
                                          lvef_baseline~"{mean} ± {sd}",
                                          #                'LV end-diastolic diameter, mm'~"{mean} ± {sd}",
                                          lad~"{mean} ± {sd}"),
                         digits = list(primary_diagnosis_age~c(1,1),
                                       echo_max_lvt~c(1,1),
                                       #'Max LV wall thickness, z-score'~c(1,1),
                                       lvef_baseline~c(1,1),
                                       #            'LV end-diastolic diameter, mm'~c(1,1),
                                       lad~c(1,1))) %>% #gtsummary::add_overall() %>% 
  gtsummary::bold_labels() %>% 
   gtsummary::as_flex_table() #%>% 
  #flextable::save_as_docx(path="Ugliest_table_ever.docx") #Code below saves to image



#####
df_dcc %>%  
  select(sex, primary_diagnosis_age, first_encounter_age, is_proband, nyha, race, 
         genetics,genes_table1,
         lvef_baseline, #event_lvoto50, 
         lvot_obs,event_lvef50,
         fu_lvsd_table1, lad, echo_max_lvt) %>% 
  mutate(genetics2 = case_when(genetics== "Not tested"~NA_character_,
                               genetics== "Genonegative"~"Non-Sarcomeric HCM",
                               genetics== "Non-Sarcomeric HCM"~"Non-Sarcomeric HCM",
                               genetics== "VUS"~"VUS",
                               T~"Sarcomeric HCM"),
         genetics2 = fct_infreq(genetics2),
         race = case_when(race=="Black"~"Black",
                          race=="Asian"~"Asian",
                          race=="White"~"White",
                          T~"Other or Not Reported"),
         race= fct_relevel(race, "White", "Black", "Asian"),
         fu_lvsd_table1 = case_when(fu_lvsd_table1 == "prevalent lvsd"& primary_diagnosis_age>18~"Prevalent LVSD, adult-onset",
                                    fu_lvsd_table1 == "prevalent lvsd"& primary_diagnosis_age<=18~"Prevalent LVSD, pediatric-onset",
                                    primary_diagnosis_age>18 ~"No LVSD, adult-onset",
                                    primary_diagnosis_age<= 18~"No LVSD, pediatric-onset",
                                    T~NA_character_)) %>% 
  select(sex, primary_diagnosis_age, first_encounter_age, is_proband, nyha, race, 
         genetics2,genes_table1,
                lvef_baseline, #event_lvoto50, 
                lvot_obs,
                fu_lvsd_table1, lad, echo_max_lvt) %>% 
  gtsummary::tbl_summary(by = fu_lvsd_table1,
                         value = list(sex~"Female"),
                         label = list(primary_diagnosis_age~"Age at HCM diagnosis, years",
                                      first_encounter_age~"Age at first SHaRe evaluation, years",
                                      is_proband~"Proband",
                                      nyha~"NYHA functional class",
                                      race~"Self-Reported Race",
                                      lvef_baseline~"Initial LVEF (%)",
                                      lvot_obs~"Obstruction Present",
                                      lad~"LA diameter, mm",
                                      echo_max_lvt~"Maximal LV wall thickness, mm",
                                      genetics2~"Genetic Status"),
                         statistic = list(primary_diagnosis_age~"{median} ({p25} to {p75})",
                                          #'Total follow-up, yrs'~"{median} ({p25} to {p75})",
                                          echo_max_lvt~"{median} ({p25} to {p75})",
                                          #  'Max LV wall thickness, z-score'~"{median} ({p25} to {p75})",
                                          lvef_baseline~"{mean} ± {sd}",
                                          #                'LV end-diastolic diameter, mm'~"{mean} ± {sd}",
                                          lad~"{mean} ± {sd}"),
                         digits = list(primary_diagnosis_age~c(1,1),
                                       echo_max_lvt~c(1,1),
                                       #'Max LV wall thickness, z-score'~c(1,1),
                                       lvef_baseline~c(1,1),
                                       #            'LV end-diastolic diameter, mm'~c(1,1),
                                       lad~c(1,1))) #%>% #gtsummary::add_overall() %>% 
  #gtsummary::bold_labels() %>% as_gt() %>% gtsave("not_as_ugly.docx") 
  #gtsummary::as_flex_table() %>% 
  #flextable::save_as_docx(path="Not_as_ugly_table.docx") #Code below saves to image
























df %>%  
  select(sex, primary_diagnosis_age, first_encounter_age, is_proband, nyha, race, 
         genetics,genes_table1,
         lvef_baseline, #event_lvoto50, 
         lvot_obs,
         fu_lvsd_table1, lad, echo_max_lvt) %>% 
  mutate(genetics2 = case_when(genetics== "Not tested"~NA_character_,
                               genetics== "Genonegative"~"Non-Sarcomeric HCM",
                               genetics== "Non-Sarcomeric HCM"~"Non-Sarcomeric HCM",
                               genetics== "VUS"~"VUS",
                               T~"Sarcomeric HCM"),
         genetics2 = fct_infreq(genetics2),
  ) %>% 
  gtsummary::tbl_summary(by = fu_lvsd_table1,
                         value = list(),
                         statistic = list(primary_diagnosis_age~"{median} ({p25} to {p75})",
                                          #'Total follow-up, yrs'~"{median} ({p25} to {p75})",
                                          echo_max_lvt~"{median} ({p25} to {p75})",
                                          #  'Max LV wall thickness, z-score'~"{median} ({p25} to {p75})",
                                          lvef_baseline~"{mean} ± {sd}",
                                          #                'LV end-diastolic diameter, mm'~"{mean} ± {sd}",
                                          lad~"{mean} ± {sd}"),
                         digits = list(primary_diagnosis_age~c(1,1),
                                       echo_max_lvt~c(1,1),
                                       #'Max LV wall thickness, z-score'~c(1,1),
                                       lvef_baseline~c(1,1),
                                       #            'LV end-diastolic diameter, mm'~c(1,1),
                                       lad~c(1,1))) %>% #gtsummary::add_overall() %>% 
  gtsummary::bold_labels() ##%>% 
  #gtsummary::as_flex_table() %>% 
  #flextable::save_as_docx(path="filament info.docx") #Code below saves to image









######################


df %>%  filter(fu_lvsd_cmp<0) %>% select(event_lvef50)


df %>%  #filter(fu_lvsd_cmp>=0) %>% 
  select(event_lvef50, 
         primary_diagnosis_age,
         event_srt,event_vad,
         event_transplant, event_death, 
         death_cause
         ) %>% 
  mutate(event_lvef50 = if_else(event_lvef50==1, "LVSD", "No LVSD"),
         death_cause = case_when(str_detect(death_cause, "Non-Card|Infect|Unkn|Proce")~"Non-cardiovascular death",
                                 str_detect(death_cause, "SCD")~"Sudden cardiac death",
                                 str_detect(death_cause, "Myoca|Cardio|Strok")~"Other cardiovascular death",
                                 str_detect(death_cause, "Fail")~"Heart failure",
                                 T~death_cause)) %>% 
  gtsummary::tbl_summary(by = event_lvef50,
                         label = list(event_srt~"Septal reduction therapy",
                                      event_vad~"Left ventricular assist device",
                                      event_transplant~"Cardiac transplantation",
                                      event_death~"All-cause mortality",
                                      death_cause~"Causes of death",
                                      primary_diagnosis_age~"Age at HCM diagnosis")) %>% 
  gtsummary::bold_labels() #%>%
  #gtsummary::add_p() %>% 
  #gtsummary::as_flex_table() %>% 
  #flextable::save_as_docx(path="outcomes.docx") #Code below saves to image




##################



df_dcc %>%  
  select(sex, primary_diagnosis_age, first_encounter_age, is_proband, nyha, race, 
         genetics,genes_table1,
         lvef_baseline, #event_lvoto50, 
         lvot_obs,event_lvef50,
         fu_lvsd_table1, lad, echo_max_lvt,hcm_ped) %>% 
  mutate(genetics2 = case_when(genetics== "Not tested"~NA_character_,
                               genetics== "Genonegative"~"Non-Sarcomeric HCM",
                               genetics== "Non-Sarcomeric HCM"~"Non-Sarcomeric HCM",
                               genetics== "VUS"~"VUS",
                               T~"Sarcomeric HCM"),
         race = case_when(race=="White"~ "White",
                          race== "Black"~"Black",
                          race=="Asian"~"Asian",
                            T~"Other or Not Reported"),
         genetics2 = fct_infreq(genetics2),
         fu_lvsd_table1 = case_when(fu_lvsd_table1 == "prevalent lvsd"& primary_diagnosis_age>18~"Prevalent LVSD, adult-onset",
                                    fu_lvsd_table1 == "prevalent lvsd"& primary_diagnosis_age<=18~"Prevalent LVSD, pediatric-onset",
                                    primary_diagnosis_age>18 ~"No LVSD, adult-onset",
                                    primary_diagnosis_age<= 18~"No LVSD, pediatric-onset",
                                    T~NA_character_)) %>% 
  select(!event_lvef50) %>% 
  gtsummary::tbl_summary(by = hcm_ped,
                         value = list(),
                         statistic = list(primary_diagnosis_age~"{median} ({p25} to {p75})",
                                          #'Total follow-up, yrs'~"{median} ({p25} to {p75})",
                                          echo_max_lvt~"{median} ({p25} to {p75})",
                                          #  'Max LV wall thickness, z-score'~"{median} ({p25} to {p75})",
                                          lvef_baseline~"{mean} ± {sd}",
                                          #                'LV end-diastolic diameter, mm'~"{mean} ± {sd}",
                                          lad~"{mean} ± {sd}"),
                         digits = list(primary_diagnosis_age~c(1,1),
                                       echo_max_lvt~c(1,1),
                                       #'Max LV wall thickness, z-score'~c(1,1),
                                       lvef_baseline~c(1,1),
                                       #            'LV end-diastolic diameter, mm'~c(1,1),
                                       lad~c(1,1))) %>% #gtsummary::add_overall() %>% 
  gtsummary::bold_labels() #%>% #as_gt() %>% gtsave("not_as_ugly.docx") 
#gtsummary::as_flex_table() %>% 
#  flextable::save_as_docx(path="Not_ugly_table.docx") #Code below saves to image
