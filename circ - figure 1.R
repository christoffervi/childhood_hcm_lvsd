##########
#Diagnosis

# Use the 'npsurv' function from the 'rms' package to estimate the cumulative incidence of LVSD from diagnosis 
# in individuals aged 18 and below
f1 <- rms::npsurv(Surv(t2_lvef50 - primary_diagnosis_age, factor(cmpr_event)) ~ 1, data = filter(cr_df_18, primary_diagnosis_age <= 18))

# Filter the 'f1' data and create a new column 'fu' which represents the follow-up time intervals
# Then group the data by 'fu' and 'state', create a new column 'r' with row numbers, and filter the data 
# to keep only the first row for each group
# This creates the 'p1_annotate' data frame that is used for annotating the plot
f1<- rms::npsurv(Surv(t2_lvef50-primary_diagnosis_age, factor(cmpr_event))~1, data = filter(cr_df_18, primary_diagnosis_age<=18))

p1_annotate<-
  f1 %>% broom::tidy() %>%filter(state %in% c("(s0)" ,"1")) %>%  
  mutate(fu = case_when(time<5~0,time<10~5,time<15~10,time<20~15,
                        time<25~20,time<30~25,T~30)) %>% 
  group_by(fu, state) %>% mutate(r=row_number()) %>% ungroup() %>% 
  filter(r==1) %>%# select(fu, n.risk, estimate, state) %>%
  mutate(estimate = estimate)

# Create the plot 'p1' using the 'ggplot' function from the 'ggplot2' package
# Use the 'geom_line' and 'geom_ribbon' functions to create the main plot, 
# representing the cumulative incidence of LVSD from diagnosis over time
# Use the 'geom_segment' function to add vertical and horizontal lines to the plot indicating the 
# time intervals and estimated probabilities
p1<-
  f1 %>% broom::tidy() %>%filter(state==1) %>% 
  ggplot(aes(x=time, y=estimate))+
  geom_line(color = scico(1,palette = "lapaz",begin = .4),
  )+
  geom_ribbon(aes(x=time, y=estimate, ymin=conf.low, ymax=conf.high), alpha = .5,
              fill = scico(1,palette = "lapaz",begin = .4)
  )+
  geom_line(aes(x=time, y=estimate), alpha = 1,
            #color = scico(1,palette = "lapaz",,begin = .4),
  )+
  ggthemes::theme_clean()+
  coord_cartesian(xlim= c(0,30), ylim = c(0,.3), expand = F)+
  scale_x_continuous(breaks = seq(0,30,5))+
  scale_y_continuous(breaks = seq(0,1,.05))+
  scale_fill_scico_d(palette = "berlin", end = .9)+
  scale_color_scico_d(palette = "berlin", end = .9)+
  labs(x= "Years from diagnosis of HCM",
       y= "Cumulative incidence of LVSD")+ 
  theme(axis.title.x = element_text(family = "Roboto", size = 12, hjust = 1),
        axis.title.y = element_text(family = "Roboto", size = 12, hjust = 1),
        legend.position = c(0.9,.2),
        legend.text = element_text(family = "Roboto", size = 10, hjust = 1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(family = "Roboto", size = 10),
        axis.text.y = element_text(family = "Roboto", size = 10),
        plot.background = element_rect(color = "white"),
        panel.grid.major.y = element_blank())+
  geom_segment(data = p1_annotate %>%  filter(fu>0,fu<30, state=="1"), 
               aes(x=fu, xend=fu, y=0, yend= estimate), color = "black", linetype = 3)+
  geom_segment(data = p1_annotate %>%  filter(fu>0,fu<30, state=="1"), 
               aes(x=fu, xend=0, y=estimate, yend= estimate), color = "black", linetype = 3)+
  annotate("text", 
           x= #filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(fu)%>% pull()-4.5, 
             c(0.5,4,8,13,18),
           y= filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(estimate) %>% pull()+.01, 
           label = paste(sep="",  
                         #filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(fu) %>% pull, 
                         #" year incidence\n ", 
                         round(filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(estimate) %>% pull,3)*100,
                         "% ",
                         "[CI: ",
                         round(filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(conf.low) %>% pull,3)*100,
                         "-",
                         round(filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(conf.high) %>% pull,3)*100,
                         "]"
           ),
           hjust=0, vjust=0, family = "Roboto", fontface="bold", size= 3)
p2<-
  f1 %>% broom::tidy() %>%filter(state=="(s0)") %>%  
  mutate(fu = case_when(time<5~0,
                        time<10~5,
                        time<15~10,
                        time<20~15,
                        time<25~20,
                        time<30~25,
                        T~30)) %>% 
  group_by(fu) %>% mutate(r=row_number()) %>% 
  ungroup() %>% 
  filter(r==1) %>% select(fu, n.risk, estimate) %>%
  mutate(estimate = 1-estimate) %>% 
  ggplot(aes(x=fu, y= 1, 
             label = n.risk))+
  geom_text()+ theme_void()+
  geom_segment(aes(x=-1.6, xend = -1, y=1, yend=1), color = scico(1,palette = "lapaz",begin = .4))+
  coord_cartesian(xlim= c(0,30), ylim = c(0.4,1.6), expand = F, clip = "off")+
  annotate("text", x= -1.5, y=1.6, label ="Numbers at risk", family="Roboto", fontface = "bold", hjust =0)

p1/p2+plot_layout(heights = c(5,1))





##############
# AGE 
# in the following the same is done using age as the time-scale and  time since first SHaRe visit
q1<- rms::npsurv(Surv(t2_lvef50, factor(cmpr_event))~1, data = filter(cr_df_18, primary_diagnosis_age<=18))

q1_annotate<-
  q1 %>% broom::tidy() %>%filter(state %in% c("(s0)" ,"1")) %>%  
  mutate(fu = case_when(time<5~0,time<10~5,time<15~10,time<20~15,
                        time<25~20,time<30~25,time<35~30,time<40~35,time<45~40,time<50~45,T~50)) %>% 
  group_by(fu, state) %>% mutate(r=row_number()) %>% ungroup() %>% 
  filter(r==1) %>%# select(fu, n.risk, estimate, state) %>%
  mutate(estimate = estimate)


q1<-
  q1 %>% broom::tidy() %>%filter(state==1) %>% 
  ggplot(aes(x=time, y=estimate))+
  geom_line(color = scico(1,palette = "lapaz",begin = .4),
  )+
  geom_ribbon(aes(x=time, y=estimate, ymin=conf.low, ymax=conf.high), alpha = .5,
              fill = scico(1,palette = "lapaz",begin = .4)
  )+
  geom_line(aes(x=time, y=estimate), alpha = 1,
            #color = scico(1,palette = "lapaz",,begin = .4),
  )+
  ggthemes::theme_clean()+
  coord_cartesian(xlim= c(0,50), ylim = c(0,.4), expand = F)+
  scale_x_continuous(breaks = seq(0,50,5))+
  scale_y_continuous(breaks = seq(0,1,.05))+
  scale_fill_scico_d(palette = "berlin", end = .9)+
  scale_color_scico_d(palette = "berlin", end = .9)+
  labs(x= "Age in Years",
       y= "Cumulative incidence of LVSD")+ 
  theme(axis.title.x = element_text(family = "Roboto", size = 12, hjust = 1),
        axis.title.y = element_text(family = "Roboto", size = 12, hjust = 1),
        legend.position = c(0.9,.2),
        legend.text = element_text(family = "Roboto", size = 10, hjust = 1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(family = "Roboto", size = 10),
        axis.text.y = element_text(family = "Roboto", size = 10),
        plot.background = element_rect(color = "white"),
        panel.grid.major.y = element_blank())+
  geom_segment(data = q1_annotate %>%  filter(fu %in% c(15,25,35,45), state=="1"), 
               aes(x=fu, xend=fu, y=0, yend= estimate), color = "black", linetype = 3)+
  geom_segment(data = q1_annotate %>%  filter(fu %in% c(15,25,35,45), state=="1"), 
               aes(x=fu, xend=0, y=estimate, yend= estimate), color = "black", linetype = 3)+
  annotate("text", 
           x= #filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(fu)%>% pull()-4.5, 
             c(5,13,23,33),
           # 1,
           y= filter(q1_annotate, state=="1",  fu %in% c(15,25,35,45)) %>% select(estimate) %>% pull()+.01, 
           label = paste(sep="",  
                         #filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(fu) %>% pull, 
                         #" year incidence\n ", 
                         round(filter(q1_annotate, state=="1",   fu %in% c(15,25,35,45)) %>% select(estimate) %>% pull,3)*100,
                         "% ",
                         "[CI: ",
                         round(filter(q1_annotate, state=="1",   fu %in% c(15,25,35,45)) %>% select(conf.low) %>% pull,3)*100,
                         "-",
                         round(filter(q1_annotate, state=="1",   fu %in% c(15,25,35,45)) %>% select(conf.high) %>% pull,3)*100,
                         "]"
           ),
           hjust=0, vjust=0, family = "Roboto", fontface="bold", size= 3)
q2<-
  f1 %>% broom::tidy() %>%filter(state=="(s0)") %>%  
  mutate(fu = case_when(time<5~0,
                        time<10~5,
                        time<15~10,
                        time<20~15,
                        time<25~20,
                        time<30~25,
                        time<35~30,
                        time<40~35,
                        time<45~40,
                        time<50~45,
                        T~50)) %>% 
  group_by(fu) %>% mutate(r=row_number()) %>% 
  ungroup() %>% 
  filter(r==1) %>% select(fu, n.risk, estimate) %>%
  mutate(estimate = 1-estimate) %>% 
  ggplot(aes(x=fu, y= 1, 
             label = n.risk))+
  geom_text()+ theme_void()+
  geom_segment(aes(x=-1.6*1.67, xend = -1*1.67, y=1, yend=1), color = scico(1,palette = "lapaz",begin = .4))+
  coord_cartesian(xlim= c(0,50), ylim = c(0.4,1.6), expand = F, clip = "off")+
  annotate("text", x= -1.5, y=1.6, label ="Numbers at risk", family="Roboto", fontface = "bold", hjust =0)

q1/q2+plot_layout(heights = c(5,1))


############
####original
f1<- rms::npsurv(Surv(t2_lvef50-echo_age0, factor(cmpr_event))~1, data = filter(cr_df_18, primary_diagnosis_age<=18))

s1_annotate<-
  f1 %>% broom::tidy() %>%filter(state %in% c("(s0)" ,"1")) %>%  
  mutate(fu = case_when(time<5~0,time<10~5,time<15~10,time<20~15,
                        time<25~20,time<30~25,T~30)) %>% 
  group_by(fu, state) %>% mutate(r=row_number()) %>% ungroup() %>% 
  filter(r==1) %>%# select(fu, n.risk, estimate, state) %>%
  mutate(estimate = estimate)


s1<-
  f1 %>% broom::tidy() %>%filter(state==1) %>% 
  ggplot(aes(x=time, y=estimate))+
  geom_line(color = scico(1,palette = "lapaz",begin = .4),
  )+
  geom_ribbon(aes(x=time, y=estimate, ymin=conf.low, ymax=conf.high), alpha = .5,
              fill = scico(1,palette = "lapaz",begin = .4)
  )+
  geom_line(aes(x=time, y=estimate), alpha = 1,
            #color = scico(1,palette = "lapaz",,begin = .4),
  )+
  ggthemes::theme_clean()+
  coord_cartesian(xlim= c(0,30), ylim = c(0,.5), expand = F)+
  scale_x_continuous(breaks = seq(0,30,5))+
  scale_y_continuous(breaks = seq(0,1,.05))+
  scale_fill_scico_d(palette = "berlin", end = .9)+
  scale_color_scico_d(palette = "berlin", end = .9)+
  labs(x= "Years from first SHaRe visit",
       y= "Cumulative incidence of LVSD")+ 
  theme(axis.title.x = element_text(family = "Roboto", size = 12, hjust = 1),
        axis.title.y = element_text(family = "Roboto", size = 12, hjust = 1),
        legend.position = c(0.9,.2),
        legend.text = element_text(family = "Roboto", size = 10, hjust = 1),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(family = "Roboto", size = 10),
        axis.text.y = element_text(family = "Roboto", size = 10),
        plot.background = element_rect(color = "white"),
        panel.grid.major.y = element_blank())+
  geom_segment(data = s1_annotate %>%  filter(fu>0,fu<30, state=="1"), 
               aes(x=fu, xend=fu, y=0, yend= estimate), color = "black", linetype = 3)+
  geom_segment(data = s1_annotate %>%  filter(fu>0,fu<30, state=="1"), 
               aes(x=fu, xend=0, y=estimate, yend= estimate), color = "black", linetype = 3)+
  annotate("text", 
           x= #filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(fu)%>% pull()-4.5, 
             c(0.5,4,8,13,18),
           y= filter(s1_annotate, state=="1", fu!=0, fu!=30) %>% select(estimate) %>% pull()+.01, 
           label = paste(sep="",  
                         #filter(p1_annotate, state=="1", fu!=0, fu!=30) %>% select(fu) %>% pull, 
                         #" year incidence\n ", 
                         round(filter(s1_annotate, state=="1", fu!=0, fu!=30) %>% select(estimate) %>% pull,3)*100,
                         "% ",
                         "[CI: ",
                         round(filter(s1_annotate, state=="1", fu!=0, fu!=30) %>% select(conf.low) %>% pull,3)*100,
                         "-",
                         round(filter(s1_annotate, state=="1", fu!=0, fu!=30) %>% select(conf.high) %>% pull,3)*100,
                         "]"
           ),
           hjust=0, vjust=0, family = "Roboto", fontface="bold", size= 3)
s2<-
  f1 %>% broom::tidy() %>%filter(state=="(s0)") %>%  
  mutate(fu = case_when(time<5~0,
                        time<10~5,
                        time<15~10,
                        time<20~15,
                        time<25~20,
                        time<30~25,
                        T~30)) %>% 
  group_by(fu) %>% mutate(r=row_number()) %>% 
  ungroup() %>% 
  filter(r==1) %>% select(fu, n.risk, estimate) %>%
  mutate(estimate = 1-estimate) %>% 
  ggplot(aes(x=fu, y= 1, 
             label = n.risk))+
  geom_text()+ theme_void()+
  geom_segment(aes(x=-1.6, xend = -1, y=1, yend=1), color = scico(1,palette = "lapaz",begin = .4))+
  coord_cartesian(xlim= c(0,30), ylim = c(0.4,1.6), expand = F, clip = "off")+
  annotate("text", x= -1.5, y=1.6, label ="Numbers at risk", family="Roboto", fontface = "bold", hjust =0)



s1

s1+s2+p1+p2+q1+q2+plot_layout(ncol = 1, heights = c(5,1,5,1,5,1))+ plot_annotation(tag_levels = list(c("A","", "B","","C","")))

ggsave("figure 1_update.tiff", compression = "lzw", units = "cm", height = 32, width = 18, dpi = 900)