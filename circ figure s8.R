library(tidyverse);library(ggtext);library(scico); library(survival);library(epiR);
library(ggpubr);library(ggthemes);library(survminer)
df1 <- df %>% filter(fu_lvsd_cmp>0)
fit1 <- Surv(df1$echo_age0, df1$t2_lvef50+.001, event = df1$event_lvef50)~1
spl <- as_tibble(survSplit(fit1, data = df1, cut = c(2,8,12,20,30,40), episode = "agegroup")) %>% 
  mutate(time = tstop-tstart) %>% 
  group_by(agegroup) %>% 
  summarise(nyear = sum(time,na.rm = T), ncas = sum(event,na.rm = T)) %>% 
  ungroup() 

dat.df<- spl %>% select(ncas,nyear) %>% 
  as.matrix() %>% 
  epiR::epi.conf(ctype = "inc.rate", method = "byar") %>%  
  bind_cols(spl) %>% 
  mutate(sur = 1-est,
         #Agegroup = factor(Age.group, labels = c("<30", "30-44", "45-59", "\u226560"))
  )
dat.df

p2<-
dat.df %>% tibble() %>%  
  ggplot(aes(x=agegroup, y= 1, 
             label = round(nyear,0)))+geom_text()+
  geom_text(size =3.5, family = "Helvetica")+ 
    theme_void()+
  #geom_segment(aes(x=-.3, xend = -.6, y=1, yend=1, color = "dodgerblue"), show.legend = F,
  #             linewidth = 1.5#color = scico(4,palette = "batlow")
#  )+
  coord_cartesian(xlim= c(.4,7.6), ylim = c(0.8,1.2), expand = F, clip = "off")+
  annotate("text", x= 1, y=1.2, label ="Total years at risk", family="Helvetica", fontface = "bold", hjust =.5)



p1<-
  dat.df %>% 
  mutate(agegroups = factor(agegroup, labels = c("<2", "2-8", "9-12","13-20", "21-30","31-40", ">40")),
         across(.cols = 1:3, ~.x*100)) %>%
  #filter(agegroup<4) %>% 
  ggplot(aes(x= agegroups, y= est, ymin = lower, ymax= upper))+
  geom_errorbar(width = .1)+
  geom_line(aes(x= c(1,2,3,4,5,6,7)))+
  geom_point(shape = 21, size = 4, show.legend = T, fill = "dodgerblue")+
  scale_fill_scico_d(palette = "berlin")+
  scale_y_continuous(breaks = seq(0,20,1))+
  theme_pubclean()+
  labs(y ="Incidence rate of LVSD per 100 person-year",
       x = "Age groups")+
  theme(axis.title.x = element_text(family = "Helvetica", size = 10, hjust = 1),
        axis.title.y = element_text(family = "Helvetica", size = 10, hjust = 1),
        legend.position = "none",
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(family = "Helvetica", size = 10),
        axis.text.y = element_text(family = "Helvetica", size = 10))
p1/p2+plot_layout(heights = c(1,.2))


ggsave("circ s8A.pdf", units = "cm", width = 20, height = 14, dpi = 1200)


fitt1_echo<-
  df %>% mutate(agegroup = case_when(primary_diagnosis_age <=2~1,
                                     #  primary_diagnosis_age <8~2,
                                     primary_diagnosis_age <12~3,
                                     primary_diagnosis_age <=18~4,
                                     #  primary_diagnosis_age <30~5,
                                     #   primary_diagnosis_age <45~6,
                                     #  primary_diagnosis_age >=45~7,
                                     T~NA_real_),
                agegroup = factor(agegroup, labels = c("<2","2-12", "12-18"))) %>% 
  group_by(agegroup) %>% #summarise(n=n()) %>% 
  filter(fu_lvsd>0) %>% 
  surv_fit(Surv(fu_lvsd, event = event_lvef50)~agegroup, data = .) 

s1<-
  ggsurvplot(fitt1_echo,
             xlim = c(0,10), break.x.by = 1,
             censor = F,
             #legend.title = "",
             palette = scico(3, palette = "batlow"),
             fun = "event",
             ylim = c(0,.2),
             break.y.by = .02,
             pval = T,
             legend = "none")
s1

s1s<-
  s1$plot+ theme_pubclean()+ 
  theme(legend.position = "none",
        axis.title.x = element_text(family = "Helvetica", hjust = 1),
        axis.title.y = element_text(family = "Helvetica", hjust = .8),
        axis.text.x = element_text(family = "Helvetica", color = "black"),
        axis.text.y = element_text(family = "Helvetica", color = "black"))+
  annotate("text", x = 3, y= c(.19,.17,.15), label = c("<2", "2-12", "12-18"),
           color = scico(3, palette = "batlow"), family = "Helvetica", fontface = "bold" )+
  
  coord_cartesian(xlim= c(0,10.4), ylim = c(0,0.21), expand = F, clip = "off")+
  labs(#title = "Cumulative incidence of LVSD",
    x= "Years from first SHaRe visit",
    y = "Cumulative incidence of LVSD")
s1s
s2<-
fitt1_echo %>% broom::tidy() %>%
  mutate(fu = case_when(time<1~0,
                        time<2~1,
                        time<3~2,
                        time<4~3,
                        time<5~4,
                        time<6~5,
                        time<7~6,
                        time<8~7,
                        time<9~8,
                        time<10~9,
                        T~10)) %>% 
  group_by(fu, strata) %>% mutate(r=row_number()) %>% 
  ungroup() %>% 
  filter(r==1) %>% select(fu, n.risk, estimate, strata) %>%
  mutate(estimate = 1-estimate) %>% 
  ggplot(aes(x=fu, y= strata, 
             label = n.risk))+
  geom_text(size =3.5, family = "Helvetica")+ theme_void()+scale_color_scico_d(direction=-1, palette = "batlow", end = 1)+
  geom_segment(aes(x=-.3, xend = -.6, y=strata, yend=strata, color = strata), show.legend = F,
               linewidth = 1.5#color = scico(4,palette = "batlow")
  )+
  coord_cartesian(xlim= c(0,10.4), ylim = c(0.4,4.6), expand = F, clip = "on")+
  annotate("text", x= 0, y=4.9, label ="Numbers at risk", family="Helvetica", fontface = "bold", hjust =.5)

s1s/s2+plot_layout(heights = c(1,.2))
ggsave("circ s8B.pdf", units = "cm", width = 20, height = 14, dpi = 1200)

(p1+s1s)/(p2+s2)+plot_layout(heights = c(1,.2,1,.2), widths = c(1,1,1,1),nrow = 2)

(s1s/s2)+plot_layout(heights = c(1,.2,1,.2))
