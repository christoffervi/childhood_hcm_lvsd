

df1<- df_dcc %>% filter(t2_d_htx_vad>first_encounter_age) %>% 
  select(pid, t2_d_htx_vad, d_htx_vad_c, echo_age0,first_encounter_age,
         event_lvef50,t2_lvef50,hcm_ped,
         lvsd_timing,sex,genetics, primary_diagnosis_age) %>% 
  mutate(fu_comp_event = t2_d_htx_vad,
  )

# Merge data to create time-dependent covariates
# This section uses the 'tmerge' function from the survival package to merge the 'df1' data frame with itself and create time-dependent covariates
df1_new<- 
  survival::tmerge(df1, df1, id= pid, 
                   lvef50 = tdc(t2_lvef50, event_lvef50,0),
                   endpt = event(t2_d_htx_vad, d_htx_vad_c),
                   
                   tstart = first_encounter_age, tstop = fu_comp_event) %>% 
  mutate(grouping = case_when(lvef50==1 & hcm_ped==1~"Adult-Dx HCM, LVSD",
                              lvef50==0 & hcm_ped==1~"Adult-Dx HCM, no LVSD",
                              lvef50==1 & hcm_ped==0~"Childhood-Dx HCM, LVSD",
                              lvef50==0 & hcm_ped==0~"Childhood-Dx HCM, no LVSD",
                              T~NA_character_),
         time = tstop-tstart)


f1<- surv_fit(Surv( time,endpt)~grouping, data = df1_new)
x1<-
  ggsurvplot(f1, fun = "event",
             xlim = c(0,10),
             conf.int = FALSE,
             risk.table = T,
             break.x.by = 1,
             censor=F,
             legend.labs = c("Adult-Dx HCM, LVSD","Adult-Dx HCM, no LVSD","Childhood-Dx HCM, LVSD","Childhood-Dx HCM, no LVSD")
  )
x1p<-
  x1$plot+scale_x_continuous(breaks = seq(0,20,1))+
  scale_y_continuous(breaks = seq(0,1,.05))+
  scale_fill_scico_d(direction=-1, palette = "batlow", end = 1)+
  scale_color_scico_d(direction=-1, palette = "batlow", end = 1)+
  labs(x= "Years from initial SHaRe visit or LVSD",
       y= "Cumulative Incidence of Death, Txp or LVAD")+ 
  theme(axis.title.x = element_text(family = "Helvetica", size = 12, hjust = 1),
        axis.title.y = element_text(family = "Helvetica", size = 12, hjust = 1),
        legend.position = c(0.8,.3),
        legend.text = element_text(family = "Helvetica", size = 8, hjust = 0),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(family = "Helvetica", size = 10),
        axis.text.y = element_text(family = "Helvetica", size = 10),
        plot.background = element_rect(color = "white"))+
  coord_cartesian(xlim = c(0,10.1),
                  ylim =c(0,.62),
                  expand = F)


x2<-
  f1 %>% broom::tidy() %>%
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
  coord_cartesian(xlim= c(0,10.1), ylim = c(0.4,4.6), expand = F, clip = "off")+
  annotate("text", x= 0, y=4.9, label ="Numbers at risk", family="Helvetica", fontface = "bold", hjust =.5)
x1p/x2+plot_layout(heights = c(1,.25))

ggsave("circ figure 5.tiff", compression = "lzw",units = "cm", width = 20, height = 14, dpi = 1200)
ggsave("circ figure 5.pdf", units = "cm", width = 20, height = 14, dpi = 1200)
