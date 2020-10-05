library(tidyverse); library(readxl); library(Distance); library(lubridate)
library(cowplot)

#load in the raw data
source('Umbrella/U_load_data.r')

# Figure 2 -----
early_obs<-DCERP_filt %>%
  group_by(Species,Year) %>%
  arrange(Date1)%>%
  slice(1) %>%
  mutate(julianday=case_when(Year=='2009'~round(julian(Date1, origin='2009-01-01'),0),
                             Year=='2010'~round(julian(Date1, origin='2010-01-01'),0))) %>%
  dplyr::select(Species, Date1, Year, julianday) 

ex_migrants <- early_obs %>%
  filter(julianday > 140)
early_obs %>% filter(Species %in% ex_migrants$Species)

rich_df<-DCERP_filt %>% 
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)]) %>%
  filter(!(Species %in% c("LAGU", "TUVU"))) %>%
  group_by(`Pt name`, Year, USFS_hab_index, Region.Label) %>%
  summarize(ntaxa=length(unique(Species)),
            taxa=paste(unique(Species), collapse=', ')) %>%
  filter(!is.na(Region.Label))
ggplot(rich_df, aes(y=ntaxa))+
  stat_summary(aes(x=as.character(Region.Label)))+
  scale_x_discrete(name='RCW Habitat Score')+
  scale_y_continuous(name='Species Richness')+
  theme_cowplot()
ggsave('Umbrella/results/Figures/Fig2_site_richness_taxa.tiff', width=3, height=3)
rich_df %>% filter(is.na(Region.Label))
DCERP_filt %>% filter(`Pt name`=='PHOTO63')


early_obs<-DCERP_filt %>%
  group_by(Species,Year) %>%
  arrange(Date1)%>%
  slice(1) %>%
  dplyr::select(Species, Date1, Year)
ggplot()+
  geom_text(data=early_obs, aes(x=1, y=Date1, label=Species),
            size=3, position='jitter')+
  facet_wrap(~Year, scales='free')+
  theme_cowplot()
ggsave('bird_first_obs_plot.jpg')
library(ggrepel)
ggplot()+
  geom_density(data=DCERP_filt, aes(x=Date1, fill=Species),
               alpha=0.3)+
  geom_text(data=early_obs, aes(x=Date1, y=.000001, label=Species),
            size=3, position=position_jitter(height = .0000003))+
  scale_fill_discrete(guide=F)+
  facet_wrap(~Year, scales='free')+
  theme_cowplot()+
  coord_flip()
ggsave('bird_sighting_timing.jpg', height=8, width=6)
write.csv(early_obs,
          'earliest_observation_each_year.csv')

ggplot()+
  geom_density(data=early_obs, aes(x=Date1), fill='goldenrod2', alpha=0.3)+
  geom_text(data=early_obs, aes(x=Date1, y=0.0000006, label=Species),
            position=position_jitter(height=.000002),
            size=3)+
  theme_cowplot()+
  facet_wrap(~Year, scales='free')+
  theme(axis.text.x = element_text(angle=30, hjust=0.9))
ggsave('bird_early_observation_plot.jpg', height=6, width=6)


# Figure 3 ----
DCERP_filt %>% 
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)]) %>%
  filter(Species=='RCWO') %>%
  group_by(`Pt name`, Replicate, Region.Label) %>%
  tally() %>% arrange(desc(Region.Label))

rcwo_summary<-summary(RCWO.mod47)
rcwo_raw_ab<-rcwo_summary$dht$individuals$bysample %>%
  mutate(`Pt name`=gsub('\\..','', Sample))%>%
  left_join(point_sample_empty %>% dplyr::select(`Pt name`, USFS_hab_index),
            by='Pt name') %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)]) 
rcwo_raw_ab %>%
  group_by(Region.Label) %>%
  summarize(mean(Nhat),
            mean(Nchat),
            mean(Dhat))

ggplot(data=rcwo_raw_ab)+
  #  geom_point(aes(x=USFS_hab_index, y=Nchat, color=as.character(Region.Label)), 
  #             alpha=0.3)+
  stat_summary(aes(x=as.character(Region.Label),
                   y=Nchat))+
  scale_y_continuous(name=expression("RCW per km"^-2),
                     breaks=0:5)+
  scale_x_discrete('USFWS Habitat index')+
  scale_color_viridis_d(guide=F, option='plasma')+
  theme_cowplot()
ggsave('Umbrella/results/fig3_rcwabun.tiff',
       width=3, height=3)


# Figure 4 ----
all_models %>% filter(species=="BACS") %>%
  arrange(AICc) %>% slice(1)
all_models %>% filter(species=="BACS") %>%
  arrange(desc(gof_cvm_p)) %>% slice(1)

bacs_summary<-summary(BACS.mod256)
bacs_raw_ab<-bacs_summary$dht$individuals$bysample %>%
  mutate(`Pt name`=gsub('\\..','', Sample))%>%
  left_join(point_sample_empty %>% dplyr::select(`Pt name`, USFS_hab_index),
            by='Pt name') %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)]) 
bacs_raw_ab %>%
  group_by(Region.Label) %>%
  summarize(mean(Nhat),
            mean(Nchat),
            mean(Dhat))

bacs_ab_plot<-ggplot(data=bacs_raw_ab)+
  #  geom_point(aes(x=USFS_hab_index, y=Nchat, color=as.character(Region.Label)), 
  #             alpha=0.3)+
  stat_summary(aes(x=as.character(Region.Label),
                   y=Nchat))+
  scale_y_continuous(name=expression("Birds per km"^-2),
                     breaks=0:5)+
  scale_x_discrete('USFWS Habitat index')+
  scale_color_viridis_d(guide=F, option='plasma')+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=30, hjust=.9))

### Presence Bachman's -----
bacs_pres<-bind_rows(read.csv('Umbrella/presence_extracted_from_figure/BACS_2009.csv',
                              header=F) %>%
                       mutate(Year='2009'),
                     read.csv('Umbrella/presence_extracted_from_figure/BACS_2010.csv',
                              header=F) %>%
                       mutate(Year='2010')) %>%
  dplyr::select(-V3, -V4) %>%
  filter(!is.na(V1)) %>%
  rename(USFS_hab_index='V1',
         prob='V2') %>%
  mutate(USFS_hab_index=round(USFS_hab_index, 2),
         prob=round(prob, 5))
bacs_pres_se<-bind_rows(read.csv('Umbrella/presence_extracted_from_figure/BACS_2009.csv',
                              header=F) %>%
                       mutate(Year='2009'),
                     read.csv('Umbrella/presence_extracted_from_figure/BACS_2010.csv',
                              header=F) %>%
                       mutate(Year='2010')) %>%
  dplyr::select(V3,V4, Year) %>%
  filter(!is.na(V3)) %>%
  mutate(USFS_hab_index=round(V3, 1)) %>%
  group_by(USFS_hab_index, Year) %>%
  summarize(ll=min(V4),
            ul=max(V4),
            difff=ul-ll) %>%
  arrange(Year, USFS_hab_index) %>%
  filter(difff!=0) %>%
  group_by(Year) %>%
  mutate(lag.ll1=lag(ll),
         lead.ll1=lead(ll),
         lag.ul1=lag(ul),
         lead.ul1=lead(ul),
         lag.ll2=lag(ll,2),
         lead.ll2=lead(ll,2),
         lag.ul2=lag(ul,2),
         lead.ul2=lead(ul,2)) %>%
  rowwise() %>%
  mutate(new.ll=mean(c(ll,lag.ll1, lead.ll1), na.rm=T),
         new.ul=mean(c(ul,lag.ul1, lead.ul1), na.rm=T)) %>%
  ungroup() %>%
  slice(seq(1, n(), 2))

head(bacs_pres_se)
bacs_pres_plot<-ggplot()+
  geom_ribbon(data=bacs_pres_se,
              aes(x=USFS_hab_index, ymin=new.ll, ymax=new.ul, group=Year),
              fill='lightgrey', alpha=0.7)+
  geom_smooth(data=bacs_pres,
              aes(x=USFS_hab_index, y=prob, color=Year),
            size=1.5, level=0)+
  scale_color_manual(values=c('black', 'darkgrey'))+
  scale_y_continuous('Proportion\nof sites occupied',
                     breaks=seq(0,1,.2),
                     limits=c(0,1))+
  scale_x_continuous('USFWS Habitat index', 
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76))+
  theme_cowplot()+
  theme(legend.position=c(.1, .8),
        axis.text.x = element_text(angle=30, hjust=.9))
        
plot_grid(bacs_pres_plot, bacs_ab_plot, 
          labels="AUTO")
ggsave('Umbrella/results/Figures/Fig3_bacs.tiff', width=6, height=3)
