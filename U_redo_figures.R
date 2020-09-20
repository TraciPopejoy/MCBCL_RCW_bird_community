
library(tidyverse); library(readxl); library(Distance); library(lubridate)
library(cowplot)

#load in the raw data
source('Umbrella/U_load_data.r')

rich_df<-DCERP_filt %>% 
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)]) %>%
  group_by(`Pt name`, Year, USFS_hab_index, Region.Label) %>%
  summarize(ntaxa=length(unique(Species)),
            taxa=paste(unique(Species), collapse=', ')) 
ggplot(rich_df, aes(y=ntaxa))+
  stat_summary(aes(x=Region.Label))+
  scale_x_continuous(name='RCWO Habitat Matrix Score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76))+
  scale_y_continuous(name='Mean species richness')+
  theme_cowplot()
ggsave('Umbrella/results/site_richness_taxa_Fig2.tiff', width=3, height=3)
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
