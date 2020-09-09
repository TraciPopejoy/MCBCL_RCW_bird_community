
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
