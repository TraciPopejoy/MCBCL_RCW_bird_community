library(tidyverse); library(cowplot)
abundance_summary <- read.csv('Umbrella/results/abundance summary 95 confidence.csv')

# PDF of bird abundance trend & relationship with habitat ----
abun_gdf <- abundance_summary %>% 
  left_join(bird.assem, by=c('species'='Species code')) %>%
  left_join(reg_sig) %>%
  left_join(booted_reg_sum)
View(abun_gdf)

pdf('Umbrella/results/spp_abundance_by_habitat.pdf')
for(j in unique(bird.assem$Habitat.group)){
  ba_plot <- abun_gdf %>%
    filter(Habitat.group==j) %>%
    ggplot()+
    geom_linerange(aes(x=HabCatfix, ymin=x2.5, ymax=x97.5))+
    geom_point(aes(x=HabCatfix, y=x50))+
    geom_text(data=reg_sig %>% filter(Habitat.group==j), 
              aes(x=1, y=-0.3, label=round(qt0,3), color=new_sig),
              hjust=0) +
    facet_wrap(~species, scales='free_y')+
    scale_y_continuous(name=expression('Birds per km'^2))+
    scale_x_discrete('RCW Habitat Quality Index')+
    scale_color_manual(values=c('black','red'), guide=F)+
    ggtitle(paste(j, 'associated species'))+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle=20, size=10, hjust=.9))
  print(ba_plot)
}
dev.off()

# Values to help write results -----
# info on RCW and Bachmann's Sparrow
abun_gdf %>%
  left_join(reg_sig) %>%
  filter(species %in% c('RCWO', 'BACS'))%>%
  select(`Mean abundance birds/km`, Hab.coeff, phi, qt0)
booted_reg %>% filter(species %in% c('RCWO', 'BACS'))

abun_gdf %>%
  filter(Habitat.group=='longleaf',
         !(species %in% c('BHNU','CHSP','PIWA','PRAW')),
         HabCatfix==2.5) %>%
  ungroup()%>%
  summarize(median(x50))

abun_gdf %>%
  filter(Habitat.group=='shrub',
         HabCatfix %in% c(2.5, 4.76)) %>%
  select(species, HabCatfix, x50, `Common name`, phi) %>%
  View()
  
abun_gdf %>% filter(Habitat.group=='hardwood') %>%
  select(phi)

abun_gdf %>% filter(Habitat.group=='shrub') 

# fuzzy ordination to assess alignment with RCW habitat index ----
library(fso); library(vegan); library(cowplot)
abundance_good<-read.csv('Umbrella/results/abundance_used_quant.csv')

#build a site x community matrix based on bird densities
ordination_df<-abundance_good %>% 
  left_join(sampled_area_km) %>%
  # hand calculate Nchat and Density
  mutate(PtName=gsub('\\..*','', old.label),
         nchat=tob/phi,
         NEWnhat=nchat/samp.area,
         raw_dens_no_phi=tob/samp.area, 
         HabCatfix=case_when(substr(HabCat,4,4)==5~sub('5','6',HabCat),
                             T~as.character(HabCat))) %>%
  replace_na(list(NEWnhat=0)) %>%
  # get the 50% quantile for each species at each point
  group_by(species, PtName, HabCatfix, USFS_hab_index) %>%
  summarize(bdens_m=quantile(NEWnhat, probs=.5)) %>%
  # create a community matrix
  pivot_wider(names_from='species', values_from='bdens_m',
              values_fill=0) %>%
  dplyr::select(PtName, HabCatfix, USFS_hab_index, #reordering
                RCWO, 
                BACS, EAWP, MODO, # might want to exclude these because models are poor fit for reality
                everything()) 
View(ordination_df)

site.com<-ordination_df[,-c(2:4)] %>% #excluding RCW as we know they should be correlated with them 
  bind_cols(., total.dens=rowSums(.[,-1])) %>% #select(total.dens) %>% arrange(desc(total.dens))
  filter(total.dens != 0) %>%
  dplyr::select(-total.dens) %>%
  ungroup()
nrow(site.com); head(site.com)
dis<-vegdist(site.com[,-1], method="bray") #create a distance matrix for the bird communities

site.env<-ordination_df[,c(1:4)] 
which(site.env$PtName!=site.com$PtName)

#run the fuzzy ordination
fuzzy_ord<-fso(~USFS_hab_index, dis, site.env, permute=10000,
               scaling=3)
#mu - fuzzy membership values for each site
#r - correlation between original vector and fuzzy set
#d - correlation of pair-wise distances among each fuzzy set
summary(fuzzy_ord)
plot(fuzzy_ord)
ord_gdf<- ordination_df %>%
  bind_cols(mu=fuzzy_ord$mu[,1]) %>%
  dplyr::select(PtName, mu, everything())

# Figure 7 fit of fuzzy ordination plot ----
ggplot()+
  geom_point(data=ord_gdf, aes(x=USFS_hab_index, y=mu), alpha=0.6)+
  geom_text(aes(x=2.5, y=0.545, label="r = 0.491\np < 0.001"))+
  #ggtitle('Fit of USFS Habitat Index ordination')+
  scale_y_continuous('mu')+
  scale_x_continuous('RCW habitat score')+
  theme_cowplot()
ggsave('Umbrella/results/Figures/Fig7_fuzzy_ord_fit.tiff',
       width=3, height=3)

# Figure 8 ordination plot ----
hab_breaks<-ord_gdf %>% arrange(USFS_hab_index) %>%
  pull(USFS_hab_index) %>% round(2)
ord_gdf_long <-ord_gdf %>%
  arrange(mu)%>%
  rowid_to_column('mu.order') %>%
  arrange(USFS_hab_index) %>%
  rowid_to_column('hab.order') %>%
  pivot_longer(cols=c(BACS:YTWA), 
               names_to="species", values_to="Birds/sqkm")%>%
  left_join(bird.assem, by=c('species'='Species code')) %>%
  mutate(HabitatF=factor(Habitat.group, 
                         levels=c('generalist', 'shrub','hardwood','longleaf'))) %>%
  group_by(PtName, USFS_hab_index) %>%
  mutate(rel.abundance=(`Birds/sqkm`/sum(`Birds/sqkm`)))

ord_gdf_long %>% summarize(tper=sum(rel.abundance)) %>% 
  arrange(tper) %>% ungroup() %>% slice(1:5, 213:218)

#plot with relative abundance
ggplot(ord_gdf_long)+
  geom_col(aes(x=hab.order, y=rel.abundance, 
               fill=HabitatF))+
  scale_fill_viridis_d(name="Habitat\nassociation")+
  scale_x_continuous(name='RCW habitat score',
                     breaks=c(seq(1,218,25),218),
                     labels=hab_breaks[c(seq(1,218,25),218)],
                     expand = c(0,0))+
  scale_y_continuous(name='Relative Abundance',
                     expand=c(0,0),
                     labels=scales::percent)+
  theme(axis.text.x = element_text(angle=90, hjust=0.9),
        axis.text.y = element_text(angle=90, hjust=0.5))+
  coord_flip()+theme_cowplot()
ggsave('Umbrella/results/Figures/Fig8_habitat_sitecom_ra.tiff',
       width=5, height=7.5)

#plot with bird density
ggplot(ord_gdf_long)+
  geom_col(aes(x=hab.order, y=`Birds/sqkm`, 
               fill=HabitatF))+
  scale_fill_viridis_d(name="Habitat\nassociation")+
  scale_x_continuous(name='RCW habitat score',
                     breaks=c(seq(1,218,25),218),
                     labels=hab_breaks[c(seq(1,218,25),218)],
                     expand = c(0,0))+
  scale_y_continuous(name=expression('Birds / km'^2),
                     expand=c(0,0))+
  theme(axis.text.x = element_text(angle=90, hjust=0.9),
        axis.text.y = element_text(angle=90, hjust=0.5))+
  coord_flip()+theme_cowplot()
ggsave('Umbrella/results/Figures/Fig8_habitat_sitecom_birddens.tiff',
       width=5, height=7.5)

# Figure 5 Total Abundance Box Plots ----
hab_abun_plot<-ord_gdf_long %>%
  group_by(PtName, HabCatfix, USFS_hab_index, HabitatF) %>%
  summarize(Total.Birds.sqkm=sum(`Birds/sqkm`), .groups='keep') %>%
  ggplot()+
  geom_boxplot(aes(x=HabCatfix, y=Total.Birds.sqkm, group=HabCatfix),
               fill='lightgrey')+
  scale_x_discrete(name='RCW habitat score')+
  theme_cowplot()+facet_wrap(~HabitatF, scales='free_y')+
  theme(axis.title.y=element_text(size=0),
        axis.text.x = element_text(size=9, angle=30, hjust=.9))
abund_plot<-ord_gdf_long %>% 
  group_by(PtName, HabCatfix, USFS_hab_index,) %>%
  summarize(Total.Birds.sqkm=sum(`Birds/sqkm`), .groups='keep') %>%
  ggplot(aes(x=HabCatfix,  
                   y=Total.Birds.sqkm, group=HabCatfix))+
  geom_jitter(alpha=0.2)+
  geom_boxplot(outlier.alpha = 0, alpha=0.8, fill='lightgrey')+
  scale_x_discrete(name='RCW habitat score')+
  scale_y_continuous(expression("Total Birds km"^-2))+
  theme_cowplot()+
  theme(axis.text.x = element_text(size=11, angle=30, hjust=.9))
plot_grid(abund_plot, hab_abun_plot, labels=c("       All Birds", ""),
          rel_widths = c(0.43,0.57))
ggsave('Umbrella/results/Figures/Fig5_bird_abundance.tiff', width=6, height=4)

write.csv(ord_gdf_long, 'Umbrella/results/plotting_df.csv')

# Table S4 and Figure S2 detection information from bootstraps -----
abundance_good<-read.csv('Umbrella/results/abundance_used_quant.csv')
# want a table of species x habitat category detection
names(abundance_good)
unique(abundance_good$HabCat)
phi.table<-abundance_good %>% 
  mutate(HabCatfix=case_when(substr(HabCat,4,4)==5~sub('5','6',HabCat),
                              T~as.character(HabCat))) %>%
  group_by(species, HabCatfix) %>%
  summarize(m.phi=mean(phi)) %>%
  pivot_wider(names_from=HabCatfix, values_from=m.phi) %>%
  mutate(w.m.phi=weighted.mean(c(`2.5`,`3.26`,`3.76`,`4.26`,`4.76`),
                               c(32,34,51,55,46)))

phi.table %>% 
  left_join(bird.assem, by=c('species'='Species code')) %>%
  arrange(Habitat.group) %>%
  select(-`Scientific name`, -species) %>%
  write.csv('Umbrella/results/TableS4_dist_detect_habcat.csv')

model.appendix.all<-read.csv('Umbrella/model.appendix.20201116.csv') # from U_
sp_det_hab<-model.appendix.all %>% group_by(species) %>%
  arrange(AIC) %>% slice(1) %>%
  filter(grepl('Hab', formula)) %>% pull(species)

phi.table %>%
  filter(species %in% sp_det_hab) %>%
  select(-w.m.phi) %>%
  pivot_longer(-species) %>%
  mutate(rcw_score=as.numeric(name)) %>%
  ggplot()+
  geom_line(aes(x=rcw_score, y=value, group=species,
                color=species), size=2, alpha=0.6)+
  geom_text(data=. %>% filter(name=='4.76'), 
                  aes(x=4.95, y=value, label=species))+
  scale_y_continuous('Detection Probability')+
  scale_x_continuous('RCW habitat score', 
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76),
                     expand = expansion(mult=c(0.1, .17)))+
  scale_color_viridis_d(option = 'inferno')+
  theme_cowplot()+
  theme(legend.position = 'none')
ggsave('Umbrella/results/Figures/FS2_distance_phi_habitat.tiff', width=6, height=4)

# Figure S3 or Figure 9 ----
ord_gdf_long %>%
  group_by(PtName, HabCatfix, USFS_hab_index, mu, Habitat.group) %>%
  summarize(sRA=sum(rel.abundance)) %>%
  ggplot()+
  geom_point(aes(x=mu, y=sRA, color=HabCatfix),
             alpha=0.5)+
  facet_wrap(~Habitat.group)+
  scale_y_continuous(name='Relative Abundance',
                     labels=scales::percent)+
  scale_color_viridis_d('RCW habitat score')+
  theme_cowplot()+
  theme(legend.position='bottom')

ggsave('Umbrella/results/Figures/Fig9_relabun_mu_habgroup.tiff',
       width=6, height=6)

#check out PRAW 
ord_gdf_long %>%
  filter(species == 'PRAW',
         rel.abundance!=0) %>%
  ggplot()+
  geom_point(aes(x=USFS_hab_index, y=`Birds/sqkm`),
             size=1.5, alpha=.7)+
  scale_x_continuous('RCW habitat score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76))+
  scale_y_continuous(expression('Birds / km'^2))+
  theme_cowplot()+
  theme(axis.title=element_text(size=10),
        axis.text=element_text(size=9))

ggsave('Umbrella/results/Figures/FigS3_praw_abun.tiff',
       width=3, height=3)  
  
# look at which species driving this
pdf('Umbrella/results/fuzz_ord_spp_rel_ab.pdf')
for(j in unique(bird.assem$Habitat.group)){
  ra_ord_plot <- ord_gdf_long %>%
    filter(rel.abundance != 0,
           Habitat.group == j) %>%
    ggplot()+
    geom_point(aes(x=mu, y=rel.abundance, color=species),
               alpha=0.5)+
    facet_wrap(~species) +
    scale_y_continuous(name='Relative Abundance',
                       labels=scales::percent)+
    scale_color_viridis_d('', option='magma')+
    ggtitle(j)+
    theme_cowplot()+
    theme(legend.position='bottom',
          axis.text.x=element_text(size=9, angle=20, hjust=.9))
  print(ra_ord_plot)
}
dev.off()
