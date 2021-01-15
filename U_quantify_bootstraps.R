library(tidyverse); library(cowplot); library(grid); library(readxl)

source('Umbrella/U_load_data.r') #loads bird.assem dataframe
rm(rare_birds, point_sample_empty, plot_types, bird.codes, files, firesup, generalist, longleaf, pocosin)

#pull in the model information and abundance information from the bootstrap -----
all_mod_run<-NULL
for(u in c(list.files('C:/Users/Owner/Downloads/bootstrap_results (3)/', 
                    pattern='model', full.names = T),
           list.files('C:/Users/Owner/Downloads/bootstrap_results_scaled/', 
                      pattern='model', full.names = T))){
  ulll<-read.csv(u) %>%
    mutate(species=ifelse(substr(u, 48, 50)=='ed/', 
                      substr(u, 51,54), substr(u, 48, 51)),
           type=ifelse(grepl('scaled', u), 'scaled', 'uns'),
           batch=ifelse(type=='scaled', substr(u, 68, 69), substr(u, 63, 63)))
  if(ulll %>% count(runN) %>% pull(n) %>% max()!=1){
    ulll <- ulll %>% mutate(runN=1:nrow(.))}
  all_mod_run<-bind_rows(all_mod_run, ulll)
}
str(all_mod_run)
unique(all_mod_run$batch)
rm(u, ulll)
all_mod_run %>% mutate(sp.run.bat=paste(species, runN, batch, type)) %>%
  count(sp.run.bat) %>% arrange(desc(n)) %>% slice(1:4)

#checking scaled == not scaled
all_mod_run %>%
  filter(species %in% c('PRAW','RCWO', 'SUTA', 'GCFL')) %>% 
  ggplot()+
  #geom_density(aes(x=Hab.Est), color='purple', size=4)+
  geom_density(aes(x=Hab.Est, fill=type), alpha=.5)+
  facet_wrap(~species)

#loading the abundance data
# each row is an abundance estimate at one sample (visit to a site)
all_abun_est<-NULL
for(u in c(list.files('C:/Users/Owner/Downloads/bootstrap_results_scaled/', 
                    pattern='abun', full.names = T),
           list.files('C:/Users/Owner/Downloads/bootstrap_results (3)/', 
                      pattern='abun', full.names = T))){
  a <- read.csv(u) %>%
    mutate(species=ifelse(substr(u, 48, 50)=='ed/', 
                          substr(u, 51,54), substr(u, 48, 51)),
           type=ifelse(grepl('scaled', u), 'scaled', 'uns'),
           batch=ifelse(type=='scaled', substr(u, 63, 64), substr(u, 59, 59))) %>%
    dplyr::select(Species, Sample.Label, old.label,
                  Year, HabCat, phi, tob, USFS_hab_index,
                  run, nchat, nhat, dens, NEWnhat, raw_dens_no_phi,
                  species, batch, type)
    if(a %>% count(run) %>% pull(n) %>% max()!=296){
    a$run <- rep(1:(nrow(a)/296), each=296)}
  all_abun_est<-bind_rows(all_abun_est, a)
}
str(all_abun_est)
rm(a, u)
unique(all_abun_est$batch)

all_abun_est %>% mutate(sp.run.bat=paste(species, run, batch, type)) %>%
  count(sp.run.bat) %>% arrange(desc(n)) %>% slice(1:4)
all_abun_est %>% mutate(sp.run.bat=paste(species, run, batch, type)) %>%
  pull(sp.run.bat) %>% unique() %>% length()
nrow(all_mod_run)

head(all_abun_est)
#checking scaled == not scaled
all_abun_est %>%
  filter(species %in% c('PRAW','RCWO', 'SUTA', 'GCFL')) %>% 
  ggplot()+
  geom_density(aes(x=NEWnhat, fill=type), alpha=.3)+
  facet_wrap(~species, scales='free')

# identify the area sampled at each point, used to quantify density later?
sampled_area_km<-read.csv('Umbrella/results/models_20201016.csv') %>%
  group_by(species) %>%
  summarize(samp.area=mean(trun.dist*trun.dist*pi/1e6)) 

# summarizing the data from the bootstraps -----
#bootstraps
(nbootstraps <- all_mod_run %>% filter(phi !=10) %>% count(species) %>% arrange(n))


#models I will use
set.seed(122015) #repeatability
booted_mod<-all_mod_run %>% 
  filter(phi < 1.1) %>%
  mutate(sp.run.bat=paste(species, runN, batch, type)) %>%
  group_by(species) %>%
  slice_sample(n=500)
#check this worked
booted_mod %>% count(species) %>% arrange(n)

booted_reg<-booted_mod %>%
  dplyr::select(-X, -runN, -n_obs, -sp.run.bat, -type, -batch) %>%
  group_by(species) %>%
  summarize_all(quantile, prob=c(.5,.025,.975)) %>%
  mutate(value_type=c('median','ll','ul')) %>%
  dplyr::select(species, value_type, everything())
write.csv(booted_reg,'Umbrella/results/abundance_model_summary.csv')

booted_reg_sum<-booted_reg %>%
  dplyr::select(species, Hab.Est, Hab.pvalue, phi, value_type) %>%
  pivot_wider(names_from='value_type', 
              values_from=c('Hab.Est', 'Hab.pvalue', 'phi')) %>%
  summarize(Hab.coeff=paste0(round(Hab.Est_median,3),
                             ' (', round(Hab.Est_ll,2),
                             ', ', round(Hab.Est_ul,2), ')'),
            Hab.pval=paste0(round(Hab.pvalue_median,3),
                            ' (', round(Hab.pvalue_ll,2),
                            ', ', round(Hab.pvalue_ul,2), ')'),
            Sig=ifelse(Hab.pvalue_median < 0.051, "*",""),
            phi=paste0(round(phi_median,2),
                       ' (', round(phi_ll,1),
                       ', ', round(phi_ul,1), ')'))
booted_reg_sum

reg_sig<-booted_mod %>% 
  group_by(species) %>%
  summarize(meanHab=median(Hab.Est),
            meanHabll=quantile(Hab.Est, .025),
            meanHabul=quantile(Hab.Est, .975),
            qt0=ecdf(Hab.Est)(0)) %>%
  left_join(bird.assem, by=c('species'='Species code')) %>%
  arrange(Habitat.group) %>%
  mutate(new_sig=ifelse(qt0 < 0.026 | qt0 > 0.974, '*',''))

abundance_good<-all_abun_est %>% as_tibble() %>% 
  mutate(sp.run.bat=paste(species, run, batch, type)) %>%
  filter(sp.run.bat %in% booted_mod$sp.run.bat)
write.csv(abundance_good, 'Umbrella/results/abundance_used_quant.csv')
str(abundance_good)
abundance_good %>% count(species) %>% mutate(n/296) %>%
  arrange(`n/296`)

abundance_summary<-abundance_good %>%
  left_join(sampled_area_km) %>%
  #calculate abundance info at site level
  dplyr::mutate(nchat=tob/phi, 
                nhat=nchat/(samp.area*1178),
                dens=(1/samp.area)*nchat,
                NEWnhat=nchat/samp.area,
                raw_dens_no_phi=tob/samp.area, 
                PtName=gsub("\\..*","", Sample.Label),
                HabCatfix=case_when(substr(HabCat,4,4)==5~sub('5','6',HabCat),
                                    T~as.character(HabCat))) %>%
  dplyr::select(-Sample.Label, -old.label, -run) %>%
  replace_na(list(nchat=0, nhat=0, dens=0, NEWnhat=0))%>%
  # within each bootstrap, calculate mean abundance
  group_by(species, sp.run.bat, HabCatfix) %>%
  summarise(mNh=mean(NEWnhat)) %>%
  # calculate quantiles of mean of each habcat in each bootstrap
  group_by(species, HabCatfix) %>%
  select(-sp.run.bat) %>%
  summarize_all(quantile, prob=c(.5,.025,.975))%>%
  mutate(var_type=c('x50','x2.5', 'x97.5')) %>%
  pivot_wider(names_from = 'var_type', values_from=mNh) %>%
  mutate(`Mean abundance birds/km`=paste0(round(x50,1),
                                          ' (', round(x2.5,1),
                                          ', ', round(x97.5,1), ')')) 
abundance_summary 
write.csv(abundance_summary, 'Umbrella/results/abundance summary 95 confidence.csv')

npts_hab<- DCERP_filt %>% group_by(`Pt name`, USFS_hab_index, HabCat) %>% slice(1) %>% 
  ungroup() %>% count(HabCat) %>% filter(!is.na(HabCat))  

abundance_summary %>% 
  left_join(npts_hab, by=c('HabCatfix'='HabCat')) %>%
  group_by(species) %>%
  summarize(tmean=weighted.mean(x50, n)) %>%
  left_join(bird.assem, by=c('species'='Species code')) %>%
  filter(Habitat.group=='generalist') %>%
  arrange(desc(tmean))

# Figure 6 -----
# csv from U_occupancy_models.R
occ_gdf<-read.csv('Umbrella/OccupancyModel/occupancy_hab_coef.csv')
top_occ_mod<-read.csv('Umbrella/OccupancyModel/top_occupancy_models.csv')
occ_table<-occ_gdf %>% filter(rowname=='USFS_s') %>%
  left_join(top_occ_mod %>% select(spp, rowname, delta), 
            by=c('species'='spp', 'model_n'='rowname')) %>%
  mutate(topM=ifelse(delta==0, '*',''))

oh<-booted_mod %>%
  select(-X, -runN, -n_obs, -sp.run.bat) %>%
  group_by(species) %>% 
  summarize(Hab.Est_mean=median(Hab.Est),
            Hab.Est_SE=median(Hab.SE),
            Hab.Est_p=median(Hab.pvalue)) %>%
  left_join(occ_table %>% select(species, Estimate, SE, topM), 
            by='species') %>%
  pivot_longer(cols=c(-species, -Hab.Est_p, -topM))%>% 
  mutate(type= case_when(grepl('Hab',name)~'Density',
                         T~'Occupancy'),
         shit=case_when(grepl('SE', name)~'se',
                        grepl('se', name)~'se',
                        T~'mean'))%>%
  pivot_wider(names_from = shit,
              values_from = value) %>%
  group_by(species, type, topM) %>%
  summarize(mean=mean(mean, na.rm=T),
            se=mean(se, na.rm=T),
            meanp=median(Hab.Est_p)) %>%
  left_join(reg_sig %>% select(qt0, species)) %>%
  mutate(sig=ifelse(qt0>= 0.975 | qt0<=0.025, '@', ''),
         allsig=case_when(topM=='*' & sig=='@'~'B',
                          topM=='*' & sig!='@'~'o',
                          topM!='*' & sig=='@'~'D',
                          T~'')) %>%
  left_join(bird.assem, by=c('species'='Species code')) %>%
  group_by(Habitat.group) %>%
  mutate(sigx=min(mean-se),
         HabF=factor(Habitat.group, levels=c('longleaf', 'generalist','hardwood', 'shrub'))) %>%
  filter(species != 'RCWO') %>%
  ggplot()+
  geom_rect(xmin=0,xmax=5.7, ymin=-Inf, ymax=Inf, 
            fill='grey85')+
  geom_vline(xintercept=0, linetype='dashed', color='black')+
  geom_linerange(aes(y=`Common name`, color=type,
                     xmin=mean-se, xmax=mean+se),
                 position=position_dodge(width=0.4),
                 alpha=0.5)+
  geom_point(aes(y=`Common name`, x=mean, color=type), 
             size=2, position=position_dodge(width=.4),
             alpha=0.5)+
  geom_text(aes(y=`Common name`, x=min(sigx), label=allsig))+
  facet_wrap(~HabF, scales='free_y', ncol=1, 
             strip.position = 'left')+
  scale_x_continuous('Coefficient', trans=scales::pseudo_log_trans(sigma=.2),
                     expand=c(0.04,0),
                     breaks = c(-1, 0,1,3,6))+
  scale_y_discrete('')+
  scale_color_viridis_d('',end=.6)+
  theme_cowplot() +
  theme(strip.placement='outside',
        strip.background = element_rect(fill=NA),
        strip.text=element_text(size=10),
        #panel.grid.major.y = element_line(color="black"),
        legend.position = 'top',
        #legend.justification = 'center',
        axis.text.y=element_text(size=8.5))
#tutorial https://stackoverflow.com/questions/52341385/how-to-automatically-adjust-the-width-of-each-facet-for-facet-wrap
ohp = ggplot_gtable(ggplot_build(oh))
#gtable::gtable_show_layout(ohp)
# get the number of unique x-axis values per facet (1 & 3, in this case)
y.var <- sapply(ggplot_build(oh)$layout$panel_scales_y,
                function(l) length(l$range$range))
# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
ohp$heights[ohp$layout$t[grepl("panel", ohp$layout$name)]] <- ohp$heights[ohp$layout$t[grepl("panel", ohp$layout$name)]] * y.var
grid.draw(ohp)
ggsave('Umbrella/results/Figures/Fig6_coeffs.jpg', grid.draw(ohp),
       width=4.5, height=7)

# Table 3 ----
#number of points
DCERP_filt %>%
  group_by(`Pt name`) %>%
  slice(1) %>%
  group_by(HabCat_nf) %>% tally() %>%
  left_join(hab_bins, by=c('HabCat_nf'='mid')) %>%
  filter(!is.na(HabCat_nf)) 
#number of unique samples
DCERP_filt %>%
  group_by(`Pt name`, USFS_hab_index) %>%
  slice(1) %>%
  group_by(HabCat_nf) %>% tally() %>%
  filter(!is.na(HabCat_nf)) 


# Table 4 ----
# t4_occ created in U_occupancy_models.R
abundance_summary %>%
  group_by(species) %>%
  select(species, HabCatfix, x50) %>%
  left_join(npts_hab %>% mutate(HabCatfix=as.numeric(paste(HabCat))) %>%
            select(-HabCat))%>%
  mutate(w.mean=weighted.mean(x50, n)) %>%
  filter(HabCatfix %in% c(2.5, 4.76)) %>%
  select(-n) %>%
  pivot_wider(names_from=HabCatfix, values_from=x50) %>%
  mutate(l_to_h=paste0(round(`2.5`,1), ' to ', round(`4.76`,1))) %>%
  left_join(booted_reg_sum) %>%
  left_join(reg_sig) %>%
  left_join(t4_occ, by=c('species'='spp')) %>%
  select(Habitat.group, `Common name`, `Occupancy Estimate`, phi, Hab.coeff, qt0, w.mean, l_to_h) %>%
  arrange(Habitat.group) %>% View()
  write.csv('Umbrella/results/Table4_20201116.csv')

# Old plots not used anymore ------
#mean abundance across all habitat type
g<-abundance_summary %>%
  group_by(species) %>% 
  summarize_if(is.numeric, ~weighted.mean(x=., w=npts_hab$n)) %>% #weight this based on sample size
  left_join(bird.assem, by=c('species'='Species code')) %>%
  ggplot()+
  geom_linerange(aes(x=species, ymin=x2.5, ymax=x97.5)) +
  geom_point(aes(x=species, y=x50), size=2.5) +
  scale_y_continuous(name=expression('Birds km '^-2))+
  facet_wrap(~Habitat.group, nrow=1, scales='free_x',
             strip.position = 'bottom') +
  theme_cowplot()+
  theme(axis.text.x = element_text(size=9, angle=30, hjust=0.9),
        strip.placement = 'outside',
        strip.text.x = element_text(size=10))
#tutorial https://stackoverflow.com/questions/52341385/how-to-automatically-adjust-the-width-of-each-facet-for-facet-wrap
gt = ggplot_gtable(ggplot_build(g))
#gtable::gtable_show_layout(gt)
# get the number of unique x-axis values per facet (1 & 3, in this case)
x.var <- sapply(ggplot_build(g)$layout$panel_scales_x,
                function(l) length(l$range$range))
# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
gt$widths[gt$layout$l[grepl("panel", gt$layout$name)]] <- gt$widths[gt$layout$l[grepl("panel", gt$layout$name)]] * x.var
grid.draw(gt)

#abundance change between high & low categories
abund.plot<-abundance_summary %>%
  filter(HabCatfix %in% c(2.5, 4.76)) %>%
  ungroup() %>%
  mutate(SpFac=as.factor(species),
         fact.level=case_when(HabCatfix==2.5~as.numeric(SpFac)-.075,
                              HabCatfix==4.76~as.numeric(SpFac)+.075))%>%
  ggplot()+
  geom_linerange(aes(x=species, ymin=x2.5, ymax=x97.5, 
                     color=as_factor(HabCatfix)),
                 position = position_dodge(width=.5),
                 size=1.25)+
  geom_line(aes(x=fact.level, y=x50, group=species),
            color='grey')+
  geom_point(aes(x=species, y=x50, fill=as_factor(HabCatfix), 
                 shape=as_factor(HabCatfix)),
             position = position_dodge(width=.5),
             size=2.5)+
  scale_color_viridis_d(name='RCW Habitat Category', 
                        end=.8)+
  scale_fill_manual(name='RCW Habitat Category', values=c('white','black'))+
  scale_shape_manual(name='RCW Habitat Category',
                     values=c(22,16))+
  scale_y_continuous(trans='log1p',
                     expression('Birds per km'^2),
                     breaks=c(.5,1,5,10,25,50,100,200),
                     labels=c(0.5, 1,5,10,25,50,100,200),
                     expand=c(0,0))+
  scale_x_discrete("", labels=rep('', 28))+
  theme_cowplot()+
  theme(#panel.grid.major.y = element_line(color='lightgrey'),
    legend.position='top',
    legend.justification = 'center',
    axis.text.x = element_text(size=-5),
    axis.title.x = element_text(size=-3))
abund.plot

# abundance regression plot
reg_plot<-booted_reg %>%
  select(species, value_type, everything()) %>%
  dplyr::select(species, value_type,Hab.Est, Hab.pvalue) %>%
  pivot_wider(names_from = 'value_type', 
              values_from=c('Hab.Est','Hab.pvalue')) %>%
  mutate(pas=ifelse(Hab.pvalue_ul <=0.5, '*', ""),
         direction=case_when(Hab.Est_ll > 0 & Hab.Est_ul > 0 ~ '+',
                             Hab.Est_ll <0 & Hab.Est_ul > 0 ~ '0',
                             Hab.Est_ll < 0 & Hab.Est_ul < 0 ~ '-')) %>%
  ggplot()+
  geom_linerange(aes(x=species, ymin=Hab.Est_ll, ymax=Hab.Est_ul))+
  geom_point(aes(x=species, y=Hab.Est_median),
             size=2)+
  geom_hline(yintercept=0, linetype='dashed', color='grey')+
  geom_text(aes(x=species, y= -.85, label=direction, color=pas), 
            size=4)+
  #geom_text(data=nbootstraps,
  #          aes(x=species, y= 1.9, label=n), size=3)+
  scale_y_continuous('Habitat Category\nCoefficient')+
  scale_x_discrete('')+
  scale_color_manual(name="Statistical Sig.",
                     labels=c('ns','p < 0.5'), values=c('black','blue'))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=30, hjust=.9, size=9),
        axis.title.x = element_text(size=-2))
plot_grid(abund.plot, reg_plot+theme(legend.position = 'none'), 
          rel_heights=c(.8, 1), ncol=1, labels='AUTO')
#density & abundance beta coeff together
ggsave('Umbrella/results/Figures/abundance_estimates_reg.jpg',
       width=6, height=6)
reg.plot<-reg_plot +coord_flip()+
  scale_y_continuous('Habitat Category Coefficient',
                     expand=c(0.02,0.02))+
  theme(axis.text.y=element_text(size=10),
        axis.text.x = element_text(size=10, angle=0),
        axis.title.x = element_text(size=12),
        legend.position = 'top',
        legend.box = 'vertical')
reg.plot #abundance beta coeff together
ggsave('Umbrella/results/Figures/hab_coeff_reg.jpg',
       width=4, height=4.25)

