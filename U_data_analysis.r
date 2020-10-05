# General Data Analysis #####
library(tidyverse); library(readxl); library(Distance); library(lubridate)
library(cowplot)

#load in the raw data
source('Umbrella/U_load_data.r')
species<-c('NOCA','NOBO','EATO','PIWA','EAWP','BHCO','MODO',
           'RHWO','PRAW','WEVI','NOFL','YTWA','RBWO','INBU',
           'SUTA','CACH','COYE','BLJA','RCWO','AMCR','BACS',
           'GCFL','BGGN','TUTI','CARW','BHNU','CHSP','EABL')
for(j in species){
  assign(paste0(j, '_dat'),
         spfilter(j))
}
model_parameters<-read.csv('Umbrella/all_possible_models (1).csv') %>%
  dplyr::select(-X) %>%
  bind_rows(data.frame(formu))
old_model_param<-read.csv('Umbrella/all_possible_models_0929.csv') %>%
  mutate(model.n=paste0('mod',row_number()))

# Habitat Associate groups from Allen et al. 2006 (the Auk) -----
#bird codes from: https://www.birdpop.org/pages/birdSpeciesCodes.php 
bird.codes<-read.csv('Umbrella/results/Cornell_list19p.csv')
longleaf<-c('Pine Warbler', 'Brown-headed Nuthatch','Red-cockaded Woodpecker',
            'Prairie Warbler',"Bachman's Sparrow", 'Chipping Sparrow',
            'Eastern Wood-Pewee', 'Red-headed Woodpecker', 'Eastern Bluebird')
firesup<-c('Red-eyed Vireo', 'Acadian Flycatcher', 'Ovenbird', 
           'Black-and-white Warbler', 'Tufted Titmouse', 'Wood Thrush', 'Yellow-throated Vireo',
           'Yellow-throated Warbler', 'Blue-gray Gnatcatcher')
pocosin<-c('Eastern Towhee', 'Common Yellowthroat', 'Carolina Wren', 'Northern Cardinal', 
           'White-eyed Vireo','Hooded Warbler')
generalist<-c('Carolina Chickadee', 'Summer Tanager', 'Great Crested Flycatcher', 'Blue Jay',
              'Red-bellied Woodpecker', 'Northern Flicker', 'Brown-headed Cowbird', 'Indigo Bunting',
              'Mourning Dove', 'Northern Bobwhite', 'American Crow')
bird.assem<-data.frame(COMMONNAME=c(longleaf, firesup, pocosin, generalist),
                       `Habitat group`=c(rep('longleaf', length(longleaf)),
                                         rep('hardwood', length(firesup)),
                                         rep('shrub', length(pocosin)),
                                         rep('generalist', length(generalist)))) %>%  
  left_join(bird.codes) %>%
  filter(SPEC %in% species) %>%
  dplyr::select(-SP, -CONF, -SPEC6, -CONF6) %>%
  rename(`Common name`='COMMONNAME', `Scientific name`='SCINAME',
         `Species code`='SPEC') %>%
  dplyr::select(`Scientific name`, `Common name`, `Species code`, Habitat.group)

#birds that are not within the habitat-association table
bird.codes %>%
  filter(SPEC %in% species[!(species %in% bird.assem$`Species code`)])

#Once all the species distance sampling models are built
#pull back in the model information after it has been run -----
model.files<-list.files(path='C:/Users/Owner/Downloads/DistanceModels (1)/', #path='Umbrella/results/DistanceModels/',
                        pattern='.csv',
                        full.names=T)
all_models<-NULL
for(b in model.files){
  if(substr(b, 45,48)!='hab_'){  
    ms<-read.csv(b) %>% left_join(old_model_param, by='model.n')} else {
    ms<-read.csv(b) %>% left_join(model_parameters, by='model.n')  
    }
  all_models<-bind_rows(all_models, ms)
}

all_models<-all_models %>% dplyr::select(-model.n, -model.name) %>%
  left_join(model_parameters) %>% 
  mutate(model.name=paste(species, model.n, sep='.'))

#Check all needed models were run for all species
models_completed<- all_models %>% 
  filter(species=='NOCA') %>% pull(model.n)
ideal_model_list<-data.frame(species=rep(unique(all_models$species), 
                                         each=length(models_completed)),
                             model.n=rep(models_completed, 28)) %>%
  mutate(model.name=paste(species, model.n, sep='.'))

ideal_model_list %>% group_by(species) %>% tally() %>% arrange(n)
ideal_model_list %>% group_by(model.n) %>% tally() %>% arrange(n)

setdiff(ideal_model_list, 
        all_models %>% dplyr::select(species, model.n, model.name)) %>%
  count(species)

# identify which models explained detection/abundance well
model.appendix<-all_models %>% as_tibble() %>%
  group_by(species) %>%
  filter(gof_cvm_p >= 0.050) %>%
  mutate(deltaAICc=AICc - min(AICc),
         dAIC=AIC-min(AIC)) %>%
  filter(deltaAICc < 2) %>%
  arrange(AICc) %>%
  dplyr::select(-model.state, -dAIC)
model.appendix 
length(unique(model.appendix$species))

#species with models that do not fit reality
species[!(species %in% model.appendix$species)]
best_bad_models<-all_models %>% 
  filter(species %in% species[!(species %in% model.appendix$species)]) %>%
  group_by(species) %>%
  arrange(desc(gof_cvm_p), AICc) %>%
  slice(1:5) %>%
  filter(AICc == min(AICc)) %>% 
  slice(1) %>%
  #left_join(model_parameters, by='model.n') %>%
  dplyr::select(-model.state)

model.appendix.all<-bind_rows(model.appendix, best_bad_models)

write.csv(model.appendix.all, 'Umbrella/model.appendix.20201005.csv')

top.models<-model.appendix.all %>%
  group_by(species) %>%
  arrange(AICc) %>%
  slice(1:2)
#top.models<-read.csv('Umbrella/results/model_parameter_summary20200921.csv') %>% 
#  dplyr::select(species, model.name,nobs, formula)
top.models %>% tally() %>% arrange(n)

top.models<-top.models %>% 
  filter(species %in% c("RCWO", "BACS","NOCA"))


# create the models that have not been created yet
conversion.factor <- convert_units(distance_units = "Metre", 
                                   effort_units=NULL, 
                                   area_units = "square kilometre")
for(h in 1:nrow(top.models)){
  if(exists(top.models$model.name[h])){
    print(paste(top.models$model.name[h], 'exists'))
  }else{
    assign(top.models$model.name[h], 
           ds(get(paste0(top.models$species[h], '_dat')),
              transect="point", key="hr",
              formula=as.formula(top.models$formula[h]), 
              adjustment = NULL, order = 0,
              truncation = '5%', 
              convert.units = conversion.factor))
    print(paste(top.models$model.name[h], 'created'))
  }
}

#summary output for ds output
#mod.res$ds == distance function
#mod.res$dht == abundance/density info
#mod.res$ddf == detection function

#pull out the top models for each species, isolate strata phi and density ----
#rows = species x variable, columns = habitat strata
label.row.names<-c(`2.5`='X1',`3.26`='X2',`3.76`='X3',
                   `4.26`='X4',`4.76`='X5', Total='X6',
                   `2.5SD`='X1.1',`3.26SD`='X2.1',`3.76SD`='X3.1',
                   `4.26SD`='X4.1',`4.76SD`='X5.1', TotalSD='X6.1')
names(label.row.names)<-paste('Nc', names(label.row.names), sep='_')
#for loop to pull out summary data (phi and density) from each top model
model_sums<-NULL
all_site_comp_tm<-NULL
for(p in unique(top.models$model.name)){
  m<-get(p)
 
  if(is.null(m$dht)){
    tm_sums<-data.frame(species=as.character(top.models[top.models$model.name==p,"species"]),
                        model.name=p,
                        model.state='broken')
  }else{
    mod.sum<-summary(m)
    #mod.gof<-gof_ds(m)
    #building species summary information
    Nc.estimates<-mod.sum$dht$individuals$bysample %>% 
      group_by(Region) %>% summarize(mean(Nchat),
                                     sd(Nchat),
                                     .groups='keep') %>%
      bind_rows(mod.sum$dht$individuals$bysample %>% 
                  summarize(Region='Total', mean(Nchat),sd(Nchat), .groups='keep'))
    
    tm_sums<-data.frame(species=as.character(top.models[top.models$model.name==p,"species"]),
                        model.name=p,
                        model.state='fine',
                        AIC=AIC(m),
                        phiMEAN=mod.sum$ds$average.p,
                        phiSE=mod.sum$ds$average.p.se,
                        trun.dist=mod.sum$ds$width,
                        nobs=nrow(mod.sum$ddf$data),
                        t(Nc.estimates[,2]),
                        t(Nc.estimates[,3])) %>% 
      rename(all_of(label.row.names))
    
    #building site summary information
    site_comp<-mod.sum$dht$individuals$bysample %>% 
      mutate(Site=gsub("\\..*","",Sample),
             model.name=p) %>% 
      group_by(Site, model.name) %>%
      summarize('Mean.Est.N'=mean(Nchat), .groups='keep',#averaging the replicates here!
                species=as.character(top.models[top.models$model.name==p,"species"]))
    all_site_comp_tm<-bind_rows(all_site_comp_tm, site_comp)
    rm(mod.sum)
  }
  model_sums<-bind_rows(model_sums, tm_sums)
}
View(model_sums)

# identify the top models for each species
topmod_sums <- model_sums %>% as_tibble() %>%
  dplyr::select(-model.state) %>%
  left_join(top.models, by = c("species", "model.name", 'nobs')) 
View(topmod_sums)
#saving topmod_sums as RDS so I don't have to run them each time
for(u in unique(topmod_sums$model.name)){
  saveRDS(get(u), 
          file=paste0('Umbrella/results/DistanceModels/rds_models/',
                      u, '_', format(Sys.Date(), '%Y%m%d'),'.rds'))
}

write.csv(topmod_sums, paste0('Umbrella/results/model_parameter_summary',
                              format(Sys.Date(), "%Y%m%d"), '.csv'))
#topmod_sums <- read.csv('Umbrella/results/model_parameter_summary20200921.csv') %>%
#  dplyr::select(-X)
write.csv(topmod_sums %>% 
            dplyr::select(-ends_with('SD')) %>%
            dplyr::select(species, formula, trun.dist, nobs, AIC.df, AIC.AIC, everything()),
          'Umbrella/model_appendix_table.csv', row.names=F)


# build trend tables for species within this list -----
pop.trend<-topmod_sums %>%
  dplyr::select(-ends_with('SD')) %>%
  dplyr::select(species, starts_with('Nc')) %>%
  dplyr::select(-Nc_Total) %>%
  pivot_longer(-species) %>%
  mutate(hab_cat=as.numeric(gsub('Nc_', '', name))) %>%
  group_by(species) %>%
  summarize(slope=summary(lm(value~hab_cat))$coefficients[2,1],
            slopeSE=summary(lm(value~hab_cat))$coefficients[2,2],
            slope.ul=slope+1.96*slopeSE,
            slope.ll=slope-1.96*slopeSE,
            pval=summary(lm(value~hab_cat))$coefficients[2,4],
            .groups='keep') %>%
  mutate(trend=case_when(slope.ul > 0 & slope.ll > 0 ~ "increases with restoration",
                         slope.ul > 0 & slope.ll < 0 ~ paste('slope overlap is',as.character(round(abs(slope.ll)-0,2))),
                         slope.ul < 0 & slope.ll < 0 ~ 'decreases with restoration')) %>%
  arrange(trend) %>%
  left_join(bird.assem[,2:4], by=c('species'='Species code')) %>% 
  arrange(Habitat.group)
write.csv(pop.trend, 
          'Umbrella/results/density_trends_byspp_20200925.csv')

#also want to make a graphing table - density at each site ----

#Density
all_site_comp<-all_site_comp_tm %>% 
  filter(model.name %in% topmod_sums$model.name) %>%
  ungroup() %>% dplyr::select(-model.name)
  
dens_gdf<-all_site_comp %>%
  left_join(DCERP %>%
              group_by(`Pt name`) %>%
              dplyr::select(`Pt name`, USFS_hab_index, Base_Area) %>%
              slice(1),
            by=c('Site'='Pt name')) %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)]) %>%
  left_join(bird.assem, by=c("species"='Species code'))

pdf('Umbrella/results/spp_abundance_20200924_hab.pdf')
com_plot<-ggplot()+scale_color_viridis_d('RCWO Habitat Matrix Score',
                                option='plasma',
                                end=0.8,
                                alpha=0.4)+
  scale_y_continuous(name='Estimate Birds * sq. km')+
  theme_bw()+theme(legend.position = 'bottom')+
  facet_wrap(~species, scales='free_y')
longleaf_dens<-com_plot+
  geom_point(data=dens_gdf[dens_gdf$Habitat.group=='longleaf',], 
             aes(x=USFS_hab_index, y=Mean.Est.N,
                 color=as.factor(Region.Label)))+
  stat_summary(data=dens_gdf[dens_gdf$Habitat.group=='longleaf',], 
               aes(x=Region.Label, y=Mean.Est.N))+
  ggtitle('Longleaf Pine-associated Species')
print(longleaf_dens)
hard_dens<-com_plot+
  geom_point(data=dens_gdf[dens_gdf$Habitat.group=='hardwood',], 
             aes(x=USFS_hab_index, y=Mean.Est.N,
                 color=as.factor(Region.Label)))+
  stat_summary(data=dens_gdf[dens_gdf$Habitat.group=='hardwood',],
               aes(x=Region.Label, y=Mean.Est.N))+
  ggtitle('Hardwood-associated Species')
print(hard_dens)
shrub_dens<-com_plot+
  geom_point(data=dens_gdf[dens_gdf$Habitat.group=='shrub',], 
             aes(x=USFS_hab_index, y=Mean.Est.N,
                 color=as.factor(Region.Label)))+
  stat_summary(data=dens_gdf[dens_gdf$Habitat.group=='shrub',],
               aes(x=Region.Label, y=Mean.Est.N))+
  ggtitle('Shrub-associated Species')
print(shrub_dens)
gen_dens<-com_plot+
  geom_point(data=dens_gdf[dens_gdf$Habitat.group=='generalist',], 
             aes(x=USFS_hab_index, y=Mean.Est.N,
                 color=as.factor(Region.Label)))+
  stat_summary(data=dens_gdf[dens_gdf$Habitat.group=='generalist',],
               aes(x=Region.Label, y=Mean.Est.N))+
  ggtitle('Habitat Generalists')
print(gen_dens)
dev.off()

abund_label_df<- dens_gdf %>%
  filter(Region.Label == 4.76) %>%
  group_by(species, Region.Label, Habitat.group) %>%
  summarize(meanNc=mean(Mean.Est.N))

ggplot()+
  scale_y_continuous(name=expression('Birds km'^2))+
  scale_x_continuous('USFWS Habitat index')+
  stat_summary(data=dens_gdf, 
               aes(x=Region.Label, y=Mean.Est.N, group=species),
               geom='line', color='lightgrey')+
  stat_summary(data=dens_gdf, 
               aes(x=Region.Label, y=Mean.Est.N),
               geom='line')+
  stat_summary(data=dens_gdf, 
               aes(x=Region.Label, y=Mean.Est.N))+
  #geom_text_repel(data=abund_label_df,
  #                aes(x=Region.Label, y=meanNc, label=species),
  #                size=2.5, color='grey')+
  facet_wrap(~Habitat.group)+
  theme_cowplot()+theme(legend.position='none')
ggsave('Umbrella/results/Figures/abundance_trends.tiff',
       width=6, height=6)

#summary table for manuscript ----
bird.codes %>%
  filter(SPEC %in% species) %>%
  rename(species="SPEC") %>%
  left_join(topmod_sums) %>%
  dplyr::select(-ends_with('SD'))

# fuzzy ordination to assess alignment with USFS habitat index #####
library(fso); library(vegan)

#build a site x community matrix based on bird densities
ordination_df<-all_site_comp %>% 
  pivot_wider(names_from = 'species', values_from='Mean.Est.N') %>%
  left_join(DCERP %>%
              group_by(`Pt name`) %>%
              dplyr::select(`Pt name`, USFS_hab_index, Base_Area) %>%
              slice(1),
            by=c('Site'='Pt name')) %>% 
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)]) %>%
  dplyr::select(Site, USFS_hab_index, Base_Area, Region.Label, 
                RCWO, 
                BACS, EAWP, MODO, #want to exclude these because models are poor fit for reality
                everything()) #reordering

site.com<-ordination_df[,-c(2:5)] %>% #excluding RCW as we know they should be correlated with them 
  bind_cols(., total.dens=rowSums(.[,-1])) %>%
  filter(total.dens != 0) %>%
  dplyr::select(-total.dens) %>%
  ungroup()
nrow(site.com); head(site.com)
dis<-vegdist(site.com[,-1], method="bray") #create a distance matrix for the bird communities

site.env<-ordination_df[,c(1:5)] %>%
  filter(Site %in% site.com$Site)

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
  dplyr::select(Site, mu, everything())
ggplot()+
  geom_point(data=ord_gdf, aes(x=USFS_hab_index, y=mu), alpha=0.6)+
  geom_text(aes(x=2.8, y=0.535, label="r = 0.562\np < 0.001"))+
  #ggtitle('Fit of USFS Habitat Index ordination')+
  scale_y_continuous('Predicted habitat index')+
  scale_x_continuous('USFWS Habitat index')+
  theme_cowplot()
ggsave('Umbrella/results/fuzzy_ord_fit.tiff',
       width=3, height=3)

mean(ord_gdf$USFS_hab_index); sd(ord_gdf$USFS_hab_index)
hab_range<-max(ord_gdf$USFS_hab_index)-min(ord_gdf$USFS_hab_index)
hist((ord_gdf$USFS_hab_index-min(ord_gdf$USFS_hab_index))/hab_range)
hist(ord_gdf$mu)
hist((ord_gdf$mu*5)+2.112)
hist((ord_gdf$USFS_hab_index-2.112)/5)

site.ordered<-ord_gdf$Site[order(ord_gdf$mu)]
ord_gdf_long<-ord_gdf %>%
  arrange(mu)%>%
  rowid_to_column('mu.order') %>%
  arrange(USFS_hab_index) %>%
  rowid_to_column('hab.order') %>%
  pivot_longer(cols=c(BACS:YTWA), 
               names_to="species", values_to="Birds/sqkm")%>%
  mutate(Site.ord=factor(Site, levels=site.ordered)) %>%
  left_join(bird.assem, by=c('species'='Species code')) %>%
  mutate(HabitatF=factor(Habitat.group, 
                         levels=c('generalist', 'shrub','hardwood','longleaf'))) %>%
  group_by(Site) %>%
  mutate(rel.abundance=(`Birds/sqkm`/sum(`Birds/sqkm`)))

#ordinated communities
mu_breaks<-ord_gdf_long %>% arrange(mu) %>%
  pull(mu) %>% unique() %>% round(3)
ggplot(data=ord_gdf_long)+
  geom_col(aes(x=mu.order, y=rel.abundance, 
               fill=HabitatF))+
  scale_fill_viridis_d(name="Habitat\nassociation")+
  scale_x_continuous(name='mu\nhigher mu corr. with higher habitat index',
                     breaks=c(seq(1,148,25),148),
                     labels=mu_breaks[c(seq(1,148,25),148)],
                     expand = c(0,0))+
  scale_y_continuous(name='Relative Abundance',
                     expand=c(0,0),
                     labels=scales::percent)+
  theme(axis.text.x = element_text(angle=90, hjust=0.9),
        axis.text.y = element_text(angle=90, hjust=0.5))+
  coord_flip()+theme_cowplot()
ggsave('Umbrella/results/mu_ordinated_sitecom_ra.tiff',
       width=5, height=9)
ggplot(data=ord_gdf_long)+
  geom_col(aes(x=mu.order, y=`Birds/sqkm`, 
               fill=HabitatF))+
  scale_fill_viridis_d(name="Habitat\nassociation")+
  scale_x_continuous(name='mu\nhigher mu corr. with higher habitat index',
                     breaks=c(seq(1,148,25),148),
                     labels=mu_breaks[c(seq(1,148,25),148)],
                     expand = c(0,0))+
  scale_y_continuous(name=expression("total birds km"^-2),
                     expand=c(0,0))+
  theme(axis.text.x = element_text(angle=90, hjust=0.9),
        axis.text.y = element_text(angle=90, hjust=0.5))+
  coord_flip()+theme_cowplot()
ggsave('Umbrella/results/mu_ordinated_sitecom_total.tiff',
       width=5, height=9)

#raw data plots
hab_breaks<-ord_gdf %>% arrange(USFS_hab_index) %>%
  pull(USFS_hab_index) %>% round(2)
ggplot(data=ord_gdf_long)+
  geom_col(aes(x=hab.order, y=rel.abundance, 
               fill=HabitatF))+
  scale_fill_viridis_d(name="Habitat\nassociation")+
  scale_x_continuous(name='USFWS habitat index',
                     breaks=c(seq(1,148,25),148),
                     labels=hab_breaks[c(seq(1,148,25),148)],
                     expand = c(0,0))+
  scale_y_continuous(name='Relative Abundance',
                     expand=c(0,0),
                     labels=scales::percent)+
  theme(axis.text.x = element_text(angle=90, hjust=0.9),
        axis.text.y = element_text(angle=90, hjust=0.5))+
  coord_flip()+theme_cowplot()
ggsave('Umbrella/results/usfws_ordinated_sitecom_ra.tiff',
       width=5, height=9)

ggplot(data=ord_gdf_long)+
  geom_col(aes(x=hab.order, y=`Birds/sqkm`, 
               fill=HabitatF))+
  scale_fill_viridis_d(name="Habitat\nassociation")+
  scale_x_continuous(name='USFWS habitat index',
                     breaks=c(seq(1,148,25),148),
                     labels=hab_breaks[c(seq(1,148,25),148)],
                     expand = c(0,0))+
  scale_y_continuous(name=expression("birds km"^-2),
                     expand=c(0,0))+
  theme(axis.text.x = element_text(angle=90, hjust=0.9),
        axis.text.y = element_text(angle=90, hjust=0.5))+
  coord_flip()+theme_cowplot()
ggsave('Umbrella/results/usfws_ordinated_sitecom_dens.tiff',
       width=5, height=9)

abun_hab_df<-ord_gdf_long %>%
  group_by(Site, USFS_hab_index, Region.Label, HabitatF) %>%
  summarize(Total.Birds.sqkm=sum(`Birds/sqkm`), .groups='keep')
abun_all_df<-ord_gdf_long %>%
  group_by(Site, USFS_hab_index, Region.Label) %>%
  summarize(Total.Birds.sqkm=sum(`Birds/sqkm`), .groups='keep')
hab_abund_plot<-ggplot()+
  #geom_jitter(data=abun_hab_df,
  #           aes(x=Region.Label, y=Total.Birds.sqkm, group=Region.Label),
  #           alpha=0.3)+
  geom_boxplot(data = abun_hab_df, 
               aes(x=as.character(Region.Label), y=Total.Birds.sqkm, group=Region.Label),
               fill='lightgrey')+
  scale_x_discrete(name='RCW habitat score')+
  theme_cowplot()+facet_wrap(~HabitatF, scales='free_y')+
  theme(axis.title.y=element_text(size=0),
        axis.text.x = element_text(size=9, angle=30, hjust=.9))
abund_plot<-ggplot()+
  geom_jitter(data=abun_all_df,
              aes(x=as.character(Region.Label), 
                  y=Total.Birds.sqkm, group=Region.Label),
              alpha=0.2)+
  geom_boxplot(data = abun_all_df, 
               aes(x=as.character(Region.Label), 
                   y=Total.Birds.sqkm, group=Region.Label),
               outlier.alpha = 0, alpha=0.8, fill='lightgrey')+
  scale_x_discrete(name='RCW habitat score')+
  scale_y_continuous(expression("Total Birds km"^-2))+
  theme_cowplot()+
  theme(axis.text.x = element_text(size=11, angle=30, hjust=.9))
plot_grid(abund_plot, hab_abund_plot, labels=c("       All Birds", ""),
          rel_widths = c(0.43,0.57))
ggsave('Umbrella/results/bird_abundance.tiff', width=6, height=4)

write.csv(ord_gdf_long, 'Umbrella/plotting_df.csv')

### Table 3 ----
DCERP_filt %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)],
         Sample.Label=paste(`Pt name`, Replicate, sep="."), #this identifies unique samples within strata
         size=`Cluster size`,
         Effort=1,
         Study.Area="MCBCL") %>%
  group_by(`Pt name`) %>%
  slice(1) %>%
  group_by(Region.Label) %>% tally() %>%
  left_join(hab_bins, by=c('Region.Label'='mid')) %>%
  bind_rows(data.frame(Region.Label=0,
                       n=sum(.$n)-1))

# Figure 3 RCW density ----
ord_gdf %>%
  dplyr::select(Site, mu, USFS_hab_index, Region.Label, RCWO)

ggplot(data=ord_gdf)+
#  geom_point(aes(x=USFS_hab_index, y=RCWO, color=as.character(Region.Label)), 
#             alpha=0.3)+
  stat_summary(aes(x=Region.Label, y=RCWO))+
  scale_y_continuous(name=expression("RCW per km"^-2))+
  scale_color_viridis_d(guide=F, option='plasma')+
  theme_cowplot()
ggsave('Umbrella/results/fig3b_rcwo_abun_dist.tiff',
       width=3, height=3)

# detection at each strata ====
model_phi_long<-NULL
for(p in unique(top.models$model.name)){
  m<-get(p)
  model_phi<-bind_cols(m$ddf$data, 
                       predict.phi=predict(m,esw=FALSE)$fitted) %>%
    as_tibble() %>%
    dplyr::select(Region.Label, predict.phi) %>%
    group_by(Region.Label) %>%
    summarize(Region=as.character(Region.Label)[1],
              mean.phi=mean(predict.phi),
              .groups='drop') %>%
    bind_rows(data.frame(Region='Total',
                         mean.phi=mean(predict(m,esw=FALSE)$fitted))) %>%
    dplyr::select(-Region.Label) %>%
    mutate(species = top.models[top.models$model.name==p,]$species)
  model_phi_long<-bind_rows(model_phi_long, model_phi)
}
phi_table<-model_phi_long  %>% 
  pivot_wider(names_from = Region, values_from=mean.phi)
write.csv(phi_table,'Umbrella/results/phi_summary_table.csv')

big_delt_sp<-model_phi_long %>%
  filter(Region != 'Total') %>%
  group_by(species) %>%
  mutate(prevReg=lag(mean.phi),
         difBTWNreg=mean.phi-prevReg) %>%
  filter(difBTWNreg !=0,
    difBTWNreg < -0.05 | difBTWNreg > 0.05) %>%
  pull(species) %>% unique()

library(ggrepel)
ggplot()+
  geom_line(data=model_phi_long[model_phi_long$species %in% big_delt_sp,],
            aes(x=Region, y=mean.phi, group=species,
                color=species), size=2, alpha=0.6)+
  geom_text_repel(data=phi_table[phi_table$species %in% big_delt_sp,],
                  aes(x='Total', y=Total, label=species))+
  scale_y_continuous('Detection Probability')+
  scale_x_discrete(expand = expansion(mult=c(0.02, .17)))+
  scale_color_viridis_d(option = 'inferno')+
  theme_cowplot()+
  theme(legend.position = 'none')
ggsave('Umbrella/largedetectchanges.tiff', width=6, height=4)

# parameter barchart ----
mod_parm<-topmod_sums %>%
  dplyr::select(species, formula) %>%
  separate(formula, into=paste0('X', 1:7)) %>%
  dplyr::select(-X1) %>%
  pivot_longer(-species) %>%
  filter(!is.na(value)) %>%
  mutate(parmF=factor(value, 
                      levels=c('1','nMinAfterMid','Year','Temp',
                               'Clouds','Noise','Replicate','Wind','OBS')))
ggplot()+
  geom_histogram(data=mod_parm,
                 aes(x=parmF), stat='count')+
  coord_flip()+
  scale_x_discrete('Model Parameter')+
  scale_y_continuous(expand=c(0,0))+
  theme_cowplot()
ggsave('Umbrella/results/common model parameters.tiff', 
       width=4, height=4)
