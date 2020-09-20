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
model_parameters<-read.csv('Umbrella/all_possible_models.csv')

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

bird.codes %>%
  filter(SPEC %in% species[!(species %in% bird.assem$`Species code`)])

#Once all the species distance sampling models are built
#pull back in the model information after it has been run -----
model.files<-list.files(path='Umbrella/results/DistanceModels/', full.names=T)
all_models<-NULL
for(b in model.files){
  ms<-read.csv(b)
  all_models<-bind_rows(all_models, ms)
}
all_models<-all_models %>% dplyr::select(-formula, -model.parm)

#Check all needed models were run for all species
models_completed<- all_models %>% 
  filter(species=='BGGN') %>% pull(model.n)
ideal_model_list<-data.frame(species=rep(unique(all_models$species), 
                                         each=length(models_completed)),
                             model.n=rep(models_completed, 28)) %>%
  mutate(model.name=paste(species, model.n, sep='.'))

ideal_model_list %>% group_by(species) %>% tally() %>% arrange(desc(n))
ideal_model_list %>% group_by(model.n) %>% tally() %>% arrange(n)

# Need to calculate AICc - can convert from AIC
# need to know sample size and degrees of freedom. Sample size can be pulled from the raw data
old.model.parm<-read.csv('Umbrella/results/model_parameter_summary20200826.csv')
parmfor.aicc<-old.model.parm %>% 
 dplyr::select(species, truncation.distance) %>%
  rowwise() %>%
  mutate(nobss=nrow(get(paste0(species, '_dat'))%>% 
                      filter(distance<=truncation.distance))) %>%
  rename(tdis='truncation.distance')
#need to figure out PIWA trun.dist

model.appendix<-all_models %>% as_tibble() %>%
  group_by(species) %>%
  left_join(parmfor.aicc, by='species') %>%
  mutate(AICc= AIC + ((2*df*(df + 1))/(nobss - df - 1))) %>%
  filter(AIC >200) %>% #removing bad models (high colinearity among variables)
  mutate(deltaAICc=AICc - min(AICc),
         dAIC=AIC-min(AIC)) %>%
  filter(deltaAICc < 2) %>%
  arrange(AICc) %>%
  left_join(model_parameters, by='model.n') %>%
  dplyr::select(-X)

model.appendix 
write.csv(model.appendix, 'Umbrella/model.appendix.lowAICc.csv')
conversion.factor <- convert_units(distance_units = "Metre", 
                                   effort_units=NULL, 
                                   area_units = "square kilometre")
for(h in 1:nrow(model.appendix)){
  if(exists(model.appendix$model.name[h])){
    print(paste(model.appendix$model.name[h], 'exists'))
  }else{
    assign(model.appendix$model.name[h], 
           ds(get(paste0(model.appendix$species[h], '_dat')),
              transect="point", key="hr",
              formula=as.formula(model.appendix$formula[h]), 
              adjustment = NULL, order = 0,
              truncation = '5%', 
              convert.units = conversion.factor))
    print(paste(model.appendix$model.name[h], 'created'))
  }
}

GOF_mod.app.miss<-NULL
for(h in miss_sp_sub$model.name){
  m<-get(h)
  mod.gof<-gof_ds(m)
  #building species summary information
  gofinf<-data.frame(model.name=h,
                     gof_cvm_w=mod.gof$dsgof$CvM$W,
                     gof_cvm_p=mod.gof$dsgof$CvM$p)
  GOF_mod.app.miss<-bind_rows(GOF_mod.app.miss,gofinf)
}

run_two<-miss_sp_sub %>% left_join(GOF_mod.app.miss) %>%
  group_by(species) %>%
  filter(gof_cvm_p > 0.05)

mod.app.fin<-model.appendix %>% left_join(GOF_mod.app)
write.csv(mod.app.fin, 'Umbrella/model_appendix_wGOF.csv')

done_spp<-mod.app.fin %>% filter(gof_cvm_p > 0.05) %>%
  bind_rows(miss_sp_sub %>% left_join(GOF_mod.app.miss) %>%
              group_by(species) %>%
              filter(gof_cvm_p > 0.05)) %>%
  group_by(species) %>% 
  arrange(desc(gof_cvm_p)) %>%
  slice(1) %>% pull(species)

miss_spp<-unique(model.appendix$species)[which(!(unique(model.appendix$species) %in% done_spp))]

test<-all_models %>% as_tibble() %>%
  group_by(species) %>%
  left_join(parmfor.aicc, by='species') %>%
  mutate(AICc= AIC + ((2*df*(df + 1))/(nobss - df - 1))) %>%
  filter(AIC >200,
         species %in% miss_spp) %>%
  mutate(deltaAICc=AICc-min(AICc))

ggplot()+geom_density(data=test, aes(x=deltaAICc))+
  scale_x_continuous(lim=c(0, 10))

miss_sp_sub<-all_models %>% as_tibble() %>%
  group_by(species) %>%
  left_join(parmfor.aicc, by='species') %>%
  mutate(AICc= AIC + ((2*df*(df + 1))/(nobss - df - 1))) %>%
  filter(AIC >200,
         species %in% miss_spp) %>%
  mutate(deltaAICc=AICc-min(AICc)) %>%
  filter(deltaAICc < 7) %>%
  arrange(desc(deltaAICc)) %>%
  left_join(model_parameters, by='model.n')

for(h in 1:nrow(miss_sp_sub)){
  if(exists(miss_sp_sub$model.name[h])){
    print(paste(miss_sp_sub$model.name[h], 'exists'))
  }else{
    assign(miss_sp_sub$model.name[h], 
           ds(get(paste0(miss_sp_sub$species[h], '_dat')),
              transect="point", key="hr",
              formula=as.formula(miss_sp_sub$formula[h]), 
              adjustment = NULL, order = 0,
              truncation = '5%', 
              convert.units = conversion.factor))
    print(paste(miss_sp_sub$model.name[h], 'created'))
  }
}


# Identify the top model, run it and pull out density & detection info
top.models.spp<-model.appendix %>%
  filter(deltaAICc == 0)
top.models.spp

# create the models if they haven't been created yet
for(h in 1:nrow(top.models.spp)){
  if(exists(top.models.spp$model.name[h])){
    print(paste(top.models.spp$model.name[h], 'exists'))
  }else{
    assign(top.models.spp$model.name[h], 
           ds(get(paste0(top.models.spp$species[h], '_dat')),
              transect="point", key="hr",
              formula=as.formula(top.models.spp$formula[h]), 
              adjustment = NULL, order = 0,
              truncation = '5%', 
              convert.units = conversion.factor))
    print(paste(top.models.spp$model.name[h], 'created'))
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
                   `2.5SE`='X1.1',`3.26SE`='X2.1',`3.76SE`='X3.1',
                   `4.26SE`='X4.1',`4.76SE`='X5.1', TotalSE='X6.1')
names(label.row.names)<-paste('Density', names(label.row.names), sep='.')
#for loop to pull out summary data (phi and density) from each top model
topmod_sums<-NULL
all_site_comp<-NULL
for(p in unique(top.models.spp$species)){
  m<-get(top.models.spp[top.models.spp$species==p,]$model.name)
 
  if(is.null(m$dht)){
    tm_sums<-data.frame(species=p, model.state='broken')
    row.names(tm_sums)<-p
  }else{
    mod.sum<-summary(m)
    mod.gof<-gof_ds(m)
    #building species summary information
    tm_sums<-data.frame(species=p,
                        model.state='fine',
                        gof_cvm_w=mod.gof$dsgof$CvM$W,
                        gof_cvm_p=mod.gof$dsgof$CvM$p,
                        phiMEAN=mod.sum$ds$average.p,
                        phiSE=mod.sum$ds$average.p.se,
                        truncation.distance=mod.sum$ds$width,
                        t(mod.sum$dht$individuals$D[,2]),
                        t(mod.sum$dht$individuals$D[,3])) %>% 
      rename(all_of(label.row.names))
    row.names(tm_sums)<-p
    
    #building site summary information
    site_comp<-mod.sum$dht$individuals$bysample %>% 
      mutate(Site=gsub("\\..*","",Sample)) %>% 
      group_by(Site) %>%
      summarize('Mean.Est.N'=mean(Dhat), .groups='keep',#averaging the replicates here!
                species=p)
    all_site_comp<-bind_rows(all_site_comp, site_comp)
    rm(mod.sum)
  }
  topmod_sums<-bind_rows(topmod_sums, tm_sums)
}
View(topmod_sums)
write.csv(topmod_sums, paste0('Umbrella/results/model_parameter_summary',
                              format(Sys.Date(), "%Y%m%d"), '.csv'))

# build trend tables for species within this list -----
pop.trend<-topmod_sums %>%
  filter(model.state!='broken') %>%
  dplyr::select(-ends_with('SE')) %>%
  dplyr::select(species, starts_with('Density')) %>%
  dplyr::select(-Density.Total) %>%
  pivot_longer(-species) %>%
  mutate(hab_cat=as.numeric(gsub('Density.', '', name))) %>%
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
  arrange(trend)
write.csv(pop.trend, 'Umbrella/results/density_trends_byspp_0908.csv')

#also want to make a graphing table - density & phi at each site ----

#Density
dens_gdf<-all_site_comp %>%
  left_join(DCERP %>%
              group_by(`Pt name`) %>%
              dplyr::select(`Pt name`, USFS_hab_index, Base_Area) %>%
              slice(1),
            by=c('Site'='Pt name')) %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)]) %>%
  left_join(bird.assem, by=c("species"='Species code'))
pdf('Umbrella/results/spp_abundance_0908_hab.pdf')
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
               aes(x=Region.Label, y=Mean.Est.N))
print(longleaf_dens)
hard_dens<-com_plot+
  geom_point(data=dens_gdf[dens_gdf$Habitat.group=='hardwood',], 
             aes(x=USFS_hab_index, y=Mean.Est.N,
                 color=as.factor(Region.Label)))+
  stat_summary(data=dens_gdf[dens_gdf$Habitat.group=='hardwood',],
               aes(x=Region.Label, y=Mean.Est.N))
print(hard_dens)
shrub_dens<-com_plot+
  geom_point(data=dens_gdf[dens_gdf$Habitat.group=='shrub',], 
             aes(x=USFS_hab_index, y=Mean.Est.N,
                 color=as.factor(Region.Label)))+
  stat_summary(data=dens_gdf[dens_gdf$Habitat.group=='shrub',],
               aes(x=Region.Label, y=Mean.Est.N))
print(shrub_dens)
gen_dens<-com_plot+
  geom_point(data=dens_gdf[dens_gdf$Habitat.group=='generalist',], 
             aes(x=USFS_hab_index, y=Mean.Est.N,
                 color=as.factor(Region.Label)))+
  stat_summary(data=dens_gdf[dens_gdf$Habitat.group=='generalist',],
               aes(x=Region.Label, y=Mean.Est.N))
print(gen_dens)
dev.off()

#summary table for manuscript ----
bird.codes %>%
  filter(SPEC %in% species) %>%
  rename(species="SPEC") %>%
  left_join(topmod_sums) %>%
  dplyr::select(-ends_with('SE'))

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
  dplyr::select(Site, USFS_hab_index, Base_Area, Region.Label, everything()) #reordering

site.com<-ordination_df[,-c(2:4,6)] %>%
  bind_cols(., total.dens=rowSums(.[,-1])) %>%
  filter(total.dens != 0) %>%
  dplyr::select(-total.dens) %>%
  ungroup()
nrow(site.com); head(site.com)
library(vegan)
dis<-vegdist(site.com[,-1], method="bray") #create a distance matrix for the bird communities

site.env<-ordination_df[,c(1:4,6)] %>%
  filter(Site %in% site.com$Site)

#run the fuzzy ordination
fuzzy_ord<-fso(~USFS_hab_index, dis, site.env, permute=1000)
#mu - fuzzy membership values for each site
#r - correlation between original vector and fuzzy set
#d - correlation of pair-wise distances among each fuzzy set
summary(fuzzy_ord)
ord_gdf<- ordination_df %>%
  bind_cols(mu=fuzzy_ord$mu[,1]) %>%
  dplyr::select(Site, mu, everything())
ggplot()+
  geom_point(data=ord_gdf, aes(x=USFS_hab_index, y=mu), alpha=0.6)+
  geom_text(aes(x=2.8, y=0.545, label="r = 0.566, p = 0.001"))+
  ggtitle('Fit of USFS Habitat Index ordination')+
  theme_bw()

site.ordered<-ord_gdf$Site[order(ord_gdf$mu)]
ord_gdf_long<-ord_gdf %>%
  pivot_longer(cols=c(AMCR:NOCA), 
               names_to="species", values_to="Birds/sqkm")%>%
  mutate(Site.ord=factor(Site, levels=site.ordered)) %>%
  left_join(bird.assem, by=c('species'='Species code')) %>%
  mutate(HabitatF=factor(Habitat.group, levels=c('longleaf', 'hardwood', 'shrub','generalist'))) %>%
  group_by(Site) %>%
  mutate(rel.abundance=(`Birds/sqkm`/sum(`Birds/sqkm`))*100)

(site.comp.plot<-ggplot()+
  geom_col(data=ord_gdf_long, 
           aes(x=Site.ord, y=`Birds/sqkm`,fill=HabitatF))+
  scale_fill_viridis_d(name="Habitat\nassociation")+
  scale_x_discrete(name="Site")+
  theme(axis.text.x = element_text(angle=90, hjust=0.9),
        axis.text.y = element_text(angle=90, hjust=0.5))+
  coord_flip()+theme_cowplot()+
  theme(axis.text.y = element_text(size=5)))
ggsave('Umbrella/results/site_comp_plot_abun_0908.tiff', site.comp.plot, 
       width=5, height=9)
(site.ra.plot<-ggplot()+
    geom_col(data=ord_gdf_long, 
             aes(x=Site.ord, y=rel.abundance,fill=HabitatF))+
    scale_fill_viridis_d(name="Habitat\nassociation")+
    scale_x_discrete(name="Site")+
    theme(axis.text.x = element_text(angle=90, hjust=0.9),
          axis.text.y = element_text(angle=90, hjust=0.5))+
    coord_flip()+theme_cowplot()+
    theme(axis.text.y = element_text(size=5)))
ggsave('Umbrella/results/site_comp_plot_relabun_0908.tiff', site.ra.plot, 
       width=5, height=9)

abun_hab_df<-ord_gdf_long %>%
  group_by(Site, USFS_hab_index, Region.Label, HabitatF) %>%
  summarize(Total.Birds.sqkm=sum(`Birds/sqkm`), .groups='keep')
abun_all_df<-ord_gdf_long %>%
  group_by(Site, USFS_hab_index, Region.Label) %>%
  summarize(Total.Birds.sqkm=sum(`Birds/sqkm`), .groups='keep')
hab_abund_plot<-ggplot()+
  geom_jitter(data=abun_hab_df,
             aes(x=Region.Label, y=Total.Birds.sqkm, group=Region.Label),
             alpha=0.4)+
  geom_boxplot(data = abun_hab_df, 
               aes(x=Region.Label, y=Total.Birds.sqkm, group=Region.Label),
               outlier.alpha = 0, alpha=0.7, fill='lightgrey')+
  theme_cowplot()+facet_wrap(~HabitatF, scales='free_y')+
  theme(axis.title.y=element_text(size=0))
abund_plot<-ggplot()+
  geom_jitter(data=abun_all_df,
              aes(x=Region.Label, y=Total.Birds.sqkm, group=Region.Label),
              alpha=0.4)+
  geom_boxplot(data = abun_all_df, 
               aes(x=Region.Label, y=Total.Birds.sqkm, group=Region.Label),
               outlier.alpha = 0, alpha=0.7, fill='lightgrey')+
  theme_cowplot()
plot_grid(abund_plot, hab_abund_plot, labels=c("   All Birds", ""))
ggsave('Umbrella/results/bird_abundance.tiff', width=6, height=4)
