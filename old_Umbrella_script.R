### Analyzing Distance Sampling data ###
library(tidyverse); library(readxl); library(Distance); library(lubridate)

files<-list.files("G:/Shared drives/RCW Conservation Postdoc/Umbrella species paper/Data",
                  full.names = T)
files #what files do we have in the Data folder

#this dataframe has plot abbreviations (Pt name), plot full names, plot type, and training area
plot_types<-read_excel(files[3]) 
unique(plot_types$Type)

#dataframe of USFWS_hab_index bins
hab_bins<-data.frame(min=c(2,3.01,3.51,4.01,4.51),
                     max=c(3,3.5,4,4.5,5)) %>%
  mutate(mid=round(((max-min)/2+min),2))

#2009 point distance sampling raw data
DCERP<-bind_rows(read_excel(files[1]), #2009 data
                      read_excel(files[2])) %>% #2010 data
  select(-`Area...2`, -`...24`, -`...25`) %>% #empty columns 
  rename(Base_Area='Area...6', #cleaning up name
         # since each Pt name is visited on 4 dates with unique round, indicates replicate
         Replicate='Round (1-4)',
         #I think Label is the USFWS RCW habitat index at that spot
         USFS_hab_index='Label') %>%
  dplyr::select(`Pt name`, Replicate, Date1, everything()) %>%
  #adding information needed by Distance
  mutate(Study.Area='MCBCL',
         Start.time= strftime(`Start time`, format="%H:%M", 
                              tz="UTC", uzetz=T),
         hoursmins=hm(Start.time),        # format to 'hours:minutes:seconds'
         nMinAfterMid=hour(hoursmins)*60 + minute(hoursmins), 
         Area=1,
         Year=year(Date1)) #estimating area size in km2
names(DCERP)
#don't know what label stands for. 
#`Pt name`:notes & Replicate would be the sample data layer, 
#`Start time`:`Cluster size` would be the observation data layer
DCERP %>%
  group_by(`Pt name`) %>%
  count(USFS_hab_index)

#Avian detections beyond 300m and from flyovers were omitted 
#subsequently species that represented less than 0.1% of all bird observations were removed.  
rare_birds <- DCERP %>% 
  group_by(Species) %>% 
  tally() %>%
  mutate(percentObs=n/sum(n)*100) %>%
  arrange(percentObs) %>%
  filter(percentObs < 0.1)

DCERP_filt <- DCERP %>%
  filter(!(`Radial distance` %in% c('?','UK','FO', 'F','NA'))) %>%
  mutate(distance=as.numeric(`Radial distance`),
         Year=year(Date1)) %>%
  filter(distance<300,
         !(Species %in% rare_birds$Species))

nrow(DCERP)-nrow(DCERP_filt)
#798 observations lost... need to make sure I didn't lose any samples?

DCERP_filt %>% group_by(`Pt name`) %>%
  summarize(n.reps=length(unique(Replicate))) %>%
  arrange(n.reps)
#lost one sample of PHOTO63, 149 sites remain total
#not going to add it back in because the notes say 'logging operation no count'
DCERP %>% group_by(`Pt name`, Replicate) %>%
  filter(`Pt name`=="PHOTO63", Replicate==4) %>%
  pull(notes)

#observations by distance histogram
ggplot()+geom_histogram(data=DCERP_filt, aes(x=distance))
hist_gdf<-hist(DCERP_filt$distance, plot = F)
gdf<-data.frame(counts=hist_gdf$counts,
                radius=hist_gdf$mids,
                area=pi*(hist_gdf$mids^2)) %>%
  mutate(scaled_counts=counts/area)
ggplot()+geom_col(data=gdf, aes(x=radius, y=scaled_counts))+
  scale_y_continuous(name="frequency/area")+
  scale_x_continuous(name="radial distance")


#creating an empty 'sample layer' dataframe to account for search area
point_sample_empty <- DCERP_filt %>% 
  group_by(`Pt name`, Replicate, Year) %>%
  slice(1) %>%
  mutate(Species=NA,
         `Radial distance`=NA,
         `Type AVB`=NA,
         `0-2`=0,
         `2-4`=0,
         `4-6`=0,
         `6-8`=0,
         `Cluster size`=NA,
         distance=NA) %>%
  ungroup()

# creating a function to filter the data & supplement rows for samples with no observations
spfilter<-function(sp_abbr){
  
  #filter the data to just include the species of interest
  sp_data <- DCERP_filt %>% filter(Species==sp_abbr)
  
  #identify samples that do include the species of interest
  sampled_pts <- sp_data %>%
    group_by(`Pt name`, Replicate) %>%
    slice(1)
  
  #pull from point_sample_empty samples that are not found within sampled_pts
  needed_rows<-union(sampled_pts,
                     point_sample_empty)%>%
    group_by(`Pt name`, Replicate) %>%
    arrange(`Pt name`, Replicate)%>%
    slice(1) %>%
    filter(is.na(distance))
  
  #bind rows to capture all species observations and empty rows for sampled areas
  #minimum number of rows is 596
  output<-rbind(sp_data, needed_rows) %>%
    mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                  left.open = F)],
           Sample.Label=paste(`Pt name`, Replicate, sep="."), #this identifies unique samples within strata
           size=`Cluster size`,
           Effort=1,
           Study.Area="MCBCL") %>%
    filter(USFS_hab_index !=0) %>%
    as.data.frame()
  print("Empty Samples + Not Empty Samples = 596")
  print(paste("for", sp_abbr, nrow(needed_rows),"+", nrow(sampled_pts),"=", nrow(needed_rows)+nrow(sampled_pts)))
  return(output)
}

#building all the data sets
#species<-paste0(unique(DCERP_filt$Species), '_dat')
length(unique(DCERP_filt$Species)) #59 species in the data
species<-c('NOCA','NOBO','EATO','PIWA','EAWP','BHCO','MODO',
           'RHWO','PRAW','WEVI','NOFL','YTWA','RBWO','INBU',
           'SUTA','CACH','COYE','BLJA','RCWO','AMCR','BACS',
           'GCFL','BGGN','TUTI','CARW','BHNU','CHSP','EABL')
un_run_spp<-unique(DCERP_filt$Species)[!(unique(DCERP_filt$Species) %in% species)]
for(j in unique(rerun_these_spmod$species)){
  assign(paste0(j, '_dat'),
         spfilter(j))
}

#run the models ####
# so Distance is just a wrapper for mrds it seems
conversion.factor <- convert_units(distance_units = "Metre", 
                                   effort_units=NULL, 
                                   area_units = "square kilometre")

# temperature, cloud cover, wind, noise, start time, observer, replicate
variables<-c('Temp', 'Clouds', 'Wind', 'Noise', 'nMinAfterMid', 'OBS', 'Replicate', 'Year')
var.correl<-cor(DCERP %>%
                  dplyr::select(all_of(variables), -OBS)%>%
                  filter(!is.na(Temp)) %>%
                  data.frame())
library(corrplot)
corrplot.mixed(var.correl)

var.combos<-data.frame(data.frame(rbind(combn(variables, 2),NA, NA, NA, NA, NA, NA),
                                  rbind(combn(variables, 3),NA, NA, NA, NA, NA),
                                  rbind(combn(variables, 4),NA, NA, NA, NA),
                                  rbind(combn(variables, 5),NA, NA, NA),
                                  rbind(combn(variables, 6),NA, NA),
                                  rbind(combn(variables, 7),NA),
                                  rbind(combn(variables, 8))) %>%
             t())
models_pre<-unite(var.combos, 'formula', sep='+', na.rm=T) %>%
  summarize(formula=paste0('~', formula)) %>%
  bind_rows(data.frame(formula=paste('~', c(variables,'1')))) %>%
  mutate(model.n=paste0('mod', row_number()),
         model.parm=gsub('+', '.', gsub('~','', formula),
                                       fixed=T))
head(models_pre)
models_subset<-models_pre %>%
  filter(!(model.parm %in% DONT_run_these_models$model.parm))

nrow(models_subset) #models total


write.csv(models_subset,'Umbrella results/all_models_tested_0830.csv', row.names=F)
#need to check out COYE, EATO, and PIWA again to complete models
species_sub<- species[-c(1:17, 23,22)]
#NOTE NEED TO DO EATO models not completed
for(u in species_sub){
  for(w in 1:length(models_subset$model.n)){
    assign(paste(u, models_subset$model.n[w], sep='.'), 
           ds(get(paste0(u, '_dat')),
              transect="point", key="hr",
              formula=as.formula(models_subset$formula[w]), 
              adjustment = NULL, order = 0,
              truncation = '5%', 
              convert.units = conversion.factor,
              er.var='P3' #point encounter rate
              ))
    print(paste(u,'species number', which(species==u), 
                'model number', w))
    model.summaries<-data.frame(species=rep(u, w),
                                model.n=models_subset$model.n[1:w]) %>%
      rowwise() %>%
      mutate(model.name=paste0(species, '.', model.n),
             df=as.numeric(AIC(get(model.name))[1]),
             AIC=as.numeric(round(AIC(get(model.name))[2],3)))
    write.csv(model.summaries, 
              file=paste0('Umbrella results/DistanceModels_0828/',u,'_models.csv'), 
              row.names=F)
  }
}

#Once all the species distance sampling models are built
#pull back in the model information after it has been run
model.files<-list.files(path='Umbrella results/DistanceModels/',
           full.names=T)
all_model_aic<-NULL
for(b in model.files){
  ms<-read.csv(b)
  all_model_aic<-bind_rows(all_model_aic, ms)
}

#investigate for abnormally low aic values - something would be wrong in these models
#Error message: 'Some variance-covariance matrix elements were NA, possible numerical problems; only estimating detection function.'
(model.errors<-all_model_aic %>% filter(AIC < 300))
all_model_aic <- all_model_aic %>% 
  group_by(species) %>%
  filter(AIC > 100) %>% #remove the bad models - likely should rerun??? 
  arrange(AIC) %>%
  mutate(deltaAIC=round(AIC-min(AIC),4))

top.models.spp<-all_model_aic %>%
  filter(deltaAIC == 0) %>%
  left_join(models, by='model.n')
top.models.spp

# create the models if they haven't been created yet
for(h in 1:nrow(top.models.spp)){
  if(exists(top.models.spp$model.name[h])){
    print(paste(top.models.spp$model.name[h], 'exists'))
  }else{
    assign(paste0(top.models.spp$species[h], '_dat'), 
           spfilter(top.models.spp$species[h]))
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

#want to pull out the top models for each species, isolate strata phi and density
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
    #building species summary information
    tm_sums<-data.frame(species=p,
                        model.state='fine',
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
write.csv(topmod_sums, paste0('Umbrella results/model_parameter_summary',
                              format(Sys.Date(), "%Y%m%d"), '.csv'))


# build trend tables for species within this list
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
write.csv(pop.trend, 'Umbrella results/density_trends_byspp.csv')

#also want to make a graphing table - density & phi at each site

#Density
dens_gdf<-all_site_comp %>%
  left_join(DCERP %>%
              group_by(`Pt name`) %>%
              dplyr::select(`Pt name`, USFS_hab_index, Base_Area) %>%
              slice(1),
            by=c('Site'='Pt name')) %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)])
pdf('Umbrella results/spp_abundance.pdf')
ggplot()+
  geom_point(data=dens_gdf, aes(x=USFS_hab_index, y=Mean.Est.N,
                                   color=as.factor(Region.Label)))+
  stat_summary(data=dens_gdf, aes(x=Region.Label, y=Mean.Est.N))+
  scale_color_viridis_d('RCWO Habitat\nMatrix Score',
                        option='plasma',
                        end=0.8,
                        alpha=0.6)+
  scale_y_continuous(name='Estimate Birds * sq. km')+
  theme_bw()+theme(legend.position = 'bottom')+
  facet_wrap(~species, scales='free_y')
dev.off()

#summary table for manuscript ----
#bird codes from: https://www.birdpop.org/pages/birdSpeciesCodes.php 
bird.codes<-read.csv('Umbrella results/Cornell_list19p.csv')
bird.codes %>%
  filter(SPEC %in% species) %>%
  rename(species="SPEC") %>%
  left_join(topmod_sums) %>%
  dplyr::select(-ends_with('SE'))


# Habitat Associate groups from Allen et al. 2006 (the Auk)
longleaf<-c('Pine Warbler', 'Brown-headed Nuthatch','Red-cockaded Woodpecker',
            'Prairie Warbler',"Bachman's Sparrow", 'Chipping Sparrow',
            'Eastern Wood-Pewee', 'Red-headed Woodpecker')
firesup<-c('Red-eyed Vireo', 'Acadian Flycatcher', 'Ovenbird', 
           'Black-and-white Warbler', 'Tufted Titmouse', 'Wood Thrush', 'Yellow-throated Vireo')
pocosin<-c('Eastern Towhee', 'Common Yellowthroat', 'Carolina Wren', 'Northern Cardinal', 
           'White-eyed Vireo','Hooded Warbler')
generalist<-c('Carolina Chickadee', 'Summer Tanager', 'Great Crested Flycatcher', 'Blue Jay',
              'Red-bellied Woodpecker', 'Northern Flicker')
unclear<-c('Blue-gray Gnatcatcher', 'Brown-headed Cowbird', 'Yellow-throated Warbler', 'Indigo Bunting',
           'American Goldfinch', 'Mourning Dove', 'Eastern Bluebird')
bird.assem<-data.frame(COMMONNAME=c(longleaf, firesup, pocosin, generalist, unclear),
                       `Habitat group`=c(rep('longleaf', length(longleaf)),
                                   rep('firesup', length(firesup)),
                                   rep('pocosin', length(pocosin)),
                                   rep('generalist', length(generalist)),
                                   rep('unclear', length(unclear)))) %>%
  left_join(bird.codes) %>%
  filter(SPEC %in% species) %>%
  dplyr::select(-SP, -CONF, -SPEC6, -CONF6) %>%
  rename(`Common name`='COMMONNAME', `Scientific name`='SCINAME',
         `Species code`='SPEC') %>%
  dplyr::select(`Scientific name`, `Common name`, `Species code`, Habitat.group)
write.csv(bird.assem, 'Umbrella results/species_appendix.csv')

bird.codes %>%
  filter(SPEC %in% species[!(species %in% bird.assem$`Species code`)])

# fuzzy ordination to assess alignment with USFS habitat index ####
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

site.com<-ordination_df[,-c(2:4)] %>%
  bind_cols(., total.dens=rowSums(.[,-1])) %>%
  filter(total.dens != 0) %>%
  dplyr::select(-total.dens) %>%
  ungroup()
nrow(site.com); head(site.com)
library(vegan)
dis<-vegdist(site.com[,-1], method="bray") #create a distance matrix for the bird communities

site.env<-ordination_df[,1:4] %>%
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
  theme_bw()

site.ordered<-ord_gdf$Site[order(ord_gdf$mu)]
ord_gdf_long<-ord_gdf %>%
  pivot_longer(cols=c(AMCR:NOCA), 
               names_to="species", values_to="Birds/sqkm")%>%
  mutate(Site.ord=factor(Site, levels=site.ordered))

site.comp.plot<-ggplot()+
  geom_col(data=ord_gdf_long, 
           aes(x=Site.ord, y=`Birds/sqkm`,fill=species))+
  scale_fill_viridis_d(guide=F)+
  scale_x_discrete(name="Site")+
  theme(axis.text.x = element_text(angle=90, hjust=0.9),
        axis.text.y = element_text(angle=90, hjust=0.5))+
  coord_flip()+theme_cowplot()
ggsave('Umbrella results/site_comp_plot.tiff', site.comp.plot, 
       width=5, height=9)

# Abundance plots ------
common_plot_stuff<-ggplot()+
  scale_color_viridis_d('RCWO Habitat\nMatrix Score',
                        option='plasma',
                        end=0.8,
                        alpha=0.6)+
  theme_bw()

# NOCA !has legend
noca.mod.res<-summary(NOCA_dat_obs_hab_noise_rep_mod)

#noca.mod.res$dht$individuals$D #estimates of density for each region and total
#noca.mod.res$ds$average.p #average phi

noca_graph_df<-noca.mod.res$dht$individuals$Nhat.by.sample %>% 
  mutate(Site=gsub("\\..*","",Sample.Label)) %>% 
  group_by(Site) %>%
  summarize(Mean.Est.N=mean(Nhat)) %>%
  left_join(DCERP %>%
              group_by(`Pt name`) %>%
              slice(1),
            by=c('Site'='Pt name')) %>% 
  select(Site, Mean.Est.N, USFS_hab_index) %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                         left.open = F)])
NOCA.hab.plot<-common_plot_stuff+
  geom_point(data=noca_graph_df, aes(x=USFS_hab_index, y=Mean.Est.N,
                                color=as.factor(Region.Label)))+
  stat_summary(data=noca_graph_df, aes(x=Region.Label, y=Mean.Est.N))+
  theme(legend.direction = 'horizontal')+
  ggtitle("NOCA abundance")
#NOCA.hab.plot+theme(legend.position='none')

# NOBO plot ----
nobo.mod.res<-summary(NOBO_dat_obs_hab_noise_rep_mod)

#nobo.mod.res$dht$individuals$D #estimates of density for each region and total
#nobo.mod.res$ds$average.p #average phi

nobo_graph_df<-nobo.mod.res$dht$individuals$Nhat.by.sample %>% 
  mutate(Site=gsub("\\..*","",Sample.Label)) %>% 
  group_by(Site) %>%
  summarize(Mean.Est.N=mean(Nhat)) %>%
  left_join(DCERP %>%
              group_by(`Pt name`) %>%
              slice(1),
            by=c('Site'='Pt name')) %>% 
  select(Site, Mean.Est.N, USFS_hab_index) %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)])
NOBO.hab.plot<-common_plot_stuff+
  geom_point(data=nobo_graph_df, aes(x=USFS_hab_index, y=Mean.Est.N,
                                color=as.factor(Region.Label)))+
  stat_summary(data=nobo_graph_df, aes(x=Region.Label, y=Mean.Est.N))+
  theme(legend.position='none')+
  ggtitle("NOBO abundance")

# PIWA plot ----
piwa.mod.res<-summary(PIWA_dat_obs_mod)

#piwa.mod.res$dht$individuals$D #estimates of density for each region and total
#piwa.mod.res$ds$average.p #average phi

piwa_graph_df<-piwa.mod.res$dht$individuals$Nhat.by.sample %>% 
  mutate(Site=gsub("\\..*","",Sample.Label)) %>% 
  group_by(Site) %>%
  summarize(Mean.Est.N=mean(Nhat)) %>%
  left_join(DCERP %>%
              group_by(`Pt name`) %>%
              slice(1),
            by=c('Site'='Pt name')) %>% 
  select(Site, Mean.Est.N, USFS_hab_index) %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)])
PIWA.hab.plot<-common_plot_stuff+
  geom_point(data=piwa_graph_df, aes(x=USFS_hab_index, y=Mean.Est.N,
                                     color=as.factor(Region.Label)))+
  stat_summary(data=piwa_graph_df, aes(x=Region.Label, y=Mean.Est.N))+
  theme(legend.position = 'none')+
  ggtitle("PIWA abundance")

# EATO plot ----
eato.mod.res<-summary(EATO_dat_obs_hab_noise_rep_mod)

#eato.mod.res$dht$individuals$D #estimates of density for each region and total
#eato.mod.res$ds$average.p #average phi

eato_graph_df<-eato.mod.res$dht$individuals$Nhat.by.sample %>% 
  mutate(Site=gsub("\\..*","",Sample.Label)) %>% 
  group_by(Site) %>%
  summarize(Mean.Est.N=mean(Nhat)) %>%
  left_join(DCERP %>%
              group_by(`Pt name`) %>%
              slice(1),
            by=c('Site'='Pt name')) %>% 
  select(Site, Mean.Est.N, USFS_hab_index) %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)])
EATO.hab.plot<-common_plot_stuff+
  geom_point(data=eato_graph_df, aes(x=USFS_hab_index, y=Mean.Est.N,
                                     color=as.factor(Region.Label)))+
  stat_summary(data=eato_graph_df, aes(x=Region.Label, y=Mean.Est.N))+
  theme(legend.position = 'none')+
  ggtitle("EATO abundance")

library(cowplot)
DCERP %>% group_by(Species) %>% tally() %>% arrange(desc(n))
legend <- get_legend(NOCA.hab.plot)
plot_grid(plot_grid(NOCA.hab.plot+theme(legend.position = 'none'),
                    EATO.hab.plot, 
                    NOBO.hab.plot, 
                    PIWA.hab.plot),
          legend,
          ncol=1, rel_heights = c(1,.15))
ggsave("prelim_distance_res.jpg", width=5, height=5)

#quantifying fit ----
summarize_ds_models(NOCA_dat_obs_hab_noise_rep_mod,
                    output="plain")
gof_ds(NOCA_dat_obs_hab_noise_rep_mod,
       chisq = T)

# Presence -------
# https://www.mbr-pwrc.usgs.gov/software/presence.html
# accessed 18 Aug 2020
install.packages('C:/Users/Owner/Downloads/RPresence.zip', 
                 repos = NULL, type = "win.binary")
