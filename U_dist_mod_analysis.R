# Analyze Distance Sampling Models #####
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
  dplyr::select(-X) 
old_model_param<-read.csv('Umbrella/all_possible_models_0929.csv') %>%
  mutate(model.n=paste0('mod',row_number()))

#Once all the species distance sampling models are built
#pull back in the model information after it has been run -----
model.files<-list.files(path='Umbrella/DistanceModels_1005/', #path='Umbrella/results/DistanceModels/',
                        pattern='.csv',
                        full.names=T)
all_models<-NULL
for(b in model.files){
  if(substr(b, 30,33)!='hab_'){  
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
#two models that wouldn't run

# identify which models explained detection well
model.appendix<-all_models %>% as_tibble() %>%
  group_by(species) %>%
  filter(gof_cvm_p >= 0.050,
         AIC > 0) %>%
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

model.appendix.all<-bind_rows(model.appendix, best_bad_models) %>%
  group_by(species) %>% mutate(model.rank=rank(deltaAICc))

write.csv(model.appendix.all, 'Umbrella/model.appendix.20201116.csv')

model.appendix.all<-read.csv('Umbrella/model.appendix.20201116.csv')
# Figure S2 parameter barchart ----
mod_parm<-model.appendix.all %>% filter(model.rank==1) %>%
  dplyr::select(species, formula) %>%
  separate(formula, into=paste0('X', 1:7)) %>%
  #arrange(X7)%>%
  dplyr::select(-X1) %>%
  pivot_longer(-species) %>%
  filter(!is.na(value))%>%
  mutate(value=recode(value,
                      'nMinAfterMid'='nMin',
                      'HabCat'='RCW score')) %>%
  mutate(parmF=factor(value, 
                      levels=c('1','nMin','Temp','Year','RCW score',
                               'Clouds','Noise','Replicate','Wind','OBS')))

hab_pred_detection<-mod_parm %>% filter(parmF=='RCW score') %>% 
  left_join(bird.assem, by=c('species'='Species code')) %>%
  arrange(Habitat.group) %>%
  dplyr::select(species, `Common name`, Habitat.group, value)

model.appendix.all %>%
  dplyr::select(species, formula) %>%
  separate(formula, into=paste0('X', 1:7)) %>%
  #arrange(X7)%>%
  dplyr::select(-X1) %>%
  pivot_longer(-species) %>%
  filter(!is.na(value)) %>%
  filter(value=="HabCat") %>%
  group_by(species) %>% slice(1) %>%
  filter(!(species %in% hab_pred_detection$species))

dist.det.p<-ggplot()+
  geom_histogram(data=mod_parm,
                 aes(x=parmF), stat='count')+
  coord_flip()+
  scale_x_discrete('')+
  scale_y_continuous(expand=c(0,0))+
  theme_cowplot()+
  theme(axis.title.y=element_blank())
ggsave('Umbrella/results/Figures/FigS2_onlyDist.jpg', width=2.75, height=2.75)
#read in and create occupancy model parameter graph
top_occ_mod <- read.csv('Umbrella/SimpleOcc/top_occupancy_models.csv')
occ.det.p<-new_mod_info %>% 
  filter(parameter=='psi.USFS_s') %>%
  select(spp, p.cloud.:df, AICc) %>% 
  mutate(across(starts_with('p.'), .fns=as.character)) %>%   
  as_tibble()%>%
  pivot_longer(cols=c(p.cloud.:p.wind.)) %>%
  filter(!is.na(value)) %>% 
  bind_rows(data.frame(spp=c("AMCR","NOFL","EATO"),
                       name=as.character(1))) %>%
  mutate(parm.name=recode(name, p.wind.='Wind', p.USFS_s.='RCW score',
                          p.temp.='Temp', p.rep.='Replicate', p.obs.='OBS',
                          p.noise.='Noise', p.minAmid.='nMin', p.cloud.='Clouds'),
         parmF=factor(parm.name, levels=c('1','Clouds', 'nMin', 'Wind','Noise','OBS','Temp', 
                                          'Replicate'))) %>%
  ggplot()+
  geom_histogram(aes(x=parmF), stat='count')+
  coord_flip()+
  scale_x_discrete('Model Parameter')+
  scale_y_continuous(expand=c(0,0))+
  theme_cowplot()

plot_grid(occ.det.p, dist.det.p, nrow=1, labels='AUTO')
ggsave('Umbrella/results/Figures/Figs1 common model parameters.tiff', 
       width=6, height=3)

# count how often and in which direction RCW score predicts detection 
all_hab_coef_detect<-NULL
for(i in list.files('Umbrella/results/model_rds/')[-c(16,39)]){
  if(substr(i, 1,4) %in% hab_pred_detection$species){
    moddd<-readRDS(paste0('Umbrella/results/model_rds/',i))
    test<-summary(moddd)
    hab_coef<-test$ds$coeff$key.scale[grep('HabCat', row.names(test$ds$coeff$key.scale)),] %>% 
      mutate(species=substr(i, 1,4),
             model=sub('_','',substr(i, 1,11)))
    all_hab_coef_detect<-bind_rows(hab_coef, all_hab_coef_detect)
  }
}
cat<-c("HabCat3.26","HabCat3.76","HabCat4.26","HabCat4.76")
all_hab_coef_detect$betasss<-cat
all_hab_coef_detect %>%
  group_by(species) %>% summarize(mean(estimate))
all_hab_coef_detect %>% 
  mutate(direction=ifelse(estimate <=0, 'neg', 'pos')) %>% 
  count(direction)

# Table S3 distance sampling model formulas and info ----
mod_used<-read.csv('C:/Users/Owner/Downloads/models_20201016.csv')
mod_used %>% group_by(species) %>%
  arrange(AIC) %>% slice(1) %>%
  left_join(bird.assem, by=c('species'='Species code')) %>%
  select(Habitat.group, species, formula, trun.dist, nobss,
         df, AICc, gof_cvm_w, gof_cvm_p) %>%
  arrange(Habitat.group, species) %>%
  write.csv('Umbrella/results/TS3_dist_detect_model_info.csv')

mod_used %>% filter(gof_cvm_p < 0.05)


