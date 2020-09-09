# Run the distance sampling models ####
library(tidyverse); library(readxl); library(Distance); library(lubridate)

#building all the data sets
length(unique(DCERP_filt$Species)) #59 species in the data
species<-c('NOCA','NOBO','EATO','PIWA','EAWP','BHCO','MODO',
           'RHWO','PRAW','WEVI','NOFL','YTWA','RBWO','INBU',
           'SUTA','CACH','COYE','BLJA','RCWO','AMCR','BACS',
           'GCFL','BGGN','TUTI','CARW','BHNU','CHSP','EABL')
un_run_spp<-unique(DCERP_filt$Species)[!(unique(DCERP_filt$Species) %in% species)]

for(j in species){
  assign(paste0(j, '_dat'),
         spfilter(j))
}

#run the models ####
# so Distance is just a wrapper for mrds it seems
variables<-c('Temp', 'Clouds', 'Wind', 'Noise', 'nMinAfterMid', 'OBS', 'Replicate', 'Year')
# temperature, cloud cover, wind, noise, start time, observer, replicate
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
head(models_pre) #all the models possible
write.csv(models_pre, 'all_possible_models.csv')

# removing models that commonly break ####
library(corrplot)
var.correl<-cor(DCERP %>%
                  dplyr::select(all_of(variables), -OBS)%>%
                  filter(!is.na(Temp)) %>%
                  data.frame())
corrplot.mixed(var.correl)

# To reduce the number of models run, I ran all models on four important species
# I then removed models that had artificially low AIC in half of the species

aug27.models<-read.csv('Umbrella/results/all_models_tested_0827.csv') #parameters
new.NOCA<-read.csv('Umbrella/results/DistanceModels_0827/NOCA_models.csv') %>%
  mutate(run='Aug27') %>%
  left_join(aug27.models) %>%
  dplyr::select(-model.n)
NOCA_bad_models<-filter(new.NOCA, AIC<100) #which models specifically do not run well?

new.NOBO<-read.csv('Umbrella/results/DistanceModels_0827/NOBO_models.csv') %>%
  mutate(run='Aug27') %>%
  left_join(aug27.models) %>%
  dplyr::select(-model.n)
NOBO_bad_models <- filter(new.NOBO, AIC<100)

new.BGGN<-read.csv('Umbrella/results/DistanceModels_0828/BGGN_models.csv') %>%
  mutate(run='Aug27') %>%
  left_join(aug27.models) %>%
  dplyr::select(-model.n)
BGGN_bad_models <- filter(new.BGGN, AIC<100)

new.GCFL<-read.csv('Umbrella/results/DistanceModels_0828/GCFL_models.csv') %>%
  mutate(run='Aug27') %>%
  left_join(aug27.models) %>%
  dplyr::select(-model.n)
GCFL_bad_models <- filter(new.GCFL, AIC<100)

DONT_run_these_models<-bind_rows(NOBO_bad_models, NOCA_bad_models, GCFL_bad_models, BGGN_bad_models) %>%
  count(model.parm) %>%
  filter(n>=2) %>%
  mutate(reason='variance/covariance problems')
DONT_run_these_models # these models don't work well

# models we are running today
models_subset<-models_pre %>%
  filter(!(model.parm %in% DONT_run_these_models$model.parm))

nrow(models_subset) #models total

write.csv(models_subset, 
  paste0('Umbrella/results/all_models_tested', format(Sys.Date(), '%m%d'),'.csv'), row.names=F)

conversion.factor <- convert_units(distance_units = "Metre", 
                                   effort_units=NULL, 
                                   area_units = "square kilometre")

for(u in species){
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
              file=paste0('Umbrella/results/DistanceModels_', format(Sys.Date(), '%m%d'),'/',u,'_models.csv'), 
              row.names=F)
  }
}

#investigate for abnormally low aic values - something would be wrong in these models
#Error message: 'Some variance-covariance matrix elements were NA, possible numerical problems; only estimating detection function.'
(model.errors<-all_model_aic %>% filter(AIC < 300))

