#investigate models 
#do models match between first run and second run
aug25.models<-read.csv('Umbrella results/all_models_tested.csv')
original.NOCA<-read.csv('Umbrella results/DistanceModels_0825/NOCA_models.csv') %>%
  mutate(run='Aug25') %>%
  left_join(aug25.models) %>%
  dplyr::select(-model.n, -deltaAIC)
head(original.NOCA)

aug27.models<-read.csv('Umbrella_results/all_models_tested_0827.csv')
new.NOCA<-read.csv('Umbrella_results/DistanceModels_0827/NOCA_models.csv') %>%
  mutate(run='Aug27') %>%
  left_join(aug27.models) %>%
  dplyr::select(-model.n)
nrow(new.NOCA)
head(new.NOCA)

all_NOCA_mod<-left_join(new.NOCA, original.NOCA, 
                        by=c('formula', 'model.parm','species', 'df'),
                        suffix=c('.new', '.orig'))
all_NOCA_mod %>% mutate(AIC_dif=AIC.new-AIC.orig) %>%
  filter(!is.na(run.orig))
View(all_NOCA_mod)
#SO! all AIC are 0, do not need to rerun those models

#which models specifically do not run well?
NOCA_bad_models<-filter(new.NOCA, AIC<100)

new.NOBO<-read.csv('Umbrella_results/DistanceModels_0827/NOBO_models.csv') %>%
  mutate(run='Aug27') %>%
  left_join(aug27.models) %>%
  dplyr::select(-model.n)
head(new.NOBO)
NOBO_bad_models <- filter(new.NOBO, AIC<100)
NOBO_bad_models

new.BGGN<-read.csv('Umbrella_results/DistanceModels_0828/BGGN_models.csv') %>%
  mutate(run='Aug27') %>%
  left_join(aug27.models) %>%
  dplyr::select(-model.n)
BGGN_bad_models <- filter(new.BGGN, AIC<100)

new.GCFL<-read.csv('Umbrella_results/DistanceModels_0828/GCFL_models.csv') %>%
  mutate(run='Aug27') %>%
  left_join(aug27.models) %>%
  dplyr::select(-model.n)
GCFL_bad_models <- filter(new.GCFL, AIC<100)


#peak at variable correlation again
variables<-c('Temp', 'Clouds', 'Wind', 'Noise', 'nMinAfterMid', 'OBS', 'Replicate', 'Year')
var.correl<-cor(DCERP %>%
                  dplyr::select(all_of(variables), -OBS)%>%
                  filter(!is.na(Temp)) %>%
                  data.frame())
corrplot.mixed(var.correl)

bind_rows(NOBO_bad_models, NOCA_bad_models) %>%
  dplyr::select(formula, model.parm) %>%
  separate(model.parm, into=c(paste0('x', 1:8))) %>%
  pivot_longer(-formula) %>%
  filter(!is.na(value)) %>%
  count(value) %>%
  arrange(desc(n))
bind_rows(NOBO_bad_models, NOCA_bad_models, GCFL_bad_models, BGGN_bad_models) %>%
  dplyr::select(formula, AIC, model.parm) %>%
  rowwise() %>%
  mutate(likely.cor=case_when(grepl('Year', model.parm)&grepl('OBS', model.parm)~'year.obs',
                              grepl('Year', model.parm)&grepl('Wind', model.parm)~'year.wind',
                              grepl('Replicate', model.parm)&grepl('Temp', model.parm)~'rep.temp',
                              grepl('Temp', model.parm)&grepl('nMinAfterMid', model.parm)&grepl('Wind', model.parm)~'min.wind.temp',
                              grepl('Wind', model.parm)&grepl('Temp', model.parm)~'wind.temp',
                              grepl('MinAfter', model.parm)&grepl('Temp', model.parm)~'min.temp',
                              grepl('MinAfter', model.parm)&grepl('Wind', model.parm)~'min.wind')) %>%
  count(likely.cor)

pot_probs<-models_subset %>%
  rowwise() %>%
  mutate(likely.cor=case_when(grepl('Year', model.parm)&grepl('OBS', model.parm)~'year.obs',
                              grepl('Year', model.parm)&grepl('Wind', model.parm)~'year.wind',
                              grepl('Replicate', model.parm)&grepl('Temp', model.parm)~'rep.temp',
                              grepl('Temp', model.parm)&grepl('nMinAfterMid', model.parm)&grepl('Wind', model.parm)~'min.wind.temp',
                              grepl('Wind', model.parm)&grepl('Temp', model.parm)~'wind.temp',
                              grepl('MinAfter', model.parm)&grepl('Temp', model.parm)~'min.temp',
                              grepl('MinAfter', model.parm)&grepl('Wind', model.parm)~'min.wind')) %>%
  filter(!is.na(likely.cor))


DONT_run_these_models<-bind_rows(NOBO_bad_models, NOCA_bad_models, GCFL_bad_models, BGGN_bad_models) %>%
  count(model.parm) %>%
  filter(n>=2) %>%
  mutate(reason='variance/covariance problems')
DONT_run_these_models


bind_rows(NOBO_bad_models, NOCA_bad_models, GCFL_bad_models, BGGN_bad_models) %>%
  count(model.parm) %>%
  filter(n>=2,
         model.parm %in% pot_probs$model.parm)


#not rerunning models for NOCA & NOBO
new.NOCA %>% filter(!(formula %in% modesl_subset$formula)&
                      !(model.parm %in% DONT_run_these_models$model.parm))
new.NOBO %>% filter(!(formula %in% modesl_subset$formula)&
                      !(model.parm %in% DONT_run_these_models$model.parm))

# what models have I already run?
model.files<-list.files(path='Umbrella_results/DistanceModels/',
                        full.names=T)
all_model_aic_2<-NULL
for(b in model.files){
  ms<-read.csv(b)
  all_model_aic_2<-bind_rows(all_model_aic, ms)
}

#0825 models (accidentally deleted)
simp_models<-read.csv('Umbrella_results/all_models_tested_0825.csv')
all_models_0827<-read.csv('Umbrella_results/all_models_tested_0827.csv')
all_models_0830<-read.csv('Umbrella_results/all_models_tested_0830.csv')
mod.name.overlap<-which(all_models_0827$model.n %in% all_models_0830$model.n)
ove<-all_models_0827[mod.name.overlap,]$model.n
all_models_0830 %>%
  filter(model.n %in% ove)
test<-full_join(all_models_0827, all_models_0830, by='model.n', suffix=c('27', '30'))
which(is.na(test$formula27 != test$formula30) |
      test$formula27 != test$formula30) 
which(test$formula27==test$formula30)

RERUN_orig_models<-test[which(is.na(test$formula27 != test$formula30)),] %>%
  filter(!(model.parm27 %in% DONT_run_these_models$model.parm)) %>%
  dplyr::select(formula27, model.n, model.parm27) %>%
  rename(formula='formula27', model.parm='model.parm27')

#which sp don't have all the models
all_model_aic %>% count(species) %>% 
  arrange(desc(n))

#which models aren't represented by all sp
all_model_aic %>% count(model.n) %>% 
  arrange(desc(n)) %>%
  filter(!(model.n %in% DONT_run_these_models$model.n),
         n!=28) # 28 species is the max (aka don't need to rerun these models)

# make a long list of models that need to be rerun -
models_completed<- all_model_aic %>% 
  filter(species=='BGGN') %>% 
  pull(model.n)

ideal_model_list<-data.frame(species=rep(unique(all_model_aic$species), 
                                         each=length(models_completed)),
                             model.n=rep(models_completed, 28)) %>%
  bind_rows(data.frame(species=rep(unique(all_model_aic$species),
                                   each=32),
                       model.n=rep(RERUN_orig_models$model.n, 28))) %>%
  mutate(model.name=paste(species, model.n, sep='.'))

ideal_model_list %>% group_by(species) %>% tally() %>% arrange(desc(n))
ideal_model_list %>% group_by(model.n) %>% tally() %>% arrange(desc(n))

which(RERUN_orig_models$model.n %in% models_completed)

rerun_these_spmod<-ideal_model_list %>% 
  filter(!(model.name %in% all_model_aic_2$model.name)) %>%
  left_join(all_models_0827)
rerun_these_spmod[which(rerun_these_spmod$model.parm %in% all_models_0830$model.parm),] %>%
  pull(species) %>% unique() #the three species I didn't complete runs on
head(rerun_these_spmod)

for(w in 305:1727){
    assign(rerun_these_spmod$model.name[w], 
           ds(get(paste0(rerun_these_spmod$species[w], '_dat')),
              transect="point", key="hr",
              formula=as.formula(rerun_these_spmod$formula[w]), 
              adjustment = NULL, order = 0,
              truncation = '5%', 
              convert.units = conversion.factor,
              er.var='P3' #point encounter rate
           ))
    print(paste(rerun_these_spmod$model.name[w], 'row number', w, 'of', nrow(rerun_these_spmod)))
    model.summaries<-rerun_these_spmod %>%
      slice(305:w) %>%
      rowwise() %>%
      mutate(model.name=paste0(species, '.', model.n),
             df=as.numeric(AIC(get(model.name))[1]),
             AIC=as.numeric(round(AIC(get(model.name))[2],3)))
    write.csv(model.summaries, 
              file=paste0('Umbrella_results/rerun/fill_in_400c.csv'), 
              row.names=F)
}
# errors in models
#COYE  mod216 COYE.mod216 ~Temp+Clouds+Wind+Noise+Replicate+Year
#EABL  mod135 EABL.mod135 ~Clouds+Noise+Replicate+Year


#where ''n'' denotes the sample size and ''k'' denotes the number of parameters.


# LOAD MODEL INFO =====
# what models have I already run?
model.files<-list.files(path='Umbrella_results/DistanceModels/',
                        full.names=T)
all_model_aic_2<-NULL
for(b in model.files){
  ms<-read.csv(b)
  if(grepl('fill',b)){
    ms<-ms %>% dplyr::select(-model.parm, -formula)
    }
  all_model_aic_2<-bind_rows(all_model_aic_2, ms)
}
head(ideal_model_list)
ideal_model_list %>%
  filter(!(model.name %in% all_model_aic_2$model.name))

all_model_aic_fin<-all_model_aic_2 %>%
  filter(duplicated(.)) %>%
  mutate(AICc= AIC + (2k(k + 1)/(n - k - 1))) #convert AIC to AICc
#AICc = AIC + [2k(k + 1)/(n - k - 1)]
#where ''n'' denotes the sample size and ''k'' denotes the number of parameters.
  
inv<-summary(EATO.mod100)
library(MuMIn)
AICc(EATO.mod100)
AIC(EATO.mod100)
str(inv)
inv$ds

model.specs.spp<-all_model_aic_2 %>%
  filter(model.n=='mod101') %>%
  left_join(all_models_0830, by='model.n') %>%
  rowwise() %>%
  mutate(existsq=exists(model.name),
         trunction.dist=summary(get(model.name))$ds$width)
model.specs.spp

old.model.parm<-read.csv('Umbrella_results/model_parameter_summary20200826.csv')
parmfor.aicc<-old.model.parm %>% 
 dplyr::select(species, truncation.distance) %>%
  rowwise() %>%
  mutate(nobss=nrow(get(paste0(species, '_dat'))%>% 
                      filter(distance<=truncation.distance))) %>%
  rename(tdis='truncation.distance')
#equal to the new shit

testing.models<-test %>% group_by(species) %>% 
  filter(AIC>100) %>%arrange(AIC) %>% mutate(dA=AIC-min(AIC)) %>% filter(dA<=3) %>%
  pull(model.name)
all_model_aic_2 %>% as_tibble() %>%
  group_by(species) %>%
  left_join(parmfor.aicc) %>%
  mutate(AICc= AIC + ((2*df*(df + 1))/(nobss - df - 1))) %>%
  filter(AIC >100) %>%
  mutate(deltaAICc=AICc - min(AICc),
         dAIC=AIC-min(AIC)) %>%
  filter(deltaAICc < 2) %>%
  arrange(AICc) %>%
  left_join(all_models_0830 %>% bind_rows(RERUN_orig_models))

models<-all_models_0830 %>% bind_rows(RERUN_orig_models)
### TO DO
#organize code into logical chunks: load data, run models, assess & summarize models, make tables/plots
#creat truncation distance & nobs dataframe for each species
#use that data to recalculate AICc
#identify top models

#TOMORROW! create table of obs at different sites? specifically nonpine sites from 2011
