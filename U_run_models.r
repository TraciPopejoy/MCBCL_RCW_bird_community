# Run the distance sampling models ####
#library(remotes)
#install_github("DistanceDevelopment/mrds")
#library(devtools)
#install_github("DistanceDevelopment/Distance")
library(tidyverse); library(readxl); library(Distance); library(lubridate)

#building all the data sets
source('/home/tracidubose/BIRDS_RCW/U_load_data.r')
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
#deciding which key function to use (depends on the general shape of the detection data's PDF)
all_birds<-DCERP_filt %>%
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)],
         Sample.Label=paste(`Pt name`, Replicate, sep="."), #this identifies unique samples within strata
         size=`Cluster size`,
         Effort=1,
         Study.Area="MCBCL") %>%
  filter(USFS_hab_index !=0) %>%
  as.data.frame()

conversion.factor <- Distance::convert_units(distance_units = "Metre", 
                                             effort_units=NULL, 
                                             area_units = "square kilometre")
all_b_hr<-ds(all_birds,
   transect="point", key="hr", order = 0,
   truncation = '5%', convert.units = conversion.factor, er.var='P3')
all_b_hn<-ds(all_birds,
             transect="point", key="hn", order = 0,
             truncation = '5%', convert.units = conversion.factor, er.var='P3')
summarize_ds_models(all_b_hr, all_b_hn, output = 'plain')

# building all possible models
variables<-c('Temp', 'Clouds', 'Wind', 'Noise', 'nMinAfterMid', 'OBS', 'Replicate', 'Year',
             'HabCat')
# temperature, cloud cover, wind, noise, start time, observer, replicate
var.combos<-data.frame(data.frame(rbind(combn(variables, 2),NA, NA, NA, NA, NA, NA, NA),
                                  rbind(combn(variables, 3),NA, NA, NA, NA, NA, NA),
                                  rbind(combn(variables, 4),NA, NA, NA, NA, NA),
                                  rbind(combn(variables, 5),NA, NA, NA, NA),
                                  rbind(combn(variables, 6),NA, NA, NA),
                                  rbind(combn(variables, 7),NA, NA),
                                  rbind(combn(variables, 8),NA),
                                  rbind(combn(variables, 9))) %>%
             t())
models_pre<-unite(var.combos, 'formula', sep='+', na.rm=T) %>%
  mutate(formula=paste0('~', formula)) %>%
  bind_rows(data.frame(formula=paste('~', c(variables,'1')))) %>%
  mutate(model.n=paste0('mod', row_number()),
         model.parm=gsub('+', '.', gsub('~','', formula),
                                       fixed=T))
head(models_pre) #all the models possible
write.csv(models_pre, '/home/tracidubose/BIRDS_RCW/all_possible_models.csv')
models_already_completed<-read.csv('/home/tracidubose/BIRDS_RCW/all_possible_models_0929.csv')
models_pre1<-models_pre  %>%
  filter(!(formula %in% models_already_completed$formula))
nrow(models_pre1); nrow(models_pre)
models_pre_all<-models_pre
models_pre<-models_pre1
nrow(models_pre); nrow(models_pre_all)

# removing models that commonly break ####
library(corrplot)
var.correl<-cor(DCERP %>%
                  dplyr::select(all_of(variables), -OBS, -HabCat, HabCat_nf)%>%
                  filter(!is.na(Temp)) %>%
                  data.frame())
corrplot.mixed(var.correl)

# actually run the models
conversion.factor <- Distance::convert_units(distance_units = "Metre", 
                                   effort_units=NULL, 
                                   area_units = "square kilometre")
done.spp<-substr(list.files('/home/tracidubose/BIRDS_RCW/DistanceModels/',
                            pattern='hab_cat'),9,12)
sp.to.still.run<-species[!(species %in% done.spp)]
#sp.to.still.run<-sp.to.still.run[sp.to.still.run!="COYE"] 
#model ~103 broke for COYE, run without that after checking
#model 108 broke for CARW, run without

for(u in sp.to.still.run[-1]){
  model.summaries<-NULL
  for(w in 1:length(models_pre$model.n)){
    assign(paste(u, models_pre$model.n[w], sep='.'), 
           ds(get(paste0(u, '_dat')),
              transect="point", key="hr",
              formula=as.formula(models_pre$formula[w]), 
              adjustment = NULL, order = 0,
              truncation = '5%', 
              convert.units = conversion.factor,
              er.var='P3' #point encounter rate
              ))
    print(paste(u,'species number', which(species==u), 
                'model number', w))
    #some models break because of NA variance-covariance matrix elements 
    if(is.null(get(paste(u, models_pre$model.n[w], sep='.'))$dht)){
      mod_out<-data.frame(species=u,
                          model.n=models_pre$model.n[w],
                          model.state='broken')
    }else{
      model.sum<-summary(get(paste(u, models_pre$model.n[w], sep='.')))
      
      model.aic<-AIC(get(paste(u, models_pre$model.n[w], sep='.')))
      mod.gof<-gof_ds(get(paste(u, models_pre$model.n[w], sep='.')),
                      plot=F)
      mod_out<-data.frame(species=u,
                          model.n=models_pre$model.n[w],
                          model.state='fine') %>%
        rowwise() %>%
        mutate(model.name=paste0(species, '.', model.n),
             df=as.numeric(model.aic[1]),
             AIC=as.numeric(round(model.aic[2],4)),
             trun.dist=model.sum$ds$width,
             nobss=nrow(model.sum$ddf$data),
             AICc= AIC + ((2*df*(df + 1))/(nobss - df - 1)),
             gof_cvm_w=mod.gof$dsgof$CvM$W,
             gof_cvm_p=mod.gof$dsgof$CvM$p)
    }
    model.summaries<-bind_rows(model.summaries, mod_out)
  }
  write.csv(model.summaries, 
            file=paste0('/home/tracidubose/BIRDS_RCW/DistanceModels/hab_cat_', 
                        u,'_', format(Sys.Date(), '%m%d'),'_models.csv'), 
            row.names=F)
}

#investigate for abnormally low aic values - something would be wrong in these models
#Error message: 
#'Some variance-covariance matrix elements were NA, possible numerical problems; 
#'#only estimating detection function.'


