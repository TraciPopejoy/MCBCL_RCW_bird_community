#library(remotes)
#install_github("DistanceDevelopment/mrds")
#library(devtools)
#install_github("DistanceDevelopment/Distance")
#remotes::install_version('Distance', version='1.0.1')
library(tidyverse); library(Distance)
source('/home/tracidubose/BIRDS_RCW/U_load_data.r')
sp<-c('NOCA','NOBO','EATO','PIWA','EAWP','BHCO','MODO',
           'RHWO','PRAW','WEVI','NOFL','YTWA','RBWO','INBU',
           'SUTA','CACH','COYE','BLJA','RCWO','AMCR','BACS',
           'GCFL','BGGN','TUTI','CARW','BHNU','CHSP','EABL')

test<-NULL
for(u in c(list.files('/home/tracidubose/BIRDS_RCW/bootstrap_results_scaled/', 
                      pattern='model', full.names = T),
           list.files('/home/tracidubose/BIRDS_RCW/bootstrap_results/', 
                      pattern='model', full.names = T))){
  a<-read.csv(u) %>%
    mutate(species=ifelse(substr(u, 51, 53)=='ed/', 
                          substr(u, 55,58), substr(u, 48, 51)),
           type=ifelse(grepl('scaled', u), 'scaled', 'uns'),
           batch=ifelse(type=='scaled', substr(u, 68, 68), substr(u, 59, 59)))
  test<-bind_rows(test, a)
}
str(test)

n_run<-test %>% as_tibble() %>% 
  filter(phi!=10) %>%
  count(species) %>% 
  filter(n < 500) %>%
  mutate(k=510-n) %>%
  arrange(k)
n_run

for(j in n_run$species){
  assign(paste0(j, '_dat'),
         spfilter(j))
}

model.appendix.all<-read.csv('/home/tracidubose/BIRDS_RCW/models_20201016.csv')
boot_sum<-NULL
for(i in n_run$species[4]){
  spp<-i
  dat<-get(paste0(spp, '_dat')) %>% 
    rename(old.label='Sample.Label') #ds needs Sample.Label to be unique to each sample
  print(spp)
  bootstrap_res<-NULL
  all_abun<-NULL
  for(k in 1:as.numeric(n_run[n_run$species==i, 3])){
    # identify the subsample for the bootstrap
    ss<-dat %>%
      group_by(old.label) %>%
      slice(1) %>%
      group_by(`Pt name`) %>%
      sample_n(2, replace=T) %>%
      pull(old.label)
    
    sp_data <- NULL
    for(h in 1:length(ss)){
      sub <- dat %>% 
        filter(old.label == ss[h]) %>% 
        mutate(Sample.Label=paste(old.label, h, sep='.'))
      sp_data <- bind_rows(sp_data, sub)
    }
    
    tm.formula <- model.appendix.all %>%
      filter(species==spp) %>%
      filter(AICc == min(AICc))
    parm<-sub('\\ ','', strsplit(as.character(tm.formula$model.parm), '\\.')[[1]])
    parm<-parm[!(parm %in% c('OBS', 'HabCat'))]
    
    boot.dsm <- ds(sp_data %>% mutate_at(parm, scale),
                   transect="point", key="hr",
                   formula=as.formula(as.character(tm.formula$formula)), 
                   adjustment = NULL, order = 0,
                   truncation = '5%')
    empty_pts<-NULL
    missing_sites<-ss[!(ss %in% boot.dsm$ddf$data$old.label)]
    for(y in 1:length(missing_sites)){
      pts<-point_sample_empty %>% 
        mutate(old.label = paste(`Pt name`, Replicate, Year, sep ='.'),
               phi=0, 
               tob=0) %>%
        filter(old.label == missing_sites[y]) %>%
        mutate(Sample.Label=paste(old.label, y, sep='.')) %>%
        dplyr::select(Species, Sample.Label, old.label, Year, HabCat, phi, tob, USFS_hab_index)
      empty_pts<-bind_rows(empty_pts, pts)  
    }
    samp.area<-boot.dsm$ddf$meta.data$width*boot.dsm$ddf$meta.data$width*pi
    
    peace<-bind_cols(boot.dsm$ddf$data, 
                     predict.phi=predict(boot.dsm, esw=FALSE)$fitted) %>%
      group_by(Species, Sample.Label, old.label, Year, HabCat) %>% 
      dplyr::summarize(phi=mean(predict.phi),
                       tob=sum(size),
                       USFS_hab_index=mean(USFS_hab_index)) %>%
      ungroup() %>%
      filter(!duplicated(.)) %>%
      bind_rows(empty_pts) 
    
    amod<-glm(tob~USFS_hab_index+offset(phi), data=peace, family='poisson')
    rsss<-c(summary(amod)$coefficients[1,],
            summary(amod)$coefficients[2,],
            ifelse(is.null(boot.dsm$dht$clusters$average.p), 
                   10, boot.dsm$dht$clusters$average.p),
            k, nrow(peace))
    names(rsss)<-c(paste('Int', c('Est','SE','z','pvalue'),sep='.'),
                   paste('Hab', c('Est','SE','z','pvalue'),sep='.'),
                   'phi', 'runN', 'n_obs')
    bootstrap_res<-bind_rows(bootstrap_res, rsss)
    rabun<-amod$data %>% mutate(run=k,
                                nchat=tob/phi, 
                                nhat=nchat/(samp.area*592),
                                dens=(1/samp.area)*nchat,
                                NEWnhat=nchat/samp.area,
                                raw_dens_no_phi=tob/samp.area)
    all_abun<-bind_rows(all_abun, rabun)
  }
  write_csv(bootstrap_res, paste0('/home/tracidubose/BIRDS_RCW/bootstrap_results_scaled/',
                                  spp,'_modelcoef_',
                                  format(Sys.Date(), '%m%d'),'_2.csv'),
            append=T)
  write_csv(all_abun, paste0('/home/tracidubose/BIRDS_RCW/bootstrap_results_scaled/',
                             spp,'_abund',
                             format(Sys.Date(), '%m%d'),'_2.csv'),
            append=T)
}

summary_bootstrap<- bootstrap_res %>%
  filter(phi<2) %>%
  select(-runN) %>%
  summarize_all(quantile, prob=c(.5,.025,.975)) %>%
  mutate(value_type=c('median','ll','ul'),
         species=spp,
         nrun=nrow(bootstrap_res %>% filter(phi<2))) %>%
  select(species, value_type, everything())
boot_sum<-bind_rows(boot_sum, summary_bootstrap)

write.csv(boot_sum, paste0('/home/tracidubose/BIRDS_RCW/bootstrap_results/TOTAL_',
                           format(Sys.Date(), '%Y%m%d'), '.csv'))

test<-NULL
for(u in list.files('/home/tracidubose/BIRDS_RCW/bootstrap_results/', 
                    pattern='model', full.names = T)){
  ulll<-read.csv(u) %>%
    mutate(species=substr(u, 48,51))
  test<-bind_rows(test, ulll)
}
taxa_todo3<-test %>% filter(phi !=10) %>% count(species) %>%
  filter(n<100) %>% pull(species)
