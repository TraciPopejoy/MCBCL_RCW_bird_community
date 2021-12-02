#Simple Occupancy models based on point counts
# TPD Mar 12, 20201
library(tidyverse); library(unmarked); library(MuMIn)
library(AICcmodavg); library(cowplot)

source('Umbrella/U_load_data.R')
head(DCERP_filt)

# Prepare common covariates (site, point, year) -----
occ_cov <- DCERP_filt %>% group_by(`Pt name`) %>% slice(1) %>%
  filter(!is.na(USFS_hab_index)) %>%
  summarize(meanUSFS=mean(USFS_hab_index)) %>%
  ungroup() %>%
  mutate(USFS_s=scale(meanUSFS)[,1])

scaled_phi_cov<-DCERP_filt %>% 
  group_by(`Pt name`, Replicate, Year) %>% slice(1) %>%
  filter(!is.na(USFS_hab_index)) %>%
  ungroup() %>%
  mutate(visit=paste0('v',Year, '.', Replicate),
         temp_s=scale(Temp)[,1], 
         cloud_s=scale(Clouds)[,1], 
         wind_s=scale(Wind)[,1], 
         noise_s=scale(Noise)[,1], 
         nMinAfterMid_s=scale(nMinAfterMid)[,1],
         Rep_s=scale(Replicate)[,1]) %>%
  select(`Pt name`, Replicate, Year, visit, 
         temp_s, cloud_s, wind_s, noise_s, nMinAfterMid_s, OBS, Rep_s) %>%
  arrange(visit)

phi_cov<-list(temp = scaled_phi_cov %>% select(`Pt name`, visit, temp_s) %>%
                pivot_wider(names_from=visit, values_from=temp_s) %>% select(-`Pt name`),
              cloud = scaled_phi_cov %>% select(`Pt name`, visit, cloud_s) %>%
                pivot_wider(names_from=visit, values_from=cloud_s) %>% select(-`Pt name`),
              wind = scaled_phi_cov %>% select(`Pt name`, visit, wind_s) %>%
                pivot_wider(names_from=visit, values_from=wind_s) %>% select(-`Pt name`),
              noise = scaled_phi_cov %>% select(`Pt name`, visit, noise_s) %>%
                pivot_wider(names_from=visit, values_from=noise_s) %>% select(-`Pt name`),
              minAmid = scaled_phi_cov %>% select(`Pt name`, visit, nMinAfterMid_s) %>%
                pivot_wider(names_from=visit, values_from=nMinAfterMid_s) %>% select(-`Pt name`),
              obs = scaled_phi_cov %>% select(`Pt name`, visit, OBS) %>%
                pivot_wider(names_from=visit, values_from=OBS) %>% select(-`Pt name`),
              rep = scaled_phi_cov %>% select(`Pt name`, visit, Rep_s) %>%
                pivot_wider(names_from=visit, values_from=Rep_s) %>% select(-`Pt name`))
species<-c('NOCA','NOBO','EATO','PIWA','EAWP','BHCO','MODO',
           'RHWO','PRAW','WEVI','NOFL','YTWA','RBWO','INBU',
           'SUTA','CACH','COYE','BLJA','RCWO','AMCR','BACS',
           'GCFL','BGGN','CARW','BHNU','CHSP','EABL')
#gcfl chsp has variance isses
#EATO (NONE) eawp rhwo praw wevi suta bacs messed up, run separately
# Set up for loop to run through and save each year x species combination ----
occ_gdf<-NULL
top_occ_mod<-NULL
occ_info<-NULL
chi.all<-NULL
species[21]
species[c(11:14,16:21,23:26,28)]
for(j in species){
  bird_det<-spfilter(j) %>% 
    select(Species, `Pt name`, Replicate, Year) %>%
    group_by(`Pt name`, Replicate, Year) %>% slice(1) %>%
    mutate(presence=ifelse(is.na(Species), 0, 1),
           visit=paste0('v',Year, '.', Replicate)) %>%
    group_by(`Pt name`) %>% 
    select(`Pt name`, visit, presence) %>%
    arrange(visit) %>%
    pivot_wider(names_from=visit, values_from=presence) %>%
    ungroup() %>% select(-`Pt name`) %>% as.matrix()
  
  for(y in c('2009','2010')){
    rn_df<-unmarkedFrameOccu(y=bird_det[,grep(y, colnames(bird_det))],
                             siteCovs = occ_cov,
                             obsCovs = lapply(phi_cov, function(x) x[grep(y, colnames(x))]))
    
    rn <- occu(~ temp+cloud+wind+noise+minAmid+obs+rep ~ USFS_s, rn_df)
    occ_d<-dredge(rn)
    tm<-get.models(occ_d, subset=T)
    chisq.out<-mb.gof.test(tm[[1]], plot.hist=F, nsim=50)
    v=1
    #print(paste(j,y,v, sep=" : "))
    #while(chisq.out$p.value < 0.1 ){
    #  chisq.out<-mb.gof.test(tm[[v]], nsim=50, plot.hist = F) 
    #  v=v+1
    #  print(paste(j,y,v, sep=" : "))
    #}
    ocm_sum<-summary(tm[[v]])
    chi1<-data.frame(pval=chisq.out$p.value,
                     chat=chisq.out$c.hat.est,
                     chi=chisq.out$chi.square,
                     spp=j,
                     year=y,
                     model_n=names(tm)[v],
                     mis.data=paste(".",grep("\\.",rownames(chisq.out$chisq.table))))
    chi.all<-bind_rows(chi.all, chi1)
    occ_gdf<-bind_rows(occ_gdf, 
                       ocm_sum$state %>% rownames_to_column() %>% 
                         mutate(spp=j,
                                year=y,
                                model_n=names(tm)[v]))
    saveRDS(tm, paste0('Umbrella/SimpleOcc/topmodels/',j,'_', y, '.rds'))
    
    det.ncoefs<-length(coef(tm[[v]], 'det'))-1
    det.est<-backTransform(linearComb(tm[[v]],
                                      c(1, rep(0, det.ncoefs)),
                                      "det"))
    occu.ncoefs<-length(coef(tm[[v]], 'state'))-1
    occu.est<-backTransform(linearComb(tm[[v]], c(1, rep(0, occu.ncoefs)), "state"))
    occ_info1<-data.frame(
      occu.est=occu.est@estimate,
      occu.025=confint(occu.est, level = 0.95)[1],
      occu.975=confint(occu.est, level = 0.95)[2],
      occu.det.est=det.est@estimate,
      occu.det.025=confint(det.est, level = 0.95)[1],
      occu.det.975=confint(det.est, level = 0.95)[2],
      year=y,
      spp=j,
      model_n=names(tm)[v],
      sites.occupied=sum(bup(ranef(tm[[v]]), stat="mode")))
    occ_info<-bind_rows(occ_info, occ_info1)
    
    if(is.na(occ_d$`psi(USFS_s)`[v])){
      score_tm<-get.models(occ_d[!is.na(occ_d$`psi(USFS_s)`),][1,], subset=T)
      ocm_sum<-summary(score_tm[[1]])
      
      occ_gdf<-bind_rows(occ_gdf, 
                         ocm_sum$state %>% rownames_to_column() %>% 
                           mutate(spp=j, year=y,
                                  model_n=names(score_tm)))
      saveRDS(score_tm,
              file=paste0('Umbrella/SimpleOcc/topmodels/',j,'_',y,'_rcwhab.rds'))
    }
    top_occ_mod<-bind_rows(top_occ_mod, 
                           occ_d %>% rownames_to_column() %>% mutate(spp=j, year=y))
  }
}

write_csv(top_occ_mod, 'Umbrella/SimpleOcc/top_occupancy_models.csv')
write_csv(occ_gdf, 'Umbrella/SimpleOcc/occupancy_hab_coef_est.csv')
write_csv(occ_info,'Umbrella/SimpleOcc/occupancy_estimates.csv')
write_csv(chi.all, 'Umbrella/SimpleOcc/occupancy_gof.csv')

chi.all %>%
  filter(pval<0.1) %>%
  left_join(top_occ_mod,
            by=c('model_n'='rowname',
                 'species'='spp',
                 'year')) %>% View()

#need to re-run these to identify model with lowest AIC and gof not significant
chi.all %>% group_by(species, year) %>%
  count() %>% arrange(desc(n))

top_occ_mod %>% filter(spp=='BHCO' & year == 2009)
#TUTI problem ----
#remove the one time it wasn't detected at all
bird_det_alt<-spfilter("TUTI") %>% 
  select(Species, `Pt name`, Replicate, Year) %>%
  group_by(`Pt name`, Replicate, Year) %>% slice(1) %>%
  filter(`Pt name` != 'EXP19') %>%
  mutate(presence=ifelse(is.na(Species), 0, 1),
         visit=paste0('v',Year, '.', Replicate)) %>%
  group_by(`Pt name`) %>% 
  select(`Pt name`, visit, presence) %>%
  arrange(visit) %>%
  pivot_wider(names_from=visit, values_from=presence) %>%
  ungroup() %>% select(-`Pt name`) %>% as.matrix()

occ_cov_alt <- DCERP_filt %>% group_by(`Pt name`) %>% 
  filter(`Pt name` != 'EXP19') %>% slice(1) %>%
  filter(!is.na(USFS_hab_index)) %>%
  summarize(meanUSFS=mean(USFS_hab_index)) %>%
  ungroup() %>%
  mutate(USFS_s=scale(meanUSFS)[,1])

rn_df_alt<-unmarkedFrameOccu(y=bird_det_alt[,grep('2009', colnames(bird_det_alt))],
                             siteCovs = occ_cov_alt,
                             obsCovs = lapply(list(temp = scaled_phi_cov %>% filter(`Pt name` != 'EXP19') %>% select(`Pt name`, visit, temp_s) %>%
                                                     pivot_wider(names_from=visit, values_from=temp_s) %>% select(-`Pt name`),
                                                   cloud = scaled_phi_cov %>% filter(`Pt name` != 'EXP19') %>%select(`Pt name`, visit, cloud_s) %>%
                                                     pivot_wider(names_from=visit, values_from=cloud_s) %>% select(-`Pt name`),
                                                   wind = scaled_phi_cov %>% filter(`Pt name` != 'EXP19') %>%select(`Pt name`, visit, wind_s) %>%
                                                     pivot_wider(names_from=visit, values_from=wind_s) %>% select(-`Pt name`),
                                                   noise = scaled_phi_cov %>% filter(`Pt name` != 'EXP19') %>%select(`Pt name`, visit, noise_s) %>%
                                                     pivot_wider(names_from=visit, values_from=noise_s) %>% select(-`Pt name`),
                                                   minAmid = scaled_phi_cov %>% filter(`Pt name` != 'EXP19') %>%select(`Pt name`, visit, nMinAfterMid_s) %>%
                                                     pivot_wider(names_from=visit, values_from=nMinAfterMid_s) %>% select(-`Pt name`),
                                                   obs = scaled_phi_cov %>% filter(`Pt name` != 'EXP19') %>%select(`Pt name`, visit, OBS) %>%
                                                     pivot_wider(names_from=visit, values_from=OBS) %>% select(-`Pt name`),
                                                   rep = scaled_phi_cov %>% filter(`Pt name` != 'EXP19') %>%select(`Pt name`, visit, Rep_s) %>%
                                                     pivot_wider(names_from=visit, values_from=Rep_s) %>% select(-`Pt name`)), 
                                              function(x) x[grep(y, colnames(x))]))

rn_alt <- occu(~ temp+cloud+wind+noise+minAmid+obs+rep ~ USFS_s, rn_df_alt)
occ_d_alt<-dredge(rn_alt)
tm_alt<-get.models(occ_d_alt, subset=T)

occ_d_alt %>% as_tibble() %>%
  rowid_to_column() %>%
  filter(!is.na(`psi(USFS_s)`)) %>% View()

tuti.alt.pred <- predict(tm_alt[[12]], type="state", 
                         newdata=usfs_grad,
                         appendData=TRUE)
ggplot(tuti.alt.pred)+
  geom_ribbon(aes(x=USFS_r, y=Predicted, 
                  ymin=lower, ymax=upper),
              alpha=0.5)+
  geom_line(aes(x=USFS_r, y=Predicted))+
  scale_y_continuous('Detection Probability',
                     limit = c(0,1))+
  scale_color_grey(guide=F, aesthetics=c('fill','color'))+
  scale_linetype(guide=F)+
  theme_cowplot()+
  theme(axis.title = element_text(size=11),
        axis.text=element_text(size=10))

saveRDS(tm_alt, 'Umbrella/SimpleOcc/topmodels/TUTI_2009.rds')
score_tm<-get.models(occ_d_alt[!is.na(occ_d_alt$`psi(USFS_s)`),][1,], subset=T)
saveRDS(score_tm, 'Umbrella/SimpleOcc/topmodels/TUTI_2009_rcwdet.rds')

#Quantify Odds Ratio ----
occ_coef_trans<-NULL
for(i in grep(list.files("./Umbrella/SimpleOcc/topmodels/", 
                         full.names = T),
              pattern = 'rcwhab',invert=T, value=T)){
  model<-readRDS(i)
  ocu_lengths<-sapply(model,function(x){length(coef(x, 'state'))})
  model_name<-names(ocu_lengths[ocu_lengths==2][1])
  top_rcw_mod<-model[model_name]
  new_occ_info<-confint(top_rcw_mod[[1]], type='state', level=0.95) %>%
    as_tibble() %>% 
    bind_cols(logit.estimate=coef(top_rcw_mod[[1]], 'state'),
              parameter=c('psi.Int','psi.USFS_s')) %>%
    rename(logit.025=`0.025`, logit.975=`0.975`) %>%
    mutate(across(starts_with('logit'), exp, .names="{.col}.OR")) %>%
    rename(OR.025=logit.025.OR, OR.975=logit.975.OR, OR.est=logit.estimate.OR) %>%
    mutate(spp=substr(i, 32, 35),
           year=as.numeric(substr(i, 37, 40)),
           mod_name=as.numeric(model_name))
  occ_coef_trans<-bind_rows(occ_coef_trans, new_occ_info)
}
write.csv(occ_coef_trans, 'Umbrella/SimpleOcc/occu_model_parameter_transformations.csv')

# Visualize Results ----
nd<-data.frame(USFS_s=seq(-2.6, 1.6, .1)) %>% mutate(RCW=(USFS_s*0.6950697)+3.907453)
all_pred_df<-NULL
for(i in list.files("./Umbrella/SimpleOcc/topmodels/", 
                    pattern = 'rds', full.names = T)){
  model<-readRDS(i)
  det_val<-ifelse(sum(grepl('USFS_s', names(coef(model[[1]], "state"))))==1, 
                  "habitat->occu", "no habitat\naffect")
  pred_df<-predict(model[[1]], type='state', newdata=nd) %>% 
    bind_cols(nd) %>%
    mutate(model=ifelse(substr(i, 48, 52)==".rds","top model","model w/ RCW score"),
           taxa=substr(i,32, 35),
           year=substr(i, 37, 40),
           rcw.det = det_val)
  all_pred_df<-bind_rows(all_pred_df, pred_df)
} 

all_pred_df  %>% group_by(model) %>% count() %>% mutate(n/86)

all_pred_df %>% left_join(bird.assem, by=c('taxa'='Species code')) %>%
  filter(Habitat.group %in% c('shrub','hardwood'),
         rcw.det=='habitat->occu') %>%
  ggplot()+
  geom_ribbon(aes(x=RCW, y=Predicted,
                  ymin=lower, ymax=upper,
                  fill=year), alpha=0.5)+
  geom_line(aes(x=RCW, y=Predicted, group=year))+
  scale_y_continuous("Occupancy Probability", 
                     breaks=seq(0,1,.33),
                     limits=c(0,1))+
  scale_x_continuous("RCW score")+
  facet_wrap(~taxa, nrow=2)

library(ggforce)
pdf('Umbrella/SimpleOcc/visualize_simple_occ.pdf')
for(k in 1:2){
  print(all_pred_df %>%
          filter(model=='top model') %>%
          ggplot()+
          geom_ribbon(aes(x=RCW, y=Predicted,
                          ymin=lower, ymax=upper,
                          fill=year), alpha=0.5)+
          geom_line(aes(x=RCW, y=Predicted, group=year))+
          scale_y_continuous("Occupancy Probability", 
                             breaks=seq(0,1,.33),
                             limits=c(0,1))+
          scale_x_continuous("RCW score")+
          facet_wrap_paginate(~taxa, ncol=3, nrow=5, page=k)+
          scale_fill_viridis_d()+
          theme_bw()+
          theme(legend.position='top'))
}
print(all_pred_df %>%
        filter(model=='model w/ RCW score') %>%
        ggplot()+
        geom_ribbon(aes(x=RCW, y=Predicted,
                        ymin=lower, ymax=upper,
                        fill=year), alpha=0.5)+
        geom_line(aes(x=RCW, y=Predicted, group=year))+
        scale_y_continuous("Occupancy Probability", 
                           breaks=seq(0,1,.33),
                           limits=c(0,1))+
        scale_x_continuous("RCW score")+
        facet_wrap(~taxa)+
        scale_fill_viridis_d()+
        theme_bw())
dev.off()


#### PICK UP HERE ####
occ_gdf<-read.csv('Umbrella/SimpleOcc/occupancy_hab_coef_est.csv')
top_occ_mod<-read.csv('Umbrella/SimpleOcc/top_occupancy_models.csv')
occ_coef_trans<-read.csv('Umbrella/SimpleOcc/occu_model_parameter_transformations.csv')
View(left_join(occ_gdf, top_occ_mod %>% rename(model_n=rowname)))
occ_gdf %>% distinct(spp, year, model_n) %>% count(spp, year)

# Fig 3 & Fig 4  -----
library(cowplot)
occ_cov %>% mutate(test=USFS_s*sd(meanUSFS)+mean(meanUSFS))
sd(occ_cov$meanUSFS); mean(occ_cov$meanUSFS)

usfs_grad <- data.frame(USFS_s=seq(-2.6, 1.6, length=50)) %>%
  mutate(USFS_r=USFS_s*0.6950697+3.907453)

#Figure 3 - RCWO occupancy and presence
library(unmarked)
abundance_summary <- read.csv('Umbrella/results/abundance summary 95 confidence.csv')
rcwo_top_occ_2009<-readRDS('Umbrella/SimpleOcc/topmodels/RCWO_2009.rds')
rcwo_top_occ_2010<-readRDS('Umbrella/SimpleOcc/topmodels/RCWO_2010.rds')
min(sapply(rcwo_top_occ_2009, AIC))==AIC(rcwo_top_occ_2009[[1]]) #checking

rcw.psi.09 <- predict(rcwo_top_occ_2009[[1]], type="state", 
                   newdata=usfs_grad, appendData=TRUE)
rcw.psi.10 <- predict(rcwo_top_occ_2010[[1]], type="state", 
                      newdata=usfs_grad, appendData=TRUE)
rcw.op<-bind_rows(rcw.psi.09 %>% mutate(year="2009"), 
                  rcw.psi.10 %>% mutate(year="2010")) %>%
  ggplot()+
  geom_ribbon(aes(x=USFS_r, y=Predicted, ymin=lower, ymax=upper, fill=year),
              alpha=0.5)+
  geom_line(aes(x=USFS_r, y=Predicted, linetype=year))+
  scale_x_continuous('RCW Habitat Score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76))+
  scale_y_continuous('Proportion of Sites Occupied',
                     limit = c(0,1))+
  scale_color_grey(guide=F, aesthetics=c('fill','color'))+
  scale_linetype(guide=F)+
  theme_cowplot()+
  theme(axis.title = element_text(size=11),
        axis.text=element_text(size=10))

rcwo.ap<-abundance_summary %>% 
  filter(species=='RCWO') %>%
  mutate(HabCat=as.numeric(HabCatfix)) %>%
  ggplot()+
  geom_linerange(aes(x=HabCat, ymin=x2.5, ymax=x97.5))+
  geom_point(aes(x=HabCat, y=x50), size=2)+
  scale_y_continuous(name=expression("RCW per km"^2))+
  scale_x_continuous('RCW Habitat Score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76),
                     expand=c(0.05,0.05))+
  theme_cowplot()+
  theme(axis.title = element_text(size=11),
        axis.text=element_text(size=10))
plot_grid(rcw.op, rcwo.ap, labels='AUTO')
ggsave('Umbrella/results/Figures/Fig3_rcwabun_20210323.tiff',
       width=5.5, height=3)


#Figure 4
bacs_top_occ_2009<-readRDS('Umbrella/SimpleOcc/topmodels/BACS_2009.rds')
bacs_top_occ_2010<-readRDS('Umbrella/SimpleOcc/topmodels/BACS_2010.rds')

bac.psi.09 <- predict(bacs_top_occ_2009[[1]], type="state", 
                      newdata=usfs_grad, appendData=TRUE)
bac.psi.10 <- predict(bacs_top_occ_2010[[1]], type="state", 
                      newdata=usfs_grad, appendData=TRUE)
bacs.op<-bind_rows(bac.psi.09 %>% mutate(year="2009"), 
                  bac.psi.10 %>% mutate(year="2010")) %>%
  ggplot()+
  geom_ribbon(aes(x=USFS_r, y=Predicted, ymin=lower, ymax=upper,
                  fill=year), alpha = 0.5)+
  geom_line(aes(x=USFS_r, y=Predicted, linetype=year))+
  scale_x_continuous('RCW Habitat Score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76))+
  scale_y_continuous('Proportion of Sites Occupied',
                     limit = c(0,1))+
  scale_color_grey(guide=F, aesthetics=c('fill','color'))+
  scale_linetype(guide=F)+
  theme_cowplot()+
  theme(axis.title = element_text(size=11),
        axis.text=element_text(size=10))

bacs.ap <- abundance_summary %>% 
  filter(species=='BACS') %>%
  mutate(HabCat=as.numeric(HabCatfix)) %>%
  ggplot()+
  geom_linerange(aes(x=HabCat, ymin=x2.5, ymax=x97.5))+
  geom_point(aes(x=HabCat, y=x50), size=2)+
  scale_y_continuous(name=expression("Bachman's Sparrow per km "^2))+
  scale_x_continuous('RCW Habitat Score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76),
                     expand=c(0.05,0.05))+
  theme_cowplot()+
  theme(axis.title = element_text(size=11),
        axis.text=element_text(size=10))
plot_grid(bacs.op, bacs.ap, labels='AUTO')
ggsave('Umbrella/results/Figures/Fig4_bacsabun_20210323.tiff',
       width=5.5, height=3)

# Table S2 -----
#occupancy detection model 
occ.det.models<-occ_coef_trans %>% 
  left_join(top_occ_mod, by=c('spp', 'year', 'mod_name'='rowname')) %>%
  filter(parameter=='psi.USFS_s') %>%
  #bind_rows(top_occ_mod %>% filter(!is.na(psi.USFS_s.)) %>%
  #            group_by(spp, year) %>% slice(1)) %>%
  filter(!duplicated(.)) %>% 
  select(spp, year, p.cloud.:df, AICc, delta) %>%
  arrange(spp) %>% #View()
  left_join(bird.assem, by=c('spp'='Species code')) %>%
  select(Habitat.group, `Common name`, everything()) %>% 
  select(-spp, -`Scientific name`) %>% 
  mutate(across(starts_with('p.'), .fns=as.character))%>%  
  pivot_longer(cols = starts_with('p.')) %>%
  group_by(`Common name`, year) %>% 
  filter(!is.na(value)) %>% #lose things here, make a new table
  select(-value, -delta) %>%
  pivot_wider(values_from = name, values_fill='') %>%
  unite(col='dmf', starts_with('p.'), sep='+') 
  
top_occ_mod %>% filter(delta==0) %>%
  bind_rows(top_occ_mod %>% filter(!is.na(psi.USFS_s.)) %>%
              group_by(spp, year) %>% slice(1)) %>%
  filter(!duplicated(.)) %>% 
  select(spp, year, p.cloud.:df, AICc, delta) %>%
  arrange(spp) %>% #View()
  left_join(bird.assem, by=c('spp'='Species code')) %>% 
  select(Habitat.group, `Common name`, year, psi.USFS_s., df, AICc, delta) %>% 
  group_by(`Common name`, year) %>%
  mutate(next.delta=lead(delta)) %>%
  arrange(Habitat.group, `Common name`) %>% 
  filter(delta==0) %>% select(-delta) %>%
  arrange(Habitat.group, `Common name`) %>%
  left_join(occ.det.models) %>%
  #arrange(next.delta)%>%View() #%>%
write_csv('Umbrella/results/Occupancy_mod_info_S2_20210323.csv')

# Figure 2 ----
DCERP_filt %>% 
  mutate(Region.Label=hab_bins$mid[findInterval(USFS_hab_index, hab_bins$min,
                                                left.open = F)]) %>%
  filter(!(Species %in% c("LAGU", "TUVU"))) %>%
  group_by(`Pt name`, Year, USFS_hab_index, Region.Label) %>%
  summarize(ntaxa=length(unique(Species)),
            taxa=paste(unique(Species), collapse=', ')) %>%
  filter(!is.na(Region.Label)) %>%
  ggplot(aes(y=ntaxa))+
  stat_summary(aes(x=as.character(Region.Label)))+
  scale_x_discrete(name='RCW Habitat Score')+
  scale_y_continuous(name='Species Richness')+
  theme_cowplot()
ggsave('Umbrella/results/Figures/Fig2_site_richness_taxa.tiff', width=3, height=3)

DCERP_filt %>% pull(Species) %>% unique() 

occ_info<-read_csv('Umbrella/SimpleOcc/occupancy_estimates.csv')
occ_info %>%
  left_join(bird.assem, by=c('spp'='Species code')) %>%
  mutate(`Occupancy Probability`=paste0(round(occu.est, 2)," (",round(occu.025,2),"-", round(occu.975,2),")"),
         `Detection Probability`=paste0(round(occu.det.est, 2)," (",round(occu.det.025,2),"-", round(occu.det.975,2),")")) %>%
  select(Habitat.group, spp, year, `Occupancy Probability`, `Detection Probability`, sites.occupied) %>%
  arrange(Habitat.group, spp) %>%
  write.csv('Umbrella/results/TS2_occupancy_estimates.csv')

top_occ_mod %>% filter(delta==0) %>%
  bind_rows(top_occ_mod %>% filter(!is.na(psi.USFS_s.)) %>%
              group_by(spp, year) %>% slice(1)) %>%
  distinct(.keep_all=T)  %>% 
  group_by(spp, year) %>% mutate(n=n()) %>%
  arrange(desc(n)) %>%
  mutate(n-1) %>%
  filter(n==2) %>%
  select(spp, year, n, `n - 1`, AICc, delta) %>% 
  filter(delta!=0) %>%
  left_join(bird.assem, by=c('spp'='Species code')) %>%
  arrange(desc(delta)) %>% group_by(spp) %>% slice(1) %>% ungroup() %>%
  count(Habitat.group)

# Figure S1 ----
library(grid)
bp<-top_occ_mod %>% 
  #fixing TUTI
  filter(species != 'TUTI' & year != 2009) %>% 
  bind_rows(data.frame(spp='TUTI', year=2009, psi.USFS_s.=NA, delta=0, AICc=808.28)) %>%
  filter(delta==0) %>% distinct(.) %>%
  rename(model_n=rowname) %>%
  dplyr::select(spp, year, starts_with('psi'), AICc, delta, model_n) %>%
  mutate(RCW_yn=factor(case_when(is.na(psi.USFS_s.)~'neutral',
                                 psi.USFS_s.<0~'negative',
                                 psi.USFS_s.>0~'positive'),
                       levels=c('positive','neutral','negative'))) %>%
  left_join(bird.assem, by=c('spp'='Species code')) %>%
  mutate(HabF=factor(Habitat.group, levels=c('longleaf', 'generalist','hardwood', 'shrub')))%>%
  ggplot(aes(y=spp, x=as.factor(year)))+
  geom_tile(aes(fill=RCW_yn))+
  facet_wrap(~HabF, scales='free',ncol = 1,
             strip.position = 'left')+
  scale_fill_manual('RCW\neffect', values=c('blue','grey','red'))+
  scale_x_discrete('Year')+
  scale_y_discrete('Species')+
  theme_bw()+
  theme(strip.placement='outside',
        strip.background = element_blank())+
  ggtitle('RCW habitat in the model')
bpp = ggplot_gtable(ggplot_build(bp))
#gtable::gtable_show_layout(ohp)
# get the number of unique x-axis values per facet (1 & 3, in this case)
y.var <- sapply(ggplot_build(bp)$layout$panel_scales_y,
                function(l) length(l$range$range))
# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
bpp$heights[bpp$layout$t[grepl("panel", bpp$layout$name)]] <- bpp$heights[bpp$layout$t[grepl("panel", bpp$layout$name)]] * y.var
grid.draw(bpp)

#assess coefficient strength
new_mod_info<-occ_coef_trans %>% 
  left_join(top_occ_mod, by=c('spp', 'year', 'mod_name'='rowname'))

ap<-new_mod_info %>% filter(parameter=='psi.USFS_s') %>%
  select(-starts_with('p.'),-starts_with('OR.'), -weight, -delta, -AICc, 
         -logLik, -df, -psi.Int., -psi.USFS_s.) %>%
  mutate(signif=case_when(logit.025 > 0 & logit.975 > 0 ~ 'positive',
                          logit.025 < 0 & logit.975 > 0 ~ 'neutral',
                          logit.025 < 0 & logit.975 < 9 ~ 'negative'),
         sigF=factor(signif, levels=c('positive','neutral','negative'))) %>%
  left_join(bird.assem, by=c('spp'='Species code'))%>%
  mutate(HabF=factor(Habitat.group, levels=c('longleaf', 'generalist','hardwood', 'shrub')))%>%
  ggplot(aes(y=spp, x=as.factor(year)))+
  geom_tile(aes(fill=signif))+
  facet_wrap(~HabF, scales='free',ncol = 1,
             strip.position = 'left')+
  scale_fill_manual('RCW effect', values=c('red','grey','blue'), guide=F)+
  scale_x_discrete('Year')+
  scale_y_discrete('Species')+
  theme_bw()+
  theme(strip.placement='outside',
        strip.background = element_blank())+
  ggtitle(expression('RCW habitat '*beta*' strength'))
aapp = ggplot_gtable(ggplot_build(ap))
y.var <- sapply(ggplot_build(ap)$layout$panel_scales_y,
                function(l) length(l$range$range))
# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
aapp$heights[aapp$layout$t[grepl("panel", aapp$layout$name)]] <- aapp$heights[aapp$layout$t[grepl("panel", aapp$layout$name)]] * y.var
grid.draw(aapp)

ggsave('Umbrella/results/Figures/FigS1a_occmod_RCWpres.svg', grid.draw(bpp),
       width=3, height=6.4)
ggsave('Umbrella/results/Figures/FigS1b_occmod_RCWpval.svg', grid.draw(aapp),
       width=3, height=6.4)
occ_gdf %>%
  filter(rowname=='USFS_s',
         spp=='RCWO' | spp=='BACS') %>%
  mutate(ul=exp(Estimate+SE),
         ll=exp(Estimate-SE),
         pe=exp(Estimate),
         ym=paste(year, model))