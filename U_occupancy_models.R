#Occupancy models based on point counts
# TPD Nov 11, 2020
library(tidyverse); library(unmarked); library(MuMIn)

source('Umbrella/U_load_data.R')
head(DCERP_filt)

# Prepare common covariates (site, point, year) -----
site_cov <- DCERP_filt %>% group_by(`Pt name`) %>% slice(1) %>%
  filter(!is.na(USFS_hab_index))%>%
  select(`Pt name`, Base_Area)

occ_cov <- DCERP_filt %>% group_by(`Pt name`) %>% slice(1) %>%
  filter(!is.na(USFS_hab_index)) %>%
  summarize(meanUSFS=mean(USFS_hab_index)) %>%
  ungroup() %>%
  mutate(USFS_s=scale(meanUSFS)[,1])

year_cov<-DCERP_filt %>% group_by(`Pt name`, Year) %>% slice(1) %>%
  filter(!is.na(USFS_hab_index))%>%
  ungroup() %>%
  mutate(USFS_s=scale(USFS_hab_index)[,1],
         year=paste0('v',Year),
         Year_c=as.character(Year))%>%
  select(`Pt name`, year, Year_c, USFS_s)%>%
  pivot_wider(names_from = year, values_from = c(Year_c, USFS_s))

year_covl<-list(year=year_cov[,2:3] %>% as.data.frame(),
                usfs=year_cov[,4:5] %>% as.data.frame())

scaled_phi_cov<-DCERP_filt %>% group_by(`Pt name`, Replicate, Year) %>% slice(1) %>%
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

# Run Occupancy Model for each Spp -----
top_occ_mod<-NULL
occ_gdf<-NULL

bird_tax <- DCERP %>% group_by(Species) %>%
  count() %>% 
  filter(n>200) %>% #,
         #!(Species %in% substr(list.files('Umbrella/OccupancyModel'),1,4))) %>%
  pull(Species)

for(u in bird_tax[-c(1:3)]){
  #build the detection information
  bird_det<-spfilter(u) %>% 
    select(Species, `Pt name`, Replicate, Year) %>%
    group_by(`Pt name`, Replicate, Year) %>% slice(1) %>%
    mutate(presence=ifelse(is.na(Species), 0, 1),
           visit=paste0('v',Year, '.', Replicate)) %>%
    group_by(`Pt name`) %>% select(`Pt name`, visit, presence) %>%
    arrange(visit) %>%
    pivot_wider(names_from=visit, values_from=presence) %>%
    ungroup() %>% select(-`Pt name`) %>% as.matrix()
  
  #create the unmarked MultFrame for the function
  umf <- unmarkedMultFrame(y=bird_det, #detection record
                           siteCovs=occ_cov[,3], #site covariates
                           yearlySiteCovs=year_covl,
                           obsCovs=phi_cov,
                           numPrimary=2)
  
  #build the giant model to dredge
  bigmod<-colext(psiformula = ~ USFS_s, #initial probability of occupancy at each site
                 gammaformula = ~1, #colonization probability
                 epsilonformula = ~1, #extinction probability
                 pformula = ~temp+cloud+wind+noise+minAmid+obs+rep+USFS_s, #detection probability 
                 umf)

  occ_d<-dredge(bigmod)
  tm<-get.models(occ_d, subset = delta <= 0)
  ocm_sum<-summary(tm[[1]])
  occ_gdf<-bind_rows(occ_gdf, 
                     ocm_sum$psi %>% rownames_to_column() %>% 
                       mutate(species=u,
                              model_n=names(tm)))
  saveRDS(tm,
          file=paste0('Umbrella/OccupancyModel/',u,'_occmod.rds'))
  if(is.na(occ_d$`psi(USFS_s)`[1])){
    score_tm<-get.models(occ_d[!is.na(occ_d$`psi(USFS_s)`),][1,], subset=T)
    ocm_sum<-summary(score_tm[[1]])
    occ_gdf<-bind_rows(occ_gdf, 
                       ocm_sum$psi %>% rownames_to_column() %>% 
                         mutate(species=u,
                                model_n=names(score_tm)))
    saveRDS(score_tm,
            file=paste0('Umbrella/OccupancyModel/',u,'_occmod_rcwhab.rds'))
  }
  top_occ_mod<-bind_rows(top_occ_mod, 
                         occ_d %>% rownames_to_column() %>% 
                           filter(delta <= 3) %>% mutate(spp=u))
}
write_csv(top_occ_mod, 'Umbrella/OccupancyModel/top_occupancy_models.csv')
write_csv(occ_gdf, 'Umbrella/OccupancyModel/occupancy_hab_coef.csv')

list.files('Umbrella/OccupancyModel/')
top_occ_mod %>% group_by(spp) %>% slice(1)

# PICK UP HERE
occ_gdf<-read.csv('Umbrella/OccupancyModel/occupancy_hab_coef.csv')
top_occ_mod<-read.csv('Umbrella/OccupancyModel/top_occupancy_models.csv')

# Fig 3 & Fig 4  -----
occ_cov %>% mutate(test=USFS_s*sd(meanUSFS)+mean(meanUSFS))
sd(occ_cov$meanUSFS); mean(occ_cov$meanUSFS)

usfs_grad <- data.frame(USFS_s=seq(-2.6, 1.6, length=50)) %>%
  mutate(USFS_r=USFS_s*0.6950697+3.907453)

#Figure 3 - RCWO occupancy and presence
abundance_summary <- read.csv('Umbrella/results/abundance summary 95 confidence.csv')
rcwo_top_occ<-readRDS('Umbrella/OccupancyModel/RCWO_occmod.rds')
rcw.psi <- predict(rcwo_top_occ$`242`, type="psi", 
                    newdata=usfs_grad, appendData=TRUE)
rcw.op<-ggplot(data=rcw.psi)+
  geom_ribbon(aes(x=USFS_r, y=Predicted, ymin=lower, ymax=upper),
              fill='lightgrey')+
  geom_line(aes(x=USFS_r, y=Predicted))+
  scale_x_continuous('RCW Habitat Score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76))+
  scale_y_continuous('Proportion of Sites Occupied',
                     limit = c(0,1))+
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
ggsave('Umbrella/results/Figures/Fig3_rcwabun.tiff',
       width=5.5, height=3)


#Figure 4
bac_top_occ<-readRDS('Umbrella/OccupancyModel/BACS_occmod.rds')
bacs.psi <- predict(bac_top_occ$`138`, type="psi", 
                    newdata=usfs_grad, appendData=TRUE)
bacs.op<-ggplot(data=bacs.psi)+
  geom_ribbon(aes(x=USFS_r, y=Predicted, ymin=lower, ymax=upper),
              fill='lightgrey')+
  geom_line(aes(x=USFS_r, y=Predicted))+
  scale_x_continuous('RCW Habitat Score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76))+
  scale_y_continuous('Proportion of Sites Occupied',
                     limit = c(0,1))+
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
ggsave('Umbrella/results/Figures/Fig4_bacsabun.tiff',
       width=5.5, height=3)

# Table S2 -----
#occupancy detection model ---
top_occ_mod %>% filter(delta==0) %>%
  bind_rows(top_occ_mod %>% filter(!is.na(psi.USFS_s.)) %>%
              group_by(spp) %>% slice(1)) %>%
  filter(!duplicated(.)) %>%
  select(spp, p.cloud.:df, AICc, delta) %>%
  arrange(spp) %>% 
  left_join(bird.assem, by=c('spp'='Species code')) %>%
  select(Habitat.group, `Common name`, everything()) %>%
  select(-spp, -`Scientific name`) %>% #View()
  mutate(RCW_hab_detection=round(p.USFS_s.,2),
         across(starts_with('p.'), .fns=as.character)) %>%  
  pivot_longer(cols = starts_with('p.')) %>%
  filter(!is.na(value)) %>% 
  group_by(`Common name`) %>%
  select(-value) %>%
  pivot_wider(values_from = name, values_fill='') %>%
  unite(col='dmf', starts_with('p.'), sep='+') %>%
  mutate(next.delta=lead(delta),
         dmf=sub('.p', '', dmf)) %>%
  arrange(Habitat.group, `Common name`) %>%
  filter(!duplicated(`Common name`)) %>% select(-delta) %>%
  arrange(next.delta) %>% View() #%>%
  #write_csv('Umbrella/results/Occupancy_mod_info_S2.csv')

# subset for Table 4 -----
t4_occ<-top_occ_mod %>%  
  group_by(spp) %>% filter(!is.na(psi.USFS_s.)) %>% arrange(delta) %>% slice(1) %>%
  left_join(occ_gdf, by=c('spp'='species')) %>%
  filter(rowname.y=='USFS_s') %>%
  select(spp, psi.USFS_s., SE, delta) %>%
  arrange(desc(delta)) %>%
  mutate(`Occupancy Estimate`=paste0(round(psi.USFS_s.,2), ', ', 
                                     round(SE,2), 
                                     ifelse(delta==0, '*', ''))) %>%
  arrange(delta) %>%
  select(spp, `Occupancy Estimate`)

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


# checking match older occupancy models ----
#looking at figures to check
noca_top_occ<-readRDS('Umbrella/OccupancyModel/NOCA_occmod.rds')
noca.psi <- predict(noca_top_occ$`242`, type="psi", 
                    newdata=usfs_grad, appendData=TRUE)
ggplot(data=noca.psi)+
  geom_ribbon(aes(x=USFS_r, y=Predicted, ymin=lower, ymax=upper),
              fill='lightgrey')+
  geom_line(aes(x=USFS_r, y=Predicted))+
  scale_x_continuous('RCW Habitat Score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76))+
  scale_y_continuous('Proportion of Sites Occupied',
                     limit = c(.8,1.1))

eabl_top_occ<-readRDS('Umbrella/OccupancyModel/EABL_occmod.rds')
eabl.psi <- predict(eabl_top_occ$`242`, type="psi", 
                    newdata=usfs_grad, appendData=TRUE)
ggplot(data=eabl.psi)+
  geom_ribbon(aes(x=USFS_r, y=Predicted, ymin=lower, ymax=upper),
              fill='lightgrey')+
  geom_line(aes(x=USFS_r, y=Predicted))+
  scale_x_continuous('RCW Habitat Score',
                     breaks=c(2.5, 3.26, 3.76, 4.26, 4.76))+
  scale_y_continuous('Proportion of Sites Occupied',
                     limit = c(0,1))



nrow(scaled_phi_cov)
cm <- DCERP_filt %>% group_by(`Pt name`, Replicate, Year) %>% slice(1) %>%
  filter(!is.na(USFS_hab_index)) %>%
  ungroup()%>%
  select(Temp, Clouds, Wind, Noise, nMinAfterMid, USFS_hab_index, Replicate) %>%
  colMeans()
csd <- DCERP_filt %>% group_by(`Pt name`, Replicate, Year) %>% slice(1) %>%
  filter(!is.na(USFS_hab_index)) %>%
  ungroup()%>%
  select(Temp, Clouds, Wind, Noise, nMinAfterMid, USFS_hab_index, Replicate) %>%
  apply(2,sd)



rescale.coefs <- function(beta,mu,sigma) {
  beta2 <- beta ## inherit names etc.
  beta2[-1] <- sigma[1]*beta[-1]/sigma[-1]
  beta2[1] <- sigma[1]*beta[1]+mu[1]-sum(beta2[-1]*mu[-1])
  beta2
}

cm
csd
coef(bac_top_occ$`242`)



colext(psiformula = ~ USFS_s, #initial probability of occupancy at each site
       gammaformula = ~1, #colonization probability
       epsilonformula = ~1, #extinction probability
       pformula = ~temp+cloud+wind+noise+minAmid+obs+rep+USFS_s, #detection probability 
       umf)


# OLDER CODE -------


rcwo_det<-spfilter('RCWO') %>% 
  #select(-`Radial distance`, -`Type AVB`, -`Cluster size`, -notes, -`C or I`,
  #       -Study.Area, -size, -Effort, -Sample.Label, -Region.Label, -distance, -Area,
  #       -`Survey effort`) %>%
  select(Species, `Pt name`, Replicate, Year) %>%
  group_by(`Pt name`, Replicate, Year) %>% slice(1) %>%
  mutate(presence=ifelse(is.na(Species), 0, 1),
         visit=paste0('v',Year, '.', Replicate)) %>%
  group_by(`Pt name`) %>% select(`Pt name`, visit, presence) %>%
  arrange(visit) %>%
  pivot_wider(names_from=visit, values_from=presence) %>%
  ungroup() %>% select(-`Pt name`) %>% as.matrix()
  
umf <- unmarkedMultFrame(y=rcwo_det, #detection record
                         siteCovs=occ_cov[,3],#site_cov, #site covariates
                         yearlySiteCovs=year_covl,
                         obsCovs=phi_cov,
                         numPrimary=2)

rcwo_0<-colext(psiformula = ~ USFS_s, #inital probability of occupancy at each site
               gammaformula = ~1, #colonization probability
               epsilonformula = ~1, #extinction probability
               pformula = ~temp+cloud+wind+noise+minAmid+obs+rep+USFS_s, #detection probability 
               umf)

library(MuMIn)
rcw_occ_d<-dredge(rcwo_0) #start time 7:57, end 8:11
rcw_occ_delta<-get.models(rcw_occ_d, subset = delta <= 2.5)
head(rcw_occ_d)
names(phi_cov)
rcw_occ_delta$`242`

quantile(occ_cov$USFS_s)
nd <- data.frame(USFS_s=seq(-2.6, 1.6, length=50))
E.psi <- predict(rcw_occ_delta$`242`, type="psi", newdata=nd, appendData=TRUE)
with(E.psi, {
  plot(USFS_s, Predicted, ylim=c(0,1), type="l",
       xlab="scaled RCWO score",
       ylab=expression(hat(psi)), cex.lab=0.8, cex.axis=0.8)
  lines(USFS_s, Predicted+1.96*SE, col=gray(0.7))
  lines(USFS_s, Predicted-1.96*SE, col=gray(0.7))
})

#TUTORIAL
# https://cran.r-project.org/web/packages/unmarked/vignettes/colext.pdf

data("crossbill")
head(crossbill)

# ok so in this example, each column starting with det is a visit to a site (3 visits for 9 years). 
# columns starting with date is the julian day a site was visited 
# the last three characters in column names are unique for each site visit

DATE <- as.matrix(crossbill[,32:58]) # julian days (detection covariate)
DATE

y.cross <- as.matrix(crossbill[,5:31]) # detection records
y.cross[is.na(DATE) != is.na(y.cross)] <- NA

# remember to scale your continuous variables
sd.DATE <- sd(c(DATE), na.rm=TRUE)
mean.DATE <- mean(DATE, na.rm=TRUE)
DATE <- (DATE - mean.DATE) / sd.DATE

#Before we can fit occupancy models, we need to format this data set appropriately.
years <- as.character(1999:2007)
years <- matrix(years, nrow(crossbill), 9, byrow=TRUE)
umf <- unmarkedMultFrame(y=y.cross, #detection record
                         siteCovs=crossbill[,2:3], #site covariates
                         yearlySiteCovs=list(year=years),
                         obsCovs=list(date=DATE),
                         numPrimary=9)

# A model with constant parameters
 fm0 <- colext(~1, ~1, ~1, ~1, umf)
# Like fm0, but with year-dependent detection
 fm1 <- colext(~1, ~1, ~1, ~year, umf)
# Like fm0, but with year-dependent colonization and extinction
 fm2 <- colext(~1, ~year-1, ~year-1, ~1, umf)
# A fully time-dependent model
 fm3 <- colext(~1, ~year-1, ~year-1, ~year, umf)
# Like fm3 with forest-dependence of 1st-year occupancy
 fm4 <- colext(~forest, ~year-1, ~year-1, ~year, umf)
# Like fm4 with date- and year-dependence of detection
fm5 <- colext(~forest, ~year-1, ~year-1, ~year + date + I(date^2),
              umf, starts=c(coef(fm4), 0, 0))
# Same as fm5, but with detection in addition depending on forest cover
fm6 <- colext(~forest, ~year-1, ~year-1, ~year + date + I(date^2) +
                    forest, umf)

models <- fitList('psi(.)gam(.)eps(.)p(.)' = fm0,
                  'psi(.)gam(.)eps(.)p(Y)' = fm1,
                  'psi(.)gam(Y)eps(Y)p(.)' = fm2,
                  'psi(.)gam(Y)eps(Y)p(Y)' = fm3,
                  'psi(F)gam(Y)eps(Y)p(Y)' = fm4,
                  'psi(F)gam(Y)eps(Y)p(YD2)' = fm5,
                  'psi(F)gam(Y)eps(Y)p(YD2F)' = fm6)
ms <- modSel(models)
ms

op <- par(mfrow=c(1,2), mai=c(0.8,0.8,0.1,0.1))
nd <- data.frame(forest=seq(0, 100, length=50))
E.psi <- predict(fm6, type="psi", newdata=nd, appendData=TRUE)
with(E.psi, {
  plot(forest, Predicted, ylim=c(0,1), type="l",
       xlab="Percent cover of forest",
       ylab=expression(hat(psi)), cex.lab=0.8, cex.axis=0.8)
  lines(forest, Predicted+1.96*SE, col=gray(0.7))
  lines(forest, Predicted-1.96*SE, col=gray(0.7))
})
nd <- data.frame(date=seq(-2, 2, length=50),
                   year=factor("2005", levels=c(unique(years))),
                   forest=50)
E.p <- predict(fm6, type="det", newdata=nd, appendData=TRUE)
E.p$dateOrig <- E.p$date*sd.DATE + mean.DATE
with(E.p, {
  plot(dateOrig, Predicted, ylim=c(0,1), type="l",
       xlab="Julian date", ylab=expression( italic(p) ),
       cex.lab=0.8, cex.axis=0.8)
  lines(dateOrig, Predicted+1.96*SE, col=gray(0.7))
  lines(dateOrig, Predicted-1.96*SE, col=gray(0.7))
})
par(op)