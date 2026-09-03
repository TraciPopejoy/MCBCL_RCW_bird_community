test<-YTWA.mod105$dht$individuals$bysample %>% as_tibble()
test 

phi<-bind_cols(YTWA.mod105$ddf$data, 
          predict.phi=predict(YTWA.mod105,esw=FALSE)$fitted) %>%
  dplyr::select(Sample.Label,predict.phi)
site.lev<-left_join(test, phi, by=c('Sample'='Sample.Label'))%>% 
  dplyr::select(-Region, -Effort)%>%
  group_by(Sample) %>%
  mutate(nchat=n/predict.phi, 
         nhat=nchat/(Sample.Area*592),
         dens=(1/Sample.Area)*nchat) %>%
  dplyr::select(-predict.phi, -Area,
                -n) %>%
  ungroup() %>%
  filter(!duplicated(.))
colSums(site.lev[,-1], na.rm=T)
YTWA.mod105$dht$individuals$Nhat.by.sample[1:5,]
sum(test$Nhat)

samp.area<-YTWA.mod10$dht$individuals$bysample$Sample.Area[1]
temp<-bind_cols(YTWA.mod10$ddf$data, 
          predict.phi=predict(YTWA.mod10,esw=FALSE)$fitted) %>%
  group_by(Sample.Label, Year) %>%
  dplyr::summarize(phi=predict.phi,
                   tob=sum(size), 
            nchat=tob/predict.phi, 
            nhat=nchat/(samp.area*592),
            dens=(1/samp.area)*nchat,
            .groups='drop')%>%
  filter(!duplicated(.))
sum(temp$tob); sum(YTWA.mod10$ddf$data$size)
sum(temp$nhat)
sum(temp$dens)/592
summary(YTWA.mod10)
0.529/(1*100) #birds/ha

test$dht$individuals$bysample %>%
  as_tibble() %>%
  mutate(Nhat/Sample.Area)
samp.area<-NOCA.mod31$dht$individuals$bysample[1,4]
beep_beep<-bind_cols(NOCA.mod31$ddf$data, 
            predict.phi=predict(NOCA.mod31,esw=FALSE)$fitted) %>%
  group_by(Sample.Label, Year) %>%
  dplyr::summarize(phi=predict.phi,
                   tob=sum(size), 
                   nchat=tob/predict.phi, 
                   nhat=nchat/(samp.area*592),
                   dens=(1/samp.area)*nchat,
                   NEWnhat=nchat/samp.area,
                   raw_dens_no_phi=tob/samp.area,
                   .groups='drop')%>%
  filter(!duplicated(.))
mean(c(beep_beep$NEWnhat, rep(0, 1178 - (beep_beep %>% nrow()))))
quantile(c(beep_beep$NEWnhat, rep(0, 1178 - (beep_beep %>% nrow()))))
quantile(c(beep_beep$phi))
quantile(c(beep_beep$raw_dens_no_phi, rep(0, 1178 - (beep_beep %>% nrow()))))
sum(beep_beep$dens/592)

fuck<-bind_cols(NOCA.mod31$ddf$data, 
          predict.phi=predict(NOCA.mod31,esw=FALSE)$fitted) %>%
  group_by(Sample.Label, Year) %>% 
  dplyr::summarize(Sample.Label=paste(Sample.Label, Year, sep='.'),
                   phi=predict.phi,
                   tob=sum(size),
                   HabCat=HabCat,
                   USFS_hab_index=USFS_hab_index,
                   .groups='drop') %>%
  bind_rows(point_sample_empty %>% 
              summarize(Sample.Label = paste(`Pt name`, Replicate, Year, sep ='.'),
                        phi=0, 
                        tob=0,
                        HabCat=HabCat,
                        USFS_hab_index=USFS_hab_index)) %>%
  group_by(Sample.Label) %>%
  arrange(desc(tob)) %>%  slice(1) %>%
  dplyr::mutate(nchat=tob/phi, 
                   nhat=nchat/(samp.area*592),
                   dens=(1/samp.area)*nchat,
                   NEWnhat=nchat/samp.area,
                   raw_dens_no_phi=tob/samp.area, 
                PtName=gsub("\\..*","", Sample.Label)) %>%
  replace_na(list(nchat=1e-6, nhat=1e-6, dens=1e-6, NEWnhat=1e-6)) %>%
  group_by(PtName, HabCat, USFS_hab_index) %>%
  dplyr::summarise(pt.phi=mean(phi),
                   pt.tob=mean(tob),
                   pt.NEWnhat=mean(NEWnhat))

hist(fuck$pt.NEWnhat)
summary(glm(pt.NEWnhat~HabCat, data=fuck, family=Gamma))
beep<-summary(lm(log1p(pt.NEWnhat)~USFS_hab_index, data=fuck))
beep
coef(beep)
ggplot(fuck, aes(x=USFS_hab_index, y=pt.NEWnhat)) +
  geom_point()+
  geom_smooth(method='lm')+
  scale_y_continuous(trans='log1p')+
  scale_x_continuous(limits = c(3,5))

# trying to do model averaging -----
library(MuMIn)
ffff<-model.appendix %>%
  filter(species=='BHCO',
         grepl('mod35', model.name)) %>% pull(formula)
A <- 10.5
call("round", A)        # round(10.5)
call("round", quote(A)) # round(A)
f <- "round"

BHCO.mod35TEST<-BHCO.mod35
terms(BHCO.mod35TEST$call)
update(as.formula(models.to.run$formula[h]), ffff)
BHCO.mod35TEST$call
model.avg(BHCO.mod35TEST, BHCO.mod95)

names(BHCO.mod35$ddf$ds$aux$ddfobj$scale$parameters)[-1]
class(BHCO.mod35$call)


install.packages('AICcmodavg')
library(AICcmodavg)
modavg()
ml<-list('one'=BHCO.mod35, 'two'=BHCO.mod11)
aictab(cand.set = ml)
class(BHCO.mod11)
library(RMark)
estimate <- c(BHCO.mod35$ddf$Nhat,
              BHCO.mod95$ddf$Nhat)
AIC.values <- c(BHCO.mod35$ddf$criterion,
                BHCO.mod95$ddf$criterion)
# Nhat.se is calculated in mrds:::summary.io, not in ddf(), so
# it takes a bit to pull out
std.err <- c(summary(BHCO.mod35$ddf)$Nhat.se,
             summary(BHCO.mod95$ddf)$Nhat.se)
MuMIn::model.avg(BHCO.mod35$ddf,
                 BHCO.mod95$ddf)

mmi.list=list(estimate=estimate, AIC=AIC.values, se=std.err)
RMark::model.average(mmi.list, revised=TRUE)

# checking if scale matters -----



conversion.factor <- convert_units(distance_units = "Metre", 
                                   effort_units=NULL, 
                                   area_units = "square kilometre")
newYTWA_dat<-YTWA_dat %>%
  mutate(Noise.scale=scale(Noise),
         Relicate.scale=scale(Replicate))
YTWA.mod105.scaled<-ds(newYTWA_dat,
   transect="point", key="hr",
   formula=~Noise.scale+OBS+Relicate.scale, 
   adjustment = NULL, order = 0,
   truncation = '5%', 
   convert.units = conversion.factor)
View(model.appendix)
YTWA.mod105.scaled$ddf$par
YTWA.mod105$ddf$par
AIC(YTWA.mod105, YTWA.mod105.scaled)

#AMCR.mod86
#Wind.Noise.nMinAfterMid
newAMCR_dat<-AMCR_dat %>%
  mutate(Noise.scale=scale(Noise),
         Wind.scale=scale(Wind),
         min.scale=scale(nMinAfterMid))
AMCR.mod86.scaled<-ds(newAMCR_dat,
                       transect="point", key="hr",
                       formula=~Wind.scale+Noise.scale+min.scale, 
                       adjustment = NULL, order = 0,
                       truncation = '5%', 
                       convert.units = conversion.factor)
AMCR.mod86<-ds(AMCR_dat,
                      transect="point", key="hr",
                      formula=~Wind+Noise+nMinAfterMid, 
                      adjustment = NULL, order = 0,
                      truncation = '5%', 
                      convert.units = conversion.factor)

AMCR.mod86.scaled$ddf$par
AMCR.mod86$ddf$par
AIC(AMCR.mod86.scaled, AMCR.mod86, AMCR.mod16, AMCR.mod70, AMCR.mod92) %>% arrange(AIC)

#so pretty confident scaling wouldn't affect the AIC values. Models selected are good. 
#BUT to predict on my models, I likely should scale this shit... wouldn't it cancel out?

bind_cols(AMCR.mod86$ddf$data, 
          predict.phi=predict(AMCR.mod86,esw=FALSE)$fitted) %>%
  slice(1:5) %>%
  group_by(Sample.Label, Year) %>% 
  dplyr::summarize(Sample.Label=paste(Sample.Label, Year, sep='.'),
                   phi=predict.phi,
                   tob=sum(size),
                   tob/phi,
                   HabCat=HabCat,
                   USFS_hab_index=USFS_hab_index,
                   .groups='drop') 

bind_cols(AMCR.mod86.scaled$ddf$data, 
          predict.phi=predict(AMCR.mod86.scaled,esw=FALSE)$fitted) %>%
  slice(1:5) %>%
  group_by(Sample.Label, Year) %>% 
  dplyr::summarize(Sample.Label=paste(Sample.Label, Year, sep='.'),
                   phi=predict.phi,
                   tob=sum(size),
                   tob/phi,
                   HabCat=HabCat,
                   USFS_hab_index=USFS_hab_index,
                   .groups='drop') 

