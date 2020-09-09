### Loading Distance Sampling data ###
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
