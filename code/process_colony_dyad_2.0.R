Packages <- c("dplyr", "ggplot2","stringr","fields", "tidyr",  "lubridate", "RColorBrewer", "igraph", "lubridate")

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

setwd('/local/home/katrinac/parentage/')

#download.file(url = "https://raw.githubusercontent.com/pinskylab/genomics/master/data/known_issues.csv", destfile = "~/parentage/r_data/known_issues.csv")
issues <- read.csv("~/parentage/r_data/known_issues.csv", header=T, stringsAsFactors = F)
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/fish_db.RData?raw=true", destfile = "~/parentage/r_data/dives_db.RData")
load("~/parentage/r_data/fish_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/anem_db.RData?raw=true", destfile = "~/parentage/r_data/anem_db.RData")
load("~/parentage/r_data/anem_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/dives_db.RData?raw=true", destfile = "~/parentage/r_data/dives_db.RData")
load("~/parentage/r_data/dives_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/fish_db.RData?raw=true", destfile = "~/parentage/r_data/dives_db.RData")
load("~/parentage/r_data/fish_db.RData")
load("~/parentage/r_data/gps_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/anems_tagged.RData", destfile = "~/parentage/r_data/anems_tagged.RData")
load("~/parentage/r_data/anems_tagged.RData")
#download.file(url = "https://github.com/pinskylab/genomics/blob/master/data/fish-obs.RData?raw=true", destfile = "~/parentage/r_data/fish-obs.RData")
fish_obs <- readRDS("~/parentage/r_data/fish-obs.RData") 
#download.file(url = "https://raw.githubusercontent.com/pinskylab/db_backups/master/ligation_1-7-20.csv", destfile = "~/parentage/r_data/ligation_db.csv")
lig_db <- read.csv("~/parentage/r_data/ligation_db.csv", header=T, stringsAsFactors = F)
#download.file(url = "https://raw.githubusercontent.com/pinskylab/db_backups/master/digest_1-7-20.csv", destfile = "~/parentage/r_data/digest_db.csv")
dig_db <- read.csv("~/parentage/r_data/digest_db.csv", header=T, stringsAsFactors = F)
#download.file(url = "https://raw.githubusercontent.com/pinskylab/db_backups/master/extraction_1-7-20.csv", destfile = "~/parentage/r_data/extraction_db.csv")
ext_db <- read.csv("~/parentage/r_data/extraction_db.csv", header=T, stringsAsFactors = F)


#read distance matrices in 
dist12 <- read.csv("/local/home/katrinac/parentage/kernel_fitting/1340_loci/input/distance_matrix12.csv", header=F)
dist13 <- read.csv("/local/home/katrinac/parentage/kernel_fitting/1340_loci/input/distance_matrix13.csv", header=F)
dist14 <- read.csv("/local/home/katrinac/parentage/kernel_fitting/1340_loci/input/distance_matrix14.csv", header=F)
dist15 <- read.csv("/local/home/katrinac/parentage/kernel_fitting/1340_loci/input/distance_matrix15.csv", header=F)
dist16 <- read.csv("/local/home/katrinac/parentage/kernel_fitting/1340_loci/input/distance_matrix16.csv", header=F)
dist17 <- read.csv("/local/home/katrinac/parentage/kernel_fitting/1340_loci/input/distance_matrix17.csv", header=F)
dist18 <- read.csv("/local/home/katrinac/parentage/kernel_fitting/1340_loci/input/distance_matrix18.csv", header=F)
dist_all <- read.csv("/local/home/katrinac/parentage/kernel_fitting/1340_loci/input/distance_matrix_sampled.csv", header = F)

"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0


##SEARCH AND REPLACE FOR YEAR 4 MATCHES
#par1 <- read.table(file= "/local/home/katrinac/parentage/colony2/20190523_1340loci/results/1340loci_2012.Maternity", header= TRUE, stringsAsFactors = F)	
#par2 <- read.table(file= "/local/home/katrinac/parentage/colony2/20190523_1340loci/results/1340loci_2012.Paternity", header= TRUE, stringsAsFactors = F)	
#
#
#par1 <- par1 %>% 
#    rename(offs_lig="OffspringID", par_lig = "InferredMum1", prob="ProbMum1")
#    
#par2 <- par2 %>% 
#    rename(offs_lig="OffspringID", par_lig = "InferredDad1", prob= "ProbDad1")
#
##bind together into 1 df
#par12 <- bind_rows(par1, par2) %>%
#    mutate(year="2012")
#

##bind all year parentage together
#allpar <- bind_rows(par12, par13, par14, par15, par16, par17, par18) %>% 
#    filter(par_lig !="L1091")#remove APCL13_128, the parent that moved (see later in script for explaination comment)
#nrow(allpar) #should be 
##replace NAs with zeros
##allpar$ProbMum1[is.na(allpar$ProbMum1)] <- 0
##allpar$ProbDad1[is.na(allpar$ProbDad1)] <- 0

#allpar_clean <- allpar %>%
#    group_by(offs_lig, par_lig) %>%
#    mutate(prob = sum(prob, na.rm = T)) %>% #add the probabilities when a parent was matched twice
#    ungroup()%>%
#    filter(prob >= 0.95) #remove low prob matches
#nrow(allpar_clean)
#nrow(distinct(allpar_clean, offs_lig))

#write.csv(allpar_clean, file="~/parentage/colony2/20190523_1340loci/allpar_clean.csv", row.names=FALSE, quote=FALSE)




allpar_clean <- read.csv(file="~/parentage/colony2/20190523_1340loci/allpar_clean.csv", header=T, stringsAsFactors = F)

nrow(allpar_clean)

#summary numbers
n_total_matches <- allpar_clean %>%
    group_by(offs_lig) %>%
    mutate(n_parents=ifelse(n()==2, 2, 1 )) %>%
    distinct(offs_lig, .keep_all = T)


#so 62 parentage matches, of which 48 are trios
#how many in each year?
sum_offs <- allpar_clean %>%
    group_by(year) %>%
    distinct(offs_lig, .keep_all = T) %>%
    summarise(n_offs_matched=n())

sum(sum_offs$n_offs_matched)#check, should be 62
sum(n_total_matches$n_parents[n_total_matches$n_parents ==1])#should be 38
sum(n_total_matches$n_parents[n_total_matches$n_parents ==2])#should be 48
#REMEMBER I REMOVED ONE PARENT FROM THE RESULTS BECAUSE THE SITE ATTACHED TO THAT PARENT ISN'T WHERE THE LARVAE LIKELY ORIGINATED FROM
#87 total matches, and of those 24 are repeats (trio matches). So, 86-48=38 maches to one parent. Then 24 matches to two parents, for a total of 62 rows of distinct offspring lig ids.


#from Allison, just putting all the meta data together (Constants_database_common_functions.R)
##### Match up other relevant info (site, date, fish_indiv, etc.) to fish in the clownfish table
# Pull out year and month into a separate column in dives_db
dives_db_processed <- dives_db %>%
  mutate(year = as.integer(substring(date,1,4))) %>%
  mutate(month = as.integer(substring(date,6,7))) %>%
  mutate(dive_date = date(date))

# Pull all APCL caught or otherwise in the clownfish table
allfish_fish <- fish_db %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, anem_table_id, recap, tag_id, color, sex, size, fish_obs_time, fish_notes) %>%
  filter(fish_spp == 'APCL') %>%
  mutate(size = as.numeric(size))  # make the size numeric (rather than chr) so can do means and such

# and their corresponding anemones
allfish_anems <- anem_db %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_notes, anem_obs_time) %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

# and the corresponding dive info
allfish_dives <- dives_db_processed %>%
  select(dive_table_id, dive_type, date, year, month, site, gps, dive_notes) %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id) 

# join together
allfish_caught <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish_caught <- left_join(allfish_caught, allfish_dives, by="dive_table_id")

# add in the gen_ids and fish_indiv (now in a separate gen_id table) - gen_id only comes in the time the fish was sequenced, not at all captures
allfish_caught <- left_join(allfish_caught, (fish_obs %>% select(fish_table_id, gen_id, fish_indiv)), by = "fish_table_id")

#get sample ids for ligation ids
ext_db <- ext_db %>%
    select(extraction_id, sample_id)
dig_db <- dig_db %>%
    select(extraction_id, digest_id)
lig_db <- lig_db %>%
    select(ligation_id, digest_id)
ext_dig <- left_join(ext_db, dig_db, by="extraction_id")
dig_lig <- left_join(ext_dig, lig_db, by="digest_id")
lig_samp <- dig_lig %>%
    select(ligation_id, sample_id) %>%
    filter(sample_id %in% allfish_caught$sample_id) #remove everything but sampled APCL

fish_meta <- left_join(allfish_caught, lig_samp,  by="sample_id") %>%
    select(fish_indiv, size, color, sex, gen_id, ligation_id, sample_id, site, date, anem_obs_time, gps)


#get specific metadata for offspring
allpar_samp1 <- fish_meta %>%
    filter(ligation_id %in% allpar_clean$offs_lig) 
colnames(allpar_samp1) <- paste("offs", colnames(allpar_samp1), sep = "_")

#get specific metadata for parents
allpar_samp2 <- fish_meta %>%
    filter(ligation_id %in% allpar_clean$par_lig) 
colnames(allpar_samp2) <- paste("par", colnames(allpar_samp2), sep = "_")

#join back into parentage match format
meta_with_offs <- left_join(allpar_clean, allpar_samp1, by=c(offs_lig="offs_ligation_id"))

par_res <- left_join(meta_with_offs, allpar_samp2, by=c(par_lig="par_ligation_id"))

#the long haul of adding lat lons
results_meta <- par_res %>%
    mutate(offs_time_date =as.character(str_c(offs_date, offs_anem_obs_time, sep = " "))) %>%
    mutate(offs_time_date = ymd_hms(offs_time_date))%>%
    mutate(offs_time_date = force_tz(offs_time_date, tzone = "Asia/Manila")) %>%
    mutate(offs_time_date = with_tz(offs_time_date, tzone = "UTC")) %>%
    mutate(offs_year = year(offs_time_date)) %>%
    mutate(offs_month = month(offs_time_date)) %>%
    mutate(offs_day = day(offs_time_date)) %>%
    mutate(offs_hour = hour(offs_time_date)) %>%
    mutate(offs_minute = minute(offs_time_date)) %>%
    mutate(par_time_date =as.character(str_c(par_date, par_anem_obs_time, sep = " "))) %>%
    mutate(par_time_date = ymd_hms(par_time_date))%>%
    mutate(par_time_date = force_tz(par_time_date, tzone = "Asia/Manila")) %>%
    mutate(par_time_date = with_tz(par_time_date, tzone = "UTC")) %>%
    mutate(par_year = year(par_time_date)) %>%
    mutate(par_month = month(par_time_date)) %>%
    mutate(par_day = day(par_time_date)) %>%
    mutate(par_hour = hour(par_time_date)) %>%
    mutate(par_minute = minute(par_time_date))

gps <- gps_db %>%
    mutate(lat=as.numeric(lat)) %>%
    mutate(lon=as.numeric(lon)) %>%
    mutate(time_date = ymd_hms(time)) %>%
    mutate(year = year(time_date)) %>%
    mutate(month = month(time_date)) %>%
    mutate(day = day(time_date)) %>%
    mutate(hour = hour(time_date)) %>%
    mutate(minute = minute(time_date)) %>%
    select(-time_date, -elev, -time) %>%
    group_by(unit, year, month, day, hour, minute) %>%
    mutate(lat =mean(lat)) %>% #create within minute averages
    mutate(lon=mean(lon)) %>%
    distinct(.keep_all = T)

add_offs_loc <- left_join(results_meta, gps, by=c(offs_gps="unit",offs_year="year", offs_month="month", offs_day="day", offs_hour="hour", offs_minute="minute")) %>%
    rename(offs_lat="lat", offs_lon="lon") 

results_dist <- left_join(add_offs_loc, gps, by=c(par_gps="unit",par_year="year", par_month="month", par_day="day", par_hour="hour", par_minute="minute")) %>%
    rename(par_lat="lat", par_lon="lon") 
dim(results_dist) #should be 86



#alldists <- rdist.earth(as.matrix(results_loc[,c('par1_lon', 'par1_lat')]), as.matrix(pairlatlon12[,c('par2_lon', 'par2_lat')]), miles=FALSE, R=6371)
#pairlatlon12$dist_trios_km <- diag(alldists)


#calculate the distance from mom then from dad, use the average between the two to get net dispersal distance
alldists <- rdist.earth(as.matrix(results_dist[,c('offs_lon', 'offs_lat')]), as.matrix(results_dist[,c('par_lon', 'par_lat')]), miles=FALSE, R=6371)
results_dist$dist_par_km <- diag(alldists)

#png('~/parentage/colony2/20190523_1340loci/results/dispersal_dist_hist.png')
hist(results_dist$dist_par_km) #cool.
#dev.off()

write.table(results_dist, file="~/parentage/kernel_fitting/1340_loci/parentage_results_allyears.csv", row.names=FALSE, quote=FALSE, col.names=T, sep=",")
write.table(results_dist, file="~/parentage/colony2/20190523_1340loci/results/parentage_results_allyears.csv", row.names=FALSE, quote=FALSE, col.names=T, sep=",")

#this cell for bode kernel estimate input
par12 <- results_dist %>%
    filter(year=="2012") %>%
    distinct(offs_sample_id, .keep_all = T) %>%
    select(year, offs_site, par_site) #sort by parent year for the table for malin
dim(par12)

par13 <- results_dist %>%
    filter(year=="2013") %>%
    distinct(offs_sample_id, .keep_all = T) %>%
    select(year, offs_site, par_site) #sort by parent year for the table for malin
dim(par13)

par14 <- results_dist %>%
    filter(year=="2014") %>%
    distinct(offs_sample_id, .keep_all = T) %>%
    select(year, offs_site, par_site) #sort by parent year for the table for malin
dim(par14)

par15 <- results_dist %>%
    filter(year=="2015") %>%
    distinct(offs_sample_id, .keep_all = T) %>%
    select(year, offs_site, par_site) #sort by parent year for the table for malin
dim(par15)

par16 <- results_dist %>%
    filter(year=="2016") %>%
    distinct(offs_sample_id, .keep_all = T) %>%
    select(year, offs_site, par_site) #sort by parent year for the table for malin
dim(par16)

par17 <- results_dist %>%
    filter(year=="2017") %>%
    distinct(offs_sample_id, .keep_all = T) %>%
    select(year, offs_site, par_site) #sort by parent year for the table for malin
dim(par17)

par18 <- results_dist %>%
    filter(year=="2018") %>%
    distinct(offs_sample_id, .keep_all = T) %>%
    select(year, offs_site, par_site) #sort by parent year for the table for malin
dim(par18)

par_all <- results_dist %>%
    select(year, offs_site, par_site)
nrow(par12)+nrow(par13)+nrow(par14)+nrow(par15)+nrow(par16)+nrow(par17)+nrow(par18)
nrow(par_all)

#write for bode dispersal kernel estimate with all matches
#write.csv(par12, file="~/parentage/kernel_fitting/1340_loci/parentage12.csv", quote=TRUE, row.names= FALSE)
#write.csv(par13, file="~/parentage/kernel_fitting/1340_loci/parentage13.csv", quote=TRUE, row.names= FALSE)
#write.csv(par14, file="~/parentage/kernel_fitting/1340_loci/parentage14.csv", quote=TRUE, row.names= FALSE)
#write.csv(par15, file="~/parentage/kernel_fitting/1340_loci/parentage15.csv", quote=TRUE, row.names= FALSE)
#write.csv(par16, file="~/parentage/kernel_fitting/1340_loci/parentage16.csv", quote=TRUE, row.names= FALSE)
#write.csv(par17, file="~/parentage/kernel_fitting/1340_loci/parentage17.csv", quote=TRUE, row.names= FALSE)
#write.csv(par18, file="~/parentage/kernel_fitting/1340_loci/parentage18.csv", quote=TRUE, row.names= FALSE)
#write.csv(par_all, file="~/parentage/kernel_fitting/1340_loci/parentage_all.csv", quote=TRUE, row.names= FALSE)


#use distance matrix to plot dispersal in each year
distall_df <- data.frame(dist=as.vector(as.matrix(dist_all)))

#get a df of all distances between sampled fish as possible dispersal distances
#both axes have to be proportions for stupid ggplot
all_sampled_fish <- fish_obs  %>%
    distinct(fish_indiv, `.keep_all` = T)
all_sampled_fish_loc <- get_latlon(all_sampled_fish$sample_id) %>%
    filter(!is.na(lat) & !is.na(lon))
all_sampled_dists <- rdist.earth(as.matrix(all_sampled_fish_loc[,c('lon', 'lat')]), as.matrix(all_sampled_fish_loc[,c('lon', 'lat')]), miles=FALSE, R=6371)

all_sampled_dists_df <- data.frame(dist=as.vector(as.matrix(all_sampled_dists))) %>%
    mutate(dist=round(dist, digits=0)) %>%
    group_by(dist) %>%
    mutate(freq=n()) %>%
    mutate(prop=freq/sum(freq))
           
obs_disp_dist_prop <- results_dist %>%
    #filter(year==2018) %>%
    distinct(offs_fish_indiv, .keep_all = T) %>% 
    mutate(dist_par_km = ifelse(dist_par_km >27, floor(dist_par_km), dist_par_km)) %>%
    mutate(dist_par_km = ifelse(dist_par_km <1, ceiling(dist_par_km), dist_par_km)) %>%
    group_by(dist_par_km) %>%
    mutate(n_freq=n()) %>%
    ungroup() %>%
    mutate(prop=n_freq/sum(n_freq))


options(scipen=999)
# use these variables to set the limits on all plots
#y1max = max(obs_disp_dist_prop$dist_par_km)

#quartzFonts(avenir = c("Avenir Book", "Avenir Black", "Avenir Book Oblique", 
#        "Avenir Black Oblique"))
par(family = "Helvetica")
#pdf("~/parentage/colony2/20190523_1340loci/results/disp_dist_pub.pdf")


par(mar=(c(5, 5, 5, 5))) 
with(all_sampled_dists_df, hist(dist, ylim = c(0, max(freq)+2000000), main=NA, axes=F, xlab=NA, ylab=NA, breaks=seq(0, 30, 1), col=scales::alpha('gray',.5), border=F))
legend("topright", legend=c("Observed dispersal", "Sampled distances"),
       fill=c(scales::alpha('orange',.5), scales::alpha('gray',.5)), cex=0.8, box.lty=0)
axis(4)
mtext(side = 4, line = 3, "Frequency")
#lines(density(per_cur), col="blue",lwd=2)

#add obs with a secondary y-axis

par(new = TRUE)
with(obs_disp_dist_prop, hist(dist_par_km, ylim = c(0, 50), main=NA, breaks=seq(0, 30, 1), xlab="Distance (km)",col=scales::alpha('orange',.5), border=F))

#dev.off()

#find means 
results_dist$year <- as.character(results_dist$year)
means <- results_dist %>%
    ungroup() %>%
    group_by(year) %>%
    distinct(offs_fish_indiv, .keep_all = T) %>%
    summarise(mean_disp_dist=mean(dist_par_km))
 

dispersal_raw_results <- left_join(sum_offs, means, by="year")
dispersal_raw_results
#write.table(dispersal_raw_results, file="~/parentage/colony2/20190523_1340loci/results/dispersal_raw_results.txt", row.names=FALSE, quote=FALSE, col.names=T, sep=" ")

#png('~/parentage/colony2/20190523_1340loci/results/dispersal_dist_hist.png')
plot(dispersal_raw_results$mean_disp_dist, dispersal_raw_results$n_offs_matched)
#dev.off()

