Packages <- c("dplyr", "ggplot2", "fields","stringr", "reshape2", "dplyr", "tidyr", "lubridate", "RColorBrewer")

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

setwd('/local/home/katrinac/parentage/kernel_fitting/')

load("~/parentage/r_data/total_sampling_across_years.RData")
load("~/parentage/r_data/sampled_area_each_year.RData")
load("~/parentage/r_data/cumulative_prop_hab_sampled_by_site.RData")
#download.file(url = "https://github.com/pinskylab/genomics/blob/master/data/fish-obs.RData?raw=true", destfile = "~/parentage/r_data/fish-obs.RData")
fish_obs <- readRDS("~/parentage/r_data/fish-obs.RData") 
load("~/parentage/r_data/site_dist_info.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/anem_db.RData?raw=true", destfile = "~/parentage/r_data/anem_db.RData")
load("~/parentage/r_data/anem_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/dives_db.RData?raw=true", destfile = "~/parentage/r_data/dives_db.RData")
load("~/parentage/r_data/dives_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/fish_db.RData?raw=true", destfile = "~/parentage/r_data/dives_db.RData")
load("~/parentage/r_data/fish_db.RData")
load("~/parentage/r_data/gps_db.RData")

"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

#read in parentage matches
par12 <- read.csv(file="~/parentage/kernel_fitting/1340_loci/input/parentage12.csv", header=TRUE)
par12$offs_site <- gsub(". ", ".", par12$offs_site, fixed=TRUE)
par12$parent_site <- gsub(". ", ".", par12$parent_site, fixed=TRUE)

par13 <- read.csv(file="~/parentage/kernel_fitting/1340_loci/input/parentage13.csv", header=TRUE)
par13$offs_site <- gsub(". ", ".", par13$offs_site, fixed=TRUE)
par13$parent_site <- gsub(". ", ".", par13$parent_site, fixed=TRUE)

par14 <- read.csv(file="~/parentage/kernel_fitting/1340_loci/input/parentage14.csv", header=TRUE)
par14$offs_site <- gsub(". ", ".", par14$offs_site, fixed=TRUE)
par14$parent_site <- gsub(". ", ".", par14$parent_site, fixed=TRUE)

par15 <- read.csv(file="~/parentage/kernel_fitting/1340_loci/input/parentage15.csv", header=TRUE)
par15$offs_site <- gsub(". ", ".", par15$offs_site, fixed=TRUE)
par15$parent_site <- gsub(". ", ".", par15$parent_site, fixed=TRUE)

par16 <- read.csv(file="~/parentage/kernel_fitting/1340_loci/input/parentage16.csv", header=TRUE)
par16$offs_site <- gsub(". ", ".", par16$offs_site, fixed=TRUE)
par16$parent_site <- gsub(". ", ".", par16$parent_site, fixed=TRUE)

par17 <- read.csv(file="~/parentage/kernel_fitting/1340_loci/input/parentage17.csv", header=TRUE)
par17$offs_site <- gsub(". ", ".", par17$offs_site, fixed=TRUE)
par17$parent_site <- gsub(". ", ".", par17$parent_site, fixed=TRUE)

par18 <- read.csv(file="~/parentage/kernel_fitting/1340_loci/input/parentage18.csv", header=TRUE)
par18$offs_site <- gsub(". ", ".", par18$offs_site, fixed=TRUE)
par18$parent_site <- gsub(". ", ".", par18$parent_site, fixed=TRUE)

#read in demography estimates
prop_samp <- cumulative_prop_hab_sampled_by_site %>%
    mutate(total_possible_sample_anems = ifelse(site=="Caridad Proper", 4, total_possible_sample_anems) ) %>%
    mutate(total_prop_hab_sampled_anems_tidied= ifelse(site=="Caridad Proper" & total_anems_sampled==4, 1, total_prop_hab_sampled_anems_tidied) ) %>%
    mutate(total_possible_sample_anems = ifelse(site=="Sitio Lonas", total_anems_sampled, total_possible_sample_anems) ) %>%
    mutate(total_prop_hab_sampled_anems_tidied= ifelse(site=="Sitio Lonas", 1, total_prop_hab_sampled_anems_tidied) )

prop_samp$site <- gsub(". ", ".", prop_samp$site, fixed=TRUE)


##read in the sites that we sampled each year
N_gen_par <- read.table(file="~/parentage/colony2/20190523_1340loci/input/all_parents_corrected.txt", header = TRUE, stringsAsFactors = F) %>%#not sure that I need the parents here
    mutate(fish_indiv=as.character(fish_indiv))
N_gen_offs <- read.table(file="~/parentage/colony2/20190523_1340loci/input/all_offspring_corrected.txt", header=T, stringsAsFactors = F) %>%
    mutate(fish_indiv=as.character(fish_indiv))

#read in site geography
#dist <- read.csv("~/parentage/kernel_fitting/site_centroids.csv.csv", header=TRUE)
area <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_area_header_nonsurveyed.csv", header=TRUE))
centroids <- read.csv("~/parentage/kernel_fitting/1340_loci/input/site_centroids.csv", header=TRUE)
sampled_sites <- prop_samp %>% select(site) %>% distinct(site)
sampled_sites$site <- gsub(". ", ".", sampled_sites$site, fixed=TRUE)
site_widths <- read.table("/local/home/katrinac/parentage/text_file/site_widths.txt", header=T, sep=",")
##read in site names
#sites12 <- read.table(file="~/migest/annual/2012/20181017_mag_cab_corr/input_sites_2012.txt", header= TRUE)
#sites12$site <- gsub("_", " ", sites12$site, fixed=TRUE)
#sites12$site <- gsub(". ", ".", sites12$site, fixed=TRUE)
#
#
#sites13 <- read.table(file="~/migest/annual/2013/20181017_mag_cab_corr/input_sites_2013.txt", header=TRUE)
#sites13$site <- gsub("_", " ", sites13$site, fixed=TRUE)
#sites13$site <- gsub(". ", ".", sites13$site, fixed=TRUE)
#
#
#sites14 <- read.table(file="~/migest/annual/2014/20181017_mag_cab_corr/input_sites_2014.txt", header=TRUE)
#sites14$site <- gsub("_", " ", sites14$site, fixed=TRUE)
#sites14$site <- gsub(". ", ".", sites14$site, fixed=TRUE)
#
#sites15 <- read.table(file="~/migest/annual/2015/20181017_mag_cab_corr/input_sites_2015.txt", header=TRUE)
#sites15$site <- gsub("_", " ", sites15$site, fixed=TRUE)
#sites15$site <- gsub(". ", ".", sites15$site, fixed=TRUE)
#
#
#


#code to get summary stats of across site proportion of habitat sampled
summary(prop_samp$total_prop_hab_sampled_anems_tidied)
#prop_samp %>% filter(total_prop_hab_sampled_anems_tidied < 0.25)
mean(prop_samp$total_prop_hab_sampled_anems_tidied, na.rm = T)
sd(prop_samp$total_prop_hab_sampled_anems_tidied, na.rm = T)

#gather the summary of total offspring sampled
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
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_notes) %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

# and the corresponding dive info
allfish_dives <- dives_db_processed %>%
  select(dive_table_id, dive_type, date, year, month, site, gps, dive_notes) %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id) 

# join together
allfish_caught <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish_caught <- left_join(allfish_caught, allfish_dives, by="dive_table_id")

# add in the gen_ids and fish_indiv (now in a separate gen_id table) - gen_id only comes in the time the fish was sequenced, not at all captures
allfish_caught <- left_join(allfish_caught, (fish_obs %>% select(fish_table_id, gen_id, fish_indiv)), by = "fish_table_id") %>%
    select(fish_indiv, sample_id, site) %>%
    mutate()

N_gen_offs_annual  <- left_join(N_gen_offs, allfish_caught, by=c("fish_indiv", "sample_id")) %>% 
    group_by(year, site) %>%
    summarise(n_offs_gen=n()) %>%
    ungroup()

N_gen_offs_annual$site <- gsub(". ", ".", N_gen_offs_annual$site, fixed=TRUE)

##for all years
N_gen_offs_all <- N_gen_offs_annual %>% 
    group_by(site) %>% 
    summarise(sampled_fish=sum(n_offs_gen, na.rm=T))

sum(N_gen_offs_all$sampled_fish)#should be 792


head(area)

#index sites to get pop numbers
all_sites <- centroids %>%
    select(site) %>%
    arrange(site)
nrow(all_sites) #should be 35x1
all_sites$index <- seq(from=1, to=35, by=1)


all_sampled_sites <- inner_join(sampled_sites, all_sites) 
all_sampled_sites_index <- t(all_sampled_sites %>% select(index))

#for all years
#write.table(all_sampled_sites_index, file="~/parentage/kernel_fitting/1340_loci/site_index_all.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")



##FOR ANNUAL
#generate the correct proportion sampled table

#correct any NaN or infinities with true values based on whether or not we went to the site

for(i in 1:nrow(prop_samp)){
    
    
    if(is.nan(prop_samp$total_prop_hab_sampled_anems_tidied[i])){prop_samp$total_prop_hab_sampled_anems_tidied[i] <- 0} 
    if(is.infinite(prop_samp$total_prop_hab_sampled_anems_tidied[i])){prop_samp$total_prop_hab_sampled_anems_tidied[i] <- 1}
    
    
}


prop_samp16 <- prop_samp %>%
    filter(end_year=="2016" & total_prop_hab_sampled_anems_tidied >0) %>%
    select(site, total_prop_hab_sampled_anems_tidied)


##generate list of site indices
sites16 <- prop_samp16 %>%
    select(site)
sites16_2 <- semi_join(all_sites, sites16, by="site") %>%
    select(index)
sites16t <- t(sites16_2)
##write the site index file
#write.table(sites16t, file="~/parentage/kernel_fitting/1640_loci/site_index16.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")

##write the corrected proportion sampled file
prop_samp16 <- prop_samp16 %>%
    select(total_prop_hab_sampled_anems_tidied)
#write.table(prop_samp16, file="~/parentage/kernel_fitting/1640_loci/prop_samp16.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")




##FOR ALL YEARS
#generate the correct proportion sampled table
prop_samp_all <- prop_samp %>%
    filter(time_frame=="2012-2018") %>%
    arrange(site) %>%
    select(total_prop_hab_sampled_anems_tidied)
#write.table(prop_samp_all, file="~/parentage/kernel_fitting/1340_loci/prop_samp_all.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")

    



#make tallies for the parentage matches by offspring site 
total_par16 <- par16 %>%
    group_by(year, offs_site, parent_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() %>%
    select(-year)

#add in sites that were sampled, but there wasn't a match made 08/22/2016 wtf was I trying to do here...
sites16_beta <- sites16 %>%
    #select(-pop) %>%
    rename(parent_site="site")

sites16_beta$offs_site <- sites16_beta$parent_site

allsites_parentage16 <- full_join(sites16_beta, total_par16, by=c("parent_site", "offs_site")) %>%
    group_by(offs_site, parent_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, parent_site) #be very careful to keep everything in alphabetical order, consitent with the "pop" indexing used in MigEst.

#check for the correct number of matches
sum(allsites_parentage16$n_matches, na.rm=TRUE)

##find sites that aren't represented in both parent and offspring groups. This could happen if we didn't sample a parent site in the year we found the offspring
#test2 <- allsites_parentage16 %>%
#    ungroup() %>%
#    select(offs_site) 
#add2par <- allsites_parentage16 %>%
#    filter(offs_site %!in% test1$parent_site)    
#
#nrow(add2par) #proceed with below is nrow =<1

#
#test1 <- allsites_parentage16 %>%
#    ungroup() %>%
#    select(parent_site) 
#
#add2offs <- allsites_parentage16 %>%
#    filter(parent_site %!in% test2$offs_site)
#nrow(add2offs) #proceed with below is nobs >=1

#############parent sites to add to offspring
#        add2offs$offs_site <- add2offs$parent_site
#        add2offs$n_matches <- 0
#
############add these missing pairwise comparisons
#        allsites_parentage16 <- bind_rows(allsites_parentage16, add2offs)
#

################offspring sites to add to parent
#        add2par$parent_site <- add2par$offs_site
#        add2par$n_matches <- 0
#
##add these missing pairwise comparisons
#        allsites_parentage16 <- bind_rows(allsites_parentage16, add2par)
#
#

#check for the correct number of matches
sum(allsites_parentage16$n_matches, na.rm=TRUE)

#create the matrix with rows (parent pop) and columns (offs pop)
parmat16 <- allsites_parentage16 %>% 
    filter(!is.na(offs_site) & !is.na(parent_site)) %>%
    group_by(offs_site, parent_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

##change NAs to 0
parmat16[is.na(parmat16)] <- 0
rownames(parmat16) <- parmat16$parent_site
parmat16$parent_site <- NULL

#check for correct number, no NAs
sum(parmat16)

#check dimensions, should be symmetrical and number of populations sampled (9, 16, 16, 16)
dim(parmat16)

#select the sites visited each year and filter for year of sampling
N_gen_offs_annual$year <- as.numeric(N_gen_offs_annual$year)

N_gen_offs2 <- N_gen_offs_annual %>%
    filter(year==2016) %>%
    group_by(site, year) %>% #add group by year if doing annual
    summarise(sampled_fish=sum(n_offs_gen)) 



#make dataframe with column for unassigned juveniles in each column (offspring site)
offs_matched_site <- allsites_parentage16 %>%
    group_by(offs_site) %>%
    summarise(n_offs=sum(n_matches, na.rm=T))

unassigned_beta <- semi_join(offs_matched_site, sites16, by=c(offs_site="site")) %>%
    ungroup() %>%
    select(offs_site, n_offs)

#unassigned_beta[is.na(unassigned_beta)] <- 0


net_unassigned <- left_join(unassigned_beta, N_gen_offs2, by=c(offs_site="site")) %>%
    mutate(n_unassigned=sampled_fish-n_offs)
#transpose to create a row, bind to the matrix

row2add <- net_unassigned %>%
    arrange(offs_site) %>%
    select(n_unassigned, offs_site) %>%
    group_by(offs_site) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    select(n_unassigned)
row2addT <- as.data.frame(t(row2add))
colnames(row2addT) <- colnames(parmat16)

parmat16 <- ungroup(parmat16)
parmat_16_full <- bind_rows(parmat16, row2addT)
parmat_16_full[is.na(parmat_16_full)] <- 0


dim(parmat_16_full) #should have an extra row

#checking totals
sum(parmat16)
sum(row2addT, na.rm = T)

#write the matrix
#write.table(parmat_16_full, file="~/parentage/kernel_fitting/1340_loci/input/parentage_matrix16_corrected_for_missing_unassigned_recruits.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")

#make an all years parentage matrix to get average kernel
all_years_par <- rbind(par12, par13, par14, par15, par16, par17, par18)
#make tallies for the parentage matches by offspring site 
total_all_years_par <- all_years_par %>%
    group_by(offs_site, parent_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() 

sum(total_all_years_par$n_matches)

#test2 <- total_all_years_all %>%
#    ungroup() %>%
#    select(offs_site) 
#
#add2par <- total_all_years_all %>%
#    filter(offs_site %!in% test1$par_site)    
#
#nrow(add2par) #proceed with below if nrow =<1

##find sites that aren't represented in both parent and offspring groups
#
#test1 <- total_all_years_all %>%
#    ungroup() %>%
#    select(par_site) 
#
#add2offs <- total_all_years_all %>%
#    filter(par_site %!in% test2$offs_site)
#nrow(add2offs) #proceed with below if nrow >=1

#############parent sites to add to offspring
#        add2offs$offs_site <- add2offs$par_site
#        add2offs$n_matches <- 0
#
############add these missing pairwise comparisons
#        total_all_years_all <- bind_rows(total_all_years_all, add2offs)
#

################offspring sites to add to parent
#        add2par$parent_site <- add2par$offs_site
#        add2par$n_matches <- 0
#
##add these missing pairwise comparisons
#        total_all_years_par <- bind_rows(total_all_years_all, add2par)
#
#

sum(total_all_years_par$n_matches, na.rm=T)

#for all years averages
#add in sites that were sampled, but there wasn't a match made 08/22/2018 wtf was I trying to do here...
sites_all_beta <- sampled_sites %>%
    rename(par_site = "site") 
#undo here if there's a mistake 06/04/2019
sites_all_beta$offs_site <- sites_all_beta$par_site


allsites_parentage_all <- full_join(sites_all_beta, total_all_years_par, by=c(par_site="parent_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site) %>%#be very careful to keep everything in alphabetical order, consitent with the "pop" indexing used in MigEst.
    ungroup()
#allsites_parentage_all <- allsites_parentage_all %>%
#    group_by(offs_site, parent_site) %>%
#    summarise(n_matches=n()) %>%
#    ungroup() 

#check for the correct number of matches
sum(allsites_parentage_all$n_matches, na.rm=TRUE) #should be 62

#for all years average
#create the matrix with rows (parent pop) and columns (offs pop)
parmat_all <- allsites_parentage_all %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    spread(offs_site, n_matches)

###change NAs to 0
parmat_all[is.na(parmat_all)] <- 0
rownames(parmat_all) <- sites_all_beta$par_site
parmat_all$par_site <- NULL

#check for correct number, no NAs
sum(parmat_all)

#check dimensions, should be symmetrical and number of populations sampled (9, 17, 15, 15)
dim(parmat_all)

#make dataframe with column for unassigned juveniles in each column (offspring site)
unassigned_beta <- left_join(sampled_sites, total_all_years_par, by=c(site="offs_site")) %>%
    select(site, n_matches) %>%
    group_by(site) %>%
    summarise(n_matches=sum(n_matches))
unassigned_beta[is.na(unassigned_beta)] <- 0
                             

net_unassigned <- left_join(unassigned_beta, N_gen_offs_all, by="site") %>%
    mutate(n_unassigned=sampled_fish-n_matches)
#transpose to create a row, bind to the matrix

row2add <- net_unassigned %>%
    arrange(site) %>%
    select(n_unassigned, site) %>%
    group_by(site) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    select(n_unassigned)
row2addT <- as.data.frame(t(row2add))
colnames(row2addT) <- colnames(parmat_all)

parmat_all <- ungroup(parmat_all)
parmat_all_full <- bind_rows(parmat_all, row2addT)
parmat_all_full[is.na(parmat_all_full)] <- 0


dim(parmat_all_full) #should have an extra row

sum(parmat_all_full) #should be total number of fish in offspring file (for 1340 loci, 792)

#write.table(parmat_all_full, file="~/parentage/kernel_fitting/1340_loci/input/20200602_parentage_matrix_all.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")



#write the matrix
#write.table(parmat_all_full, file="~/parentage/kernel_fitting/1340_loci/parentage_matrix_allyears.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")


#create distance matrix from site centroids
library(geosphere)
### List of source locations
sites_source <- centroids

### List of destination locations
sites_dest <- centroids

###
#dist_mat_km <- rdist.earth(sites_source[,c('lon','lat')], sites_dest[,c('lon','lat')], miles=FALSE, R=6371) #this formula is apparentally less accurate than the Vincenty formula. It gives different distances on the order of about 0.1 decimal places, but that's enough to change the kernel fit estimates SLIGHTLY. I am going to use Vincenty because that's probbly the most accurate and I already did everything that way anyhow. https://stackoverflow.com/questions/38248046/is-the-haversine-formula-or-the-vincentys-formula-better-for-calculating-distan
dist_mat_m <- distm(sites_source[,c('lon','lat')], sites_source[,c('lon','lat')], fun=distVincentyEllipsoid)
dist_mat_old_km <- dist_mat_m*10^-3

#write.table(dist_mat_km, file="~/parentage/kernel_fitting/1340_loci/input/distance_matrix_unsurveyed.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")


#site area DOESN'T NEED TO CHANGE FOR EACH YEAR
area_10per <- area %>%
    arrange(site) %>%
    filter(site %!in% c("near_north_full1", "near_north_full2", "near_north_full3", "near_south_full1", "near_south_full2", "near_south_full3")) %>%
    mutate(kmsq=msq*10^-6) %>%
    select(kmsq)
#DOESN'T NEED TO CHANGE FOR EACH YEAR
#write.table(area_10per, file="~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")

#centroids DOESN'T NEED TO CHANGE FOR EACH YEAR
centroids <- centroids %>%
    arrange(site) %>%
    select(lon, lat)
#write the corrected centroids file
#write.table(centroids, file="~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", row.names=FALSE, quote=FALSE, col.names=FALSE, sep=",")


