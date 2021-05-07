
#iteratively downsample the input to the kernelt fitting parameters, generate a null disribution of kernels 
Packages <- c("dplyr","nleqslv","cubature", "stringr","pracma","data.table", "gridExtra","viridis", "ggsignif", "broom", "ggpubr", "caret","cowplot","ggplot2","fields","bbmle", "dplyr", "tidyr", "lubridate", "RColorBrewer")

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

setwd('/local/home/katrinac/parentage/kernel_fitting/')


load("~/parentage/r_data/total_sampling_across_years.RData")
load("~/parentage/r_data/sampled_area_each_year.RData")
load("~/parentage/r_data/cumulative_prop_hab_sampled_by_site.RData")
load("~/parentage/r_data/anem_obs_db.RData")


#download.file(url = "https://github.com/pinskylab/genomics/blob/master/data/fish-obs.RData?raw=true", destfile = "~/parentage/r_data/fish-obs.RData")
fish_obs <- readRDS("~/parentage/r_data/fish-obs_april2020.RData") 
#download.file(url="https://github.com/pinskylab/genomics/blob/ca6ce13310385e498a8bee54f48511ce3d1557f9/data/fish-obs.RData", destfile = "~/parentage/r_data/old-fish-obs.RData")
old_fish_obs <- readRDS("~/parentage/r_data/fish-obs.RData") 
load("~/parentage/r_data/site_dist_info.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/anem_db.RData?raw=true", destfile = "~/parentage/r_data/anem_db.RData")
load("~/parentage/r_data/anem_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/dives_db.RData?raw=true", destfile = "~/parentage/r_data/dives_db.RData")
load("~/parentage/r_data/dives_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/Data_from_database/fish_db.RData?raw=true", destfile = "~/parentage/r_data/dives_db.RData")
load("~/parentage/r_data/fish_db.RData")
load("~/parentage/r_data/gps_db.RData")
#download.file(url = "https://github.com/pinskylab/Clownfish_persistence/blob/master/Data/anems_tagged.RData", destfile = "~/parentage/r_data/anems_tagged.RData")
load("~/parentage/r_data/anems_tagged.RData")
load("~/parentage/r_data/encounters_list_anem.RData")



source("~/parentage/kernel_fitting/1340_loci/functions/ll_kt_both_bbmle.R")
source("~/parentage/kernel_fitting/1340_loci/functions/ll_kt_both_optim.R")
source("~/parentage/kernel_fitting/1340_loci/functions/GenGausKernInt_sum0.5.R") #integrate_kernel_sum1
source("~/parentage/kernel_fitting/1340_loci/functions/GenGausKernInt_sum1.R") #integrate_kernel_sum0.5
source("~/parentage/kernel_fitting/1340_loci/functions/cdf_solve.R") 
source("~/parentage/kernel_fitting/1340_loci/functions/cdf_solve90.R") 


"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

#read in true prop sampled data
prop_samp <- cumulative_prop_hab_sampled_by_site %>%
    mutate(total_possible_sample_anems = ifelse(site=="Caridad Proper", 4, total_possible_sample_anems) ) %>%
    mutate(total_prop_hab_sampled_anems_tidied= ifelse(site=="Caridad Proper" & total_anems_sampled==4, 1, total_prop_hab_sampled_anems_tidied) ) %>%
    mutate(total_possible_sample_anems = ifelse(site=="Sitio Lonas", total_anems_sampled, total_possible_sample_anems) ) %>%
    mutate(total_prop_hab_sampled_anems_tidied= ifelse(site=="Sitio Lonas", 1, total_prop_hab_sampled_anems_tidied) )# %>%
    #mutate(total_prop_hab_sampled_anems_tidied= ifelse(is.nan(total_prop_hab_sampled_anems_tidied), 0.1, total_prop_hab_sampled_anems_tidied ))#make NANs really small instead
#correct space in North/South Magbangon names
prop_samp$site <- gsub(". ", ".", prop_samp$site, fixed=TRUE)
#read in information necessary for kernel fitting
centroids <- read.csv("~/parentage/kernel_fitting/1340_loci/input/site_centroids.csv", header=TRUE)
sampled_sites <- prop_samp %>% select(site) %>% distinct(site)
sampled_sites$site <- gsub(". ", ".", sampled_sites$site, fixed=TRUE)
site_widths <- read.table("/local/home/katrinac/parentage/text_file/site_widths.txt", header=T, sep=",")

#index sites to get pop numbers that replace site names
all_sites <- centroids %>%
    select(site) %>%
    arrange(site)
nrow(all_sites) #should be 35x1
all_sites$index <- seq(from=1, to=35, by=1)

#from Allison, just putting all the meta data together (Constants_database_common_functions.R)
##### Match up other relevant info (site, date, fish_indiv, etc.) to fish in the clownfish table
# Pull out year and month into a separate column in dives_db
dives_db_processed <- dives_db %>%
  mutate(year = as.integer(substring(date,1,4))) %>%
  mutate(month = as.integer(substring(date,6,7))) %>%
  mutate(dive_date = date(date))

# Pull all APCL caught or otherwise in the clownfish table
allfish_fish <- fish_db %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, anem_table_id, recap, tag_id, color, sex, size, fish_obs_time) %>%
  filter(fish_spp == 'APCL') %>%
  mutate(size = as.numeric(size))  # make the size numeric (rather than chr) so can do means and such

# and their corresponding anemones
allfish_anems <- anem_db %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_obs_time) #%>%
  #filter(anem_table_id %in% allfish_fish$anem_table_id)

# and the corresponding dive info
allfish_dives <- dives_db_processed %>%
  select(dive_table_id, dive_type, date, year, month, site, gps, dive_type) %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id) 

# join together
allfish_caught1 <- left_join(allfish_anems, allfish_fish, by="anem_table_id")
allfish_caught <- left_join(allfish_caught1, allfish_dives, by="dive_table_id")

# add in the gen_ids and fish_indiv
allfish_caught <- left_join(allfish_caught, (fish_obs %>% select(fish_table_id, gen_id, fish_indiv)), by = "fish_table_id")

fish_meta <- allfish_caught %>%
    select(fish_indiv, size, color, sex, gen_id, sample_id, site, date, anem_table_id, anem_obs_time,anem_obs, anem_id, gps, dive_type) %>%
    mutate(date=ymd(date)) %>%
    mutate(year=year(date)) #%>%
    #filter(dive_type != "R") %>% #don't count recapture dives
    #distinct(anem_obs, fish_indiv, year, .keep_all = T) #to avoid double counting anemones that were visited twice in one year, keep only one observationz


#read in fish sample data
#parentage results
par_res <- read.csv(file="~/parentage/colony2/20200605_1340loci/results/parentage_results_allyears.csv", header= T) %>%
    distinct(offs_fish_indiv, .keep_all = T) %>% #just need the sites for parent and offspring, will join to fish metadata with fish_indiv
    select(offs_fish_indiv, par_fish_indiv, offs_site, par_site, year) %>% 
    mutate(offs_fish_indiv=as.character(offs_fish_indiv)) %>%
    rename(year_match="year") #change this year to be the year of the parentage match

#genotyped potential offspring
N_gen_offs <- read.table(file="~/parentage/colony2/20190523_1340loci/input/all_offspring_corrected.txt", header=T, stringsAsFactors = F) %>%
    mutate(input="offspring") %>%
    mutate(fish_indiv=as.character(fish_indiv)) %>%
    select(fish_indiv, input, year)

SampledRecruitsPar <- left_join(N_gen_offs, par_res, by=c(fish_indiv="offs_fish_indiv")) 

AllFishObsWithPar <- right_join(fish_meta, SampledRecruitsPar, by=c("fish_indiv", "year")) %>%
    ungroup() %>%
    filter(sample_id != "APCL13_626" & sample_id != "APCL16_774") %>%
    mutate(matched_offs=ifelse(!is.na(year_match), "Y", "N")) %>%
    distinct(fish_indiv, .keep_all = T)
    #mutate(year_match = ifelse(sample_id=="APCL13_626", NA, year_match)) %>% #one offspring was sampled twice in the same year, but on two different anemones. Keep the observation because both anemones were sampled, but take away the match year to make sure it doesn't end up double counted
    #mutate(input = ifelse(sample_id=="APCL13_626", "NA", input)) #fish_indiv == 2400 

#correct site names
AllFishObsWithPar$site <- gsub(". ", ".", AllFishObsWithPar$site, fixed=TRUE)
AllFishObsWithPar$offs_site <- gsub(". ", ".", AllFishObsWithPar$offs_site, fixed=TRUE)
AllFishObsWithPar$par_site <- gsub(". ", ".", AllFishObsWithPar$par_site, fixed=TRUE)

nrow(AllFishObsWithPar %>% filter(!is.na(year_match))) #should be 71
nrow(AllFishObsWithPar %>% filter(input=="offspring")) #should be 791
nrow(AllFishObsWithPar %>% filter(matched_offs == "Y")) #should be 71

#for future use
#write.csv(AllFishObsWithPar, file="~/parentage/kernel_fitting/1340_loci/final_results/tables/AllFishObsWithPar.csv", row.names=FALSE)

#make an empty data frame to hold annual simulated values
col <- c("year", "k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
SimulatedKernels <- as.data.frame(matrix(nrow=0, ncol=7), stringsAsFactors = FALSE)
colnames(SimulatedKernels) <- col


#don't print warnings for this loop, they are only for setting row names in a tibble which is depracated. But it works. 
options(warn=-1)

#set a progress bar to monitor
pb <- txtProgressBar(min = 0, max = 10000, style = 3)#7 years in each interation

StartTime <- Sys.time()

for(n in 1:10000){
#get a one column data frame of the years associated
Years <- AllFishObsWithPar %>%
    select(year)

#shuffle the years of this data frame
ShuffleYears <- data.frame(Years[sample(1:nrow(Years), replace=FALSE),])
colnames(ShuffleYears) <- "year"

#remove true years
AllFishObsWithParBeta <- AllFishObsWithPar %>%
    select(-year)

#join back in shuffled years
AllFishObsWithParSim <- bind_cols(AllFishObsWithParBeta, ShuffleYears) 

#EVERYTHING THAT FOLLOWS formats data for kernel fitting, fits the kernels, and stores the results into a data frame called SimulatedKernels
#now break up into parentage per year

par12 <- AllFishObsWithParSim %>%
    filter(year=="2012" & matched_offs=="Y") %>%
    distinct(offs_fish_indiv, .keep_all = T) %>%
    select(year, offs_site, par_site)

par13 <- AllFishObsWithParSim %>%
    filter(year=="2013" & matched_offs=="Y")  %>%
    distinct(offs_fish_indiv, .keep_all = T) %>%
    select(year, offs_site, par_site) 

par14 <- AllFishObsWithParSim %>%
    filter(year=="2014" & matched_offs=="Y")  %>%
    distinct(offs_fish_indiv, .keep_all = T) %>%
    select(year, offs_site, par_site) 

par15 <- AllFishObsWithParSim %>%
    filter(year=="2015" & matched_offs=="Y")  %>%
    distinct(offs_fish_indiv, .keep_all = T) %>%
    select(year, offs_site, par_site) 

par16 <- AllFishObsWithParSim %>%
    filter(year=="2016"& matched_offs=="Y")  %>%
    distinct(offs_fish_indiv, .keep_all = T) %>%
    select(year, offs_site, par_site) 

par17 <- AllFishObsWithParSim %>%
    filter(year=="2017" & matched_offs=="Y")  %>%
    distinct(offs_fish_indiv, .keep_all = T) %>%
    select(year, offs_site, par_site)

par18 <- AllFishObsWithParSim %>%
    filter(year=="2018"& matched_offs=="Y")  %>%
    distinct(offs_fish_indiv, .keep_all = T) %>%
    select(year, offs_site, par_site) 

total_par12 <- par12 %>%
    group_by(year, offs_site, par_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() %>%
    select(-year)

total_par13 <- par13 %>%
    group_by(year, offs_site, par_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() %>%
    select(-year)

total_par14 <- par14 %>%
    group_by(year, offs_site, par_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() %>%
    select(-year)

total_par15 <- par15 %>%
    group_by(year, offs_site, par_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() %>%
    select(-year)

total_par16 <- par16 %>%
    group_by(year, offs_site, par_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() %>%
    select(-year)

total_par17 <- par17 %>%
    group_by(year, offs_site, par_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() %>%
    select(-year)

total_par18 <- par18 %>%
    group_by(year, offs_site, par_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() %>%
    select(-year)

#prop_vec <- seq(0.1, 1, 0.05), sample(prop_vec, 1) instead of static average

#add in sites that were sampled but there was no match
for(i in 1:nrow(prop_samp)){
    
    
    if(is.nan(prop_samp$total_prop_hab_sampled_anems_tidied[i])){prop_samp$total_prop_hab_sampled_anems_tidied[i] <- 0.1} 
    if(is.infinite(prop_samp$total_prop_hab_sampled_anems_tidied[i])){prop_samp$total_prop_hab_sampled_anems_tidied[i] <- 1}
    ifelse(prop_samp$total_prop_hab_sampled_anems_tidied[i] == 0, 0.48 , prop_samp$total_prop_hab_sampled_anems_tidied[i])
    
}


prop_samp12 <- prop_samp %>%
    filter(end_year=="2012") %>%
    select(site, total_prop_hab_sampled_anems_tidied)

#use nrow(sites) to get the dimensions that I need to trim the matrix to

sites12 <- prop_samp12 %>%
    select(site)
sites12_2 <- suppressWarnings(semi_join(all_sites, sites12, by="site")) %>%
    select(index)
sites12t <- t(sites12_2)



sites12_beta <- sites12 %>%
    #select(-pop) %>%
    rename(par_site="site")

sites12_beta$offs_site <- sites12_beta$par_site

allsites_parentage12 <- full_join(sites12_beta, total_par12, by=c("par_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site)


prop_samp13 <- prop_samp %>%
    filter(end_year=="2013" ) %>%
    select(site, total_prop_hab_sampled_anems_tidied)


##generate list of site indices
sites13 <- prop_samp13 %>%
    select(site)
sites13_2 <- suppressWarnings(semi_join(all_sites, sites13, by="site")) %>%
    select(index)
sites13t <- t(sites13_2)

sites13_beta <- sites13 %>%
    #select(-pop) %>%
    rename(par_site="site")

sites13_beta$offs_site <- sites13_beta$par_site

allsites_parentage13 <- full_join(sites13_beta, total_par13, by=c("par_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site)


prop_samp14 <- prop_samp %>%
    filter(end_year=="2014" ) %>%
    select(site, total_prop_hab_sampled_anems_tidied)


##generate list of site indices
sites14 <- prop_samp14 %>%
    select(site)
sites14_2 <- suppressWarnings(semi_join(all_sites, sites14, by="site")) %>%
    select(index)
sites14t <- t(sites14_2)

sites14_beta <- sites14 %>%
    #select(-pop) %>%
    rename(par_site="site")

sites14_beta$offs_site <- sites14_beta$par_site

allsites_parentage14 <- full_join(sites14_beta, total_par14, by=c("par_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site)

prop_samp15 <- prop_samp %>%
    filter(end_year=="2015" ) %>%
    select(site, total_prop_hab_sampled_anems_tidied)


##generate list of site indices
sites15 <- prop_samp15 %>%
    select(site)
sites15_2 <- suppressWarnings(semi_join(all_sites, sites15, by="site")) %>%
    select(index)
sites15t <- t(sites15_2)

sites15_beta <- sites15 %>%
    #select(-pop) %>%
    rename(par_site="site")

sites15_beta$offs_site <- sites15_beta$par_site

allsites_parentage15 <- full_join(sites15_beta, total_par15, by=c("par_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site)

prop_samp16 <- prop_samp %>%
    filter(end_year=="2016" ) %>%
    select(site, total_prop_hab_sampled_anems_tidied)


##generate list of site indices
sites16 <- prop_samp16 %>%
    select(site)
sites16_2 <- suppressWarnings(semi_join(all_sites, sites16, by="site")) %>%
    select(index)
sites16t <- t(sites16_2)

sites16_beta <- sites16 %>%
    #select(-pop) %>%
    rename(par_site="site")

sites16_beta$offs_site <- sites16_beta$par_site

allsites_parentage16 <- full_join(sites16_beta, total_par16, by=c("par_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site)

prop_samp17 <- prop_samp %>%
    filter(end_year=="2017") %>%
    select(site, total_prop_hab_sampled_anems_tidied)


##generate list of site indices
sites17 <- prop_samp17 %>%
    select(site)
sites17_2 <- suppressWarnings(semi_join(all_sites, sites17, by="site")) %>%
    select(index)
sites17t <- t(sites17_2)

sites17_beta <- sites17 %>%
    #select(-pop) %>%
    rename(par_site="site")

sites17_beta$offs_site <- sites17_beta$par_site

allsites_parentage17 <- full_join(sites17_beta, total_par17, by=c("par_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site)

prop_samp18 <- prop_samp %>%
    filter(end_year=="2018" ) %>%
    select(site, total_prop_hab_sampled_anems_tidied)


##generate list of site indices
sites18 <- prop_samp18 %>%
    select(site)
sites18_2 <- suppressWarnings(semi_join(all_sites, sites18, by="site")) %>%
    select(index)
sites18t <- t(sites18_2)

sites18_beta <- sites18 %>%
    #select(-pop) %>%
    rename(par_site="site")

sites18_beta$offs_site <- sites18_beta$par_site

allsites_parentage18 <- full_join(sites18_beta, total_par18, by=c("par_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site)


#turn into a matrix of the correct format
parmat12 <- allsites_parentage12 %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

##change NAs to 0
parmat12[is.na(parmat12)] <- 0
rownames(parmat12) <- suppressWarnings(parmat12$par_site)
parmat12$par_site <- NULL

parmat12 <- (parmat12[1:nrow(sites12_beta), 1:nrow(sites12_beta)])


parmat13 <- allsites_parentage13 %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

##change NAs to 0
parmat13[is.na(parmat13)] <- 0
rownames(parmat13) <- suppressWarnings(parmat13$par_site)
parmat13$par_site <- NULL

parmat13 <- (parmat13[1:nrow(sites13_beta), 1:nrow(sites13_beta)])

parmat14 <- allsites_parentage14 %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

##change NAs to 0
parmat14[is.na(parmat14)] <- 0
rownames(parmat14) <- suppressWarnings(parmat14$par_site)
parmat14$par_site <- NULL

parmat14 <- (parmat14[1:nrow(sites14_beta), 1:nrow(sites14_beta)])

parmat15 <- allsites_parentage15 %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

##change NAs to 0
parmat15[is.na(parmat15)] <- 0
rownames(parmat15) <- suppressWarnings(parmat15$par_site)
parmat15$par_site <- NULL

parmat15 <- (parmat15[1:nrow(sites15_beta), 1:nrow(sites15_beta)])

parmat16 <- allsites_parentage16 %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

##change NAs to 0
parmat16[is.na(parmat16)] <- 0
rownames(parmat16) <- suppressWarnings(parmat16$par_site)
parmat16$par_site <- NULL

parmat16 <- (parmat16[1:nrow(sites16_beta), 1:nrow(sites16_beta)])

parmat17 <- allsites_parentage17 %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

##change NAs to 0
parmat17[is.na(parmat17)] <- 0
rownames(parmat17) <- suppressWarnings(parmat17$par_site)
parmat17$par_site <- NULL

parmat17 <- (parmat17[1:nrow(sites17_beta), 1:nrow(sites17_beta)])


#create the matrix with rows (parent pop) and columns (offs pop)
parmat18 <- allsites_parentage18 %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

##change NAs to 0
parmat18[is.na(parmat18)] <- 0
rownames(parmat18) <- suppressWarnings(parmat18$par_site)
parmat18$par_site <- NULL

parmat18 <- (parmat18[1:nrow(sites18_beta), 1:nrow(sites18_beta)])


#add in the unassigned juveniles
N_gen_offs_all <- AllFishObsWithParSim %>% 
    filter(input=="offspring") %>%
    group_by(year, site) %>%
    summarise(n_offs_gen=n())


N_gen_offs2 <- N_gen_offs_all %>%
    filter(year==2012 ) %>%
    group_by(site, year) %>% #add group by year if doing annual
    summarise(sampled_fish=sum(n_offs_gen)) 


#make dataframe with column for unassigned juveniles in each column (offspring site)
offs_matched_site <- allsites_parentage12 %>%
    group_by(offs_site) %>%
    summarise(n_offs=sum(n_matches, na.rm=T))

unassigned_beta <- semi_join(offs_matched_site, sites12, by=c(offs_site="site")) %>%
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
colnames(row2addT) <- colnames(parmat12)

parmat12 <- ungroup(parmat12)
parmat_12_full <- bind_rows(parmat12, row2addT)
parmat_12_full[is.na(parmat_12_full)] <- 0


N_gen_offs2 <- N_gen_offs_all %>%
    filter(year==2013 ) %>%
    group_by(site, year) %>% #add group by year if doing annual
    summarise(sampled_fish=sum(n_offs_gen)) 


#make dataframe with column for unassigned juveniles in each column (offspring site)
offs_matched_site <- allsites_parentage13 %>%
    group_by(offs_site) %>%
    summarise(n_offs=sum(n_matches, na.rm=T))

unassigned_beta <- semi_join(offs_matched_site, sites13, by=c(offs_site="site")) %>%
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
colnames(row2addT) <- colnames(parmat13)

parmat13 <- ungroup(parmat13)
parmat_13_full <- bind_rows(parmat13, row2addT)
parmat_13_full[is.na(parmat_13_full)] <- 0

N_gen_offs2 <- N_gen_offs_all %>%
    filter(year==2014 ) %>%
    group_by(site, year) %>% #add group by year if doing annual
    summarise(sampled_fish=sum(n_offs_gen)) 


#make dataframe with column for unassigned juveniles in each column (offspring site)
offs_matched_site <- allsites_parentage14 %>%
    group_by(offs_site) %>%
    summarise(n_offs=sum(n_matches, na.rm=T))

unassigned_beta <- semi_join(offs_matched_site, sites14, by=c(offs_site="site")) %>%
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
colnames(row2addT) <- colnames(parmat14)

parmat14 <- ungroup(parmat14)
parmat_14_full <- bind_rows(parmat14, row2addT)
parmat_14_full[is.na(parmat_14_full)] <- 0


N_gen_offs2 <- N_gen_offs_all %>%
    filter(year==2015 ) %>%
    group_by(site, year) %>% #add group by year if doing annual
    summarise(sampled_fish=sum(n_offs_gen)) 


#make dataframe with column for unassigned juveniles in each column (offspring site)
offs_matched_site <- allsites_parentage15 %>%
    group_by(offs_site) %>%
    summarise(n_offs=sum(n_matches, na.rm=T))

unassigned_beta <- semi_join(offs_matched_site, sites15, by=c(offs_site="site")) %>%
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
colnames(row2addT) <- colnames(parmat15)

parmat15 <- ungroup(parmat15)
parmat_15_full <- bind_rows(parmat15, row2addT)
parmat_15_full[is.na(parmat_15_full)] <- 0

N_gen_offs2 <- N_gen_offs_all %>%
    filter(year==2016 ) %>%
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


N_gen_offs2 <- N_gen_offs_all %>%
    filter(year==2017 ) %>%
    group_by(site, year) %>% #add group by year if doing annual
    summarise(sampled_fish=sum(n_offs_gen)) 


#make dataframe with column for unassigned juveniles in each column (offspring site)
offs_matched_site <- allsites_parentage17 %>%
    group_by(offs_site) %>%
    summarise(n_offs=sum(n_matches, na.rm=T))

unassigned_beta <- semi_join(offs_matched_site, sites17, by=c(offs_site="site")) %>%
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
colnames(row2addT) <- colnames(parmat17)

parmat17 <- ungroup(parmat17)
parmat_17_full <- bind_rows(parmat17, row2addT)
parmat_17_full[is.na(parmat_17_full)] <- 0


N_gen_offs2 <- N_gen_offs_all %>%
    filter(year==2018 ) %>%
    group_by(site, year) %>% #add group by year if doing annual
    summarise(sampled_fish=sum(n_offs_gen)) 


#make dataframe with column for unassigned juveniles in each column (offspring site)
offs_matched_site <- allsites_parentage18 %>%
    group_by(offs_site) %>%
    summarise(n_offs=sum(n_matches, na.rm=T))

unassigned_beta <- semi_join(offs_matched_site, sites18, by=c(offs_site="site")) %>%
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
colnames(row2addT) <- colnames(parmat18)

parmat18 <- ungroup(parmat18)
parmat_18_full <- bind_rows(parmat18, row2addT)
parmat_18_full[is.na(parmat_18_full)] <- 0

#end file prep

#fit kernels 
Assignments <- parmat_12_full
Adult_sample_proportions <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/prop_samp12.csv", header=FALSE))
Sampled_reefs <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_index12.csv", header=FALSE))
Distances <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/distance_matrix_unsurveyed.csv", header=FALSE))
Reef_sizes <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", header=FALSE))
Centroids <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", header=T))

a=-10
b=10

x <- list(Distances=Distances, Assignments=Assignments, Sampled_reefs=Sampled_reefs, Reef_sizes=Reef_sizes, Adult_sample_proportions=Adult_sample_proportions) #put inputs into a list because that's the bbmle format

fit_both12 <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=1, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=x, control=list(maxit=500)))
k12_b <- coef(fit_both12)[1]
theta12_b <- coef(fit_both12)[2]
MDD12_b <- cubintegrate(integrate_kernel_sum1, lower = 0, upper = Inf, k=k12_b, theta=theta12_b, method = "pcubature")$integral
theta_eval <- theta12_b
k_eval <- k12_b
MedianDispDist12 <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2)
Dist90Retained12 <- round(nleqslv(x = 7, fn = cdf_solve90)$x, 2) 

col <- c("year", "k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
BootstrapKernels12 <- as.data.frame(matrix(nrow=1, ncol=7), stringsAsFactors = FALSE)
colnames(BootstrapKernels12) <- col

BootstrapKernels12$year <- "2012"
BootstrapKernels12$k <- k12_b
BootstrapKernels12$theta <- theta12_b
BootstrapKernels12$MDD <- MDD12_b
BootstrapKernels12$MedianDispDist <- MedianDispDist12
BootstrapKernels12$Dist90Retained <- Dist90Retained12
BootstrapKernels12$iteration <- as.character(n)



Assignments <- parmat_13_full
Adult_sample_proportions <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/prop_samp13.csv", header=FALSE))
Sampled_reefs <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_index13.csv", header=FALSE))
Distances <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/distance_matrix_unsurveyed.csv", header=FALSE))
Reef_sizes <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", header=FALSE))
Centroids <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", header=T))

a=-10
b=10

x <- list(Distances=Distances, Assignments=Assignments, Sampled_reefs=Sampled_reefs, Reef_sizes=Reef_sizes, Adult_sample_proportions=Adult_sample_proportions) #put inputs into a list because that's the bbmle format

fit_both13 <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=1, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=x, control=list(maxit=500)))
k13_b <- coef(fit_both13)[1]
theta13_b <- coef(fit_both13)[2]
MDD13_b <- cubintegrate(integrate_kernel_sum1, lower = 0, upper = Inf, k=k13_b, theta=theta13_b, method = "pcubature")$integral
theta_eval <- theta13_b
k_eval <- k13_b
MedianDispDist13 <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2)
Dist90Retained13 <- round(nleqslv(x = 7, fn = cdf_solve90)$x, 2) 

col <- c("year", "k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
BootstrapKernels13 <- as.data.frame(matrix(nrow=1, ncol=7), stringsAsFactors = FALSE)
colnames(BootstrapKernels13) <- col

BootstrapKernels13$year <- "2013"
BootstrapKernels13$k <- k13_b
BootstrapKernels13$theta <- theta13_b
BootstrapKernels13$MDD <- MDD13_b
BootstrapKernels13$MedianDispDist <- MedianDispDist13
BootstrapKernels13$Dist90Retained <- Dist90Retained13
BootstrapKernels13$iteration <- as.character(n)




Assignments <- parmat_14_full
Adult_sample_proportions <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/prop_samp14.csv", header=FALSE))
Sampled_reefs <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_index14.csv", header=FALSE))
Distances <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/distance_matrix_unsurveyed.csv", header=FALSE))
Reef_sizes <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", header=FALSE))
Centroids <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", header=T))

a=-10
b=10

x <- list(Distances=Distances, Assignments=Assignments, Sampled_reefs=Sampled_reefs, Reef_sizes=Reef_sizes, Adult_sample_proportions=Adult_sample_proportions) #put inputs into a list because that's the bbmle format

fit_both14 <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=1, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=x, control=list(maxit=500)))
k14_b <- coef(fit_both14)[1]
theta14_b <- coef(fit_both14)[2]
MDD14_b <- cubintegrate(integrate_kernel_sum1, lower = 0, upper = Inf, k=k14_b, theta=theta14_b, method = "pcubature")$integral
theta_eval <- theta14_b
k_eval <- k14_b
MedianDispDist14 <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2)
Dist90Retained14 <- round(nleqslv(x = 7, fn = cdf_solve90)$x, 2) 

col <- c("year", "k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
BootstrapKernels14 <- as.data.frame(matrix(nrow=1, ncol=7), stringsAsFactors = FALSE)
colnames(BootstrapKernels14) <- col

BootstrapKernels14$year <- "2014"
BootstrapKernels14$k <- k14_b
BootstrapKernels14$theta <- theta14_b
BootstrapKernels14$MDD <- MDD14_b
BootstrapKernels14$MedianDispDist <- MedianDispDist14
BootstrapKernels14$Dist90Retained <- Dist90Retained14
BootstrapKernels14$iteration <- as.character(n)
    

Assignments <- parmat_15_full
Adult_sample_proportions <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/prop_samp15.csv", header=FALSE))
Sampled_reefs <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_index15.csv", header=FALSE))
Distances <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/distance_matrix_unsurveyed.csv", header=FALSE))
Reef_sizes <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", header=FALSE))
Centroids <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", header=T))

a=-10
b=10

x <- list(Distances=Distances, Assignments=Assignments, Sampled_reefs=Sampled_reefs, Reef_sizes=Reef_sizes, Adult_sample_proportions=Adult_sample_proportions) #put inputs into a list because that's the bbmle format

fit_both15 <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=1, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=x, control=list(maxit=500)))
k15_b <- coef(fit_both15)[1]
theta15_b <- coef(fit_both15)[2]
MDD15_b <- cubintegrate(integrate_kernel_sum1, lower = 0, upper = Inf, k=k15_b, theta=theta15_b, method = "pcubature")$integral
theta_eval <- theta15_b
k_eval <- k15_b
MedianDispDist15 <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2)
Dist90Retained15 <- round(nleqslv(x = 7, fn = cdf_solve90)$x, 2) 

col <- c("year", "k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
BootstrapKernels15 <- as.data.frame(matrix(nrow=1, ncol=7), stringsAsFactors = FALSE)
colnames(BootstrapKernels15) <- col

BootstrapKernels15$year <- "2015"
BootstrapKernels15$k <- k15_b
BootstrapKernels15$theta <- theta15_b
BootstrapKernels15$MDD <- MDD15_b
BootstrapKernels15$MedianDispDist <- MedianDispDist15
BootstrapKernels15$Dist90Retained <- Dist90Retained15
BootstrapKernels15$iteration <- as.character(n)


Assignments <- parmat_16_full
Adult_sample_proportions <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/prop_samp16.csv", header=FALSE))
Sampled_reefs <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_index16.csv", header=FALSE))
Distances <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/distance_matrix_unsurveyed.csv", header=FALSE))
Reef_sizes <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", header=FALSE))
Centroids <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", header=T))

a=-10
b=10

x <- list(Distances=Distances, Assignments=Assignments, Sampled_reefs=Sampled_reefs, Reef_sizes=Reef_sizes, Adult_sample_proportions=Adult_sample_proportions) #put inputs into a list because that's the bbmle format

fit_both16 <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=1, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=x, control=list(maxit=500)))
k16_b <- coef(fit_both16)[1]
theta16_b <- coef(fit_both16)[2]
theta_eval <- theta16_b
k_eval <- k16_b
MedianDispDist16 <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2)
Dist90Retained16 <- round(nleqslv(x = 7, fn = cdf_solve90)$x, 2) 

col <- c("year", "k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
BootstrapKernels16 <- as.data.frame(matrix(nrow=1, ncol=7), stringsAsFactors = FALSE)
colnames(BootstrapKernels16) <- col

BootstrapKernels16$year <- "2016"
BootstrapKernels16$k <- k16_b
BootstrapKernels16$theta <- theta16_b
BootstrapKernels16$MDD <- MDD16_b
BootstrapKernels16$MedianDispDist <- MedianDispDist16
BootstrapKernels16$Dist90Retained <- Dist90Retained16
BootstrapKernels16$iteration <- as.character(n)


Assignments <- parmat_17_full
Adult_sample_proportions <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/prop_samp17.csv", header=FALSE))
Sampled_reefs <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_index17.csv", header=FALSE))
Distances <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/distance_matrix_unsurveyed.csv", header=FALSE))
Reef_sizes <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", header=FALSE))
Centroids <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", header=T))

a=-10
b=10

x <- list(Distances=Distances, Assignments=Assignments, Sampled_reefs=Sampled_reefs, Reef_sizes=Reef_sizes, Adult_sample_proportions=Adult_sample_proportions) #put inputs into a list because that's the bbmle format

fit_both17 <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=1, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=x, control=list(maxit=500)))
k17_b <- coef(fit_both17)[1]
theta17_b <- coef(fit_both17)[2]
MDD17_b <- cubintegrate(integrate_kernel_sum1, lower = 0, upper = Inf, k=k17_b, theta=theta17_b, method = "pcubature")$integral
theta_eval <- theta17_b
k_eval <- k17_b
MedianDispDist17 <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2)
Dist90Retained17 <- round(nleqslv(x = 7, fn = cdf_solve90)$x, 2) 

col <- c("year", "k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
BootstrapKernels17 <- as.data.frame(matrix(nrow=1, ncol=7), stringsAsFactors = FALSE)
colnames(BootstrapKernels17) <- col

BootstrapKernels17$year <- "2017"
BootstrapKernels17$k <- k17_b
BootstrapKernels17$theta <- theta17_b
BootstrapKernels17$MDD <- MDD17_b
BootstrapKernels17$MedianDispDist <- MedianDispDist17
BootstrapKernels17$Dist90Retained <- Dist90Retained17
BootstrapKernels17$iteration <- as.character(n)
    

Assignments <- parmat_18_full
Adult_sample_proportions <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/prop_samp18.csv", header=FALSE))
Sampled_reefs <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_index18.csv", header=FALSE))
Distances <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/distance_matrix_unsurveyed.csv", header=FALSE))
Reef_sizes <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", header=FALSE))
Centroids <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", header=T))

a=-10
b=10

x <- list(Distances=Distances, Assignments=Assignments, Sampled_reefs=Sampled_reefs, Reef_sizes=Reef_sizes, Adult_sample_proportions=Adult_sample_proportions) #put inputs into a list because that's the bbmle format

fit_both18 <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=1, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=x, control=list(maxit=500)))
k18_b <- coef(fit_both18)[1]
theta18_b <- coef(fit_both18)[2]
MDD18_b <- cubintegrate(integrate_kernel_sum1, lower = 0, upper = Inf, k=k18_b, theta=theta18_b, method = "pcubature")$integral
theta_eval <- theta18_b
k_eval <- k18_b
MedianDispDist18 <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2)
Dist90Retained18 <- round(nleqslv(x = 7, fn = cdf_solve90)$x, 2) 
    
    
col <- c("year", "k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
BootstrapKernels18 <- as.data.frame(matrix(nrow=1, ncol=7), stringsAsFactors = FALSE)
colnames(BootstrapKernels18) <- col

BootstrapKernels18$year <- "2018"
BootstrapKernels18$k <- k18_b
BootstrapKernels18$theta <- theta18_b
BootstrapKernels18$MDD <- MDD18_b
BootstrapKernels18$MedianDispDist <- MedianDispDist18
BootstrapKernels18$Dist90Retained <- Dist90Retained18
BootstrapKernels18$iteration <- as.character(n)


#bind together all simulated kernel fits
#they aren't bootstraps technically but whatever it's a name
SimulatedKernelsBeta <- bind_rows(BootstrapKernels12, BootstrapKernels13, BootstrapKernels14, BootstrapKernels15, BootstrapKernels16, BootstrapKernels17, BootstrapKernels18)

SimulatedKernels <- bind_rows(SimulatedKernels, SimulatedKernelsBeta)

setTxtProgressBar(pb, n)


}
close(pb)
EndTime <- Sys.time()
EndTime-StartTime
options(warn=0) #turn warnings back on
write.csv(SimulatedKernels, file="~/parentage/kernel_fitting/1340_loci/final_results/simulations/SimulatedAnnual.csv", row.names=FALSE, quote=FALSE)


nrow(SimulatedKernels %>% filter(theta==5 | theta==0.15))/nrow(SimulatedKernels) 

SimulatedKernels <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/simulations/SimulatedAnnual.csv", header=T) %>% #load simulations
    mutate(iteration=as.character(iteration))

kernels <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/tables/kernel_fitting_summary.csv", header=T, stringsAsFactors = F) %>%
        filter(Year!="2012-18" )  #load empirical results, test removing 2016 to see if that's exaggerating CV


#compare simulations and empirical
SimulatedKernelsVar <- SimulatedKernels %>% #Med
    group_by(iteration) %>%
    mutate(cvtheta=sd(theta)/mean(theta, na.rm=T)) %>%
    mutate(sdk=sd(k)) %>%
    mutate(MDD=ifelse(MDD <0, 0,MDD)) %>% #if MDD is less than zero because of negative k and low theta (near 0.10), replace with 0
    mutate(cvMDD=sd(MDD)/mean(MDD, na.rm=T)) %>%
    mutate(cvMed=sd(MedianDispDist)/mean(MedianDispDist, na.rm=T)) %>%
    mutate(cv90=sd(Dist90Retained)/mean(Dist90Retained, na.rm=T)) %>%
    group_by(iteration) %>%
    distinct(iteration, .keep_all = T) %>%
    select(-year) %>%
    ungroup() %>%
    sample_n(1000, replace=F) #get back to 1000 simulations

RealKernelVar <- kernels %>%
    mutate(cvMDD=sd(MeanDispDist)/mean(MeanDispDist, na.rm=T))%>%
    mutate(cvMed=sd(MedianDispDist)/mean(MedianDispDist, na.rm=T))%>%
    mutate(cvtheta=sd(best_theta)/mean(best_theta, na.rm=T)) %>%
    mutate(sdk=sd(best_k)) %>%
    mutate(cv90=sd(Dist90Retained)/mean(Dist90Retained, na.rm=T)) %>%
    distinct(cvMDD, cvtheta, sdk, cvMed, cv90)



nrow(SimulatedKernelsVar)

RealKernelVar

(nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(sdk > RealKernelVar$sdk)))+1

#what percentile are our observations in?

((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(sdk > RealKernelVar$sdk))))/1001
((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(cvtheta > RealKernelVar$cvtheta))))/1001
((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(cvMDD > RealKernelVar$cvMDD))))/1001
((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(cvMed > RealKernelVar$cvMed))))/1001
((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(cv90 > RealKernelVar$cv90))))/1001


pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/AnnualMDDCV.pdf")
hist(SimulatedKernelsVar$cvMDD, col="grey", breaks=seq(0, 1.5, 0.01), main=NA, xlab= "Mean dispersal distance CV")
abline(v = RealKernelVar$cvMDD, col="red", lwd=3, lty=2)
dev.off()

pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/AnnualMedCV.pdf")
hist(SimulatedKernelsVar$cvMed,col="grey", breaks=seq(0, 1.5, 0.01), main=NA, xlab= "Median dispersal distance CV")
abline(v = RealKernelVar$cvMed, col="red", lwd=3, lty=2)
dev.off()

pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/Annual90CV.pdf")
hist(SimulatedKernelsVar$cv90,col="grey", breaks=seq(0, 1.5, 0.01), main=NA, xlab= "0.90 dispersal distance CV")
abline(v = RealKernelVar$cv90, col="red", lwd=3, lty=2)
dev.off()

pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/AnnualThetaCV.pdf")
hist(SimulatedKernelsVar$cvtheta, col="grey", breaks=seq(0, 1.5, 0.01), main=NA, xlab="theta CV")
abline(v = RealKernelVar$cvtheta, col="red", lwd=3, lty=2)
dev.off()

pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/AnnualKSD.pdf")
hist(SimulatedKernelsVar$sdk, col="grey",main=NA, breaks=seq(0, 3, 0.01), xlab= "k SD")
abline(v = RealKernelVar$sdk, col="red", lwd=3, lty=2)
dev.off()



SWM_recruits <- AllFishObsWithPar %>%
    filter(size >=4.5 & size < 6) %>%
    mutate(season="SWM")
NEM_recruits <- AllFishObsWithPar %>%
    filter(size <= 3.5) %>%
    mutate(season="NEM")

nrow(SWM_recruits) #should be 428
nrow(NEM_recruits) #should be 132

nrow(SWM_recruits %>% filter(matched_offs=="Y")) #should be 35
nrow(NEM_recruits %>% filter(matched_offs=="Y")) #should be 11
AllFishObsWithParSeason <- bind_rows(SWM_recruits, NEM_recruits)
nrow(AllFishObsWithParSeason) #should be 560

#simulate seasonal kernels

#don't print warnings for this loop, they are just about row names for tibbles being deprecated
options(warn=-1)

pb <- txtProgressBar(min = 0, max = 10000, style = 3)

StartTime <- Sys.time()

#make an empty data frame to hold simulated values
col <- c("season", "k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
SimulatedKernelsSeas <- as.data.frame(matrix(nrow=0, ncol=7), stringsAsFactors = FALSE)
colnames(SimulatedKernelsSeas) <- col

for(n in 1:10000){

    
#shuffle seasons only
shuff_seasons <- AllFishObsWithParSeason %>%
    select(season) %>%
    sample_n(nrow(AllFishObsWithParSeason))

#drop and re-add year column
AllFishObsWithParSeasonSimBeta <- AllFishObsWithParSeason %>%
    select(-season)

AllFishObsWithParSeasonSim  <- bind_cols(AllFishObsWithParSeasonSimBeta , shuff_seasons)

    
SWM_recruits <- AllFishObsWithParSeasonSim %>%
    filter(season=="SWM")%>%
    select(season, offs_site, par_site, matched_offs)
NEM_recruits <- AllFishObsWithParSeasonSim %>%
    filter(season=="NEM")%>%
    select(season, offs_site, par_site, matched_offs)
    

total_parSWM <- SWM_recruits %>%
    filter(matched_offs=="Y") %>%
    group_by(offs_site, par_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() 

total_parNEM <- NEM_recruits %>%
    filter(matched_offs=="Y") %>%
    group_by(offs_site, par_site) %>%
    summarise(n_matches=n()) %>%
    ungroup() 

for(i in 1:nrow(prop_samp)){
    
    
    if(is.nan(prop_samp$total_prop_hab_sampled_anems_tidied[i])){prop_samp$total_prop_hab_sampled_anems_tidied[i] <- 0.01} 
    if(is.infinite(prop_samp$total_prop_hab_sampled_anems_tidied[i])){prop_samp$total_prop_hab_sampled_anems_tidied[i] <- 1}
    ifelse(prop_samp$total_prop_hab_sampled_anems_tidied[i] == 0, 0.48, prop_samp$total_prop_hab_sampled_anems_tidied[i])

    
}

prop_sampSeas <- prop_samp %>%
    filter(time_frame=="2012-2018" ) %>%
    select(site, total_prop_hab_sampled_anems_tidied)

sitesSWM <- prop_sampSeas %>%
    select(site)
sitesSWM_2 <- suppressWarnings(semi_join(all_sites, sitesSWM, by="site")) %>%
    select(index)
sitesSWMt <- t(sitesSWM_2)


sitesSWM_beta <- sitesSWM %>%
    #select(-pop) %>%
    rename(par_site="site")

sitesSWM_beta$offs_site <- sitesSWM_beta$par_site

allsites_parentageSWM <- full_join(sitesSWM_beta, total_parSWM, by=c("par_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site)

##generate list of site indices
sitesNEM <- prop_sampSeas %>%
    select(site)
sitesNEM_2 <- suppressWarnings(semi_join(all_sites, sitesNEM, by="site")) %>%
    select(index)
sitesNEMt <- t(sitesNEM_2)

sitesNEM_beta <- sitesNEM %>%
    #select(-pop) %>%
    rename(par_site="site")

sitesNEM_beta$offs_site <- sitesNEM_beta$par_site

allsites_parentageNEM <- full_join(sitesNEM_beta, total_parNEM, by=c("par_site", "offs_site")) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    arrange(offs_site, par_site)

#turn into a matrix of the correct format
parmatSWM <- allsites_parentageSWM %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

###change NAs to 0
parmatSWM[is.na(parmatSWM)] <- 0
rownames(parmatSWM) <- suppressWarnings(parmatSWM$par_site)
parmatSWM$par_site <- NULL

#parmatSWM <- (parmatSWM[1:nrow(sitesSWM_beta), 1:nrow(sitesSWM_beta)])
#
#
parmatNEM <- allsites_parentageNEM %>% 
    filter(!is.na(offs_site) & !is.na(par_site)) %>%
    group_by(offs_site, par_site) %>%
    filter(row_number()==1) %>%
    spread(offs_site, n_matches)

##change NAs to 0
parmatNEM[is.na(parmatNEM)] <- 0
rownames(parmatNEM) <- suppressWarnings(parmatNEM$par_site)
parmatNEM$par_site <- NULL

N_gen_offs_all <- AllFishObsWithParSeasonSim %>% 
    filter(input=="offspring") %>%
    group_by(season, site) %>%
    summarise(n_offs_gen=n())
sum(N_gen_offs_all$n_offs_gen) #missing 5 fish, because they are marked J based on tail color but we don't have their size
#nrow(N_gen_offs %>% filter(fish_indiv %!in% all_seasons$fish_indiv))

N_gen_offs2 <- N_gen_offs_all %>%
    group_by(site, season) %>% #add group by year if doing annual
    summarise(sampled_fish=sum(n_offs_gen)) 


#make dataframe with column for unassigned juveniles in each column (offspring site)
offs_matched_site <- allsites_parentageSWM %>%
    group_by(offs_site) %>%
    summarise(n_offs=sum(n_matches, na.rm=T))

unassigned_beta <- semi_join(offs_matched_site, sitesSWM, by=c(offs_site="site")) %>%
    ungroup() %>%
    select(offs_site, n_offs)

#unassigned_beta[is.na(unassigned_beta)] <- 0


net_unassigned <- left_join(unassigned_beta, (N_gen_offs2 %>% filter(season=="SWM")), by=c(offs_site="site")) %>%
    mutate(n_unassigned=sampled_fish-n_offs)

net_unassigned$sampled_fish[is.na(net_unassigned$sampled_fish)] <- 0
net_unassigned$n_unassigned[is.na(net_unassigned$n_unassigned)] <- 0

#transpose to create a row, bind to the matrix

row2add <- net_unassigned %>%
    arrange(offs_site) %>%
    select(n_unassigned, offs_site) %>%
    group_by(offs_site) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    select(n_unassigned)
row2addT <- as.data.frame(t(row2add))
colnames(row2addT) <- colnames(parmatSWM)

parmatSWM <- ungroup(parmatSWM)
parmat_SWM_full <- bind_rows(parmatSWM, row2addT)
parmat_SWM_full[is.na(parmat_SWM_full)] <- 0
parmat_SWM_full[parmat_SWM_full <0 ] <- 0

#make dataframe with column for unassigned juveniles in each column (offspring site)
offs_matched_site <- allsites_parentageNEM %>%
    group_by(offs_site) %>%
    summarise(n_offs=sum(n_matches, na.rm=T))

unassigned_beta <- semi_join(offs_matched_site, sitesNEM, by=c(offs_site="site")) %>%
    ungroup() %>%
    select(offs_site, n_offs)

#unassigned_beta[is.na(unassigned_beta)] <- 0
net_unassigned <- left_join(unassigned_beta, (N_gen_offs2 %>% filter(season=="NEM")), by=c(offs_site="site")) %>%
    mutate(n_unassigned=sampled_fish-n_offs)

net_unassigned$sampled_fish[is.na(net_unassigned$sampled_fish)] <- 0
net_unassigned$n_unassigned[is.na(net_unassigned$n_unassigned)] <- 0

#transpose to create a row, bind to the matrix

row2add <- net_unassigned %>%
    arrange(offs_site) %>%
    select(n_unassigned, offs_site) %>%
    group_by(offs_site) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    select(n_unassigned)
row2addT <- as.data.frame(t(row2add))
colnames(row2addT) <- colnames(parmatNEM)

parmatNEM <- ungroup(parmatNEM)
parmat_NEM_full <- bind_rows(parmatNEM, row2addT)
parmat_NEM_full[is.na(parmat_NEM_full)] <- 0
parmat_NEM_full[parmat_NEM_full < 0] <- 0
    
    
#fit kernels 
Assignments <- parmat_SWM_full
Adult_sample_proportions <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/prop_samp18.csv", header=FALSE))
Sampled_reefs <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_index18.csv", header=FALSE))
Distances <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/distance_matrix_unsurveyed.csv", header=FALSE))
Reef_sizes <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", header=FALSE))
Centroids <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", header=T))

a=-10
b=10

x <- list(Distances=Distances, Assignments=Assignments, Sampled_reefs=Sampled_reefs, Reef_sizes=Reef_sizes, Adult_sample_proportions=Adult_sample_proportions) #put inputs into a list because that's the bbmle format

fit_bothSWM <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=1, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=x, control=list(maxit=500)))
kSWM_b <- coef(fit_bothSWM)[1]
thetaSWM_b <- coef(fit_bothSWM)[2]
MDDSWM_b <- cubintegrate(integrate_kernel_sum1, lower = 0, upper = Inf, k=kSWM_b, theta=thetaSWM_b, method = "pcubature")$integral
theta_eval <- thetaSWM_b
k_eval <- kSWM_b
MedianDispDistSWM <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2)
Dist90RetainedSWM <- round(nleqslv(x = 7, fn = cdf_solve90)$x, 2) 


col <- c("season","k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
BootstrapKernelsSWM <- as.data.frame(matrix(nrow=1, ncol=7), stringsAsFactors = FALSE)
colnames(BootstrapKernelsSWM) <- col

BootstrapKernelsSWM$season <- "SWM"
BootstrapKernelsSWM$k <- kSWM_b
BootstrapKernelsSWM$theta <- thetaSWM_b
BootstrapKernelsSWM$MDD <- MDDSWM_b
BootstrapKernelsSWM$MedianDispDist <- MedianDispDistSWM
BootstrapKernelsSWM$Dist90Retained <- Dist90RetainedSWM
BootstrapKernelsSWM$iteration <- as.character(n)



Assignments <- parmat_NEM_full
Adult_sample_proportions <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/prop_samp18.csv", header=FALSE))
Sampled_reefs <- as.matrix(read.csv("~/parentage/kernel_fitting/1340_loci/input/site_index18.csv", header=FALSE))
Distances <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/distance_matrix_unsurveyed.csv", header=FALSE))
Reef_sizes <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/area_unsurveyed.csv", header=FALSE))
Centroids <- as.matrix(read.csv("~/parentage/kernel_fitting/894_loci/centroids_unsurveyed.csv", header=T))

a=-10
b=10

x <- list(Distances=Distances, Assignments=Assignments, Sampled_reefs=Sampled_reefs, Reef_sizes=Reef_sizes, Adult_sample_proportions=Adult_sample_proportions) #put inputs into a list because that's the bbmle format

fit_bothNEM <- suppressWarnings(mle2(LL_kt_bbmle, start=list(k=1, theta=1), lower=c(-10, 0.15), upper=c(10, 5), method="L-BFGS-B", data=x, control=list(maxit=500)))
kNEM_b <- coef(fit_bothNEM)[1]
thetaNEM_b <- coef(fit_bothNEM)[2]
MDDNEM_b <- cubintegrate(integrate_kernel_sum1, lower = 0, upper = Inf, k=kNEM_b, theta=thetaNEM_b, method = "pcubature")$integral
theta_eval <- thetaNEM_b
k_eval <- kNEM_b
MedianDispDistNEM <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2)
Dist90RetainedNEM <- round(nleqslv(x = 7, fn = cdf_solve90)$x, 2) 

col <- c("season","k", "theta", "MDD", "MedianDispDist","Dist90Retained", "iteration")
BootstrapKernelsNEM <- as.data.frame(matrix(nrow=1, ncol=7), stringsAsFactors = FALSE)
colnames(BootstrapKernelsNEM) <- col

BootstrapKernelsNEM$season <- "NEM"
BootstrapKernelsNEM$k <- kNEM_b
BootstrapKernelsNEM$theta <- thetaNEM_b
BootstrapKernelsNEM$MDD <- MDDNEM_b
BootstrapKernelsNEM$MedianDispDist <- MedianDispDistNEM
BootstrapKernelsNEM$Dist90Retained <- Dist90RetainedNEM
BootstrapKernelsNEM$iteration <- as.character(n)




#bind together all simulated kernel fits

SimulatedKernelsBetaSeas <- bind_rows(BootstrapKernelsNEM, BootstrapKernelsSWM)

SimulatedKernelsSeas <- bind_rows(SimulatedKernelsSeas, SimulatedKernelsBetaSeas)

setTxtProgressBar(pb, n)


}
close(pb)
EndTime <- Sys.time()
EndTime-StartTime
options(warn=0)
write.csv(SimulatedKernelsSeas, file="~/parentage/kernel_fitting/1340_loci/final_results/simulations/SimulatedKernelsSeasonal.csv", row.names=T, quote=FALSE)


head(SimulatedKernelsSeas)

SeasonSimMed <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/simulations/SimulatedKernelsSeasonal.csv", header=T, , stringsAsFactors = F) #%>%#load simulations 
#    select(-X)
seasonal_kernels <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/tables/RecruitSizeAsSeasonlity_summary.csv", header=T) #load empirical
head(seasonal_kernels)
#summary(SeasonSimMed %>% filter(season=="NEM"))
#summary(SeasonSimMed %>% filter(season=="SWM"))


#compare simulation results to empirical
SimulatedKernelsVar <- SimulatedKernelsSeas %>%
    group_by(iteration) %>%
    mutate(MDD=ifelse(MDD <0, 0,MDD)) %>% #if MDD is less than zero because of negative k and low theta (near 0.10), replace with 0
    mutate(cv_k=sd(k)/mean(k, na.rm=T)) %>%
    mutate(cvtheta=sd(theta)/mean(theta, na.rm=T)) %>%
    mutate(sdk=sd(k)) %>%
    mutate(cvMDD=sd(MDD)/mean(MDD, na.rm=T)) %>%
    #mutate(meanMDD=mean(MDD)) %>%
    mutate(cvMed=sd(MedianDispDist)/mean(MedianDispDist, na.rm=T))%>%
    #mutate(sdMDD=sd(MDD)) %>%
    mutate(sdMed=sd(MedianDispDist)) %>%
    mutate(cv90=sd(Dist90Retained)/mean(Dist90Retained, na.rm=T)) %>%
    group_by(iteration) %>%
    distinct(iteration, .keep_all = T) %>%
    ungroup() %>%
    sample_n(1000, replace=F) #get back to 1000 simulations

seasonal_kernels <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/tables/RecruitSizeAsSeasonlity_summary.csv", header=T) 
RealKernelVar <- seasonal_kernels %>%
    mutate(meanMDD=mean(MeanDispersalDistance)) %>%
    mutate(sdMDD=sd(MeanDispersalDistance)) %>%
    mutate(cvMDD=sd(MeanDispersalDistance)/mean(MeanDispersalDistance, na.rm=T))%>%
    mutate(cvMed=sd(MedianDispersalDistance)/mean(MedianDispersalDistance, na.rm=T))%>%
    mutate(cvtheta=sd(best_theta)/mean(best_theta, na.rm=T)) %>%
    mutate(sdk=sd(best_k)) %>%
    mutate(sdMed=sd(MedianDispersalDistance)) %>%
    mutate(cv90=sd(Dist90Retained)/mean(Dist90Retained, na.rm=T))%>%
    distinct(cvMDD, cvtheta, sdk, cvMed, cv90)

    

RealKernelVar #empirical values

#what percentile are our observations in?

((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(sdk > RealKernelVar$sdk)))+1)/1001
((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(cvtheta > RealKernelVar$cvtheta)))+1)/1001
((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(cvMDD > RealKernelVar$cvMDD)))+1)/1001
((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(cvMed > RealKernelVar$cvMed)))+1)/1001
((nrow(SimulatedKernelsVar)-nrow(SimulatedKernelsVar %>% filter(cv90 > RealKernelVar$cv90)))+1)/1001


pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/SeasonalMedCV.pdf")
hist(SimulatedKernelsVar$cvMed, main=NA, col="grey", breaks=seq(0, 1, 0.01), xlab= "Median dispersal distance CV")
abline(v = RealKernelVar$cvMed, col="red", lwd=3, lty=2)
dev.off()

pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/SeasonalMDDCV.pdf")
hist(SimulatedKernelsVar$cvMDD, col="grey", breaks=seq(0, .5, 0.01), main=NA, xlab= "Mean dispersal distance CV")
abline(v = RealKernelVar$cvMDD, col="red", lwd=3, lty=2) #from cv calculated in testing_seasonal.ipynb 02/20/2020
dev.off()

pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/SeasonalThetaCV.pdf")
hist(SimulatedKernelsVar$cvtheta, col="grey", breaks=seq(0, 1, 0.01), main=NA, xlab="theta CV")
abline(v = RealKernelVar$cvtheta, col="red", lwd=3, lty=2)
dev.off()

pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/SeasonalKSD.pdf")
hist(SimulatedKernelsVar$sdk, col="grey", main=NA, breaks=seq(0, 1.5, 0.01), xlab= "k SD")
abline(v = RealKernelVar$sdk, col="red", lwd=3, lty=2)
dev.off()

pdf("~/parentage/kernel_fitting/1340_loci/final_results/simulations/Seasonal90CV.pdf")
hist(SimulatedKernelsVar$cv90,col="grey", breaks=seq(0, 1.5, 0.01), main=NA, xlab= "0.90 dispersal distance CV")
abline(v = RealKernelVar$cv90, col="red", lwd=3, lty=2)
dev.off()

hist(SimulatedKernelsVar$cvMDD, breaks=seq(0, 3, 0.01), main=NA, col="black", xlab= "Mean dispersal distance CV")
abline(v = 1.233, col="red", lwd=3, lty=2) #from cv calculated in testing_seasonal.ipynb 02/20/2020













#it looks like the bimodal distribution of CV is driven by the large right skew in the MDD data
adj_max <- 100 #max(SimulatedKernels$MDD)
exp <- matrix(nrow=adj_max, ncol=1)
x <- matrix(seq(1,adj_max, 1))


for(i in 1:nrow(exp)) {
    
exp[i,] = exp(x[i,])

}

cv <- matrix(nrow=max(SimulatedKernels$MDD), ncol=1)

pb <- txtProgressBar(min = 0, max = nrow(exp), style = 3)

for(i in 1:nrow(exp)){
    
    exp_eval <- exp[sample(nrow(exp),size=7,replace=TRUE),]
    exp_eval[is.infinite(exp_eval)] <- NA
    cv[i,] = sd(exp_eval, na.rm=T)/mean(exp_eval, na.rm=T)

setTxtProgressBar(pb, i)    
    
}
close(pb)
hist(cv)

hist(SimulatedKernels$MDD)
