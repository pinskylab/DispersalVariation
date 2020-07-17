Packages <- c("dplyr",  "RMySQL",  "tidyr", "lubridate")
invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))
setwd('/local/home/katrinac/parentage')

load("~/parentage/r_data/sampled_area_each_year.RData")
load("~/parentage/r_data/anems_visited_by_year.RData")
load("~/parentage/r_data/total_sampling_across_years.RData")

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
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0


nrow(fish_obs)


#dat_gen <- read.table("/local/home/katrinac/parentage/colony2/20190422_894loci/input/seq33-03_norecap.gen", colClasses= "character",  quote="", header=TRUE, row.names=NULL, sep=" ")
#loci_to_remove <- read.table("/local/home/katrinac/parentage/colony2/20190422_894loci/input/20190523_894_plink.prune.out", quote="", header=T, stringsAsFactors = F)

#dim(dat_gen)
#dim(loci_to_remove)
#218/2

#dim(dat_gen) #697 loci total
##dim(loci_to_keep)
#
#
##do this after filtering for only loci not in LD and maf 0.05
##names(dat_gen) <- gsub("dDocent_Contig_" , "",  names(dat_gen))  
##names(dat_gen) <- gsub("\\d_" , "",  names(dat_gen))
#
#
#
##find and drop the loci that have three alleles
#
##flip the data frame so you can use filter
#
##dat_gen$ligation_id
#
#row.names(dat_gen) <-dat_gen$ligation_id
#
##dat_gen$ligation_id <- NULL
#
#
###could use this code to split columns rather than doing it as in terminal20180517 --> cSplit(mydf, grep("Allele", names(mydf)), "", stripWhite = FALSE)
#
#
#dat_gen_t <- as.data.frame(t(dat_gen), stringsAsFactors = F)#look for loci appended with 2 for the 3 allele loci
##dat_gen_t$ligation_id <- rownames(dat_gen_t)
#
#dat_gen_t$locus <- rownames(dat_gen_t)
##head(dat_gen_t)
#
#ligation_ids <- dat_gen %>%
#    select(ligation_id)
#
#no_ld <- dat_gen_t %>%
#  filter(locus=="ligation_id" | locus %!in% loci_to_remove$locus)
##
#dim(no_ld)
##
##
###drop the columns with a 03 value
##
#dat_gen_t_beta <- no_ld %>% mutate(third_allele = rowSums(no_ld=="03"))
##
#sum(dat_gen_t_beta$third_allele != 0) #one locus of the 697 is in LD
##
##head(dat_gen_t_beta)
##
#dat_gen_beta2 <- dat_gen_t_beta %>% filter(third_allele == 0)
##
###head(dat_gen_beta2)
##
#col <- dat_gen_beta2$locus
#dat_gen_beta2$locus <- NULL
#dat_gen_beta2$third_allele <- NULL
##
#dat_gen <- t(dat_gen_beta2)
##
##
#dat_gen <- as.data.frame(dat_gen, stringsAsFactors = F, row.names = NULL, col.names = col)
##
#names(dat_gen) <- col
##
###dat_gen$ligation_id <- row.names(dat_gen)
#row.names(dat_gen) <- NULL
##
#dim(dat_gen)
##
##dat_gen2 <- dat_gen %>%
##    filter(ligation_id %in% dat_gen_FORIDS$ligation_id)

#2895-2680
#2680/2
#215/2

#dim(dat_gen)

#write.table(dat_gen, file="/local/home/katrinac/parentage/colony2/20190523_1340loci/input/20190523_1340_loci_seq33-03_noLD.txt", col.names=T, quote=FALSE, row.names=FALSE)



#load the sequencing data already modified
dat_gen <- read.table("/local/home/katrinac/parentage/colony2/20190523_1340loci/input/20190523_1340_loci_seq33-03_noLD.txt", colClasses= "character",  quote="", header=TRUE, row.names=NULL)



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

fish_in_gen_lig <- dat_gen %>% #pull of the ligation_ids of fish in dat_gen, and remove those in the issues table
    select(ligation_id) %>%
    filter(ligation_id %!in% issues$ligation_id) #5 fish in issues table, these are dropped
nrow(dat_gen) -nrow(fish_in_gen_lig) 

fish_in_gen_samp <- lig_samp %>% 
    filter(ligation_id %in%dat_gen$ligation_id) %>%
    distinct(ligation_id, sample_id, .keep_all = T) #all of the ligation_ids in fish_in_gen_samp are in fish_in_gen_lig. Are the extra 7 sample_ids just multiples for regenotyped fish? YES see the following when this line is commented out: test %>% group_by(ligation_id, sample_id) %>% summarise(n=n()) %>% filter(n >1)

nrow(fish_in_gen_lig)- nrow(fish_in_gen_samp) #should be 0

sampled_fish <- left_join(allfish_caught, fish_in_gen_samp,  by="sample_id") %>%
    select(fish_indiv, year, size, color, sex, gen_id, ligation_id, sample_id) %>%
    filter(!is.na(sample_id))

fish_meta <- sampled_fish %>%#join with the rest of the data
    select(fish_indiv, year, size, color, sex, gen_id, ligation_id, sample_id)
###join with genetic data
dat_gen_meta <-  left_join(dat_gen, fish_meta, by="ligation_id")


dim(dat_gen_meta)

#dat_gen_meta3 <- row_count(dat_gen_meta2, count = "00", append = T)#change the "00" genotypes to NA to count loci genotyped and select the ligation with the greatest number of genotyped loci
#this step is unnecessary for the 894 loci because Michelle did this step in her protocol before giving me the genepop
#best_genos <- dat_gen_meta3 %>% 
#    group_by(gen_id) %>%
#    #filter(n()>1) %>%
#    slice(which.min(rowcount)) %>% #select the ligation_id that have the fewest NA loci
#    select(-rowcount) %>%
#    ungroup()
#dim(best_genos)
#
#why are there 7 fish not included?only with 697 loci
#test <- dat_gen_meta3 %>% filter(ligation_id %!in% best_genos$ligation_id)
#dim(test)
#do they have multiple ocurrences of their gen_id? YES 
#dim(test %>% filter(gen_id %in% dat_gen_meta3$gen_id))

#dat_gen_meta4 <- dat_gen_meta2 %>%
#    select(year, gen_id, color, size, sex, ligation_id) #now we just need the metadata for sorting
#dim(dat_gen_meta2)#should be the same number of rows, just fewer columns
#dim(dat_gen_meta4)#should be the same number of rows, just fewer columns

years <- as.matrix(c(2012, 2013, 2014, 2015, 2016, 2017, 2018)) #years to evaluate
indiv_ids <- as.matrix(dat_gen_meta %>% #gen_ids to sort into parent/offspring files
                     select(fish_indiv) %>% 
                     distinct(fish_indiv))

parents <- as.data.frame(matrix(nrow=0, ncol=8))  #copy previous year parent file to current parent file by starting this file out of the loop and building on it in the loop. mutate a column with the year it will be a part of, and sort through that as "for 2013, grab everything in the parent file >=2013" after the loop
colnames(parents) <- names(fish_meta)

offs_annual <- as.data.frame(matrix(nrow=0, ncol=8))
colnames(offs_annual) <- names(fish_meta)
    
for(i in 1:nrow(years)){
    
    year_eval <- years[i]
    
    for(k in 1:nrow(indiv_ids)){ 
        
        indiv_id_eval <- indiv_ids[k]
        
        meta_eval <- fish_meta %>% #this will be the parents
            filter(fish_indiv == indiv_id_eval & year == year_eval) %>% #year == year_eval, taking this out in testing with AD, add it back in later if you want and sort for year later in the loop
            filter(size >= 6 | sex=="M" | sex=="F")  %>%
            #filter(sex != "J" ) %>%
            #filter(fish_indiv == indiv_id_eval) %>% #& ligation_id %in% best_genos$ligation_id) #commented out 06/01/2020
            distinct(fish_indiv, .keep_all = T) #06/01/2020, shouldn't need this unless I fish was sampled multiple times in a year, which could happen if they were caught on a recap dive, but then there wouldn't be a sample id for them. and even if there was it wouldn't matter because the observation would be the same size/year probably so distinct is works to get rid of that if it's happening
        
        if(indiv_id_eval %in% meta_eval$fish_indiv){#then it's a parent
       
            if(indiv_id_eval %!in% parents$fish_indiv){ #if it's not already in the parent file, add it #maybe change this to elseif
                
                parents <-  bind_rows(meta_eval, parents)  

                
            }
            
#then it's an offspring
        }else{
            
          
             if(indiv_id_eval %!in% offs_annual$fish_indiv){ #if it's not already in the offspring file, add it
                    
                    add_to_offspring <- fish_meta %>%
                            filter(fish_indiv == indiv_id_eval & year==year_eval)%>% #& ligation_id %in% best_genos$ligation_id)
                            distinct(fish_indiv, .keep_all = T)
                 offs_annual <- bind_rows(offs_annual, add_to_offspring) 
                 
        
            }
        }

    }
}
parents <- parents #%>%
    #filter(!is.na(gen_id))
offs_annual <- offs_annual #%>%
    #filter(!is.na(gen_id))


nrow(parents)
nrow(offs_annual) #any fish with a weird tail color here (YP or O) has been corrected by Michelle, as indicated in the "fish_correction" column.
nrow(parents)+nrow(offs_annual)#should be 2625

#print and upload the histogram of parent/offspring sizes, the total fish captured in each year, total offspring in the offspring file for each year, total in the parents of each year, and the number of each tail color in the parent and offspring file so everyone can look on github, total number of parents in parentage.

#print and output summaries for Clownfish team
#double check total number
nrow(parents)+nrow(offs_annual) #should be 2483

#png(file="~/parentage/colony2/20190523_1340loci/input/size_distribution_parents.png")
hist(parents$size, main="size distribution of fish in parent files (all <6cm are YP/F or O/M)")#looks good, but make sure the fish smaller than 6 in the parents. YES, see test below
#dev.off()
#png(file="~/parentage/colony2/20190523_1340loci/input/size_distribution_offspring.png")
hist(offs_annual$size, main="size distribution of fish in offspring files") #looks good
#dev.off()
test <- parents %>% filter(size <6)
(test$sex)

#table of summaries

#how many fish from each year?
n_samples <- fish_meta %>% filter(!is.na(sample_id)) %>% group_by(year) %>% summarise(n_observed=n()) 
sum(n_samples$n_observed)
n_samples
#write.table(n_samples, file="~/parentage/colony2/20190523_1340loci/input/annual_number_fish_geno.txt", col.names = T, quote=F, row.names = F)

#what are the number of fish of each sex in parents and offspring file?
par_sex <- parents %>% group_by(sex) %>% summarise(n=n()) 
#write.table(par_sex, file="~/parentage/colony2/20190523_1340loci/input/parent_nsex.txt", col.names = T, quote=F, row.names = F)
offs_sex <- offs_annual %>% group_by(sex) %>% summarise(n=n()) 
#write.table(offs_sex, file="~/parentage/colony2/20190523_1340loci/input/offspring_nsex.txt", col.names = T, quote=F, row.names = F)



#change NA genotypes to 0 for Colony
#parents_gen[is.na(parents_gen)] <- 00
#offs_gen[is.na(offs_gen)] <- 00



#sort into Male/Female/Offspring and annual files
#2012

all_offs12 <- offs_annual %>%
    filter(year == 2012)  
all_offs12_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_offs12$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T) %>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)

all_par12 <- parents %>%
    filter(year==2012 & fish_indiv %!in% all_offs12$fish_indiv)
all_par12_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_par12$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)



#2013

all_offs13 <- offs_annual %>%
    filter(year==2013)
all_offs13_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_offs13$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)

all_par13 <- parents %>%
    filter(year %in% c(2012, 2013) & fish_indiv %!in% all_offs13$fish_indiv)
all_par13_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_par13$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)



#2014
all_offs14 <- offs_annual %>%
    filter(year==2014)
all_offs14_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_offs14$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)


all_par14 <- parents %>%
    filter(year %in% c(2012, 2013, 2014) & fish_indiv %!in% all_offs14$fish_indiv) 
all_par14_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_par14$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)



##2015
all_offs15 <- offs_annual %>%
    filter(year==2015)
all_offs15_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_offs15$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)

all_par15 <- parents %>%
    filter(year %in% c(2012, 2013, 2014, 2015) & fish_indiv %!in% all_offs15$fish_indiv) 
all_par15_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_par15$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)



##2016
all_offs16 <- offs_annual %>%
    filter(year==2016)
all_offs16_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_offs16$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)


all_par16 <- parents %>%
    filter(year %in% c(2012, 2013, 2014, 2015, 2016) & fish_indiv %!in% all_offs16$fish_indiv) 
all_par16_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_par16$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)



#2017
all_offs17 <- offs_annual %>%
    filter(year==2017) 
all_offs17_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_offs17$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)


all_par17 <- parents %>%
    filter(year %in% c(2012, 2013, 2014, 2015, 2016, 2017) & fish_indiv %!in% all_offs17$fish_indiv)  
all_par17_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_par17$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)



#2018
all_offs18 <- offs_annual %>%
    filter(year==2018)
all_offs18_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_offs18$fish_indiv)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)


all_par18 <- parents %>%
    filter(year %in% c(2012, 2013, 2014, 2015, 2016, 2017, 2018) & fish_indiv %!in% all_offs18$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)
all_par18_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% all_par18$fish_indiv) %>%
    distinct(fish_indiv, .keep_all = T)%>% 
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id)



#pooled analysis
pooled <- bind_rows(parents, offs_annual) %>%
    distinct(fish_indiv, .keep_all = T)
pooled_gen <- dat_gen_meta %>%
    filter(fish_indiv %in% pooled$fish_indiv)%>% 
    distinct(fish_indiv, .keep_all = T) %>%
    select(-fish_indiv, -year, -size, -color, -sex, -gen_id, -sample_id) 




dim(pooled_gen)



#fish_obs %>% 
#    filter(!is.na(gen_id)) %>% 
#    filter(!is.na(tag_id))%>% 
#    group_by(fish_indiv) %>% 
#    filter(n()>1) %>% 
#    ungroup() %>% 
#    #distinct(gen_id, .keep_all = T) %>% 
#    group_by(gen_id) %>% 
#    filter(n()>1)
#
#

input_summary <- as.data.frame(matrix(nrow=8, ncol=3))
col <- c("year", "n_offspring", "n_parents")
colnames(input_summary) <- col
#
input_summary$year <- c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "all_years")
input_summary$n_offspring <- c(nrow(all_offs12_gen), nrow(all_offs13_gen), nrow(all_offs14_gen), nrow(all_offs15_gen), nrow(all_offs16_gen), nrow(all_offs17_gen), nrow(all_offs18_gen), nrow(offs_annual))
input_summary$n_parents <- c(nrow(all_par12_gen), nrow(all_par13_gen), nrow(all_par14_gen), nrow(all_par15_gen), nrow(all_par16_gen), nrow(all_par17_gen), nrow(all_par18_gen), nrow(parents))
write.table(input_summary, file="~/parentage/colony2/20190523_1340loci/input/input_summary.txt", col.names=T, quote=FALSE, row.names=FALSE )




input_summary
nrow(pooled_gen)#should be 2405
sum(input_summary$n_offspring)-nrow(offs_annual)#should be 903
nrow(parents)


#write the offspring and parent files for the N_captured values in kernel_fit_file_prep.ipynb
#write.table(parents, file="~/parentage/colony2/20190523_1340loci/input/all_parents_corrected.txt", col.names=T, quote=FALSE, row.names=FALSE )
#write.table(offs_annual, file="~/parentage/colony2/20190523_1340loci/input/all_offspring_corrected.txt", col.names=T, quote=FALSE, row.names=FALSE )


#write.table(pooled_gen, file="~/parentage/colony2/20190523_1340loci/input/full_genetics_pooled.txt", col.names=FALSE, quote=FALSE, row.names=FALSE )



## make a markers file

dup <- parents_gen %>% select(-ligation_id,-year, -contains(".02"))

#make loci only number for colony
names(dup) <- gsub("dDocent_Contig_" , "",  names(dup))  
names(dup) <- gsub("\\d_" , "",  names(dup))
names(dup) <- gsub(".01" , "",  names(dup))  


locus_names <- as.matrix(colnames(dup))
length(locus_names) #should be 696
error_matrix <- matrix(nrow=4, ncol=894)
error_matrix[1,] <- locus_names[1:894]
error_matrix[2,] <- 0
error_matrix[3,] <- 0
error_matrix[4,] <- 0.026
head(error_matrix)
#write.table(error_matrix, file="/local/home/katrinac/parentage/colony2/20190422_894loci/input/894_loci_error_mat.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)



#proportion sampled in each time frame for colony "prob of parent in parent file" input
total_sampling_across_years$total_area_method <- as.character(total_sampling_across_years$total_area_method)
total_sampling_across_years$total_anems_method <- as.character(total_sampling_across_years$total_anems_method)

total_sampling_across_years %>% filter(total_area_method=="anems tidied" &total_anems_method=="metal tags")

input_summary #2018 n_parents should be the full number of parents
sum(input_summary$n_offspring) #792

