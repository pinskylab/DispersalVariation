Packages <- c("dplyr",  "RMySQL",  "tidyr", "lubridate")
invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))
setwd('/local/home/katrinac/parentage')

#source("~/scripts/conleyte.R")
#source("~/scripts/conlabor.R")
#source("~/scripts/sample_latlon.R")
#source("~/scripts/get_date.R")
#source("~/scripts/get_samp.R")
#source("~/scripts/anem_meta.R")
#source("~/scripts/get_fish.R")
#source("~/scripts/get_lig.R")
#source("~/scripts/create_text.R")
load("~/parentage/r_data/sampled_area_each_year.RData")
load("~/parentage/r_data/anems_visited_by_year.RData")
load("~/parentage/r_data/total_sampling_across_years.RData")

# download the fish_Obsfile
download.file(url = "https://github.com/pinskylab/genomics/blob/master/data/fish-obs.RData?raw=true", destfile = "~/parentage/r_data/fish-obs.RData")

# read in the data
fish_obs <- readRDS("~/parentage/r_data/fish-obs.RData") 

#fish_obs= readRDS("~/parentage/r_data/fish-obs.RData")


#labor <- conlabor()
#leyte <- conleyte()

"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0
#issues <- as.data.frame(labor %>% tbl("known_issues") %>% select(ligation_id))
#issues <- get_samp(issues$ligation_id)

nrow(fish_obs)

dat_gen <- read.table("/local/home/katrinac/parentage/colony2/20190422_894loci/input/seq33-03_norecap.gen", colClasses= "character",  quote="", header=TRUE, row.names=NULL, sep=" ")
loci_to_remove <- read.table("/local/home/katrinac/parentage/colony2/20190422_894loci/input/20190523_894_plink.prune.out", quote="", header=T, stringsAsFactors = F)

dim(dat_gen)
dim(loci_to_remove)
#218/2

dim(dat_gen) #697 loci total
#dim(loci_to_keep)
#
#
#do this after filtering for only loci not in LD and maf 0.05
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
row.names(dat_gen) <-dat_gen$ligation_id
#
#dat_gen$ligation_id <- NULL
#
#
###could use this code to split columns rather than doing it as in terminal20180517 --> cSplit(mydf, grep("Allele", names(mydf)), "", stripWhite = FALSE)
#
#
dat_gen_t <- as.data.frame(t(dat_gen), stringsAsFactors = F)#look #for loci appended with 2 for the 3 allele loci
dat_gen_t$ligation_id <- rownames(dat_gen_t)

dat_gen_t$locus <- rownames(dat_gen_t)

ligation_ids <- dat_gen %>%
   select(ligation_id)

no_ld <- dat_gen_t %>%
  filter(locus=="ligation_id" | locus %!in% loci_to_remove$locus)

dim(no_ld)

###drop the columns with a 03 value
dat_gen_t_beta <- no_ld %>% mutate(third_allele = rowSums(no_ld=="03"))

sum(dat_gen_t_beta$third_allele != 0) #one locus of the 697 is in LD

dat_gen_beta2 <- dat_gen_t_beta %>% filter(third_allele == 0)

head(dat_gen_beta2)

col <- dat_gen_beta2$locus
dat_gen_beta2$locus <- NULL
dat_gen_beta2$third_allele <- NULL
dat_gen <- t(dat_gen_beta2)

dat_gen <- as.data.frame(dat_gen, stringsAsFactors = F, row.names = NULL, col.names = col)

names(dat_gen) <- col

dat_gen$ligation_id <- row.names(dat_gen)
row.names(dat_gen) <- NULL

dat_gen2 <- dat_gen %>%
    filter(ligation_id %in% dat_gen_FORIDS$ligation_id)

dat_gen <- dat_gen2


#check
sum(dat_gen=="00")

test <- dat_gen$ligation_id
head(fish_obs)

fish_in_gen_lig <- dat_gen %>% #pull of the ligation_ids of fish in dat_gen, and remove those in the issues table
    select(ligation_id) %>%
    filter(ligation_id %!in% issues$ligation_id) #5 fish in issues table, these are dropped
nrow(dat_gen) -nrow(fish_in_gen_lig) 

fish_in_gen_samp <- get_samp(fish_in_gen_lig$ligation_id) %>% 
    distinct(ligation_id, sample_id, .keep_all = T) -
nrow(fish_in_gen_lig)- nrow(fish_in_gen_samp) #should be 0

fish <- get_fish() %>% #get all identifying info for a sample
    filter(sample_id %in% fish_obs$sample_id) %>%
    distinct(sample_id, .keep_all = T)

#join all meta-data to know who an individual fish is, including all of its potential ids
fish_ids <- left_join(fish_obs, fish, by=c("sample_id", "tag_id", "fish_table_id")) 
fish_ids <- left_join(fish_ids, fish_in_gen_samp, by="sample_id")
fish_date <- get_date(fish_obs$sample_id)%>%# add the date to the sample_id
    distinct(sample_id, .keep_all = T) %>% 
    mutate(date=ymd(date)) %>% #format for lubridate
    mutate(year=year(date)) %>%
    select(-date)
fish_meta <- left_join(fish_ids, fish_date, by="sample_id") %>%#join with the rest of the data
    select(fish_indiv, year, size, color, sex, gen_id, ligation_id, sample_id)
###join with genetic data
dat_gen_meta <-  left_join(dat_gen, fish_meta, by="ligation_id")

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
        
        meta_eval <- fish_meta %>%
            filter(fish_indiv == indiv_id_eval & year == year_eval) %>% #year == year_eval, taking this out in testing with AD, add it back in later if you want and sort for year later in the loop
            filter(size >= 6 | sex=="M" | sex=="F")  %>%
            #filter(sex != "J" ) %>%
            filter(fish_indiv == indiv_id_eval) %>% #& ligation_id %in% best_genos$ligation_id)
            distinct(fish_indiv, .keep_all = T)
        
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


hist(parents$size, main="size distribution of fish in parent files (all <6cm are YP/F or O/M)")#looks good, but make sure the fish smaller than 6 in the parents. YES, see test below

hist(offs_annual$size, main="size distribution of fish in offspring files") #looks good

test <- parents %>% filter(size <6)
(test$sex)

#table of summaries

#how many fish from each year?
n_samples <- fish_meta %>% filter(!is.na(sample_id)) %>% group_by(year) %>% summarise(n_observed=n()) 
sum(n_samples$n_observed)
n_samples
#what are the number of fish of each sex in parents and offspring file?
par_sex <- parents %>% group_by(sex) %>% summarise(n=n()) 

offs_sex <- offs_annual %>% group_by(sex) %>% summarise(n=n()) 

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

input_summary$year <- c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "all_years")
input_summary$n_offspring <- c(nrow(all_offs12_gen), nrow(all_offs13_gen), nrow(all_offs14_gen), nrow(all_offs15_gen), nrow(all_offs16_gen), nrow(all_offs17_gen), nrow(all_offs18_gen), nrow(offs_annual))
input_summary$n_parents <- c(nrow(all_par12_gen), nrow(all_par13_gen), nrow(all_par14_gen), nrow(all_par15_gen), nrow(all_par16_gen), nrow(all_par17_gen), nrow(all_par18_gen), nrow(parents))


input_summary


## make a Colony markers file

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

#proportion sampled in each time frame for colony "prob of parent in parent file" input
total_sampling_across_years$total_area_method <- as.character(total_sampling_across_years$total_area_method)
total_sampling_across_years$total_anems_method <- as.character(total_sampling_across_years$total_anems_method)

total_sampling_across_years %>% filter(total_area_method=="anems tidied" &total_anems_method=="metal tags")

input_summary 
sum(input_summary$n_offspring)


