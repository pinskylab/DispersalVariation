
Packages <- c("dplyr", "ggplot2", "fields", "stringr", "reshape2", "splitstackshape", "readr",  "dplyr", "tidyr", "stringr", "tidyverse", "tibble", "ggfortify", "lubridate", "RColorBrewer", "vegan", "vcfR")

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

setwd('/local/home/katrinac/migest')

source("~/scripts/conleyte.R")
source("~/scripts/conlabor.R")
labor <- conlabor()
leyte <- conleyte()


#read in parentage matches
mums <- read.csv(file="~/parentage/colony2/20180712_2188loci/results/20180713colony_migest_mums_allyears.csv", header=TRUE)
dads <- read.csv(file="~/parentage/colony2/20180712_2188loci/results/20180713colony_migest_dads_allyears.csv", header=TRUE)

#read in demography estimates
prop_samp <- read.csv(file="~/migest/input/20180808prop_habitat_samp.csv", header=TRUE)

N_cap <- read.csv(file="~/migest/input/20180713_Nsampled.csv", header = TRUE)

N_cap_offs <- read.csv(file="~/migest/input/20180713_Nsampled_offs.csv")

N_trios <- read.csv(file="~/parentage/colony2/20180712_2188loci/results/20180713colony_migest_trios_allyears.csv", header=TRUE)

#iterate the input generation over all sites

#generate site list
prop_samp$site <- as.character(prop_samp$site)
prop_samp <- prop_samp %>% 
    arrange(site)
#    filter(site != "Sitio Tugas" & site!= "Sitio Lonas" & site!="Hicgop") %>%



#prop_samp

prop_samp$prop_habitat_sampled <- round(prop_samp$prop_habitat_sampled, digits=4)

###SENSITIVITY TO PROPORTION SAMPLED, COMMENT OUT THIS CELL OTHERWISE. What does the connectivity matrix look like when we add 10% or subtract 10% from prop_samp?
#prop_samp <- prop_samp %>%
   #mutate(prop_habitat_sampled_double=prop_habitat_sampled*2) #>%
   #mutate(prop_habitat_sampled_minus.1=prop_habitat_sampled-.1)
    

#prop_samp$prop_habitat_sampled_double

###SENSITIVITY TO N_offs, COMMENT OUT THIS CELL OTHERWISE. What does the connectivity matrix look like when we add 15% or subtract 15% from n_cap_offs?
#N_cap_offs <- N_cap_offs %>%
 # mutate(N_cap_offs_plus.15=n_captured_fish+(.15*n_captured_fish)) %>%
 # mutate(N_cap_offs_minus.15=n_captured_fish-(.15*n_captured_fish))
    

###SENSITIVITY TO N_cap, COMMENT OUT THIS CELL OTHERWISE. What does the connectivity matrix look like when we add 15% or subtract 15% from n_cap?
#N_cap <- N_cap %>%
 # mutate(N_cap_plus.15=n_captured_fish+(.15*n_captured_fish)) %>%
 # mutate(N_cap_minus.15=n_captured_fish-(.15*n_captured_fish))
    

#there has to be a more elegant way to make an iput file per year, but for now, just manually filter N_cap, N_cap_offs, N_trios, and mums/dads to evaluate for in the for loop
mums <- mums %>%
    filter(year==2013)
dads <- dads %>%
    filter(year==2013)
N_trios <- N_trios %>%
    filter(year== 2013)
N_cap_offs <- N_cap_offs %>%
    filter(year==2013)
N_cap <- N_cap %>%
    filter(year==2013)
prop_samp <- prop_samp %>%
    filter(year==2013) %>%
    select(site, prop_habitat_sampled)

#table of number of adults sampled

N_cap$year <- as.numeric(N_cap$year)
N_cap2 <- N_cap %>%
filter(year >= 2012 & year <= 2015) %>%
group_by(site) %>%
summarise(sampled_fish=sum(n_captured_fish))
    
N_cap3 <- N_cap2 %>%
mutate(n_male=(sampled_fish/2)) %>%
mutate(n_female=(sampled_fish/2))

#table of number of juveniles sampled
N_cap_offs2 <- N_cap_offs %>%
filter(year >= 2012 & year <= 2015) %>%
group_by(site) %>%
summarise(sampled_fish=sum(n_captured_fish)) 


N_cap3$n_female <- round(N_cap3$n_female, digits=0)
N_cap3$n_male <- round(N_cap3$n_male, digits=0)

N_cap_offs2$sampled_fish <- round(N_cap_offs2$sampled_fish, digits=0)


#use Ncap to generate sites list
sites1 <- N_cap_offs2 %>%
    ungroup() %>%
    select(site) %>%
    arrange(site)

sites2 <- N_cap3 %>%
    ungroup() %>%
    select(site) %>%
    arrange(site)
sites <- bind_rows(sites1, sites2)
sites <- sites %>%
    group_by(site) %>%
    filter(row_number()==1) %>%
    ungroup()
sites$site <- as.character(sites$site)
#migest won't take site names, just pop numbers. The sites are in alphabetical order, here is their order with corresponding pop numbers
sites$pop <- c(seq(from=1, to=nobs(sites), by=1))

#having a problem here-- migest is using our sampling in that year... but some years, we find a parentage match where the parent isn't in a population sampled that year. I think that I will drop these from the migest analysis, because including them means we have to include that population as "sampled" in that particular year, when it wasn't. And we would have to input parameters about it, which would largely be arbitrary. I think the error introduced by including these far outweighs the information gained. Malin agrees.



#(prop_samp) #we are way underestimating proportion sampled

#create empty data frame to add individual site tables to 
empt_dat2 <- as.data.frame(matrix(nrow=0, ncol=(nobs(sites))))
col <- c(sites$pop)
colnames(empt_dat2) <- col

#create an empty data frame to store site names
site_inc <- as.data.frame(matrix(nrow=0, ncol=1))
col <- "site"
colnames(site_inc) <- col



#check for correcct parentage sum in these matrices at the end
sum_par <- as.data.frame(matrix(nrow=nobs(sites), ncol=1))
sum_trios <- as.data.frame(matrix(nrow=nobs(sites), ncol=1))


for(i in 1:length(sites$site)){
    #year <- 2012 #comment one of these two out for either annual (this row) or pooled (next line)
    #year <- c(2012, 2013, 2014, 2015)
    site_eval <- sites$site[i]
    #prepare the parentage matrix input
  
     if(site_eval %in% mums$offs_site){
    site_matF <- mums %>%
        filter(offs_site == site_eval) %>%
        select(offs_site, par1_site, n_mum) %>%
        group_by(offs_site, par1_site) %>%
        summarise(n_mum=(sum(n_mum))) }
    else {site_matF <- as.data.frame(matrix(nrow=1, ncol=3))
          site_matF[1,1] <- site_eval
          site_matF[1,2] <- site_eval
          site_matF[1,3] <- 0
          col <- c("offs_site", "par1_site", "n_mum")
        colnames(site_matF) <- col
         }
    site_matF2 <- suppressWarnings(left_join(sites, site_matF, by=c(site ="par1_site"))) 
        site_matF3 <- site_matF2 %>%
        select(n_mum)

    if(site_eval %in% dads$offs_site){
    site_matM <- dads %>%
        filter(offs_site == site_eval) %>%
        select(offs_site, par2_site, n_dad) %>%
        group_by(offs_site, par2_site) %>%
        summarise(n_dad=(sum(n_dad))) }
    else {site_matM <- as.data.frame(matrix(nrow=1, ncol=3))
          site_matM[1,1] <- site_eval
          site_matM[1,2] <- site_eval
          site_matM[1,3] <- 0
          col <- c("offs_site", "par2_site", "n_dad")
        colnames(site_matM) <- col
         }
    site_matM2 <- suppressWarnings(left_join(sites, site_matM, by=c(site ="par2_site"))) 
        site_matM3 <- site_matM2 %>%
        select(n_dad)
    
#now bind_cols
       #order matters here--> mums need to be in the righthand column, join so sites without matches appear as zero values rather than absent
   site_mat_fin <- bind_cols(site_matM3, site_matF3)
        
    #need to transpose so that it's formatted for migest
    site_mat_T <- as.data.frame(t(site_mat_fin))
    site_mat_T[is.na(site_mat_T)] <- 0
    #generate header
    site_head <- as.data.frame(matrix(nrow=2, ncol=2))
    col <- c("male", "female")
    colnames(site_head) <- col
    
    #proportion sampled
    if(site_eval %in% N_cap3$site){
    site_head_prop <- prop_samp %>% 
    filter(site== site [i]) %>%
    rename(n_or_p="prop_habitat_sampled") }
    else {site_head_prop <- as.data.frame(matrix(nrow=1, ncol=1))
          site_head_prop[1,1] <- 0
          col <- "n_or_p"
        colnames(site_head_prop) <- col
        
    }
    
    
    #number of captured adults
    if(site_eval %in% N_cap3$site){
        site_head_ncap <- N_cap3 %>%
        filter(site== site_eval) %>%
        select(n_male, n_female) %>%
        rename(n_or_p_M="n_male", n_or_p_F="n_female")} 
    else {site_head_ncap <- as.data.frame(matrix(nrow=1, ncol=2))
          site_head_ncap[1,1] <- 0
          site_head_ncap[1,2] <- 0
          col <- c("n_or_p_M", "n_or_p_F")
        colnames(site_head_ncap) <- col
         }
    
    #make a male/female column of proportion sampled
    site_head_prop2 <- site_head_prop %>%
    select(n_or_p) %>%
    mutate(n_or_p_F=n_or_p) %>%
    rename(n_or_p_M=n_or_p)
    
    #bind header rows
    site_head <- bind_rows(site_head_prop2, site_head_ncap)
    
    #add 16 extra columns so the table writes properly
    empt_dat1 <- as.data.frame(matrix(nrow=2, ncol=(nobs(sites)-2)))
    site_head2 <- bind_cols(site_head, empt_dat1)
    
    #rename columns consistently for binding
    col <- c(sites$pop)
    colnames(site_head2) <- col
    colnames(site_mat_T) <- col
    
    #add the parentage matrix
    site_with_parmat <- bind_rows(site_head2, site_mat_T)
    
    #add trio counts
    site_trios <- N_trios %>%
    filter(offs_site==site_eval)
    #filter(year %in% year)


    
    site_trios_with_sites <- suppressWarnings(left_join(sites, site_trios, by=c(site="par1_site")))
    site_trios_with_sites <- site_trios_with_sites %>% select(n_trios)
    site_trios_with_sites_T <- as.data.frame(t(site_trios_with_sites))
    site_trios_with_sites_T[is.na(site_trios_with_sites_T)] <- 0
    
    #bind trios to parentage matrix
    colnames(site_trios_with_sites_T) <- col
    site_trios_with_sites_T[is.na(site_trios_with_sites_T)] <- 0
    site_with_trios <- bind_rows(site_with_parmat, site_trios_with_sites_T)
    
    #add the number of offspring sampled
    if(site_eval %in% N_cap_offs2$site){
    site_offs <- N_cap_offs2 %>%
    filter(site== site_eval) %>%
    select(sampled_fish)}
    else {site_offs <- as.data.frame(matrix(nrow=1, ncol=1))
          site_offs[1,1] <- 0
           }
    
    #build data frame for offspring sampled count
    offs_dat <- as.data.frame(matrix(nrow=1, ncol=nobs(sites)))
    offs_dat[1,1] <- site_offs
    colnames(offs_dat) <- col
    
    #pull together whole data table
    site_with_offs <- bind_rows(site_with_trios, offs_dat) #site_with_offs
    
    #add space line so populations are separeated in the input table
    space <- as.data.frame(matrix(nrow=1, ncol=(nobs(sites))))
    colnames(space) <- col
    site_fin <- bind_rows(site_with_offs, space)
    colnames(site_fin) <- col
    
    #bind together complete tables for each site
    if(rowSums(site_fin[1,], na.rm=TRUE)>0) {empt_dat2 <- bind_rows(empt_dat2, site_fin)}
    
    #check sums of parentage matches
    sum_par[i,] <-(sum(site_mat_T, na.rm=TRUE))
    sum_trios[i,] <-(sum(site_trios_with_sites_T, na.rm=TRUE))
    
    #output sites included
    site_eval <- as.data.frame(site_eval)
    col <- "site"
    colnames(site_eval) <- col
    site_eval$site <- as.character(site_eval$site)
    if(rowSums(site_fin[1,], na.rm=TRUE)>0) {site_inc <- bind_rows(site_inc, site_eval)}

}
       

dim(empt_dat2) #should be npop*7 X npop

sum(dads$n_dad, na.rm= TRUE)
sum(mums$n_mum, na.rm=TRUE)
#should add to 76, plus 14 trios and we have 90

sum(sum_par, na.rm=TRUE) #should match the number of singular parentage matches


sum(sum_trios, na.rm=TRUE) #should match the number of trios

n_pop <- nobs(empt_dat2)/7

n_pop

dim(site_inc)
dim(sites)

#write the site lists for use in processing output
sites_all <- site_inc %>%
    select(site) %>%
    arrange(site) 
sites_all$pop <- seq(from=1, to=n_pop, by=1)
dim(sites_all)
sites_all$site <- gsub(" ", "_", sites_all$site, fixed=TRUE)


#prepare the fille table header
header <- as.data.frame(matrix(nrow=4, ncol=(nobs(site_inc))))
col <- c(sites_all$pop)
colnames(header) <- col
header[1,1] <- (nobs(empt_dat2)/7) #7 rows of input for 1 population 
if((nobs(empt_dat2)/7) < 18){header[2,1] <- 0} else {header[2,1] <- 1} #1 yes all populations sampled, 0 no
header[3,1] <- 0.0 #no prior weight for now
header[4,1] <- "Migest2015.out"
header_fin <- bind_rows(header, space)

migest_in <- rbind(header_fin, empt_dat2)


write.table(sites_all, file="~/migest/annual/2015/20180809_propsamp_corr/input_sites_2015.txt", row.names = FALSE, col.names = TRUE, quote=FALSE)

write.table(migest_in, file="~/migest/annual/2015/20180809_propsamp_corr/MigEst.in", quote=FALSE, na=" ", sep=" ", col.names=FALSE, row.names=FALSE,)
