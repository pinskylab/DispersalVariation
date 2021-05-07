
Packages <- c("MASS","dplyr", "GGally", "broom","bbmle", "cowplot","ggplot2","stringr","fields","tidyr","lubridate", "RColorBrewer", "igraph", "lubridate", "lmtest", "car", "caret", "ROCR",  "lme4", "glmmTMB")  

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

setwd('/local/home/katrinac/parentage/')

load("~/parentage/r_data/sampled_area_each_year.RData")
load("~/parentage/r_data/anems_visited_by_year.RData")
load("~/parentage/r_data/total_sampling_across_years.RData")
load("~/parentage/r_data/cumulative_prop_hab_sampled_by_site.RData")
fish_meta <- readRDS("~/parentage/r_data/fish_meta.rds")



#
#leyte <- conleyte()

"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

# download the fish_Obsfile
#download.file(url = "https://github.com/pinskylab/genomics/blob/master/data/fish-obs.RData?raw=true", destfile = "~/parentage/r_data/fish-obs.RData")

# read in the data
#fish_obs <- readRDS("~/parentage/r_data/fish-obs.RData") 




#three possible ways to look at predictors of dispersal
#A: all sampled site combinations, (19x19) 0/1 disp occurred (binomial), doesn't include info about number of 
#disp events of that trajectory
#B: same as above, but with count data rather than 0/1 *** 07/23/2019, I think this is the best option
#C: all sampled fish combos (nrow=sampled offs, ncol=sampled "adults"). Makes some weird assumptions, like that
#all adults are equally likely to have produced the offspring, better might be all offs anem/par anem combos? 


head(disp_dist)

##read in dispersal distance data
#disp_dist <- read.csv(file="~/parentage/colony2/20200605_1340loci/results//20200624colony_dispersaldirection.csv", header=T, stringsAsFactors=F) #%>%
##fix magbangon spaces
#disp_dist$offs_site <- gsub(". ", ".", disp_dist$offs_site, fixed=TRUE)
#disp_dist$par_site <- gsub(". ", ".", disp_dist$par_site, fixed=TRUE)



#read in parent and offspring files for all potential disp trajectories
parents <- read.table(file="~/parentage/colony2/20190523_1340loci/input/all_parents_corrected.txt", header=T, stringsAsFactors = F)
offs <-read.table(file="~/parentage/colony2/20190523_1340loci/input/all_offspring_corrected.txt", header=T, stringsAsFactors = F)


#attach sample id to the directionality 


##for seasonal
disp_dist <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/tables/RecruitsByMonsoon.csv", header=T, stringsAsFactors=F) %>%
    filter(matched_offs=="Y")
disp_dist$par_site <- gsub(". ", ".", disp_dist$par_site, fixed=TRUE)
disp_dist$offs_site <- gsub(". ", ".", disp_dist$offs_site, fixed=TRUE)

#read in directionality data, don't need the below for annual
sites_ns <- read.table(file="~/parentage/text_file/sites_NS.txt", header=T, sep=",", stringsAsFactors = F)
#fix magbangon spaces
sites_ns$site <- gsub(". ", ".", sites_ns$site, fixed=TRUE)

offs_index <- left_join(disp_dist, sites_ns, by=c(offs_site="site")) %>%
    rename(offs_index="index")
#don't need below join for interannual
par_directionality <- left_join(offs_index, sites_ns, by=c(par_site="site")) %>%
    rename(par_index="index") %>%
    mutate(direction=ifelse(par_index > offs_index, "south", ifelse(par_index < offs_index, "north", "self"))) #each offspring and parent site has an index number 1-19, North to South. Add in a column for dirrectionality



head(disp_dist)

disp_dist_sum <- par_directionality %>%
    group_by(direction) %>%
    summarise(n_obs=n())
disp_dist_sum

disp_dist_sum <- par_directionality %>%
    group_by(season) %>%
    mutate(n_obs=n()) %>%
    mutate(North = sum(direction=="north")/n_obs) %>%
    mutate(South = sum(direction=="south")/n_obs) %>%
    mutate(Self = sum(direction=="self")/n_obs) %>%
    distinct(season, North, South, Self) %>%
    gather(direction, proportion, 2:4) %>%
    ungroup()

disp_dist_sum

disp_dist_sum <- par_directionality %>%
    group_by(season) %>%
    mutate(n_obs=n()) %>%
    mutate(North = sum(direction=="north")/n_obs) %>%
    mutate(South = sum(direction=="south")/n_obs) %>%
    mutate(Self = sum(direction=="self")/n_obs) %>%
    distinct(season, North, South, Self) %>%
    gather(direction, proportion, 2:4) %>%
    ungroup()

directionality_plot <- ggplot(data=disp_dist_sum, aes(x=season, y=proportion, fill=direction, color=direction)) +
    #geom_point(size=8) +
    geom_bar(stat = "identity") +
    annotate("text", x ="NEM", y = .95, vjust=1, hjust=.5, label = "n=11", size=2, family="Helvetica") +
    annotate("text", x ="SWM", y = .95, vjust=1, hjust=.5, label = "n=35", size=2, family="Helvetica") +
    #coord_flip() +
    #scale_shape_manual(name="Direction", labels=c("North", "Self", "South"), values=c(1, 6, 11, 16)) +
    scale_color_manual(name="Direction", labels=c("North", "Self", "South"), values=c("darkseagreen", "orange", "deepskyblue2"))+
    scale_fill_manual(name="Direction", labels=c("North", "Self", "South"), values=c("darkseagreen", "orange", "deepskyblue2"))+
    ylab("Proportion dispersed") +
    xlab("Season") +
    theme(axis.title.y=element_text(size=8, color="black", family = "Helvetica"),
    axis.title.x=element_text(size=8, color="black", family = "Helvetica"),
    axis.text.x = element_text(size=6, color="black",  hjust=1, vjust=1, family = "Helvetica"),
    axis.text.y = element_text(size=6, color="black", family = "Helvetica")) +
    theme(panel.border = element_rect(fill=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.line = element_line(colour = "black"), 
    plot.margin=unit(c(0.2,0.2,0.2,2),"cm"))+
    #theme_linedraw()+
    #labs(shape='Direction') +
    expand_limits(y = c(0, 1.02)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.position="bottom", legend.justification = c(.45, 1), legend.direction="horizontal", legend.key = element_blank(), legend.title= element_blank(), legend.text= element_text(size=6, family = "Helvetica")) #+
 #   scale_colour_manual(values = c("2012" = "darkgoldenrod1", "2013" = "darkseagreen4", "2014"="darkorange2", "2015"="deepskyblue4"))


directionality_plot
ggsave(filename="disp_direction_pub_season.pdf", plot=directionality_plot,width=30, height=70, units="mm", path="~/parentage/colony2/20200605_1340loci/results/")


##make a figure with pie charts indicating proportion sr/unassigned
kernels <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/tables/kernel_fitting_summary.csv", header=T, as.is=T)
annualSR_unassigned <- kernels %>%
    filter(Year !="2012-2018") %>%
    mutate(n_offs_unassigned=NumOffsSampled-NumParentageMatches) %>%
    select(Year, NumParentageMatches, n_offs_unassigned, NumOffsSampled) %>%
    gather(2:3, key=type, value=number) 


kernels

(disp_dist %>%
    group_by(year) %>%
    summarise(n_matches=n()))

#simple plots of proportion dispersed in each direction
disp_dist_sum <- disp_dist %>%
    mutate(year=as.character(year)) %>%
    group_by(year) %>%
    #group_by(season) %>%
    mutate(n_obs=n()) %>%
    mutate(North = sum(direction=="north")/n_obs) %>%
    mutate(South = sum(direction=="south")/n_obs) %>%
    mutate(Self = sum(direction=="self")/n_obs) %>%
    distinct(year, North, South, Self) %>%
    gather(direction, proportion, 2:4) %>%
    ungroup()

disp_dist_sum$year <-factor(disp_dist_sum$year, levels=c("2012","2013", "2014", "2015", "2016", "2017", "2018"))


#plot the proportions
directionality_plot <- ggplot(data=disp_dist_sum, aes(x=year, y=proportion, fill=direction, color=direction)) +
    #geom_point(size=8) +
    geom_bar(stat = "identity") +
    annotate("text", x ="2012", y = .95, vjust=1, hjust=.5, label = "n=3", size=2, family="Helvetica") +
    annotate("text", x ="2013", y = .95, vjust=1, hjust=.5, label = "n=21", size=2, family="Helvetica") +
    annotate("text", x ="2014", y = .95, vjust=1, hjust=.5, label = "n=13", size=2, family="Helvetica") +
    annotate("text", x ="2015", y = .95, vjust=1, hjust=.5, label = "n=11", size=2, family="Helvetica") +
    annotate("text", x ="2016", y = .95, vjust=1, hjust=.5, label = "n=6", size=2, family="Helvetica") +
    annotate("text", x ="2017", y = .95, vjust=1, hjust=.5, label = "n=13", size=2, family="Helvetica") +
    annotate("text", x ="2018", y = .95, vjust=1, hjust=.5, label = "n=4", size=2, family="Helvetica") +
    #coord_flip() +
    #scale_shape_manual(name="Direction", labels=c("North", "Self", "South"), values=c(1, 6, 11, 16)) +
    scale_color_manual(name="Direction", labels=c("North", "Self", "South"), values=c("darkseagreen", "orange", "deepskyblue2"))+
    scale_fill_manual(name="Direction", labels=c("North", "Self", "South"), values=c("darkseagreen", "orange", "deepskyblue2"))+
    ylab("Proportion dispersed") +
    xlab("Year") +
    theme(axis.title.y=element_text(size=8, color="black", family = "Helvetica"),
    axis.title.x=element_text(size=8, color="black", family = "Helvetica"),
    axis.text.x = element_text(size=6, color="black",  hjust=.5, vjust=1, family = "Helvetica"),
    axis.text.y = element_text(size=6, color="black", family = "Helvetica")) +
    theme(panel.border = element_rect(fill=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.line = element_line(colour = "black"), 
    plot.margin=unit(c(0.2,0.2,0.2,2),"cm"))+
    #theme_linedraw()+
    #labs(shape='Direction') +
    expand_limits(y = c(0, 1.02)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.position="bottom", legend.justification = c(.45, 1), legend.direction="horizontal", legend.key = element_blank(), legend.title= element_blank(), legend.text= element_text(size=6, family = "Helvetica")) #+
 #   scale_colour_manual(values = c("2012" = "darkgoldenrod1", "2013" = "darkseagreen4", "2014"="darkorange2", "2015"="deepskyblue4"))


directionality_plot
ggsave(filename="disp_direction_pub.eps", plot=directionality_plot,width=83, height=70, units="mm", path="~/parentage/colony2/20200605_1340loci/results/")




mycols <- c("darkgray", "black")

pie18 <- ggplot((annualSR_unassigned %>%
    filter(Year=="2018")), aes(x = "", y = number, fill=type), colors=type) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0)+
    #geom_text(aes(y = number-(number/2), label = number), color="white", size=12, family="Helvetica")+
    scale_fill_manual(name=NULL, labels=c("Unassigned juveniles", "Assigned juveniles"), values = mycols) +
    scale_color_manual(values = mycols) +
    theme_void()+
    theme(legend.position="top", legend.margin=margin(0,0,0,0), legend.key = element_blank(), legend.spacing.x = unit(1.0, 'cm'), legend.spacing. = unit(.5, 'cm'),
    legend.box.margin=margin(c(0,0,0,0)), legend.title= element_text(size=12, family = "Helvetica"), legend.text= element_text(size=12, family = "Helvetica")) 
    #theme(legend.position="none")
#only a legend the bottom pie

pie18

pierow <- plot_grid( pie12 + theme(legend.position="none"),
        pie13 + theme(legend.position="none"),
        pie14 + theme(legend.position="none"),
        pie15 + theme(legend.position="none"),
        pie16 + theme(legend.position="none"),
        pie17 + theme(legend.position="none"),
        pie18 + theme(legend.position="none"),
        align = 'h',
        hjust = -1,
        nrow = 1
           )

#ggsave(filename="pie18.png", plot=pie18, path="~/parentage/colony2/20190523_1340loci/results/")

legend <- get_legend(pie12 + theme(legend.position="top",legend.margin=margin(0,0,0,0)))

# add the legend to the row we made earlier. Give it one-third of the width
# of one plot (via rel_widths).
pie_fin <- suppressWarnings(plot_grid(legend, pierow, directionality_plot, ncol=1,  rel_heights = c(.3, .5, 3), rel_widths = c(.3, .5, 2)))
suppressWarnings(pie_fin)
ggsave(filename="directionality_pie.pdf", plot=pie_fin, path="~/parentage/colony2/20200605_1340loci/results/")



dim(annualSR_unassigned)
dim(disp_dist)
dim(disp_dist_obs_beta)
dim(disp_dist_obs_beta2)

#read in directionality data
sites_ns <- read.table(file="~/parentage/text_file/sites_NS.txt", header=T, sep=",", stringsAsFactors = F)
#fix magbangon spaces
sites_ns$site <- gsub(". ", ".", sites_ns$site, fixed=TRUE)

#read in centroids for distance calculation
centroids <- read.csv("~/parentage/kernel_fitting/1340_loci/input/site_centroids.csv", header=TRUE, stringsAsFactors = F) %>%
    filter(site %in% sites_ns$site)
#account for sites sampled in each year

#read in demography estimates, correct for NaN/Inf values
prop_samp <- cumulative_prop_hab_sampled_by_site %>%
    filter(time_frame=="2012-2018") %>% #for seasons
    mutate(total_possible_sample_anems = ifelse(site=="Caridad Proper", 4, total_possible_sample_anems) ) %>%
    mutate(total_prop_hab_sampled_anems_tidied= ifelse(site=="Caridad Proper" & total_anems_sampled==4, 1, total_prop_hab_sampled_anems_tidied) ) %>%
    mutate(total_possible_sample_anems = ifelse(site=="Sitio Lonas", total_anems_sampled, total_possible_sample_anems) ) %>%
    mutate(total_prop_hab_sampled_anems_tidied= ifelse(site=="Sitio Lonas", 1, total_prop_hab_sampled_anems_tidied) ) %>%
    select(site, end_year, total_prop_hab_sampled_anems_tidied, total_possible_sample_anems) %>%
    rename(prop_habitat_sampled = "total_prop_hab_sampled_anems_tidied", year="end_year", total_anems= "total_possible_sample_anems") 

prop_samp$prop_habitat_sampled[is.nan(prop_samp$prop_habitat_sampled)] <- 0
prop_samp$site <- gsub(". ", ".", prop_samp$site, fixed=TRUE)


#join with centroids to get all distance

annual_sampled_dists <- left_join(centroids, prop_samp, by="site")


#
##calculate the distance from all potential parents and all potential offspring
all_possible_dists <- as.data.frame(rdist.earth(as.matrix(centroids[,c('lon', 'lat')]), as.matrix(centroids[,c('lon', 'lat')]), miles=FALSE, R=6371))
#
##attach the sample_ids to each distance, so you can also get site and year
colnames(all_possible_dists) <- centroids$site
all_possible_dists$site_i <- centroids$site
#
##gather into tidy 
all_possible_dists_tidy <- all_possible_dists %>%
    select(site_i, everything()) %>%
    gather(2:20, key=site_j, value=dist) 
nrow(all_possible_dists_tidy)
19*19 #check that these two are the same

#add back in the data on year and proportion sampled
dists_meta <- left_join(all_possible_dists_tidy, prop_samp, by=c(site_i="site")) %>%
                    rename(prop_samp_i="prop_habitat_sampled", total_anems_i = "total_anems") %>%
                    select(-year) #only for seasons
dists_meta2 <- left_join(dists_meta, prop_samp, by=c(site_j="site")) %>% #join by year too for annual
                    rename(prop_samp_j="prop_habitat_sampled", total_anems_j = "total_anems") %>% #for seasons get rid of the years
                    distinct(site_i, site_j, .keep_all=T)  %>%#comment out for year
                    select(-year) #only for seasons

#add in directionality

#site_i==origin, site_j==destination
site_i_index <- left_join(dists_meta2, sites_ns, by=c(site_i ="site")) %>%
    rename(site_i_index="index")

site_directionality <- left_join(site_i_index, sites_ns, by=c(site_j="site")) %>%
    rename(site_j_index="index") %>%
    mutate(direction=ifelse(site_j_index > site_i_index, "south", ifelse(site_j_index < site_i_index, "north", "self"))) #%>%  #each offspring and parent site has an index number 1-19, North to South. Add in a column for dirrectionality
    #filter(prop_samp_j !=0)#remove any site-site rows where we didn't sample site_j (offspring site) we couldn't observe these dispersal trajectories
#for seasons, need 2 site_directionality, one for NEM, one for SWM
site_directionality1 <- site_directionality %>%
    mutate(season="NEM")
site_directionality2 <- site_directionality %>%
    mutate(season="SWM")
site_directionality <- bind_rows(site_directionality1, site_directionality2)
#finally, join in the data on observed trajectories
disp_dist_obs_beta <- disp_dist %>%
    #distinct(fish_indiv, .keep_all=T) %>%
    select(offs_site, par_site, season) %>%
    #mutate(year=as.character(year)) %>%
    group_by(season, offs_site, par_site) %>%
    summarise(n_obs=n()) %>%
    ungroup()

#disp_dist_obs_beta2 <- left_join(disp_dist_obs_beta, annualSR_unassigned, by=c(year="Year")) %>%
#    select(-type) %>%
#    #group_by(par_site, offs_site, year) %>%
#    distinct(season, par_site,offs_site, .keep_all = T)

disp_dist_obs <- full_join(site_directionality, disp_dist_obs_beta, by=c(site_j="offs_site", site_i="par_site", "season")) %>%     #join by year for annual
    mutate(n_obs=ifelse(is.na(n_obs), 0, n_obs)) %>% 
    distinct(site_i, site_j, season, .keep_all = T) %>%
    #group_by(season) %>%
    ##mutate(prop_disp = n_obs/NumOffsSampled) %>%
    ##mutate(prop_disp=ifelse(is.na(prop_disp), 0, prop_disp)) %>%
    mutate(bi_disp = ifelse(n_obs>0, 1, 0)) %>%
    mutate(bi_disp=as.factor(bi_disp)) %>%
    #ungroup() %>%
    #mutate(dist_direction=ifelse(direction=="north", 1*dist, ifelse(direction=="south", -1*dist, 0))) %>% # code direction as a negative or positive sign on the distance value, self is 0
    mutate(bi_direction=ifelse(direction=="north", 1, ifelse(direction=="south", -1, 0))) #%>%
    #filter(!is.na(season))#for seasonal, omit larvae who weren't assigned to a season

sum(disp_dist_obs$n_obs) #should be 71 for annual, 46 for seasonal
## add "season" for seasonal
nrow(disp_dist_obs)

361*7
361*2

head(disp_dist_obs)

#annual
full <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_i+ total_anems_j +dist+ bi_direction*year, data=disp_dist_obs, family="binomial")
no_int <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_i+ total_anems_j + dist+ bi_direction +year, data=disp_dist_obs, family="binomial")
no_year <- glm(bi_disp ~ prop_samp_j + prop_samp_i +dist+bi_direction, data=disp_dist_obs, family="binomial")
no_dist_direction <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_i+ total_anems_j + bi_direction, data=disp_dist_obs, family="binomial")
no_prop_sampi <- glm(bi_disp ~ prop_samp_j + + total_anems_i+ total_anems_j + dist+bi_direction +year, data=disp_dist_obs, family="binomial")
no_prop_sampj <- glm(bi_disp ~ prop_samp_i + total_anems_i+ total_anems_j + dist+bi_direction +year, data=disp_dist_obs, family="binomial")
no_anemsi <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_j + dist+bi_direction +year, data=disp_dist_obs, family="binomial")
no_anemsj <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_i + dist+bi_direction +year, data=disp_dist_obs, family="binomial")
no_prop_samp <- glm(bi_disp ~dist+bi_direction + total_anems_i+ total_anems_j +year, data=disp_dist_obs, family="binomial")
no_prop_samp_with_int <- glm(bi_disp ~ total_anems_i+ total_anems_j + dist+bi_direction*year, data=disp_dist_obs, family="binomial")

null <- glm(bi_disp ~dist, data=disp_dist_obs, family="binomial")




summary(no_prop_samp)

directionality_results <- tidy(no_prop_samp)
#write.csv(directionality_results, file="~/parentage/colony2/20200605_1340loci/results/DispLogRegressionModSummaryYear.csv",row.names=F,  quote=F)


aic <- as.data.frame(AIC(null, no_year, no_dist_direction, no_prop_samp, no_prop_sampi, no_prop_sampj, no_anemsi, no_anemsj, no_int, no_prop_samp_with_int, full))
aic$model <- row.names(aic)
aic <- aic %>% arrange(AIC)

#write.csv(aic, file="~/parentage/colony2/20200605_1340loci/results/DispLogRegressionModAICAnnual.csv",row.names=F,  quote=F)


aic %>% arrange(AIC)


summary(no_prop_samp)

directionality_results <- tidy(no_prop_samp)
#write.csv(directionality_results, file="~/parentage/colony2/20200605_1340loci/results/DispLogRegressionModSummaryAnnual.csv",row.names=F,  quote=F)






nrow(disp_dist_obs)

#seasonal
full <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_i+ total_anems_j +dist+ bi_direction*season, data=disp_dist_obs, family="binomial")
no_int <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_i+ total_anems_j + dist+ bi_direction +season, data=disp_dist_obs, family="binomial")
no_season <- glm(bi_disp ~ prop_samp_j + prop_samp_i +dist+bi_direction, data=disp_dist_obs, family="binomial")
no_dist_direction <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_i+ total_anems_j + bi_direction, data=disp_dist_obs, family="binomial")
no_prop_sampi <- glm(bi_disp ~ prop_samp_j + + total_anems_i+ total_anems_j + dist+bi_direction +season, data=disp_dist_obs, family="binomial")
no_prop_sampj <- glm(bi_disp ~ prop_samp_i + total_anems_i+ total_anems_j + dist+bi_direction +season, data=disp_dist_obs, family="binomial")
no_anemsi <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_j + dist+bi_direction +season, data=disp_dist_obs, family="binomial")
no_anemsj <- glm(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_i + dist+bi_direction +season, data=disp_dist_obs, family="binomial")
no_prop_samp <- glm(bi_disp ~dist+bi_direction + total_anems_i+ total_anems_j +season, data=disp_dist_obs, family="binomial")
no_prop_samp_with_int <- glm(bi_disp ~ total_anems_i+ total_anems_j + dist+bi_direction*season, data=disp_dist_obs, family="binomial")
### best model is bi_disp ~dist_direction + total_anems_i+ total_anems_j + season

null <- glm(bi_disp ~dist, data=disp_dist_obs, family="binomial")





aic <- as.data.frame(AIC(null, no_season, no_dist_direction, no_prop_samp, no_prop_sampi, no_prop_sampj, no_anemsi, no_anemsj, no_int, no_prop_samp_with_int, full))
aic$model <- row.names(aic)
aic <- aic %>% arrange(AIC)

#write.csv(aic, file="~/parentage/colony2/20200605_1340loci/results/DispLogRegressionModAICSeason.csv",row.names=F,  quote=F)


aic %>% arrange(AIC)


summary(no_prop_samp)

directionality_results <- tidy(no_prop_samp)
#write.csv(directionality_results, file="~/parentage/colony2/20200605_1340loci/results/DispLogRegressionModSummarySeason.csv",row.names=F,  quote=F)




model.names <- c("null", "no_season", "no_dist_direction", "no_prop_samp", "no_prop_sampi", "no_prop_sampj", "no_anemsi", "no_anemsj", "no_int", "no_prop_samp_with_int", "full")

summ.table <- do.call(rbind, lapply(list(null, no_season, no_dist_direction, no_prop_samp, no_prop_sampi, no_prop_sampj, no_anemsi, no_anemsj, no_int, no_prop_samp_with_int, full), broom::glance))

table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid. Df", "Resid. Dev", "AIC")
reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
reported.table$model <- model.names

full_mod_results <- left_join(aic, reported.table, by="model") %>%
    select(model, everything())
full_mod_results %>% filter(AIC==min(AIC))
full_mod_results <- full_mod_results %>% 
        arrange(AIC)
full_mod_results

write.csv(full_mod_results, file="~/parentage/colony2/20200605_1340loci/results/DispLogRegressionModAIC_fulltableSeason.csv",row.names=F,  quote=F)


disp_dist_obs <- disp_dist_obs %>%
    mutate(total_anems_i= ifelse(total_anems_i==0, 1, total_anems_i)) %>%
    mutate(total_anems_j= ifelse(total_anems_j==0, 1, total_anems_j)) %>%
    mutate(total_anems_i= log(total_anems_i)) %>%
    mutate(total_anems_j= log(total_anems_j))

summary(disp_dist_obs)

#mixed model
full_mix <- glmer(bi_disp ~ prop_samp_j + prop_samp_i+ total_anems_i+ total_anems_j +dist+ (1+direction|year), data=disp_dist_obs, family="binomial")




summary(full_mix)

ranef(full_mix)$year



lrt <- as.data.frame(lrtest(no_year, no_dist_direction, no_prop_samp, no_prop_sampi, no_prop_sampj, no_anemsi, no_anemsj, no_int, no_prop_samp_with_int, full))

#write.csv(lrt, file="~/parentage/colony2/20190523_1340loci/results/final_tables/DispLogRegressionModLRT.csv",row.names=F,  quote=F)

##what if it's binomial
#
#all_pred_disp_bi <- all_pred_disp %>%
#    mutate(prob_disp= ifelse(n_obs>0, 1, 0))
#
#all_pred_mod_bi <- glm(prob_disp ~ dist + year+ prop_samp_j + prop_samp_i  +direction+ direction * year, data=all_pred_disp_bi, family=binomial(link="logit"))
#
#mod_sig_bi <- as.data.frame(anova(all_pred_mod_bi, test="LRT")) 
#mod_sig_bi
#mod_sig

#https://stats.stackexchange.com/questions/86351/interpretation-of-rs-output-for-binomial-regression
summary(no_prop_samp)


#plot the data
ggplot() +
    geom_histogram(data=all_pred_disp %>%
                   filter(year=="2012"), aes(x=n_obs, fill=direction), alpha=0.5, position="dodge") +
    ggtitle("Dispersal Directionality in 2012") +
    ylab("frequency of direction") +
    xlab("number of dispersal events observed")

ggplot() +
    geom_histogram(data=all_pred_disp %>%
                   filter(year=="2013" & n_obs >0), aes(x=n_obs, fill=direction),alpha=0.5, position="dodge") +
    ggtitle("Dispersal Directionality in 2013") +
    ylab("frequency of direction") +
    xlab("number of dispersal events observed")

ggplot() +
    geom_histogram(data=all_pred_disp %>%
                   filter(year=="2014"), aes(x=n_obs, fill=direction), alpha=0.5, position="dodge") +
    ggtitle("Dispersal Directionality in 2014") +
    ylab("frequency of direction") +
    xlab("number of dispersal events observed")

ggplot() +
    geom_histogram(data=all_pred_disp %>%
                   filter(year=="2015" & n_obs >0), aes(x=n_obs, fill=direction), alpha=0.5, position="dodge") +
    ggtitle("Dispersal Directionality in 2015") +
    ylab("frequency of direction") +
    xlab("number of dispersal events observed")

ggplot() +
    geom_histogram(data=all_pred_disp %>%
                   filter(year=="2016" & n_obs >0), aes(x=n_obs, fill=direction), alpha=0.5, position="dodge") +
    ggtitle("Dispersal Directionality in 2016") +
    ylab("frequency of direction") +
    xlab("number of dispersal events observed")

ggplot() +
    geom_histogram(data=all_pred_disp %>%
                   filter(year=="2017" & n_obs >0), aes(x=n_obs, fill=direction), alpha=0.5, position="dodge") +
    ggtitle("Dispersal Directionality in 2017") +
    ylab("frequency of direction") +
    xlab("number of dispersal events observed")

ggplot() +
    geom_histogram(data=all_pred_disp %>%
                   filter(year=="2018" & n_obs >0), aes(x=n_obs, fill=direction), alpha=0.5, position="dodge") +
    ggtitle("Dispersal Directionality in 2018") +
    ylab("frequency of direction") +
    xlab("number of dispersal events observed")

###using option C, all combos between potential offspring and parents input into colony


disp_dist <- read.csv(file="~/parentage/colony2/20190523_1340loci/results/20190620colony_dispersaldirection.csv", header=T, stringsAsFactors=F) %>%
    mutate(disp=1) %>%
    select(offs_samp, par_samp, disp)

#some offs and parents don't have lat/lon, fill them in manually 
centroids <- read.csv("~/parentage/kernel_fitting/site_centroids.csv", header=TRUE)
centroids %>% filter(site=="Poroc San Flower")


#add lat/lon 
offs_loc <- get_latlon(offs$sample_id) %>%
    rename(offs_lat="lat", offs_lon="lon", offs_samp ="sample_id", offs_site="site") %>%
    mutate(offs_lat=ifelse(is.na(offs_lat), 10.7641, offs_lat)) %>%
    mutate(offs_lon=ifelse(is.na(offs_lon), 124.7853, offs_lon))

offs_loc <- left_join(offs, offs_loc, by=c(sample_id="offs_samp"))

p_need_latlon <- offs_loc %>%
    filter(is.na(offs_lat) | is.na(offs_lon))

par_loc <- get_latlon(parents$sample_id) %>%
    rename(par_lat="lat", par_lon="lon", par_samp ="sample_id", par_site="site") %>%
    mutate(par_lat=ifelse(is.na(par_lat), 10.7641, par_lat)) %>%
    mutate(par_lon=ifelse(is.na(par_lon), 124.7853, par_lon))



#calculate the distance from all potential parents and all potential offspring
all_possible_dists <- as.data.frame(rdist.earth(as.matrix(offs_loc[,c('offs_lon', 'offs_lat')]), as.matrix(par_loc[,c('par_lon', 'par_lat')]), miles=FALSE, R=6371))

#attach the sample_ids to each distance, so you can also get site and year
colnames(all_possible_dists) <- parents$sample_id
all_possible_dists$offs_samp <- offs$sample_id

#gather into tidy df and add meta data
all_possible_dists_tidy <- all_possible_dists %>%
    select(offs_samp, everything()) %>%
    gather(2:1720, key=par_samp, value=dist) 
nrow(all_possible_dists_tidy)
#dang that's a lot of possible disp trajectories

#join in date to get year, and
meta <- get_fish() %>%
    mutate(time_date = as.character(str_c(date, 
    anem_obs_time, sep = " "))) %>% mutate(time_date = ymd_hms(time_date)) %>% 
    mutate(time_date = force_tz(time_date, tzone = "Asia/Manila")) %>% 
    mutate(time_date = with_tz(time_date, tzone = "UTC")) %>% 
    mutate(year = year(time_date)) %>% mutate(month = month(time_date)) %>% 
    mutate(day = day(time_date)) %>% mutate(hour = hour(time_date)) %>% 
    mutate(minute = minute(time_date)) %>% select(-time_date, 
    -date)
meta <- left_join(fish_obs, meta, by="sample_id")

all_possible_dists_meta <- inner_join(all_possible_dists_tidy, meta, by=c(offs_samp="sample_id")) %>%
    rename(offs_year="year", offs_site="site", offs_fish_indiv="fish_indiv") %>%
    select(offs_samp, par_samp, offs_year, offs_site, dist, offs_fish_indiv) %>%
    mutate(offs_year=as.factor(offs_year))


all_possible_dists_meta2 <- inner_join(all_possible_dists_meta, meta, by=c(par_samp="sample_id")) %>%
    rename(par_year="year", par_site="site", par_fish_indiv="fish_indiv" ) %>%
    select(offs_samp, par_samp, offs_year, par_year, offs_site, par_site, dist, par_fish_indiv, offs_fish_indiv) %>%
    distinct(offs_fish_indiv, par_fish_indiv, .keep_all = T) %>%
    mutate(par_year=as.factor(par_year))




### Assesing directionality
#read in sites NS
sites_ns <- read.table(file="~/parentage/text_file/sites_NS.txt", header=T, sep=",", stringsAsFactors = F)

#join to all possible distances
offs_index <- left_join(all_possible_dists_meta2, sites_ns, by=c(offs_site ="site")) %>%
    rename(offs_index="index")
par_directionality <- left_join(offs_index, sites_ns, by=c(par_site="site")) %>%
    rename(par_index="index") %>%
    filter(offs_index!=par_index) %>% 
    mutate(direction=ifelse(par_index > offs_index, "south", "north")) #each offspring and parent site has an index number 1-19, North to South. Add in a column for dirrectionality

#add in the self-recruitment values
sr <- left_join(offs_index, sites_ns, by=c(par_site="site")) %>%
    rename(par_index="index") %>%
    filter(offs_index==par_index) %>% 
    mutate(direction="self")
all_possible_disp <- bind_rows(par_directionality, sr) %>% #subset rows because there are way too many points here
    sample_n(100000)

#add fish_obs data to disp dist
disp_dist_offs_obs <- left_join(disp_dist, fish_obs, by=c(offs_samp="sample_id")) %>%
    select(-fish_table_id, -tag_id, -gen_id) %>%
    rename(offs_fish_indiv="fish_indiv")
disp_dist_obs <- left_join(disp_dist_offs_obs, fish_obs, by=c(par_samp="sample_id")) %>%
    select(-fish_table_id, -tag_id, -gen_id) %>%
    rename(par_fish_indiv="fish_indiv")


disp_df <- full_join(all_possible_disp, disp_dist_obs) %>% #add a binary disp 1/no disp 0 score column
     distinct(offs_fish_indiv, par_fish_indiv, .keep_all=T)%>%
     select(-offs_index, -par_index)

disp_df$disp[is.na(disp_df$disp)] <- 0

#add in prop_habitat sampled to be used as a predictor variable
#read in demography estimates
prop_samp <- cumulative_prop_hab_sampled_by_site %>%
    mutate(total_possible_sample_anems = ifelse(site=="Caridad Proper", 4, total_possible_sample_anems) ) %>%
    mutate(total_prop_hab_sampled_anems_tidied= ifelse(site=="Caridad Proper" & total_anems_sampled==4, 1, total_prop_hab_sampled_anems_tidied) ) %>%
    mutate(total_possible_sample_anems = ifelse(site=="Sitio Lonas", total_anems_sampled, total_possible_sample_anems) ) %>%
    mutate(total_prop_hab_sampled_anems_tidied= ifelse(site=="Sitio Lonas", 1, total_prop_hab_sampled_anems_tidied) ) %>%
    select(site, end_year, total_prop_hab_sampled_anems_tidied) %>%
    rename(prop_habitat_sampled = "total_prop_hab_sampled_anems_tidied", year="end_year") %>%
    mutate(year=as.factor(year))


prop_samp$prop_habitat_sampled[is.nan(prop_samp$prop_habitat_sampled)] <- 0

#join to potential offspring
disp_df_beta <- left_join(disp_df, prop_samp, by=c(offs_site="site", offs_year="year")) %>%
    rename(offs_prop_samp="prop_habitat_sampled")

disp_df <- left_join(disp_df_beta, prop_samp, by=c(par_site="site", par_year="year")) %>%
    rename(par_prop_samp="prop_habitat_sampled") %>%
    filter(!is.na(dist) & !is.na(par_site) & !is.na(offs_site) & !is.na(par_year) & !is.na(offs_year)) %>%
    mutate(dist=ifelse(dist==0, 0.001, dist))


#run logistic regressions
dist_mod <- glm(disp ~ dist, data=disp_df, family="binomial")
dist_par_samp_mod <- glm(disp ~ dist + par_prop_samp, data=disp_df, family=binomial(link="logit"))
dist_offs_samp_mod <- glm(disp ~ dist + offs_prop_samp, data=disp_df, family=binomial(link="logit"))
dist_year_mod <- glm(disp ~ dist + offs_year, data=disp_df, family=binomial(link="logit"))
dist_direction_mod <- glm(disp ~ dist + direction, data=disp_df, family=binomial(link="logit"))
dist_offspar_samp_mod <- glm(disp ~ dist + offs_prop_samp+par_prop_samp, data=disp_df, family=binomial(link="logit"))
dist_yeardirection_mod <- glm(disp ~ dist + direction * offs_year, data=disp_df, family=binomial(link="logit"))



#dist_offs_samp_year_mod <- glm(disp ~ dist + offs_prop_samp, data=disp_df, family=binomial(link="logit"))

all_pred_mod <- glm(disp ~ dist + offs_year+ offs_prop_samp +par_prop_samp +offs_year +direction+ direction * offs_year, data=disp_df, family=binomial(link="logit"))

#https://datascienceplus.com/perform-logistic-regression-in-r/

anova(all_pred_mod, test="Chisq")
