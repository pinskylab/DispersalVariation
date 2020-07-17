library(dplyr)
library(nleqslv)
library(pracma)
library(cubature)
library(ggplot2)
source("~/parentage/kernel_fitting/1340_loci/functions/GenGausKernInt_sum0.5.R") 
source("~/parentage/kernel_fitting/1340_loci/functions/cdf_solve.R") 

kernels <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/tables/kernel_fitting_summary.csv", stringsAsFactors = F, header=T)


for(j in 1:nrow(kernels)){
    
theta_eval <- kernels$best_theta[j]
k_eval <- kernels$best_k[j]

kernels$MedianDispDist[j] <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2) # answer is stored in $x

}

profile12_95 <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/likelihood_profiles_grid_search/profile95CI_2012.csv", header=T) %>%
        filter(k_eval <=10 & k_eval >=-10 & theta_eval <=5)

min12 <- min(profile12_95$log_like)

profile12_weighted <- profile12_95 %>%
        mutate(deviation=min12/log_like) %>%
        sample_n(size=1000, weight=deviation, replacement=T)%>%
        mutate(Year="2012")

profile12_weighted$median <- NA

pb <- txtProgressBar(min = 0, max =1000, style = 3)

for(i in 1:nrow(profile12_weighted)){
    
theta_eval <- profile12_weighted$theta_eval[i]
k_eval <- profile12_weighted$k_eval[i]

profile12_weighted$median[i] <- round(nleqslv(x = 7, fn = cdf_solve)$x, 2) 

setTxtProgressBar(pb, i)

}

close(pb)

#write.csv(profile12_weighted, file="~/parentage/kernel_fitting/1240_loci/final_results/mean_disp_dist/MedPdf2012.csv", row.names=F) 



profileall_weighted$Year <- "2012-18" #to make the year label fit on the graph axis
all_med <- bind_rows(profile12_weighted, profile13_weighted, profile14_weighted, profile15_weighted, profile16_weighted, profile17_weighted, profile18_weighted, profileall_weighted)
kernels$MedianDispDist_CI95_lower <- c(round(min(profile12_weighted$median), 2), round(min(profile13_weighted$median), 2), round(min(profile14_weighted$median), 2), round(min(profile15_weighted$median), 2), round(min(profile16_weighted$median), 2), round(min(profile17_weighted$median), 2), round(min(profile18_weighted$median),2), round(min(profileall_weighted$median), 2))
kernels$MedianDispDist_CI95_upper <- c(round(max(profile12_weighted$median), 2), round(max(profile13_weighted$median), 2), round(max(profile14_weighted$median), 2), round(max(profile15_weighted$median), 2), round(max(profile16_weighted$median), 2), round(max(profile17_weighted$median), 2), round(max(profile18_weighted$median), 2), round(max(profileall_weighted$median), 2))
write.csv(kernels, file="~/parentage/kernel_fitting/1340_loci/final_results/tables/kernel_fitting_summary.csv", row.names=F) 


# plot the median distributions
all_med <- all_med %>%
    arrange(Year)%>%
    mutate(Year=factor(Year, levels=c("2012-18", "2018", "2017", "2016", "2015", "2014", "2013", "2012"))) %>%
    arrange(Year)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette_rev <- rev(cbbPalette)

med_violin <- ggplot(data=all_med, aes(x=Year, y=median, color=Year, fill=Year), alpha=0.5) +
    geom_violin() +
    geom_point(data=kernels, aes(x=Year, y=MedianDispDist), fill="snow",color="darkgray", shape=21) +
    coord_flip() +
    #scale_x_continuous(limits = c(0, 100), expand=c(0,0)) +
    scale_y_continuous(limits = c(0, 83), expand = c(0,0))+
    theme(panel.grid.major = element_blank(),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), #,
    axis.line = element_line(colour = "black")) +
    xlab("Year") + 
    ylab("Median dispersal distance (km)") +
    theme(axis.text.x = element_text(size=8, color="black", family="Helvetica"),#15 for publication, 20 for presentation #element_text(size=15, color="black", family="Helvetica"),
    axis.text.y =  element_text(size=8, color="black", family="Helvetica"),
    axis.title.y =  element_text(size=10, color="black", family="Helvetica"), 
    axis.title.x =  element_text(size=10, color="black", family="Helvetica"),    
    legend.position = "none") + 
    scale_colour_manual(values=cbbPalette_rev)+
    scale_fill_manual(values=cbbPalette_rev)# +  
    #scale_y_continuous(breaks = 1:6, labels = c(0:200,"break",:8))
    

med_violin
ggplot2::ggsave(filename="All95CImedianDispDistanceViolin_pub.pdf",  plot=med_violin, width=83, height=70, units="mm", path="~/parentage/kernel_fitting/1340_loci/final_results/mean_disp_dist/")



