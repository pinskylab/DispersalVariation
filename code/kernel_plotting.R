Packages <- c("dplyr","pracma","cubature","data.table", "gridExtra","viridis", "ggsignif", "broom", "ggpubr", "caret","cowplot","ggplot2","fields","bbmle", "dplyr", "tidyr", "lubridate", "RColorBrewer")

invisible(suppressPackageStartupMessages(lapply(Packages, library, character.only = TRUE)))

setwd('/local/home/katrinac/parentage/kernel_fitting/')

load("~/parentage/r_data/site_dist_info.RData")

#final scripts to use
source("~/parentage/kernel_fitting/1340_loci/functions/ll_kt_both_bbmle.R")
source("~/parentage/kernel_fitting/1340_loci/functions/ll_kt_both_grid_search.R")
source("~/parentage/kernel_fitting/1340_loci/functions/ll_kt_both_optim.R")
source("~/parentage/kernel_fitting/1340_loci/functions/GenGausKernInt_sum0.5.R")
source("~/parentage/kernel_fitting/1340_loci/functions/GenGausKernInt_sum1.R")
source("~/parentage/kernel_fitting/1340_loci/functions/GenGausKernPDF.R")

"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

# read in kernel results for plotting
kernels <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/tables/kernel_fitting_summary.csv", header=T)
kernels

profile_16_df <- read.csv(file="~/parentage/kernel_fitting/1340_loci/final_results/likelihood_profiles_grid_search/profile95CI_2016.csv", header=T) 

profile_16_df95 <- profile_16_df %>% 
        #filter(log_like < cutoff_16) %>%
        sample_n(500)


k_16 <- as.numeric(kernels %>%
                     filter(Year=="2016") %>%
                     select(best_k))
theta_16 <- as.numeric(kernels %>%
                     filter(Year=="2016") %>%
                     select(best_theta))

MDD_16 <- as.numeric(kernels %>%
                     filter(Year=="2016") %>%
                     select(MeanDispDist))


max_dist <- 60
possible_dist <- matrix(seq(0,max_dist,0.1))
to_plot <- matrix(nrow=nrow(profile_16_df95), ncol=nrow(possible_dist))

for(i in 1:nrow(profile_16_df95)) {
     
    k <- profile_16_df95$k_eval[i]
    theta <- profile_16_df95$theta_eval[i]
    
    for(j in 1:nrow(possible_dist)) {
    
    dist <- possible_dist[j]
    to_plot[i,j] <- predicted_disp(k=k, theta=theta, d=dist)
    
    }
}

to_plot_df <- as.data.frame(to_plot) 

col <- as.character(possible_dist)
colnames(to_plot_df) <- col
to_plot_df <- bind_cols(profile_16_df95, to_plot_df) %>%    
    group_by(k_eval, theta_eval) 

to_plot_df2 <- to_plot_df %>%
    gather(4:604, key=dist, value=pdf) %>%
    mutate(dist=as.numeric(dist))
dim(to_plot_df2)
to_plot_df2$iter <- paste0(to_plot_df2$k_eval, '_', to_plot_df2$theta_eval) # create an identifier for each sample from the LL surface

##make a df to project the kernel
kernel_16 <- data.frame(dist=seq(0,max_dist,0.1))
k= k_16
theta=theta_16

kernel_16 <- kernel_16 %>%
    mutate(pdf= predicted_disp(k=k, theta=theta, d=dist))




disp_16_pub <- ggplot(data=kernel_16, aes(x=dist, y=pdf))+ 
    geom_line(data=to_plot_df2, aes(x=dist, y=pdf, group=iter), color="gray60", size=.2, alpha=0.2)+
    #geom_area(data=kernel15, aes(x=dist, y=norm_pdf, fill=as.factor(median)), alpha=0.1)+
    annotate("text", x =Inf, y = 0.17, vjust=1, hjust=1, label = "2016", size=6, family="Helvetica")+#8 for publication, 15 for presentation, fontface="bold", col="deepskyblue2") + #for 15 year  fontface="bold", col="deepskyblue2"
    #annotate("text", x =MDD15, y = 0.01, vjust=1, hjust=1, label = "mu", color="deepskyblue2", size=15, family="Helvetica", parse=T)+#, fontface="bold", col="deepskyblue2") + #for 15 year  fontface="bold", col="deepskyblue2"
    #scale_fill_manual(values= c("whitesmoke", "cadetblue3")) +
    #scale_color_manual(values= c("whitesmoke", "cadetblue3")) +
    geom_line(color="black", size=.6) +  
    #geom_vline(xintercept = MDD15, linetype="dotted", color = "deepskyblue2", size=1) +
    theme(axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    #plot.background = element_rect(colour = "black", size = 1),
    axis.text.x = element_text(size=15, color="black", family="Helvetica"),#15 for publication, 20 for presentation #element_text(size=15, color="black", family="Helvetica"),
    axis.text.y =  element_text(size=15, color="black", family="Helvetica"), 
    panel.grid.major = element_blank(),
    plot.margin=unit(c(0,0,0.2,0.2),"cm"),#,
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), #,
    axis.line = element_line(colour = "black")) +
    xlab("Distance (km)") + #ylab("Probability density") + 
    guides(fill=FALSE)+
    scale_x_continuous(limits = c(0,52), expand = c(0, 0), breaks=c(seq(0, 50, 10))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,.18)) 


disp_16_pub
#ggsave(filename="disp12_pub.pdf", plot=disp12_pub, path="~/parentage/kernel_fitting/1240_loci/final_results/kernel_plots" )


#add to plot
all <- plot_grid(disp_12_pub, disp_13_pub, disp_14_pub, disp_15_pub, disp_16_pub, disp_17_pub, disp_18_pub, disp_all_pub, ncol=4, scale=0.8) #, disp_15_pub,disp_16_pub, disp_17_pub, disp_18_pub, 
                 #disp_all_pub,  ncol=1) 
all_ann <- grid.arrange(arrangeGrob(all, bottom=grid::textGrob(label= "Distance (km)",
                    gp= gpar(fontsize=25, family="Helvetica", col="black", vjust=0)), 
                    left=grid::textGrob(label="Dispersal strength", rot=90, 
                    gp= gpar(fontsize=25, family="Helvetica", col="black", hjust=1))))
#png
#ggplot2::ggsave(filename="disp_panel_pub_4x4.png", width=40, height=20, units="cm", plot=all_ann, path="~/parentage/kernel_fitting/1340_loci/final_results/kernel_plots")

