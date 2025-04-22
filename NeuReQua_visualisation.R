# VISUALISATION OF NeuReQua METRICS
#april 2025, Zoé Darrasse (zoe.darrasse@cnrs.fr)


#library to load:  (if they need to be installed: install.packages("xx"))

library(ggplot2)
library(viridis)
library(dplyr)
library(stringr)
library(rstatix)
library(crayon)
library(grid)
library(plotly)


#import table from NeuReQua in the .csv format 
tbl <- read.csv("C:/Users/darrasse/Documents/NeuReQua/All_database-tracking-SU.csv", header = TRUE, sep = ",")

#path where figures will be saved 
savepath <- "C:/Users/darrasse/Documents/NeuReQua/stats/"

#------------ D A T A   V I S U A L I S A T I O N --------------

#change the scale of toe filtered RMS from V to µV
tbl$RMS_filter <- tbl$RMS_filter * (10^6)

#sort the electrodes so all 'd' (droite = right) electrode come before the 'v' (vänster = left)
#this will also order vr1 ..vr9 vr10 instead of vr1 vr10 vr2..vr9
tbl$electrodes <- factor(tbl$electrodes, levels = unique(tbl$electrodes[order(grepl("^v", tbl$electrodes))]))

#adds a column with the name of the corresponding macro electrode for each micro filament
#this can be use to group electrodes together based on the macro electrode they are on
tbl$gp_macro <- str_extract(tbl$electrodes, "^[a-zA-Z]+")

#extract the first letter so that we can group by left and right switch V to left and d to right 
tbl$lat <- substr(tbl$electrode, 1, 1)
tbl <- tbl %>% mutate(lat = case_when(lat == "v" ~ "left", lat == "d" ~ "right"))

#adds a colums thta specify the tetrode the microfilment is one (ex tetrode 1)
for (subject in unique(tbl$sub)) {
  sub_subset <- subset(tbl, sub == subject)
  for(task in unique(sub_subset$run)) {
    run_subset <- subset(sub_subset, run == task)
    for (macro in unique(run_subset$gp_macro)) {
        len <- nrow(subset(run_subset, gp_macro == macro))
        if (len %% 4 == 0) {
          nb_tetrode <- len/4
          tet <- rep(1:nb_tetrode, each=4)
          for (i in 1:len){
            tbl[tbl$sub == subject & tbl$run == task & tbl$gp_macro == macro, "gp_tetrode"] <- paste('tetrode', as.character(tet), sep = ' ')} 
        }else{
          for (i in 1:len){
            tbl[tbl$sub == subject & tbl$run == task & tbl$gp_macro == macro, "gp_tetrode"] <- paste('missing tetrode')} 
}}}}
 
#define palette for "run" colours
colours_run <- c("#9D1642","#D14E72","#E7A04B","#9CCC65","#66C2A5","#3288BD","#01579B","#29258B","#6F00A8","#B6308B","#D14E72","#E7A04B","#E9C13D","#F0F921")

#We will use the "plasma" palette prom viridis
#this crops the values too yellow out of the plasma palette
full_palette <- viridis::plasma(length(unique(tbl$sub))+3)
cropped_palette <- full_palette[1:length(unique(tbl$sub))]

#Check the colours : scales::show_col(cropped_palette)



# A. visualisation of all metrics from NeuReQua -------------------

#extract the name of all metrics (we remove the column names that are not metrics)
metrics <- setdiff(colnames(tbl), c("sub", "run", "electrodes","Nb_SU","SU","region",
                                    "gp_macro", "gp_tetrode","Day","lat","tot_sub_su"))

#define the name of varianble and its units
units <- c("RMS (V)","Filtered RMS (300–3000 Hz) (µV)","Variance (V^2)","normalised variance","Microwires correlation",
           "Deviation (V)","Kurtosis","SNR","% of Artefacts in signal","Hurst componant")



# 1. plot each metricsacross all participants ----

#set up titles for each plot
tt <-"from each microwire & task, across subjects (from a random 5-min segment per task)"
titles_bp <- c(paste("RMS",tt),
               paste("Filtered RMS (300–3000 Hz)",tt),
               paste("Variance",tt),
               paste("Normalised variance",tt),
               "Microwires correlation within a tetrodes for each microwire & task, across subjects (from a random 5-min segment per task)",
               paste("Deviation",tt),
               paste("Kurtosis",tt),
               paste( "Signal Noise ratio (SNR)",tt),
               "Percentage of signal containing artefacts for each task, across subjects (from a random 5-min segment per task)",
               paste("Hurst component" ))

tt <- "averaged per microwire, across subjects (from a random 5-min segment per task)"
titles_avg <- c(paste("RMS",tt),
                 paste("Filtered RMS (300–3000 Hz)",tt),
                 paste("Variance",tt),
                 paste("Normalised variance",tt),
                 "Microwires correlation within a tetrodes averaged per microwire, across subjects (from a random 5-min segment per task)",
                 paste("Deviation",tt),
                 paste("Kurtosis",tt),
                 paste( "Signal Noise ratio (SNR)",tt),
                 "Pourcentage of signal containing artefacts across subjects from a random 5 minute segment per task",
                 paste("Hurst component" ))


#this loop will create all plot across participant automatically for each metrics 
for (i in 1:length(metrics)) {
  
  # calculate inter quartil space to remove outliers from visualization 
  iq <- tbl %>%
        group_by(sub) %>%
        summarize(IQR_value = IQR(.data[[metrics[i]]], na.rm = TRUE))
  
  max.iq <- 5*max(iq$IQR_value, na.rm = TRUE)
  
  #calculate the mean
  average.metric <- mean(tbl[[metrics[i]]],na.rm = TRUE)
  
  #set the ylimits (used to crop at average +- 5 x the biggest interquartil space )
  
  if (metrics[i] == "tetrode_cor"){  #for correlation values
    ylimits <- c(-1.1,  1.1)
  }else if (metrics[i] == "Artefact") {  #for % values
    ylimits <- c(0,  100)
  }else if (min(tbl[[metrics[i]]], na.rm = TRUE) >=0) { #for positive values
    ylimits <- c(0, max(tbl[[metrics[i]]]))
    if (abs(average.metric + max.iq) > abs(max(tbl[[metrics[i]]], na.rm = TRUE))){
      ylimits[2] <- max(tbl[[metrics[i]]])
    }else{ylimits[2] <- average.metric + max.iq}
  }else { ylimits <- c(0,0) #chooses either the +- interquatil definie eearlier or the min/max value 
  if (abs(average.metric + max.iq) > abs(max(tbl[[metrics[i]]], na.rm = TRUE))){
    ylimits[2] <- max(tbl[[metrics[i]]])
  }else{ylimits[2] <- average.metric + max.iq}
  if (abs(average.metric - max.iq) > abs(min(tbl[[metrics[i]]], na.rm = TRUE))){ 
    ylimits[1] <- min(tbl[[metrics[i]]])
  }else{ylimits[1] <- average.metric - max.iq}
  }
 
  #calculate average metrics[i] for each filament (used in sc_microavg)
  avg_tbl <- tbl %>%
    group_by(sub, electrodes) %>%
    summarize(avg_metric = mean(.data[[metrics[i]]], na.rm = TRUE), .groups = "drop")
  
  #in avg tbl, add a column for left right (for colouring purposes) and factor by electrodes
  avg_tbl$lat <- substr(avg_tbl$electrodes, 1, 1)
  avg_tbl <- avg_tbl %>% mutate(lat = case_when(lat == "v" ~ "left", lat == "d" ~ "right"))
  avg_tbl$electrodes <- factor(avg_tbl$electrodes, levels = unique(avg_tbl$electrodes[order(grepl("^v", avg_tbl$electrodes))]))

 
  if (metrics[i] != 'RMS_filter'& metrics[i] != 'Artefact'){
  #a. boxplot of the metric across all participants
  
  bp_overall <- ggplot(tbl, aes(x = factor(sub), y = .data[[metrics[i]]],fill= factor(sub), color = factor(sub))) +
                geom_boxplot(alpha = 0.45, size = 1) + 
                geom_smooth(method = "lm", aes(x = as.numeric(factor(sub)), group = 1), 
                            linetype = "dashed", alpha = 0.2, se = TRUE, color = "#5D4037") +  
                scale_fill_manual(values = cropped_palette) +
                scale_color_manual(values = cropped_palette) +
                labs(title = titles_bp[i], x = "Subject", y = units[i]) +
                theme_minimal() +
                theme(plot.title = element_text(size = 15, face = "bold.italic")) + 
                theme(legend.position = "none")+
                coord_cartesian(ylim = ylimits) 

  
  #b. scattered plot for each metric for every micro filament for each sub colored by laterality of the electrode
  sp_micro <- ggplot(tbl, aes(x = factor(sub), y = .data[[metrics[i]]],color = factor(lat))) +
              geom_point(size = 3,alpha = 0.3) +
              #scale_color_viridis_d(option = "plasma") +
              scale_color_manual(values = c("#DA5A6AFF","#FDC229FF"),name = NULL,
                                 guide = guide_legend(override.aes = list(alpha = 1))) +
              geom_smooth(method = "lm", aes(x = as.numeric(factor(sub)), group = 1),  #add the lm 
                          linetype = "dashed", alpha = 0.2, se = TRUE, color = "#5D4037") +  
              labs(title = titles_bp[i], x = "Subject", y =  units[i]) +
              theme_minimal() +
              theme(plot.title = element_text(size = 15, face = "bold.italic")) + 
              coord_cartesian(ylim = ylimits) 

  
  #c. scattered plot for  metric averaged per micro filament for each 
  sp_microavg <- ggplot(avg_tbl, aes(x = factor(sub), y = avg_metric,color = factor(lat))) +
               geom_point(size = 3, alpha = 0.2) +
               #scale_color_viridis_d(option = "plasma") +
               scale_color_manual(values = c("#DA5A6AFF","#FDC229FF"),name = NULL,
                                  guide = guide_legend(override.aes = list(alpha = 1))) +
               geom_smooth(method = "lm", aes(x = as.numeric(factor(sub)), group = 1),  #add the lm 
                           linetype = "dashed", alpha = 0.2, se = TRUE, color = "#5D4037") +  
               labs(title = titles_avg[i], x = "Subject", y =  units[i]) +
               theme_minimal() +
               theme(plot.title = element_text(size = 15, face = "bold.italic")) + 
               coord_cartesian(ylim = ylimits) 
  }
  
  #for metrics = RMS filter we will use a log a scale on the Y axis for display 
  if (metrics[i] == 'RMS_filter'){
    #adapt the limits for the log scale 
    if (ylimits[1] == 0){ylimits[2]<-log(ylimits[2])} 
    bp_overall <- ggplot(tbl, aes(x = factor(sub), y = log(.data[[metrics[i]]]),fill= factor(sub), color = factor(sub))) +
                  geom_boxplot(alpha = 0.45, size = 1) + 
                  geom_smooth(method = "lm", aes(x = as.numeric(factor(sub)), group = 1), 
                              linetype = "dashed", alpha = 0.2, se = TRUE, color = "#5D4037") +  
                  scale_fill_manual(values = cropped_palette) +
                  scale_color_manual(values = cropped_palette) +
                  #labs(title = titles_bp[i], x = "Subject", y = paste("log(",units[i],")")) +
                  labs(title = titles_bp[i], x = "Subject", y = units[i]) +
                  theme_minimal() +
                  theme(plot.title = element_text(size = 15, face = "bold.italic")) + 
                  theme(legend.position = "none") +
                  coord_cartesian(ylim = ylimits) 
    
    sp_micro <- ggplot(tbl, aes(x = factor(sub), y = log(.data[[metrics[i]]]),color = factor(lat))) +
                geom_point(size = 3,alpha = 0.3) +
                scale_color_manual(values = c("#DA5A6AFF","#FDC229FF"),name = NULL,
                                   guide = guide_legend(override.aes = list(alpha = 1))) +
                geom_smooth(method = "lm", aes(x = as.numeric(factor(sub)), group = 1),  #add the lm 
                            linetype = "dashed", alpha = 0.2, se = TRUE, color = "#5D4037") +  
                labs(title = titles_bp[i], x = "Subject", y =  units[i]) +
                theme_minimal()+
                theme(plot.title = element_text(size = 15, face = "bold.italic"))+ 
                coord_cartesian(ylim = ylimits)   
                
    sp_microavg <- ggplot(avg_tbl, aes(x = factor(sub), y = log(avg_metric),color = factor(lat))) +
                   geom_point(size = 3, alpha = 0.2) +
                   scale_color_manual(values = c("#DA5A6AFF","#FDC229FF"),name = NULL,
                                      guide = guide_legend(override.aes = list(alpha = 1))) +
                   geom_smooth(method = "lm", aes(x = as.numeric(factor(sub)), group = 1),  #add the lm 
                              linetype = "dashed", alpha = 0.2, se = TRUE, color = "#5D4037") +  
                   labs(title = titles_avg[i], x = "Subject", y =  units[i]) +
                   theme_minimal()+
                   theme(plot.title = element_text(size = 15, face = "bold.italic"))+ 
                   coord_cartesian(ylim = ylimits) 
    
    #add a line for filter RMS at 8µV
    if (metrics[i] == 'RMS_filter'){
      bp_overall <- bp_overall + geom_hline(yintercept = log(8), color = "#3288BD", linetype = "dashed") +
                                 scale_y_continuous(breaks = log(c(1,5,8, 10, 25, 50, 75, 100, 125, 150, 175)),
                                                    labels = c("0","5", expression(bold("8")),"10", "25", "50", "", "100", "", "150", ""))

      sp_micro <- sp_micro + geom_hline(yintercept = log(8), color = "#3288BD", linetype = "dashed") +
                             scale_y_continuous(breaks = log(c(1,5,8, 10, 25, 50, 75, 100, 125, 150, 175)),
                                                labels = c("0","5", expression(bold("8")),"10", "25", "50", "", "100", "", "150", ""))
                          
      sp_microavg <- sp_microavg + geom_hline(yintercept = log(8), color = "#3288BD", linetype = "dashed") +
                                   scale_y_continuous(breaks = log(c(1,5,8, 10, 25, 50, 75, 100, 125, 150, 175)),
                                                      labels = c("0","5", expression(bold("8")),"10", "25", "50", "", "100", "", "150", ""))
    }
    
   
    }
  
  if (metrics[i] != 'Artefact'){ 
    #saving
    #if the file already exist it will remove it before it is saved 
    filepath <- paste0(savepath,"overall_", metrics[i],"_boxplot.png")
    if (file.exists(filepath)) {file.remove(filepath) 
      message(black("Existing file deleted: ", filepath))}
    #save boxplot in the "savepath" provided at the top of the script
    ggsave(filepath, plot = bp_overall, width = 15, height = 6)
    
    filepath <- paste0(savepath,"overall_", metrics[i],"_filaments.png")
    if (file.exists(filepath)) {file.remove(filepath) 
      message(black("Existing file deleted: ", filepath))}
    ggsave(filepath, plot = sp_micro, width = 15, height = 6)
    
    filepath <- paste0(savepath,"overall_", metrics[i], "_avg-filaments.png")
    if (file.exists(filepath)) {file.remove(filepath) 
      message(black("Existing file deleted: ", filepath))}
    ggsave(filepath, plot = sp_microavg, width = 15, height = 6)
  }
  
  #Artefact metrics need to be represented differently du to being percentage 
  if (metrics[i] == 'Artefact'){
    sp_run <- ggplot(tbl, aes(x = factor(sub), y = .data[[metrics[i]]],color = factor(run))) +
                geom_point(size = 4,alpha = 0.3,shape =19) +
                scale_color_manual(values = colours_run, name = NULL,
                                   guide = guide_legend(override.aes = list(alpha = 1))) +
                geom_smooth(method = "lm", aes(x = as.numeric(factor(sub)), group = 1),  #add the lm 
                            linetype = "dashed", alpha = 0.2, se = TRUE, color = "#5D4037") +  
                labs(title = titles_bp[i], x = "Subject", y =  units[i]) +
                theme_minimal() + 
                theme(plot.title = element_text(size = 15, face = "bold.italic")) + 
                coord_cartesian(ylim = ylimits) 
    
    hist_overall <- ggplot(avg_tbl, aes(x = factor(sub), y = avg_metric,color = factor(sub),fill = factor(sub))) +
                    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
                    scale_fill_manual(values = cropped_palette) +
                    scale_color_manual(values = cropped_palette) +
                    geom_smooth(method = "lm", aes(x = as.numeric(factor(sub)), group = 1),  #add the lm 
                                linetype = "dashed", alpha = 0.2, se = TRUE, color = "#5D4037") + 
                    labs(title = "Average percentage of signal containing artefacts, across subjects (from a random 5-min segment per task)", x = "Subject", y = "% of Artefacts in signal" ) +
                    theme_minimal()+
                    theme(plot.title = element_text(size = 15, face = "bold.italic")) + 
                    theme(legend.position = "none") +
                    coord_cartesian(ylim = ylimits) 
    
    #saving
    filepath <- paste0(savepath,"overall_", metrics[i],"_scatter.png")
    if (file.exists(filepath)) {file.remove(filepath) 
        message(black("Existing file deleted: ", filepath))}
    ggsave(filepath, plot = sp_run, width = 15, height = 6)
    
    filepath <- paste0(savepath,"overall_", metrics[i],"_histogram.png")
    if (file.exists(filepath)) {file.remove(filepath) 
        message(black("Existing file deleted: ", filepath))}
    ggsave(filepath, plot = hist_overall, width = 15, height = 6)
  }
}

message(green("Overall plots were successfully saved in:\n",savepath))



# 2. for each sub plot a boxplot of the metric per run and a line plot of a the metric for each electrode per run -----


for (i in 1:length(metrics)) {
  for (subject in unique(tbl$sub)) {
    
    #create a subset for each subject, 
    #make sur the electrodes are in alphabetical order 
    sub_subset <- subset(tbl, sub == subject)
    sub_subset$electrodes <- factor(sub_subset$electrodes, levels = unique(sub_subset$electrodes[order(grepl("^v", sub_subset$electrodes))]))
    
    # group day and task together so it is presented in chronoligocal oreder even if you have different tasks 
    sub_subset$daytask <- paste(sub_subset$Day, sub_subset$run)
    sub_subset$macro_tet <- paste(sub_subset$gp_macro, sub_subset$gp_tetrode)
    
    #for graphic representation of tetrode 
    macro_tet <- sub_subset %>%
                  mutate(electrodes = factor(electrodes, levels = unique(electrodes))) %>%  
                  group_by(macro_tet) %>%
                  summarize(xmin = min(as.numeric(electrodes)) - 0.5,
                            xmax = max(as.numeric(electrodes)) + 0.5)
    
   # calculate inter quartil space to remove outliers from visualisation 
    iq <- sub_subset %>%
          group_by(sub) %>%
          summarize(IQR_value = IQR(.data[[metrics[i]]], na.rm = TRUE))
    
    max.iq <- 5*max(iq$IQR_value, na.rm = TRUE) 
    average.metric <- mean(sub_subset[[metrics[i]]],na.rm = TRUE)
    
    #set the ylimits
    if (metrics[i] == 'tretrode_cor'){
      ylimits <- c(0,  1.1)
    }else if (metrics[i] == 'Artefact') { 
      ylimits <- c(0,  100)
    } else if (min(sub_subset[[metrics[i]]], na.rm = TRUE) >=0) {          #for positive values
      ylimits <- c(0,  average.metric + max.iq)
    } else{
      ylimits <- c(average.metric - max.iq,  average.metric + max.iq)
    }
   
  if (metrics[i] != "Artefact"){
     
   #a. boxplot of the metric per run        
    bp_sub <- ggplot(sub_subset, aes(x = factor(daytask), y = .data[[metrics[i]]],fill= factor(run), color = factor(run))) +
              geom_boxplot(alpha = 0.4, size = 1.2) + 
              geom_smooth(method = "lm", aes(x = as.numeric(factor(daytask)), group = 1),  #add the lm 
                          linetype = "dashed", alpha = 0.2, se = TRUE, color = "darkgrey") +  
              scale_fill_viridis_d(option = "plasma") +  
              scale_color_viridis_d(option = "plasma") +
              labs(title = paste(metrics[i], "per task for subject", subject), x = "Run", y = paste(metrics[i], units[i])) +
              theme_minimal() +
              theme(legend.position = "none") +
              theme(plot.title = element_text(size = 15, face = "bold.italic")) +
              coord_cartesian(ylim = ylimits)
            
    
    #b. line plot of a the metric for each electrode per run
    lp_sub <- ggplot(sub_subset, aes(x = factor(electrodes), y = .data[[metrics[i]]], group = factor(run),color = factor(run))) +
              geom_rect(data = macro_tet, 
                        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = factor(macro_tet)), 
                        inherit.aes = FALSE, alpha = 0.2) +
              geom_line(alpha = 1, size = 1.2) +
              geom_point(size = 3,alpha = 0.8) +
              scale_color_manual(values = colours_run)  +
              #scale_color_brewer(palette = "Set1") +
              labs(title = paste(metrics[i],"of eachtask across microwires for subject", subject), x = "run", y = paste(metrics[i], units[i])) +
              theme_minimal()+
              theme(plot.title = element_text(size = 15, face = "bold.italic")) 
    
    
    #add a line at 8µV in filterd RMS lineplot 
    if (metrics[i] == 'RMS_filter'){
      lp_sub <- lp_sub +
                geom_hline(yintercept = 8, color = "#3288BD", linetype = "dashed") +
                scale_y_continuous(breaks = c(0,5,8, 10, 15,20, 25, 50, 75, 100, 125, 150, 175),
                                   labels = c("0", "5", expression(bold("8")),"10","15","20","25", "50", "", "100", "", "150", ""))
    }
    
    #this line plot as cut line between each tetrode for better reading 
    if (metrics[i] == 'tetrode_cor'){
      lp_sub <- ggplot(sub_subset, aes(x = factor(electrodes), y = .data[[metrics[i]]], group = interaction(run, macro_tet),color = factor(run))) +
        geom_rect(data = macro_tet, 
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = factor(macro_tet)), 
                  inherit.aes = FALSE, alpha = 0.2) +
        geom_line(alpha = 1, size = 1.2) +
        geom_point(size = 3,alpha = 0.8) +
        scale_color_manual(values = colours_run)  +
        labs(title = paste(metrics[i],"of eachtask across microwires for subject", subject), x = "run", y = paste(metrics[i], units[i])) +
        theme(plot.title = element_text(size = 15, face = "bold.italic")) +
        theme_minimal()
    }
    
    
    #saving
    filepath <- paste0(savepath,subject,"_",metrics[i],"_boxplot_", ".png")
    if (file.exists(filepath)) {file.remove(filepath) 
      message(black("Existing file deleted: ", filepath))}
    ggsave(filepath, plot = bp_sub, width = 9.21, height = 6.47)
   
    filepath <- paste0(savepath,subject,"_",metrics[i],"_lineplot_", ".png")
    if (file.exists(filepath)) {file.remove(filepath) 
      message(black("Existing file deleted: ", filepath))}
    ggsave(filepath, plot = lp_sub, width = 18, height = 6.47)
  }
    
    #artefacts need to be represented differently 
    
    if (metrics[i] == 'Artefact'){
      
      hist_sub <- ggplot(sub_subset, aes(x = factor(daytask), y = .data[[metrics[i]]],fill = factor(daytask))) +
                  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
                  scale_fill_manual(values = colours_run,name = NULL,
                                    guide = guide_legend(override.aes = list(alpha = 1))) +
                  geom_smooth(method = "lm", aes(x = as.numeric(factor(sub)), group = 1),  #add the lm 
                              linetype = "dashed", alpha = 0.2, se = TRUE, color = "#5D4037") + 
                  labs(title = "Average percentage of signal containing artefacts, across subjects (from a random 5-min segment per task)", x = "Subject", y = "% of Artefacts in signal" ) +
                  theme_minimal()+
                  theme(plot.title = element_text(size = 15, face = "bold.italic")) +
                  coord_cartesian(ylim = ylimits) 
      
      #saving
      filepath <- paste0(savepath,subject,"_",metrics[i],"_histogram.png")
      if (file.exists(filepath)) {file.remove(filepath) 
        message(black("Existing file deleted: ", filepath))}
      ggsave(filepath, plot = hist_sub, width = 15, height = 6)
    }
    
}}

message(green("All plots for each sub were successfully saved in:\n",savepath))



# B. visualisation of single units -------------------------------------------------------------

# 1. histogram of number of single units identified per subjects ---- 

#create a subset where select only the subject that have 4 screening or more
sc4_subset <- data.frame()
for (subject in unique(tbl$sub)) {
 sub_subset <- subset(tbl, sub == subject)
 if ("screening4" %in% sub_subset$run){
   sc4_subset <- rbind(sc4_subset, sub_subset)
 }}

#colour patte that fits sc4 (removing the brighter yellow)
full_palette4 <- viridis::plasma(length(unique(sc4_subset$sub))+2)
cropped_palette4 <- full_palette4[1:length(unique(sc4_subset$sub))]


#a. histogram #single units per sub:
for_avg <- data.frame(sub = unique(tbl$sub))

#in our data the number of single units is counted per tetrode, 
#so we have to  divide the total number by 4 to have the real number
for (subject in unique(tbl$sub)) {
  singu <- sum(tbl$Nb_SU[tbl$sub == subject], na.rm = TRUE)/4
  tbl$tot_sub_su[tbl$sub == subject] <- singu 
  for_avg$avg_SU[for_avg$sub == subject] <- singu
}

msu <- mean(for_avg$avg_SU,na.rm = TRUE)
sd_su <- sd(for_avg$avg_SU,na.rm = TRUE)
  

hist_su <- ggplot(tbl, aes(x = factor(sub), y = tot_sub_su,fill= factor(sub), color = factor(sub))) +
           geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
           scale_fill_manual(values = cropped_palette) +
           scale_color_manual(values = cropped_palette) +
           labs(title = "Number of single units recorded on tetrodes per subject", x = "Subject", y = "#single units") +
           geom_text(aes(label = tot_sub_su), vjust = -0.5, size = 3.5, color = "#4E342E") +  
           annotate("label", x = 1, y = max(tbl$tot_sub_su, na.rm = TRUE) * 1.05,  # Slightly above the tallest bar
                     label = paste0("Mean: ", round(msu, 1), "\nSD: ", round(sd_su, 1)), 
                     size = 4, hjust = 0, vjust = 1, fill = "white", color = "black") +
           theme_minimal() +
           theme(legend.position = "none")


#same barplot but with only 4 and more screening
for_avg <- data.frame(sub = unique(sc4_subset$sub))

for (subject in unique(sc4_subset$sub)) {
  singu <- sum(sc4_subset$Nb_SU[sc4_subset$sub == subject], na.rm = TRUE)/4
  sc4_subset$tot_sub_su[sc4_subset$sub == subject] <- singu 
  for_avg$avg_SU[for_avg$sub == subject] <- singu
}

msu <- mean(for_avg$avg_SU,na.rm = TRUE)
sd_su <- sd(for_avg$avg_SU,na.rm = TRUE)


hist_su4 <- ggplot(sc4_subset, aes(x = factor(sub), y = tot_sub_su,fill= factor(sub), color = factor(sub))) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
            scale_fill_manual(values = cropped_palette4) +
            scale_color_manual(values = cropped_palette4) +
            labs(title = "Number of single units recorded on tetrodes for subject that partook in 4 or more screening", x = "Subject", y = "#single units") +
            geom_text(aes(label = tot_sub_su), vjust = -0.5, size = 3.5, color = "#4E342E") +  
            annotate("label", x = 1, y = max(sc4_subset$tot_sub_su, na.rm = TRUE) * 1.05,  # Slightly above the tallest bar
                     label = paste0("Mean: ", round(msu, 1), "\nSD: ", round(sd_su, 1)), 
                     size = 4, hjust = 0, vjust = 1, fill = "white", color = "black") +
            theme_minimal() +
            theme(legend.position = "none")

#saving
filepath <- paste0(savepath,"overall_single-units_histogram.png")
if (file.exists(filepath)) {file.remove(filepath) 
  message(black("Existing file deleted: ", filepath))}
ggsave(filepath, plot = hist_su, width = 15, height = 6)

filepath <- paste0(savepath,"overall_single-units_over-4-screening_histogram.png")
if (file.exists(filepath)) {file.remove(filepath) 
  message(black("Existing file deleted: ", filepath))}
ggsave(filepath, plot = hist_su4, width = 15, height = 6)



# 2.for each subject, number of single units recorded on each electrode per run ----    

for (subject in unique(tbl$sub)) {
  
  #creats a subset with only the data from one suject
  sub_subset <- subset(tbl, sub == subject)
  
  #replace electode names with a the macro and tetrode name + factor
  sub_subset$electrodes <- interaction(sub_subset$gp_macro, sub_subset$gp_tetrode, sep = " ")
  sub_subset$electrodes <- factor(sub_subset$electrodes, levels = unique(sub_subset$electrodes[order(sub_subset$gp_macro, sub_subset$gp_tetrode)]))
  
  #add the total number of single units for a single tetrode
  sub_subset <- sub_subset %>%
                group_by(electrodes) %>%
                mutate(su_tet = sum(Nb_SU, na.rm = TRUE)/4)
  
  # for graphic representation of macros
  gp_macro <- sub_subset %>%
              group_by(gp_macro) %>%
              summarize(xmin = min(as.numeric(electrodes)) - 0.5,
                        xmax = max(as.numeric(electrodes)) + 0.5)
                        
  
  #a. histogram of single units per electrodes:
  
  hist_su_sub <- ggplot(sub_subset, aes(x = factor(electrodes), y = Nb_SU/4,fill= factor(run))) +
                 geom_bar(stat = "identity", position = "stack", alpha = 1) +
                 scale_fill_manual(values = colours_run)  +
                 labs(title = paste("Number of single units recorded on tetrodes in each run for",subject,"(total number of SU =", sub_subset$tot_sub_su[1],")"), 
                     x = "tetrodes", 
                     y = "#single units") +
                 geom_text(aes(label = su_tet, x=factor (electrodes)),y = sub_subset$su_tet, vjust = -0.5, size = 3.5, color = "#4E342E") +
                 theme_minimal()
  
  
  #saving
  filepath <- paste0(savepath, subject,"_single-units_histogram.png")
  if (file.exists(filepath)) {file.remove(filepath) 
    message(black("Existing file deleted: ", filepath))}
  ggsave(filepath, plot = hist_su_sub, width = 18, height = 6.47)
  
  
}
message(green("All single units plots for each sub were successfully saved in:\n",savepath))


# 3. 3D  scatter plot----
#⚠️  only saved in HTML format to keep the plot interactive , you can extract it manually too 

savepath2 <- paste(savepath,"interactive-plots/",sep = "")

for (i in 1:length(metrics)) {
  sca_3D <- plot_ly(x=factor(sc4_subset$sub), y= sc4_subset[[metrics[i]]], z=sc4_subset$tot_sub_su, type="scatter3d", mode="markers", color=factor(sc4_subset$sub))
  sca_3D <- sca_3D %>% add_markers()
  sca_3D <- sca_3D %>% layout(scene = list(xaxis = list(title = 'sub'),
                                           yaxis = list(title = units[i]),
                                           zaxis = list(title = '# signle units')))

  # Save as html
  htmlwidgets::saveWidget(sca_3D,
    file = paste(savepath2,"3D_plot_",sub,".html",sep = ""),
    libdir = paste(savepath2,"supporting_files/","3D_plot_",sub,sep = ""),
    selfcontained = FALSE)
  
}

#------------ S T A T I S T I C S --------------------------------------------

# A. correlations ------------------------------------------------------------
#select all numerical value that I want to test in the correlation matrix  
cor.data <- tbl %>%
  select(RMS, RMS_filter, variance, variance_norm, tetrode_cor, deviation, kurtosis, SNR,Artefact,Hurst,Nb_SU)

cor.mat <- cor_mat(cor.data)
print(cor.mat)
cor.pval <- cor_get_pval(cor.mat)

#plot correlation
cor.mat %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)




