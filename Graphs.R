# Load needed libraries
library(tidyverse)
library(purrr)
library(R.matlab)
library(plotly)
library(ggplot2)
library(GGally)

# Import data file patient sample 
#Group
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/Group_Pt.mat") 
group <- unlist( list.import.mat[1] )
group <- array( group , dim = c( 54 , 1 ))

#HR 
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/HR_Pt.mat") 
HR_Pt <- unlist( list.import.mat[1] )
HR_Pt <- array( HR_Pt , dim = c( 54 , 1 ))
HR_Erotic_Pt <- HR_Pt[group == 2]
HR_Erotic_1_Pt <- HR_Erotic_Pt[c(1,3,5,7,9,11,13,15)]
HR_Erotic_2_Pt <- HR_Erotic_Pt[c(2,4,6,8,10,12,14,16)]
HR_Neutral_Pt <- HR_Pt[group == 3]
HR_Neutral_1_Pt <-HR_Neutral_Pt[c(1,3,5,7,9,11,13,15)]
HR_Neutral_2_Pt <-HR_Neutral_Pt[c(2,4,6,8,10,12,14,16)]

#HRV 
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/HRV_Pt.mat") 
HRV_Pt <- unlist( list.import.mat[1] )
HRV_Pt <- array( HRV_Pt , dim = c( 54 , 899))
HRV_Erotic_Pt <- HRV_Pt[group == 2,]
HRV_Erotic_1_Pt <- HRV_Erotic_Pt[c(1,3,5,7,9,11,13,15),]
HRV_Erotic_2_Pt <- HRV_Erotic_Pt[c(2,4,6,8,10,12,14,16),]
HRV_Neutral_Pt <- HRV_Pt[group == 3,]
HRV_Neutral_1_Pt <- HRV_Neutral_Pt[c(1,3,5,7,9,11,13,15),]
HRV_Neutral_2_Pt <- HRV_Neutral_Pt[c(2,4,6,8,10,12,14,16),]

list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/SD_HRV_Pt.mat") 
SD_HRV_Pt <- unlist( list.import.mat[1] )
SD_HRV_Pt <- array( SD_HRV_Pt , dim = c( 54 , 1))
SD_HRV_Erotic_Pt <- SD_HRV_Pt[group == 2]
SD_HRV_Erotic_1_Pt <- SD_HRV_Erotic_Pt[c(1,3,5,7,9,11,13,15)]
SD_HRV_Erotic_2_Pt <- SD_HRV_Erotic_Pt[c(2,4,6,8,10,12,14,16)]
SD_HRV_Neutral_Pt <- SD_HRV_Pt[group == 3]
SD_HRV_Neutral_1_Pt <- SD_HRV_Neutral_Pt[c(1,3,5,7,9,11,13,15)]
SD_HRV_Neutral_2_Pt <- SD_HRV_Neutral_Pt[c(2,4,6,8,10,12,14,16)]


#ISC
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/ISC_Pt.mat") 
ISC_Pt <- unlist( list.import.mat[1] )
ISC_Pt <- array( ISC_Pt , dim = c( 61 , 1 ))
ISC_Pt <- na.omit(ISC_Pt) # first item is NAN, exclude it
ISC_Erotic_Pt <- ISC_Pt[group == 2]
ISC_Erotic_1_Pt <-ISC_Erotic_Pt[c(1,3,5,7,9,11,13,15)]
ISC_Erotic_2_Pt <- ISC_Erotic_Pt[c(2,4,6,8,10,12,14,16)]
ISC_Neutral_Pt <- ISC_Pt[group == 3]
ISC_Neutral_1_Pt <- ISC_Neutral_Pt[c(1,3,5,7,9,11,13,15)]
ISC_Neutral_2_Pt <- ISC_Neutral_Pt[c(2,4,6,8,10,12,14,16)]

list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/Pval_Erotic1_Pt.mat") 
Pval_Erotic1_Pt <- unlist( list.import.mat[1] )
Pval_Erotic1_Pt <- array( Pval_Erotic1_Pt , dim = c( 8 , 1 ))

list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/Pval_Erotic2_Pt.mat") 
Pval_Erotic2_Pt <- unlist( list.import.mat[1] )
Pval_Erotic2_Pt <- array( Pval_Erotic2_Pt , dim = c( 8 , 1 ))

list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/Pval_Neutral1_Pt.mat") 
Pval_Neutral1_Pt <- unlist( list.import.mat[1] )
Pval_Neutral1_Pt <- array( Pval_Neutral1_Pt , dim = c( 8 , 1 ))

list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/Pval_Neutral2_Pt.mat") 
Pval_Neutral2_Pt <- unlist( list.import.mat[1] )
Pval_Neutral2_Pt <- array( Pval_Neutral2_Pt , dim = c( 8 , 1 ))


# Graph HR
  # Parallel plot --> show how each individual changes
    #make data set for the plot
HR <- data.frame(HR_Erotic_1_Pt,HR_Neutral_1_Pt,HR_Erotic_2_Pt, HR_Neutral_2_Pt)
colnames(HR) <- c("Erotic 1", "Neutral 1", "Erotic 2", "Neutral 2")
MeanHR <- c(mean(HR_Erotic_1_Pt),mean(HR_Neutral_1_Pt), mean(HR_Erotic_2_Pt), mean(HR_Neutral_2_Pt))
HR <- rbind(HR,MeanHR)
HR$Value <- "subject"
HR$Value[9] <- "mean"

  # create the plot 
HR_plot_Pt <- ggparcoord(HR,
                columns = 1:4, scale = "globalminmax", showPoints = TRUE, groupColumn = 5 
                ) +
                labs(x = "Condition and Session", y = "HR (bpm)")
                
ggplotly(HR_plot_Pt)

  # Distribution plot --> show the differences in distribution across conditions

HR <- append(HR_Erotic_1_Pt,HR_Neutral_1_Pt)
HR <- append(HR, HR_Erotic_2_Pt)
HR <- append(HR, HR_Neutral_2_Pt)
HR <- as.data.frame(HR)
HR$condition[1:8] <- "Erotic 1"
HR$condition[9:16] <- "Neutral 1"
HR$condition[17:24] <- "Erotic 2"
HR$condition[25:32] <- "Neutral 2"

HR %>% 
  ggplot(aes( x= condition, y = HR, fill = condition)) +
  
  ggdist::stat_halfeye(
    adjust= 0.6,
    justification = -0.1,
    .width = 0, 
    point_colour = NA
  ) + 
  
  geom_boxplot(
    width = 0.12, 
    outlier.color = NA, 
    alpha = 0.5
  ) + 
  
  stat_summary(
    fun = mean, geom = "point", col = "red"
  ) +
  
  labs(x = "Condition and Session", y = "HR (bpm)")


# Graph HRV 
  # Parallel plot --> show how each individual changes 
    #make data set for the plot
HRV <- data.frame(SD_HRV_Erotic_1_Pt,SD_HRV_Neutral_1_Pt,SD_HRV_Erotic_2_Pt, SD_HRV_Neutral_2_Pt)
colnames(HRV) <- c("Erotic 1", "Neutral 1", "Erotic 2", "Neutral 2")
MeanHRV <- c(mean(SD_HRV_Erotic_1_Pt),mean(SD_HRV_Neutral_1_Pt), mean(SD_HRV_Erotic_2_Pt), mean(SD_HRV_Neutral_2_Pt))
HRV <- rbind(HRV,MeanHRV)
HRV$Value <- "subject"
HRV$Value[9] <- "mean"

   # create the plot 
HRV_plot_Pt <- ggparcoord(HRV,
                         columns = 1:4, scale = "globalminmax", showPoints = TRUE, groupColumn = 5 
) +
  labs(x = "Condition and Session", y = "HRV")

ggplotly(HRV_plot_Pt)

# Distribution plot --> show the differences in distribution across conditions
HRV <- append(SD_HRV_Erotic_1_Pt,SD_HRV_Neutral_1_Pt)
HRV <- append(HRV, SD_HRV_Erotic_2_Pt)
HRV <- append(HRV, SD_HRV_Neutral_2_Pt)
HRV <- as.data.frame(HRV)
HRV$condition[1:8] <- "Erotic 1"
HRV$condition[9:16] <- "Neutral 1"
HRV$condition[17:24] <- "Erotic 2"
HRV$condition[25:32] <- "Neutral 2"

HRV %>% 
  ggplot(aes( x= condition, y = HRV, fill = condition)) +
  
  ggdist::stat_halfeye(
    adjust= 0.6,
    justification = -0.1,
    .width = 0, 
    point_colour = NA
  ) + 
  
  geom_boxplot(
    width = 0.12, 
    outlier.color = NA, 
    alpha = 0.5
  ) + 
  
  stat_summary(
    fun = mean, geom = "point", col = "red"
  ) +
  
  labs(x = "Condition and Session", y = "HRV")

#Graph ISC 
  # Parallel plot --> show how each individual changes 
    #make data set for the plot
ISC <- data.frame(ISC_Erotic_1_Pt,ISC_Neutral_1_Pt,ISC_Erotic_2_Pt, ISC_Neutral_2_Pt)
colnames(ISC) <- c("Erotic 1", "Neutral 1", "Erotic 2", "Neutral 2")
MeanISC <- c(mean(ISC_Erotic_1_Pt),mean(ISC_Neutral_1_Pt), mean(ISC_Erotic_2_Pt), mean(ISC_Neutral_2_Pt))
ISC <- rbind(ISC,MeanISC)
ISC$Value <- "subject"
ISC$Value[9] <- "mean"

    # create the plot 
ISC_plot_Pt <- ggparcoord(ISC,
                          columns = 1:4, scale = "globalminmax", showPoints = TRUE, groupColumn = 5 
) +
  labs(x = "Condition and Session", y = "ISC-HR")

ggplotly(ISC_plot_Pt)

    # make data set for p values 
PVal <- data.frame(Pval_Erotic1_Pt,Pval_Neutral1_Pt,Pval_Erotic2_Pt,Pval_Neutral2_Pt)
PVal[PVal < 0.05] <- "Significant"
PVal[PVal != "Significant"] <- "Non-Significant"

# I STILL NEED TO MAKE A GRAPH WHOICH SHOWS THE ISC FOR EACH PARTICIPANT IN EACH CONDITION BUT ALSO SHOW
# WHICH WERE SIGNIFICANT AND WHICH WEREN'T 

########################################################################################################
########################################################################################################



# Import data files for healthy sample 
  #Group
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/Group_Ht.mat") 
group <- unlist( list.import.mat[1] )
group <- array( group , dim = c( 60 , 1 ))

  #HR 
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/HR_Ht.mat") 
HR_Ht <- unlist( list.import.mat[1] )
HR_Ht <- array( HR_Ht , dim = c( 60 , 1 ))
HR_Erotic_Ht <- HR_Ht[group == 2]
HR_Neutral_Ht <- HR_Ht[group == 3]


  #HRV 
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/HRV_Ht.mat") 
HRV_Ht <- unlist( list.import.mat[1] )
HRV_Ht <- array( HRV_Ht , dim = c( 60 , 899 ))
HRV_Erotic_Ht <- HRV_Ht[group == 2,]
HRV_Neutral_Ht <- HRV_Ht[group == 3,]

list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/SD_HRV_Ht.mat") 
SD_HRV_Ht <- unlist( list.import.mat[1] )
SD_HRV_Ht <- array( SD_HRV_Ht , dim = c( 60 , 1))
SD_HRV_Erotic_Ht <- SD_HRV_Ht[group == 2]
SD_HRV_Neutral_Ht <- SD_HRV_Ht[group == 3]

  #ISC
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/ISC_Ht.mat") 
ISC_Ht <- unlist( list.import.mat[1] )
ISC_Ht <- array( ISC_Ht , dim = c( 60 , 1 ))
ISC_Erotic_Ht <- ISC_Ht[group == 2]
ISC_Neutral_Ht <- ISC_Ht[group == 3]

list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/Pval_Erotic_Ht.mat") 
Pval_Erotic_Ht <- unlist( list.import.mat[1] )
Pval_Erotic_Ht <- array( Pval_Erotic_Ht , dim = c( 30 , 1 ))

list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/Pval_Neutral_Ht.mat") 
Pval_Neutral_Ht <- unlist( list.import.mat[1] )
Pval_Neutral_Ht <- array( Pval_Neutral_Ht , dim = c( 30 , 1 ))



# Graph HR
  # Parallel plot --> show how each individual changes 
    #make data set for the plot
HR <- data.frame(HR_Erotic_Ht,HR_Neutral_Ht)
colnames(HR) <- c("Erotic", "Neutral")
MeanHR <- c(mean(HR_Erotic_Ht),mean(HR_Neutral_Ht))
HR <- rbind(HR, MeanHR)
HR$Value <- "subject"
HR$Value[31] <- "mean"

    # create the plot 
HR_plot_Pt <- ggparcoord(HR,columns = 1:2, scale = "globalminmax", 
                         showPoints = TRUE, groupColumn = 3, alphaLines = 1
                         ) +
  
                        labs(x = "Condition", y = "HR (bpm)")

ggplotly(HR_plot_Pt)

 # Distribution plot --> show the differences in distribution across conditions
HR <- as.data.frame(HR_Ht[group == 2 | group == 3])
colnames(HR)[1] <- "HR"
HR$condition[group == 2] <- "Erotic"
HR$condition[group == 3] <- "Neutral"

HR %>% 
  ggplot(aes( x= condition, y = HR, fill = condition)) +
  
  ggdist::stat_halfeye(
    adjust= 0.5,
    justification = -0.1,
    .width = 0, 
    point_colour = NA
  ) + 
  
  geom_boxplot(
    width = 0.12, 
    outlier.color = NA, 
    alpha = 0.5
  ) + 
  
  stat_summary(
    fun = mean, geom = "point", col = "red"
  ) + 
  
  labs(x = "Condition", y = "HR (bpm)")

# Graph HRV 
  # Parallel plot --> show how each individual changes 
#make data set for the plot
HRV <- data.frame(SD_HRV_Erotic_Ht,SD_HRV_Neutral_Ht)
colnames(HRV) <- c("Erotic", "Neutral")
MeanHRV <- c(mean(SD_HRV_Erotic_Ht),mean(SD_HRV_Neutral_Ht))
HRV <- rbind(HRV,MeanHRV)
HRV$Value <- "subject"
HRV$Value[31] <- "mean"

# create the plot 
HRV_plot_Pt <- ggparcoord(HRV,
                          columns = 1:2, scale = "globalminmax", showPoints = TRUE, groupColumn = 3 
) +
  labs(x = "Condition and Session", y = "HRV")

ggplotly(HRV_plot_Pt)

  # Distribution plot --> show the differences in distribution across conditions
HRV <- as.data.frame(SD_HRV_Ht[group == 2 | group == 3])
colnames(HRV)[1] <- "HRV"
HRV$condition[group == 2] <- "Erotic"
HRV$condition[group == 3] <- "Neutral"

HRV %>% 
  ggplot(aes( x= condition, y = HRV, fill = condition)) +
  
  ggdist::stat_halfeye(
    adjust= 0.5,
    justification = -0.1,
    .width = 0, 
    point_colour = NA
  ) + 
  
  geom_boxplot(
    width = 0.12, 
    outlier.color = NA, 
    alpha = 0.5
  ) + 
  
  stat_summary(
    fun = mean, geom = "point", col = "red"
  ) + 
  
  labs(x = "Condition", y = "HRV")

#Graph ISC 
  #make data set for the plot
ISC <- data.frame(ISC_Erotic_Ht,ISC_Neutral_Ht)
colnames(ISC) <- c("Erotic", "Neutral")
MeanISC <- c(mean(ISC_Erotic_Ht),mean(ISC_Neutral_Ht))
ISC <- rbind(ISC,MeanISC)
ISC$Value <- "subject"
ISC$Value[31] <- "mean"

# create the plot 
ISC_plot_Pt <- ggparcoord(ISC,
                          columns = 1:2, scale = "globalminmax", showPoints = TRUE, groupColumn = 3 
) +
  labs(x = "Condition and Session", y = "HRV")

ggplotly(ISC_plot_Pt)

  # make data set for p values 
PVal <- data.frame(Pval_Erotic_Ht,Pval_Neutral_Ht)
PVal[PVal < 0.05] <- "Significant"
PVal[PVal != "Significant"] <- "Non-Significant"




