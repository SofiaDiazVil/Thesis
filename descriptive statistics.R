
# Libraries 
library(R.matlab)
library(lsr)

## ANALYSIS FOR PATIENT GROUP ## 

# Import data files and make data sets 
#Group
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/Group_Pt.mat") 
group <- unlist( list.import.mat[1] )
group <- array( group , dim = c( 54 , 1 ))

#HR 
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/HR_Pt.mat") 
HR_Pt <- unlist( list.import.mat[1] )
HR_Pt <- array( HR_Pt , dim = c( 54 , 1 ))
HR_Erotic_Pt <- HR_Pt[group == 2]
HR_Erotic_1_Pt <- HR_Erotic_Pt[c(1,3,5,7,9,11,13,15,17)]
HR_Erotic_2_Pt <- HR_Erotic_Pt[c(2,4,6,8,10,12,14,16,18)]
HR_Neutral_Pt <- HR_Pt[group == 3]
HR_Neutral_1_Pt <-HR_Neutral_Pt[c(1,3,5,7,9,11,13,15,17)]
HR_Neutral_2_Pt <-HR_Neutral_Pt[c(2,4,6,8,10,12,14,16,18)]

#HRV 
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/HRV_Pt.mat") 
HRV_Pt <- unlist( list.import.mat[1] )
HRV_Pt <- array( HRV_Pt , dim = c( 54 , 899))
HRV_Erotic <- HRV_Pt[group == 2]
HRV_Neutral <- HRV_Pt[group == 3]

list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/SD_HRV_Pt.mat") 
SD_HRV_Pt <- unlist( list.import.mat[1] )
SD_HRV_Pt <- array( SD_HRV_Pt , dim = c( 54 , 1))
SD_HRV_Erotic_Pt <- SD_HRV_Pt[group == 2]
SD_HRV_Neutral_Pt <- SD_HRV_Pt[group == 3]

#ISC
list.import.mat <- R.matlab::readMat("D:/Users/PERSONAL/Documents/psychologie Master thesis/Data Sets Results/ISC_Pt.mat") 
ISC_Pt <- unlist( list.import.mat[1] )
ISC_Pt <- array( ISC_Pt , dim = c( 61 , 1 ))
ISC_Pt <- na.omit(ISC_Pt) # first item is NAN, exclude it
ISC_Erotic <- ISC_ALL[group == 2]
ISC_Neutral <- ISC_ALL[group == 3]


# Decriptive statistics 

#HR
mean(HR_Erotic_1_Pt) # 76.712
sd(HR_Erotic_1_Pt) # 15.388
mean(HR_Erotic_2_Pt) # 84.069
sd(HR_Erotic_2_Pt) # 11.541
mean(HR_Neutral_1_Pt) # 79.075
sd(HR_Neutral_1_Pt) # 12.021
mean(HR_Neutral_2_Pt) # 83.777
sd(HR_Neutral_2_Pt) # 9.898



## ANALYSIS FOR HEALTHY GROUP ##

# Import data files 
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
HRV_Ht <- array( HRV_Ht , dim = c( 60 , 899))
HRV_Erotic_Ht <- HRV_Ht[group == 2]
HRV_Neutral_Ht <- HRV_Ht[group == 3]

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

# Decriptive statistics 

#HR
mean(HR_Erotic_Ht) # 75.661
sd(HR_Erotic_Ht) # 8.324
mean(HR_Neutral_Ht) # 77.543
sd(HR_Neutral_Ht) # 8.282


#HRV
mean(SD_HRV_Neutral_Ht) 
sd(SD_HRV_Neutral_Ht) 
mean(SD_HRV_Erotic_Ht) 
sd(SD_HRV_Erotic_Ht) 
# cohens d
cohensD(SD_HRV_Erotic_Ht, SD_HRV_Neutral_Ht)

#ISC 
mean(ISC_Erotic_Ht) # 0.029
sd(ISC_Erotic_Ht) # 0.036
mean(ISC_Neutral_Ht) # 0.0147
sd(ISC_Neutral_Ht) #0.047



