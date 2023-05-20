# set the correct working directory 
setwd('D:/Users/PERSONAL/Documents/psychologie Master thesis/experiment_gent/TimeStamps_exp/files')

# open the needed libraries 
library("dplyr")   # Load dplyr package
library("plyr")    # Load plyr package
library("readr")   # Load readr package
library("ggplot2") # Load ggplot2 package
library("reshape2")# Load reshape2 package
library("writexl") # Load writexl package
library ("tidyverse") # Load tidyverse package

#download all time stamp files and merge them 
df <- list.files(path= 'D:/Users/PERSONAL/Documents/psychologie Master thesis/experiment_gent/TimeStamps_exp/files',    
                 pattern = NULL , full.names = FALSE) 
df <- df[1:30] %>% #select only the files for the time stamps 
  lapply(read_csv) %>%
  bind_rows  %>% #bind them together to make one df
  as.data.frame(df)

class(df) #check it is a data frame

#fix the participant number, age, gender
i <- 1
while (i < dim(df)[1]) {
  df$participant_number[i] <- df$participant_number[i+1]
  df$gender[i]             <- df$gender[i+1]
  df$age[i]                <- df$age[i+1]
  
  i <- i+3
}

# Calculate the average valence and arousal of each video 
  #videos
df_videos <- df[!df$condition == -99,]
  #averages (0 = erotic, 1 = neutral)
aggregate(df_videos$arousal, list(df_videos$condition), FUN=mean) # 0 = 5.766/ 1 = 3.366
aggregate(df_videos$arousal, list(df_videos$condition), FUN=sd) # 0 = 1.1633/ 1 = 1.564
aggregate(df_videos$valence, list(df_videos$condition), FUN=mean) # 0 = 6.233/ 1 = 5.833
aggregate(df_videos$valence, list(df_videos$condition), FUN=sd) # 0 = 1.04/ 1 = 1.126

# t test appraisal video 
Erotic_vid <- df_videos[df_videos$condition == 0,5:6]
Neutral_vid <- df_videos[df_videos$condition == 1,5:6]
t.test(Erotic_vid$valence, Neutral_vid$valence, alternative = "greater", paired = TRUE)
  #t(29) = 1.934, p < 0.05
t.test(Erotic_vid$arousal, Neutral_vid$arousal, "greater", paired = TRUE)
  #t(29) = 9.2, p < 0.001

df_videos[df_videos$condition == "0",2] <- "Erotic"
df_videos[df_videos$condition == "1",2] <- "Neutral"

df_videos %>% 
  ggplot(aes(x = factor(condition), y = valence, fill = factor(condition))) +
  geom_dotplot(binaxis = "y", stackdir = "center") + 
  stat_summary(fun = mean, geom = "point", col = "red") +
  labs(x = "Condition", y = "Valence Rating")


df_videos %>% 
  ggplot(aes(x = factor(condition), y = arousal, fill = factor(condition))) +
  geom_dotplot(binaxis = "y", stackdir = "center") + 
  stat_summary(fun = mean, geom = "point", col = "red") +
  labs(x = "Condition", y = "Arousal Rating")

# Age 
index <- 1:dim(df)[1]
index <- index[index %% 3 == 1]
Age <- df$age[index]

mean(Age) #25.4
sd(Age) #6.262

# Gender
index <- 1:dim(df)[1]
index <- index[index %% 3 == 1]
Gender <- df$gender[index]
table(Gender)







