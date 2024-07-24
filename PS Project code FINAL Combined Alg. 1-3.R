## Reading and previewing code
library(R.matlab)
library(tidyverse)
library(dplyr)
data1 <- readMat("C:/Users/Vartotojas/Desktop/PS Project data/facebook100/Yale4.mat")

## Using detailed data instead of the sparse matrix that is also in the dataset.
data1$local.info
View(data1$local.info)

# Creating a separate dataframe of the dataset that will be used.
datafr <- as.data.frame(data1$local.info)
datafr

# Adding a new column/variable to the dataset using bernoulli randomisation.
# This variable represents the binary treatment value. 
# pi = 0.5 or 50%
datafr$W <- c(rbinom(8578,1,0.5))

# Adding/Creating candidate exposure function variable in the dataset
datafr$H<-NA
datafr[datafr$W==1,"H"] <- sum(datafr$W ==1)-1 
datafr[datafr$W==0,"H"] <- sum(datafr$W ==1) 

# Calculating rows to specify sampling: rows - 8578, rows/2=4289
nrow(datafr)

# Alg. 1.1.
# Splitting data to focal and auxiliary units.
Sample_rows <- sample(1:nrow(datafr), 4289, replace = F)

Data_foc <- datafr[Sample_rows, ]
Data_aux <- datafr[-c(Sample_rows), ]

Data1_foc<-Data_foc # not necessary
Data2_aux<-Data_aux # not necessary

# Alg. 1.2.
# Regression for T - W is excluded as it has a linear relationship with H
Obtain_T = lm(V3 ~V1 + V2 + V4 + V5 + V6 + V7 + H, data = Data_foc)
summary(Obtain_T)

# Naming and saving the coefficient for further calculations
H_0 <- Obtain_T$coefficients[8]

B=200

for (i in 1:B)
{
  Data2_aux$W <- c(rbinom(4289,1,0.5)) # 1.3.1 Regenerating treatments...
  Data1_foc$H<-NA # 1.3.2 recomputing the candidate exposure...
  Data1_foc[Data1_foc$W==1,"H"] <- sum(Data1_foc$W ==1)+sum(Data2_aux$W ==1)-1 
  Data1_foc[Data1_foc$W==0,"H"] <- sum(Data1_foc$W ==1)+sum(Data2_aux$W ==1)
  Tb = lm(V3 ~V1 + V2 + V4 + V5 + V6 + V7 + H, data = Data1_foc) # 1.3.3
  summary(Tb) # 1.3.3
  H_b <- Tb$coefficients[8]
  if (H_0<=H_b) # for sum of indicator function notating T0<=Tb
  {
    i=i+1
  }
}
# Output
P=(1/(B+1))*(1+i)
P

# Algorithm 2.
# Data for 2 experiments
dset1 <- as.data.frame(data1$local.info) 
dset2 <- as.data.frame(data1$local.info)

# Generating treatments for 2 experiments
# and calculating exposure function for each candidate
dset1$W <- c(rbinom(8578,1,0.5))
View(dset1)

dset1$H<-NA
dset1[dset1$W==1,"H"] <- sum(dset1$W ==1)-1 
dset1[dset1$W==0,"H"] <- sum(dset1$W ==1) 

dset2$W <- c(rbinom(8578,1,0.5))
View(dset2)

dset2$H<-NA
dset2[dset2$W==1,"H"] <- sum(dset2$W ==1)-1 
dset2[dset2$W==0,"H"] <- sum(dset2$W ==1)

# 2.1. Creating the set of units whose treatment didn't change over the experiments
# Creating additional variables to include k2 variables in k1 dataset for simplicity
dset1$W2 <-NA
dset1$W2 <- dset2$W
dset1$H2 <-NA
dset1$H2 <- dset2$H

# Calculating Ydiff for the dataset. It will be used later on in 2.2.
Ydiff <- dset2$V3-dset1$V3
dset1$Ydiff<-NA
dset1$Ydiff<- Ydiff_foc

# Ind_nc - Set of units whose treatment didn't change over the experiments
Ind_nc <- subset(dset1, dset1$W==dset1$W2)

Sr <- sample(1:nrow(Ind_nc), nrow(Ind_nc)/2, replace = F)

Ind_foc <- Ind_nc[Sr, ]
Ind_aux <- Ind_nc[-c(Sr), ]

# Additional variables to use in for loop (not necessary, for clarity)
Ind_foc1 <- Ind_foc
Ind_aux1 <- Ind_aux

# Creating additional Hdff=H2-H1 variable,
# to express differences of candidate exposure function between experiments.
# This variable is to be used in regression model 
# to obtain the coefficient for further calculations.
Ind_foc$Hdiff<-NA
Ind_foc$Hdiff<-Ind_foc$H2-Ind_foc$H

# 2.2. 
Alg2Obtain_T = lm(Ydiff ~V1 + V2 + V4 + V5 + V6 + V7 + H + Hdiff , data = Ind_foc)
summary(Alg2Obtain_T)
Alg2_T0 <- Alg2Obtain_T$coefficients[8]

# 2.3.
for (j in 1:B)
{
  # 2.3.1 Permutate treatments - permutate rows
  # For W
  perm1<-as.data.frame(Ind_aux1$W) 
  perm1<- perm1 %>% sample_n(nrow(.))
  Ind_aux1$W<-perm1
  # For W2
  perm2<-as.data.frame(Ind_aux1$W2)
  perm2<- perm2 %>% sample_n(nrow(.))
  Ind_aux1$W2<-perm2
  # 2.3.2
  # Recompute candidate exposure Hfoc for k1 and k2, or H and H2 in the combined dataset 
  Ind_foc1$H<-NA 
  Ind_foc1[Ind_foc1$W==1,"H"] <- sum(Ind_foc1$W ==1)+sum(Ind_aux1$W ==1)-1 
  Ind_foc1[Ind_foc1$W==0,"H"] <- sum(Ind_foc1$W ==1)+sum(Ind_aux1$W ==1)
  
  Ind_foc1$H2<-NA 
  Ind_foc1[Ind_foc1$W2==1,"H2"] <- sum(Ind_foc1$W2 ==1)+sum(Ind_aux1$W2 ==1)-1 
  Ind_foc1[Ind_foc1$W2==0,"H2"] <- sum(Ind_foc1$W2 ==1)+sum(Ind_aux1$W2 ==1)
  # 2.3.3
  Tb_Alg2 = lm(Ydiff ~V1 + V2 + V4 + V5 + V6 + V7 + H + Hdiff, data = Ind_foc1) # 1.3.3
  summary(Tb_Alg2) 
  H_b_Alg2 <- Tb_Alg2$coefficients[8]
  if (Alg2_T0<=H_b_Alg2) # for sum of indicator function notating T0<=Tb
  {
    j=j+1
  }
}
Alg2_P=(1/(B+1))*(1+j)
Alg2_P

# Algorithm 3.
# Using the same 2 experiment data used in algorithm 2
#dset1, dset2
# 3.1.
# Subsetting data from both experiments in to two groups - W=0 and W=1
Ind_0 <- subset(dset1, dset1$W==dset1$W2 & dset1$W==0)
Ind_1 <- subset(dset1, dset1$W==dset1$W2 & dset1$W==1)
# 3.2.
# Arranging data by year as a form of semi-random matching.
Ind_0 <- arrange(Ind_0,Ind_0$V6)
Ind_1 <- arrange(Ind_1,Ind_1$V6)
Ind_0 <- sample_n(Ind_0,nrow(Ind_1)) # making the data groups same row length

Ind_0l <- Ind_0 
Ind_1l <- Ind_1

# 3.3.
# Computing Ydiff_Ind1 and Test statistic
Ind_1$Ydiff_Ind1<-NA
Ind_1$Ydiff_Ind1<- Ind_1$V3-Ind_0$V3

# Computing an alternative test statistic
# Subtracting the same Ydiff as covariates between experiments were unchanged
Alg3_T0 <- mean(Ind_1$Ydiff_Ind1)-mean(Ind_1$Ydiff_Ind1)

# 3.4.
rowsl<-nrow(Ind_1l)
# For loop 3.4.
for (k in 1:B)
{
  # Nested loop for random permutation - binomial randomization with probability 0.5
  for (m in 1:rowsl)
  {
    binom_perm <- rbinom(1, size=1, prob=0.5)
    if (binom_perm==1)
    {
      # Outcomes across experiments are identical in this case,
      # thus permutating only for outcomes of matched groups/indices.
      # In real case scenario when outcomes differ,
      # an identical algorithm can be used for difference between experiments.
      # a1 & b1 as mediator variables
      # to avoid duplication of variables between experiments
      a1<-Ind_1l$V3[m] 
      b1<-Ind_0l$V3[m]
      Ind_1l$V3[m]<-b1
      Ind_0l$V3[m]<-a1
    }}
  # Recomputing Ydiff
  Ind_1l$Ydiff_Ind1<- NA
  Ind_1l$Ydiff_Ind1<- Ind_1l$V3-Ind_0l$V3
  
  # Recomputing test statistic - difference in means
  Alg3_Tb <- mean(Ind_1l$Ydiff_Ind1)-mean(Ind_1l$Ydiff_Ind1)
  if (Alg3_T0<=Alg3_Tb) # for sum of indicator function notating T0<=Tb
  {
    k=k+1
  }
}
Alg3_P=(1/(B+1))*(1+k)
Alg3_P