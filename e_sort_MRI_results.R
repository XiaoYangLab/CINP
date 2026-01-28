# sort results of multimodel MRI analysis derived from UKB pipeline
# Author: Yang Xiao, PKU, 2024
# @: xiaoyang9604@gmail.com

rm(list=ls())
library(readxl)
library(glmnet)
library(stringr)
library(data.table)

# sort data 
dat_idx = " " 
mypath = ' '

################# 1. structures ######################
struc_volume <- list()
myfile = list.files(path = str_c(mypath, dat_idx, '/structures/'),
                    pattern = "_cortex\\.csv$", full.names = TRUE)
for(i in 1:length(myfile)){
tmp <- read.csv(myfile[i], sep = "", header = TRUE)
struc_volume[[i]] <- t(tmp$volume) 

}
struc_volume_df <- do.call(rbind, struc_volume)
struc_volume_df <- data.frame(struc_volume_df)
names <- paste('volume', tmp$region, sep = "_")
colnames(struc_volume_df) <- names
rownames(struc_volume_df) <- gsub("_cortex\\.csv$", "", basename(myfile))

struc_areas <- list()
for(i in 1:length(myfile)){
  tmp <- read.csv(myfile[i], sep = "", header = TRUE)
  struc_areas[[i]] <- t(tmp$areas) 
  
}
struc_areas_df <- do.call(rbind, struc_areas)
struc_areas_df <- data.frame(struc_areas_df)
names <- paste('areas', tmp$region, sep = "_")
colnames(struc_areas_df) <- names
rownames(struc_areas_df) <- gsub("_cortex\\.csv$", "", basename(myfile))

struc_thick <- list()
for(i in 1:length(myfile)){
  tmp <- read.csv(myfile[i], sep = "", header = TRUE)
  struc_thick[[i]] <- t(tmp$thick) 
  
}
struc_thick_df <- do.call(rbind, struc_thick)
struc_thick_df <- data.frame(struc_thick_df)
names <- paste('thick', tmp$region, sep = "_")
colnames(struc_thick_df) <- names
rownames(struc_thick_df) <- gsub("_cortex\\.csv$", "", basename(myfile))

struc_SV <- list()
myfile = list.files(path = str_c(mypath, dat_idx, '/structures/'),
                    pattern = "_subcortex\\.csv$", full.names = TRUE)
for(i in 1:length(myfile)){
  tmp <- read.csv(myfile[i], sep = "", header = TRUE)
  struc_SV[[i]] <- t(tmp$volume) 
  
}
struc_SV_df <- do.call(rbind, struc_SV)
struc_SV_df <- data.frame(struc_SV_df)
names <- paste('volume', tmp$region, sep = "_")
colnames(struc_SV_df) <- names
rownames(struc_SV_df) <- gsub("_subcortex\\.csv$", "", basename(myfile))
struc_SV_df <- struc_SV_df[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,
                                        15,19,21,22,23,24,25,26,27,28,29,34)]

struc_df <- cbind(struc_volume_df, 
                       struc_areas_df, 
                       struc_thick_df,
                       struc_SV_df)
write.csv(struc_df, str_c(mypath, "struc_", dat_idx, "_df.csv"))

########################### 2. fiber tracts ##################
struc_FA <- list()
myfile = list.files(path = str_c(mypath, dat_idx, '/fibers/'),
                    full.names = TRUE)
for(i in 1:length(myfile)){
  tmp <- read.csv(myfile[i], sep = "", header = TRUE)
  struc_FA_[[i]] <- t(tmp$FA) 
  
}
struc_FA_df <- do.call(rbind, struc_FA)
struc_FA_df <- data.frame(struc_FA_df)
names <- paste('FA', tmp$Fiber, sep = "_")
colnames(struc_FA_df) <- names
rownames(struc_FA_df) <- gsub("_fibers\\.csv$", "", basename(myfile))

struc_MD <- list()
for(i in 1:length(myfile)){
  tmp <- read.csv(myfile[i], sep = "", header = TRUE)
  struc_MD[[i]] <- t(tmp$MD) 
  
}
struc_MD_df <- do.call(rbind, struc_MD)
struc_MD_df <- data.frame(struc_MD_df)
names <- paste('MD', tmp$Fiber, sep = "_")
colnames(struc_MD_df) <- names
rownames(struc_MD_df) <- gsub("_fibers\\.csv$", "", basename(myfile))

struc_MO <- list()
for(i in 1:length(myfile)){
  tmp <- read.csv(myfile[i], sep = "", header = TRUE)
  struc_MO[[i]] <- t(tmp$MO) 
  
}
struc_MO_df <- do.call(rbind, struc_MO)
struc_MO_df <- data.frame(struc_MO_df)
names <- paste('MO', tmp$Fiber, sep = "_")
colnames(struc_MO_df) <- names
rownames(struc_MO_df) <- gsub("_fibers\\.csv$", "", basename(myfile))

struc_L1 <- list()
for(i in 1:length(myfile)){
  tmp <- read.csv(myfile[i], sep = "", header = TRUE)
  struc_L1[[i]] <- t(tmp$L1) 
  
}
struc_L1_df <- do.call(rbind, struc_L1)
struc_L1_df <- data.frame(struc_L1_df)
names <- paste('L1', tmp$Fiber, sep = "_")
colnames(struc_L1_df) <- names
rownames(struc_L1_df) <- gsub("_fibers\\.csv$", "", basename(myfile))

struc_L2 <- list()
for(i in 1:length(myfile)){
  tmp <- read.csv(myfile[i], sep = "", header = TRUE)
  struc_L2[[i]] <- t(tmp$L2) 
  
}
struc_L2_df <- do.call(rbind, struc_L2)
struc_L2_df <- data.frame(struc_L2_df)
names <- paste('L2', tmp$Fiber, sep = "_")
colnames(struc_L2_df) <- names
rownames(struc_L2_df) <- gsub("_fibers\\.csv$", "", basename(myfile))

struc_L3 <- list()
for(i in 1:length(myfile)){
  tmp <- read.csv(myfile[i], sep = "", header = TRUE)
  struc_L3[[i]] <- t(tmp$L3) 
  
}
struc_L3_df <- do.call(rbind, struc_L3)
struc_L3_df <- data.frame(struc_L3_df)
names <- paste('L3', tmp$Fiber, sep = "_")
colnames(struc_L3_df) <- names
rownames(struc_L3_df) <- gsub("_fibers\\.csv$", "", basename(myfile))
fibers_df <- cbind(struc_FA_df, 
                       struc_MD_df, 
                       struc_MO_df,
                       struc_L1_df,
                       struc_L2_df,
                       struc_L3_df)
write.csv(fibers_df, str_c(mypath, "fibers_", dat_idx, "_df.csv"))

# functions
struc_IC25 <- read_excel(str_c(mypath, dat_idx, str_c('/functions/', dat_idx, "_ICs_25.xlsx")))
struc_IC25 <- as.data.frame(struc_IC25)
names<- struc_IC25[[1]]
struc_IC25 <- struc_IC25[,-1]
rownames(struc_IC25) <- names
struc_IC25 <- struc_IC25[,c('NODEamps25 1',
                                      'NODEamps25 2',
                                      'NODEamps25 6',
                                      'NODEamps25 7',
                                      'NODEamps25 9',
                                      'NODEamps25 16',
                                      'NODEamps25 18',
                                      'NET25 11',
                                      'NET25 12',
                                      'NET25 16',
                                      'NET25 35',
                                      'NET25 60',
                                      'NET25 73',
                                      'NET25 82',
                                      'NET25 135',
                                      'NET25 144','NET25 148',
                                      'NET25 150','NET25 156',
                                      'NET25 162')]

struc_IC100 <- read_excel(str_c(mypath, dat_idx, str_c('/functions/', dat_idx, "_ICs_100.xlsx")))
struc_IC100 <- as.data.frame(struc_IC100)
names<- struc_IC100[[1]]
struc_IC100 <- struc_IC100[,-1]
rownames(struc_IC100) <- names
struc_IC100 <- struc_IC100[,c('NET100 16',
                                      'NET100 31',
                                      'NET100 55',
                                      'NET100 57',
                                      'NET100 78',
                                      'NET100 99',
                                      'NET100 113','NET100 118',
                                      'NET100 152','NET100 161',
                                      'NET100 185','NET100 188',
                                      'NET100 207','NET100 208',
                                      'NET100 211','NET100 218',
                                      'NET100 231','NET100 233',
                                      'NET100 239','NET100 243',
                                      'NET100 268','NET100 283',
                                      'NET100 284','NET100 285',
                                      'NET100 292','NET100 312',
                                      'NET100 322','NET100 328',
                                      'NET100 344','NET100 357',
                                      'NET100 380','NET100 405',
                                      'NET100 412','NET100 419',
                                      'NET100 437','NET100 443',
                                      'NET100 458','NET100 487',
                                      'NET100 488','NET100 512',
                                      'NET100 569','NET100 595',
                                      'NET100 605','NET100 615',
                                      'NET100 647','NET100 653',
                                      'NET100 656','NET100 658',
                                      'NET100 663','NET100 676',
                                      'NET100 677','NET100 748',
                                      'NET100 778','NET100 795',
                                      'NET100 807','NET100 821',
                                      'NET100 828','NET100 848',
                                      'NET100 856','NET100 860',
                                      'NET100 872','NET100 891',
                                      'NET100 900','NET100 903',
                                      'NET100 929','NET100 994',
                                      'NET100 1000','NET100 1016',
                                      'NET100 1018','NET100 1030',
                                      'NET100 1060','NET100 1072',
                                      'NET100 1078','NET100 1079',
                                      'NET100 1090','NET100 1129',
                                      'NET100 1132','NET100 1165',
                                      'NET100 1176','NET100 1196',
                                      'NET100 1200','NET100 1212',
                                      'NET100 1220','NET100 1221',
                                      'NET100 1231','NET100 1287',
                                      'NET100 1293','NET100 1328',
                                      'NET100 1344','NET100 1349',
                                      'NET100 1365','NET100 1382',
                                      'NET100 1385','NET100 1391',
                                      'NET100 1394','NET100 1407',
                                      'NET100 1428','NET100 1438')]
funs_df <- cbind(struc_IC25, 
                        struc_IC100)
write.csv(funs_df, str_c(mypath, "functions_", dat_idx, "_df.csv"))
