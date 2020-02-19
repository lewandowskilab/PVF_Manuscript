
###### Notes #####
# This script generates new dataframes from the objects I received. 
# Making these allowed me to get to know accustom myself with the data and create a DF I know. 
# It makes a new dataframe for Utrecht, Ulm, Leuven and combines them into a big dataframe for combined analyses. 
# Tests are included. 
# Date last modified: 2020/01/31 # F.S.Sanders

### Load workspace and files ####
setwd("~/Documents/ALS-research/")

#For the scaled data
load(file="VASC2 ALS Ulm Leuven and Utrecht data median scaled.Rdata")

als.ma.utrecht <- als.utrecht.scale
als.ma.ulm <- als.ulm
als.ma.leuven <- als.leuven.scale

###### I run all packages I generally use in other scripts working with 
library(SBA)
library(dplyr)
library(gplots)
library(beeswarm)
library(ggbeeswarm)
library(lubridate)
library(survival)
library(survminer)
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(testthat)

#### Proteins of interest ####
proteins <- c(53,123,193)

################## Utrecht Dataframe #######################

# Most of the code to generate specific columns are adapted from original script made by Julia and Anna
# Survival calculations
# Calculating survival months for alive patients, last check-up in April 2018
year <- substr(as.character(als.ma.utrecht@sample$Date_of_onset), 1, 4)
year.onset <- as.numeric(year)
month <- substr(as.character(als.ma.utrecht@sample$Date_of_onset), 6, 7)
all.months <- c("01","02","03","04","05","06","07","08","09","10","11","12")
month.onset <- match(month,all.months)
year.check <- rep(2018,length(year.onset))
month.check <- rep(4,length(month.onset))
survival.months <- (year.check-year.onset)*12+month.check-month.onset # Calculate months, disregard days
death.year <- substr(as.character(als.ma.utrecht@sample$Date.of.Death), 1, 4)
death.year <- as.numeric(death.year)
death.month <- substr(as.character(als.ma.utrecht@sample$Date.of.Death), 6, 7)
death.month <- match(death.month,all.months)
death.survival <- (death.year-year.onset)*12+death.month-month.onset # Calculate months, disregard days
survival.months[which(is.na(death.survival)==FALSE)] <- death.survival[which(is.na(death.survival)==FALSE)] # Replace the known deaths
als.ma.utrecht@sample$survival_months <- survival.months
died <- rep(FALSE,length(year))
died[which(is.na(death.survival)==FALSE)] <- TRUE
als.ma.utrecht@sample$died <- died
als.ma.utrecht@sample$survobj <- with(als.ma.utrecht@sample, Surv(survival_months, died))

# Site on onset at Bulbar true or false
onset <- rep(FALSE,length(als.ma.utrecht@sample$Site_of_Onset))
onset[which(als.ma.utrecht@sample$Site_of_Onset=="Bulbar")] <- TRUE
als.ma.utrecht@sample$onset_bulbar <- onset

#Get protein MFI values from als.ma.utrecht@X
#Spp1, Col6a1,NEFL 
mfiv <- als.ma.utrecht@X[,proteins]
mfiv <- as.data.frame(mfiv)
colnames(mfiv) <- paste("RawB.",colnames(mfiv),sep = "")

#Creating new dataframe
#Subsetting columns of interest, adding protein values, inculde cohort signature.
als.utrecht.cntrl <- als.ma.utrecht@sample[,c("sample_name","class", "gender", "age_at_onset",
                                        "sampling_age","Date_of_diagnosis", "sampling_date",
                                        "survival_months", "died","survobj", "onset_bulbar")]
als.utrecht.cntrl <- cbind(als.utrecht.cntrl, mfiv)
als.utrecht.cntrl$cohort <- "Utrecht"

#Only cases dataframe
utrechtALS <- subset(als.utrecht.cntrl, als.utrecht.cntrl$class == "Case")

#Adding gap (sampling delay information), "immortal timeline"
utrechtALS <- mutate(utrechtALS, gap = utrechtALS$sampling_age - utrechtALS$age_at_onset)
utrechtALS$gap[utrechtALS$gap == -1] <- 0

# Following function to test if sample, MFI values and survival months are still correctly linked
# Note: Rownames ID are transferred from original dataframe, BUT the sample_id column in original dataframe contains an empty additional attribute.
# Consequently test is without checking attributes so -> check attributes is FALSE.
# The function gives an error if they are not equal
test_that("New DF patients and values are similar to original DF",{
  x <- cbind(as.character(als.ma.utrecht@sample$sample_name), als.ma.utrecht@sample$survival_months ,als.ma.utrecht@X[,53])
  y <- cbind(as.character(als.utrecht.cntrl$sample_name), als.utrecht.cntrl$survival_months ,als.utrecht.cntrl$RawB.53)
  expect_equal(x, y, check.attributes = FALSE)
})

# used this to compare attributes 
# x <- cbind(as.character(als.ma.ulm@sample$sample_name), als.ma.ulm@sample$survival_months ,als.ma.ulm@X[,53])
# y <- cbind(as.character(als.ulm.cntrl$sample_name), als.ulm.cntrl$survival_months ,als.ulm.cntrl$RawB.53)
# str(x)
# str(y)

#################### ULM dataframe #####################

# Adding survival months for alive patients, last check-up in April 2018
als.ma.ulm@sample$died <- als.ma.ulm@sample$died_update_2018
survival_months <- als.ma.ulm@sample$survival_months
died <- als.ma.ulm@sample$died
als.ma.ulm@sample$survobj <- with(als.ma.ulm@sample, Surv(survival_months, died))

# Site on onset at Bulbar true or false
onset <- als.ma.ulm@sample$init_symp_bulbar_update_2018
als.ma.ulm@sample$onset_bulbar <- onset

#G#Get protein MFI values from als.ma.ulm@X
mfiv <- als.ma.ulm@X[,proteins]
mfiv <- as.data.frame(mfiv)
colnames(mfiv) <- paste("RawB.",colnames(mfiv),sep = "")

#Creating new dataframe
#Subsetting columns of interest, adding protein values, inculde cohort signature.
als.ulm.cntrl <- als.ma.ulm@sample[,c("sample_name","class", "gender", "age_at_onset","sampling_age",
                                     "survival_months", "died","survobj" , "onset_bulbar")]
als.ulm.cntrl <- cbind(als.ulm.cntrl, mfiv)
als.ulm.cntrl$cohort <- "Ulm"

#Select patients only, into new DF
ulmALS <- subset(als.ulm.cntrl, als.ulm.cntrl$class == "Case")

#Remove NA survival months, which is present in Ulm data
ulmALS <- ulmALS[!is.na(ulmALS$survival_months),]
#20 rows are removed.

# Add gap variable to account for the delay between onset and sampling.
ulmALS <- mutate(ulmALS, gap = ulmALS$sampling_age - ulmALS$age_at_onset)
ulmALS$gap[ulmALS$gap == -1] <- 0

# Same test as before
test_that("New DF patients and values are similar to original DF",{
  x <- cbind(as.character(als.ma.ulm@sample$sample_name), als.ma.ulm@sample$survival_months ,als.ma.ulm@X[,53])
  y <- cbind(as.character(als.ulm.cntrl$sample_name), als.ulm.cntrl$survival_months ,als.ulm.cntrl$RawB.53)
  expect_equal(x, y, check.attributes = FALSE)
})

###################### Leuven dataframe ###################

# Adding survival months for alive patients, last check-up in April 2018
als.ma.leuven@sample$died <- als.ma.leuven@sample$died
survival_months <- als.ma.leuven@sample$survival_months
died <- als.ma.leuven@sample$died
als.ma.leuven@sample$survobj <- with(als.ma.leuven@sample, Surv(survival_months, died))

# Site on onset at Bulbar true or false
onset <- rep(FALSE,nrow(als.ma.leuven@sample))
onset[which(als.ma.leuven@sample$Site_of_Onset_for_analysis=="Bulbar")] <- TRUE
als.ma.leuven@sample$onset_bulbar <- onset

#Selecting MFI values from als.ma.leuven@X
mfiv <- als.ma.leuven@X[,proteins]
mfiv <- as.data.frame(mfiv)
colnames(mfiv) <- gsub("\\_.*", "", colnames(mfiv))
colnames(mfiv) <- gsub("B", "RawB", colnames(mfiv))
#Creating new dataframe
als.leuven.cntrl <- als.ma.leuven@sample[,c("sample_name","class3", "gender", "Age_of_onset","Age_at_sampling",
                                           "survival_months", "died","survobj" , "onset_bulbar")]

# als.leuven.cntrl <- cbind(als.leuven.cntrl, quartiles)
als.leuven.cntrl <- cbind(als.leuven.cntrl, mfiv)
als.leuven.cntrl <- als.leuven.cntrl %>%
  rename("age_at_onset" = "Age_of_onset", "class" = "class3","sampling_age" = "Age_at_sampling")
als.leuven.cntrl$cohort <- "Leuven"

#Select cases for patient only DF
leuvenALS <- subset(als.leuven.cntrl, als.leuven.cntrl$class == "ALS")

#Remove NA survival months, which is present in leuven data
leuvenALS <- leuvenALS[!is.na(leuvenALS$survival_months),]

#Gap
leuvenALS <- mutate(leuvenALS, gap = leuvenALS$sampling_age - leuvenALS$age_at_onset)
leuvenALS$gap[leuvenALS$gap == -1] <- 0
# Same test as before
test_that("New DF patients and values are similar to original DF",{
  x <- cbind(as.character(als.ma.leuven@sample$sample_name), als.ma.leuven@sample$survival_months ,als.ma.leuven@X[,53])
  y <- cbind(as.character(als.leuven.cntrl$sample_name), als.leuven.cntrl$survival_months ,als.leuven.cntrl$RawB.53)
  expect_equal(x, y, check.attributes = FALSE)
})

####### Combining dataframes #########
# For further analyses we want to combine datasets into one big DF 
# Remove columns not present in Ulm and Leuven
utrechtALS <- select(utrechtALS,-c(Date_of_diagnosis, sampling_date))
als.utrecht.cntrl <- select(als.utrecht.cntrl, -c(Date_of_diagnosis, sampling_date))

# Combine the three cohorts
### NOTE Survobj does not transfer after rbind. Need to recompute the column with Surv function. 
combinedALS <- rbind(utrechtALS,ulmALS,leuvenALS)
combinedALS$survobj <- with(combinedALS, Surv(survival_months, died))
cohorts <- rbind(als.utrecht.cntrl,als.ulm.cntrl,als.leuven.cntrl)

#Clean global environment
rm(als.ma.utrecht)
rm(als.ma.ulm)
rm(als.ma.leuven)
rm(als.leuven.scale)
rm(als.utrecht.scale)
rm(als.ulm)
cat("Script done:\nUtrecht, Ulm, Leuven and cohorts combined loaded in Global Environment.")

# Save essential two dataframes
# save("combinedALS","cohorts", file = "Combined_Cohort_Dataframes.Rdata")
