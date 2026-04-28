# Qiyuan Peng
# R scripts for preparing the National Assessment Grades (8th)

#### 1. Basic Set Up for Data Preparation ####

rm(list = ls())

library(tidyverse)  
library(ggplot2)   
library(dplyr)     
library(data.table) 
library(foreign)

NT <- data.table::fread("N:/durable/data/registers/SSB/01_data/data_v5.0/EDUCATION_VGS_GRS/csv/EDUCATION_NASJONALE_PROVER.CSV", 
                        na.strings = c("NA", ""))

NT <- NT %>% 
  select(w19_0634_lnr, testyear = AARGANG, schoolcity = SKOLEKOM, schooltype = SKOLETYPE, 
         attendance = DELTATTSTATUS, testtype = PROVE, scorelevel = MESTRINGSNIVAA, score = POENG)

str(NT)
head(NT) 
summary(NT)

#### 2. Select children who attended national assessment and view the data distributions ####

# When attendance == "D", it means the students attend the test
NT <- NT %>% 
  filter(attendance == "D") 

# View the number of records for each type of test
table(NT$testtype)

# Group and summarize the data to see the dataset structure
summary_by_year_type <- NT %>%
  filter(testtype %in% c("NPENG05", "NPLES05", "NPREG05", "NPENG08", "NPLES08", "NPREG08", "NPLES09", "NPREG09")) %>%
  group_by(testyear, testtype) %>%
  summarize(
    min_score = min(score, na.rm = TRUE),   
    max_score = max(score, na.rm = TRUE),   
    mean_score = mean(score, na.rm = TRUE), 
    sd_score = sd(score, na.rm = TRUE),     
    .groups = 'drop'
  )

View(summary_by_year_type) 
rm(summary_by_year_type)

table(NT$testtype[NT$testtype %in% c("NPENG05", "NPENG08", "NPLES05", "NPLES08",
                                     "NPLES09", "NPREG05", "NPREG08", "NPREG09")],
      NT$testyear[NT$testtype %in% c("NPENG05", "NPENG08", "NPLES05", "NPLES08",
                                     "NPLES09", "NPREG05", "NPREG08", "NPREG09")])

# For TrioGCTA paper, we only focus on 8th grades.

#### 3. 8th Grade National Assessment ####

#### 3.1 Select and clean the subsample: 8th Grade National Assessment ####

grade08 <- NT %>%
  filter(testtype %in% c("NPENG08", "NPLES08", "NPREG08"))

head(grade08) 
str(grade08)
summary(grade08)

#### 3.2 Check missing data, impossible data, and infinite values ####

grade08 %>% summarise_all(list(~sum(is.na(.))))  
grade08 %>% summarise_all(list(~sum(is.nan(.))))
grade08 %>% summarise_all(list(~sum(is.infinite(.))))

View(subset(grade08, is.na(score)))

grade08 %>%
  summarise(
    lessthan0 = sum(score < 0, na.rm = TRUE),   # number of scores below 0
    equalto0  = sum(score == 0, na.rm = TRUE),  # number of scores exactly 0
    infinite  = sum(is.infinite(score))         # number of infinite values in score
  )

View(subset(grade08, score == 0))
table(grade08$scorelevel[grade08$score == 0])

#### 3.3 Clean the subsample: retain only records where all 8th grade tests occurred in the same year ####

# This part mainly follows Dr. Rob Eves’s script for data cleaning.

# Order observations by test year (earliest first)
grade08 <- grade08[order(grade08$testyear), ]

### Step 1: Check the distribution of duplicates

# Identify duplicates: same child, same subject (across any year)
grade08 <- grade08 %>%
  group_by(w19_0634_lnr, testtype) %>%
  mutate(duplicates = n()) %>%
  ungroup()

table(grade08$duplicates)
table(grade08$duplicates == 2)

# Identify same-year duplicates: same child, same subject, same test year
grade08 <- grade08 %>%
  group_by(w19_0634_lnr, testtype, testyear) %>%
  mutate(duplicates_sameyear = n()) %>%
  ungroup()

table(grade08$duplicates_sameyear)
View(subset(grade08, duplicates == 2))
View(subset(grade08, duplicates_sameyear == 2))

### Step 2: Resolve duplicates

# If the same test was taken multiple times:
# - In different years: keep the first (earliest year)
# - In the same year: take the average score (Neil Davies / Rob Eves approach)

grade08 <- grade08 %>%
  group_by(w19_0634_lnr, testtype) %>%
  mutate(
    num_tests = n(),
    num_tests_ID = row_number(),
    score_alt = ifelse(
      num_tests > 1 & (testyear[1] == testyear[2]),
      mean(score, na.rm = TRUE),
      score[1]
    )
  ) %>%
  ungroup()

grade08_filtered <- subset(grade08, num_tests_ID == 1)

### Step 3: Unify the test year

# If a child took English, Norwegian (Reading), and Math in different years,
# we use the mode (most frequent year) as the unified test year.
# If there are multiple modes, the earliest year is chosen.
# This step is necessary for later Z-score calculation.

id_diff_years <- grade08_filtered %>%
  select(w19_0634_lnr, testyear) %>%
  group_by(w19_0634_lnr) %>%
  summarise(n_testyears = n_distinct(testyear)) %>%
  filter(n_testyears > 1) %>%
  pull(w19_0634_lnr)

grade08_diff_years <- grade08_filtered %>%
  filter(w19_0634_lnr %in% id_diff_years)

grade08_filtered$testyear_alt <- as.numeric(grade08_filtered$testyear)

grade08_filtered <- grade08_filtered %>%
  group_by(w19_0634_lnr) %>%
  mutate(
    testyear_alt = {
      freq <- table(testyear)
      max_freq <- max(freq)
      modes <- as.numeric(names(freq)[freq == max_freq])
      min(modes)
    }
  ) %>%
  ungroup()

grade08_filtered %>%
  select(w19_0634_lnr, testyear_alt) %>%
  group_by(w19_0634_lnr) %>%
  summarise(n_years = n_distinct(testyear_alt)) %>%
  filter(n_years > 1) %>%
  summarise(n_individuals = n())   

unique_ids <- unique(grade08_filtered$w19_0634_lnr[grade08_filtered$testyear_alt != as.numeric(grade08_filtered$testyear)])
length(unique_ids)

#### 3.4 Reshape the subsample: wide format and Z-score calculation ####

# Keep the essential variables: child ID, subject, final score, and unified test year
grade08_filtered_reduced <- subset(grade08_filtered, select = c("w19_0634_lnr", "testtype", "score_alt", "testyear_alt"))

# Reshape from long to wide format
grade08_filtered_reduced_wide <- reshape(
  as.data.frame(grade08_filtered_reduced),
  idvar = c("w19_0634_lnr", "testyear_alt"),
  timevar = "testtype",
  direction = "wide"
)

# Rename the test score columns for clarity
grade08_filtered_reduced_wide <- grade08_filtered_reduced_wide %>%
  rename(
    testyear  = testyear_alt,
    math08    = score_alt.NPREG08,
    reading08 = score_alt.NPLES08,
    english08 = score_alt.NPENG08
  )

# Standardize scores within each test year (Z-scores)
grade08_filtered_reduced_wide$math08_z <- ave(
  grade08_filtered_reduced_wide$math08,
  grade08_filtered_reduced_wide$testyear,
  FUN = scale
)

grade08_filtered_reduced_wide$reading08_z <- ave(
  grade08_filtered_reduced_wide$reading08,
  grade08_filtered_reduced_wide$testyear,
  FUN = scale
)

grade08_filtered_reduced_wide$english08_z <- ave(
  grade08_filtered_reduced_wide$english08,
  grade08_filtered_reduced_wide$testyear,
  FUN = scale
)

# Inspect observations with missing test scores
grade08_filtered_reduced_wide %>% summarise_all(list(~sum(is.na(.)))) 

#### 3.5 Save the subsample dataset ####

file_path <- "M:/p805-qiyuanp/TrioGCTA/data/"

# Grade 8 Math Data without NA.
math08_data <- grade08_filtered_reduced_wide %>%
  filter(!is.na(math08)) %>%
  select(w19_0634_lnr, testyear, math08, math08_z)
save(math08_data, file = paste0(file_path, "1_02_data_math08_normalized.rda"))

# Grade 8 Reading Data without NA.
reading08_data <- grade08_filtered_reduced_wide %>%
  filter(!is.na(reading08)) %>%
  select(w19_0634_lnr, testyear, reading08, reading08_z)
save(reading08_data, file = paste0(file_path, "1_02_data_reading08_normalized.rda"))

# Grade 8 English Data without NA.
english08_data <- grade08_filtered_reduced_wide %>%
  filter(!is.na(english08)) %>%
  select(w19_0634_lnr, testyear, english08, english08_z)
save(english08_data, file = paste0(file_path, "1_02_data_english08_normalized.rda"))

# Grade 8 Math + Reading without NA.
math08_reading08_data <- grade08_filtered_reduced_wide %>%
  filter(!is.na(math08) & !is.na(reading08)) %>%
  select(w19_0634_lnr, testyear, math08, math08_z, reading08, reading08_z)
save(math08_reading08_data, file = paste0(file_path, "1_02_data_math08_reading08_normalized.rda"))

# Grade 8 Math + Reading + English without NA.
complete08_data <- grade08_filtered_reduced_wide %>%
  filter(!is.na(math08) & !is.na(reading08) & !is.na(english08)) %>%
  select(w19_0634_lnr, testyear,
         math08, math08_z,
         reading08, reading08_z,
         english08, english08_z)
save(complete08_data, file = paste0(file_path, "1_02_data_math08_reading08_english08_normalized.rda"))