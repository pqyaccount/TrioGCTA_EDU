# Qiyuan Peng
# Merge GPA data (data from 1_01_data_preprocessing_gpa_EDUCATION_F_UTD_KURS_5.0) 
# with MoBa genetic data to get the trios dataset
# Reference: prep_mobadata_trio_data.R from Espen Eilertsen

rm(list = ls())

#### 1. Library and Set up (Obtain the relevant data file) ####

library(tidyverse)
library(data.table)
library(haven)

#### 1.1 Get individual SNPs information from fam file) ####

genedir <- "N:/durable/data/genetics/MoBaPsychGen_v1_1mil"
genefile <- "MoBaPsychGen_v1_1m"

fam <- read.table(paste0(genedir, "/", genefile, ".fam"), sep = " ", na.strings = "0",
                  col.names = c("fid", "iid", "father", "mother", "sex", "phenotype")) 
dim(fam)
glimpse(fam)
head(fam)

#### 1.2 Get genotype tech info ####

mobadir <- "N:/durable/data/moba"

ce_child <- read.csv(paste0(mobadir, "/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Child_20220829.csv"), 
                     na.strings = c("", " "))  
ce_mother <- read.csv(paste0(mobadir, "/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Mother_20220829.csv"), 
                      na.strings = c("", " "))  
ce_father <- read.csv(paste0(mobadir, "/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Father_20220829.csv"), 
                      na.strings = c("", " "))  

ce_child <- ce_child %>%
  mutate(PREG_ID_2601 = as.character(PREG_ID_2601),
         BARN_NR = as.character(BARN_NR))

#### 2. Merge the valid genotyped data(from .fam) and MoBa data to get correct list of individual ####

c <- inner_join(fam, ce_child, by = c("iid" = "SENTRIX_ID"))
m <- inner_join(fam, ce_mother, by = c("iid" = "SENTRIX_ID"))
f <- inner_join(fam, ce_father, by = c("iid" = "SENTRIX_ID"))

# Merge the data from c and m （child+mother）
cm <- inner_join(c, m, by = c("mother" = "iid"), suffix = c("", "_m"))  

# Merge the dataframe cm and f, then, we get the family dataframe
# One row, we get the information of child, his or her mother and father
cmf <- inner_join(cm, f, by = c("father" = "iid"), suffix = c("_c", "_f"))  

colnames(cmf)[colnames(cmf) == "father"] <- "father_c" 

#### 3. Link registered, cleaned and normailized gpa data with MoBa IDs ####

load("M:/p805-qiyuanp/TrioGCTA/data/1_01_data_gpa_normalized.rda")
gpa <- gpa_normalizeddata
rm(gpa_normalizeddata)

linkage <- fread(paste0(mobadir, "/linkage/merged_IDs/MoBa_SSB_IDs_20250317.csv"))  

linkage <- linkage %>%
  mutate(PREG_ID_2601 = as.character(PREG_ID_2601),
         CHILD_NR = as.character(CHILD_NR))

gpa_moba <- inner_join(linkage, gpa, by = "w19_0634_lnr")  

#### 4. Link the GPA data with the genetic data from MoBa based on the child ("PREG_ID_2601", "BARN_NR") ####

glimpse(cmf)
head(cmf)

gpa_genetic_child <- inner_join(gpa_moba, cmf,
                                by =c ("SENTRIX_ID" = "iid",
                                       "PREG_ID_2601" = "PREG_ID_2601",
                                       "CHILD_NR" = "BARN_NR")
                                )

na_counts <- colSums(is.na(gpa_genetic_child))
print(na_counts)

missing_sex_c <- gpa_genetic_child %>%
  filter(is.na(sex_c))
print(missing_sex_c)

# Load the demographics from SSB
POPULATION_FASTE_OPPLYSNINGER <- read.csv("N:/durable/data/registers/SSB/01_data/data_v5.0/CORE/csv/POPULATION_FASTE_OPPLYSNINGER_reduced.csv", 
            na = c("NA", ""))
demographics <- POPULATION_FASTE_OPPLYSNINGER
rm(POPULATION_FASTE_OPPLYSNINGER)

# Match sex and birthyear information from demographics to gpa_genetic_child
demographics <- demographics %>%
  rename(
    sex = kjoenn,           
    birthdate = foedsels_aar_mnd   
  ) %>%
  mutate(
    birthyear = substr(birthdate, 1, 4)  
  )
head(demographics)

gpa_genetic_child <- gpa_genetic_child %>%
  left_join(demographics, by = "w19_0634_lnr") %>%
  select(-sex_c)

# Recheck the missing values
na_counts <- colSums(is.na(gpa_genetic_child))
print(na_counts)

#### 5. Filter for only one child per family to fit the Trio model ####

#### 5.1 Check if there are totally duplicated individuals and delete the duplicated individuals ####
any_duplicated <- any(duplicated(gpa_genetic_child) | duplicated(gpa_genetic_child, fromLast = TRUE))
print(any_duplicated)  


#### 5.2 Only keep the families who have the fid_c=fid_m=fid_f ####

indices <- which(gpa_genetic_child$fid_c != gpa_genetic_child$fid_m | 
                   gpa_genetic_child$fid_c != gpa_genetic_child$fid_f |
                   gpa_genetic_child$fid_m != gpa_genetic_child$fid_f)

observations_with_different_values <- gpa_genetic_child[indices, ]

gpa_genetic_child_samefid <- gpa_genetic_child[-indices, ]  

#### 5.3 Check for multiple children in a single birth ####

pregnancy_morekid <- gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$PREG_ID_2601) | 
                                                 duplicated(gpa_genetic_child_samefid$PREG_ID_2601, fromLast = TRUE), ]$PREG_ID_2601

gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$PREG_ID_2601) | 
                            duplicated(gpa_genetic_child_samefid$PREG_ID_2601, fromLast = TRUE), ]

length(gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$PREG_ID_2601)|
                                   duplicated(gpa_genetic_child_samefid$PREG_ID_2601, fromLast = TRUE), ]$PREG_ID_2601)

#### 5.4 Check for siblings on the mother side ####

# M_ID_2601
mother_morekid <- gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$M_ID_2601) |
                                              duplicated(gpa_genetic_child_samefid$M_ID_2601, fromLast = TRUE), ]$M_ID_2601

gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$M_ID_2601) | 
                            duplicated(gpa_genetic_child_samefid$M_ID_2601, fromLast = TRUE), ]

length(gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$M_ID_2601)|
                                   duplicated(gpa_genetic_child_samefid$M_ID_2601, fromLast = TRUE), ]$M_ID_2601)

# mother_c
mother_morekid_fam <- gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$mother_c) |
                                                  duplicated(gpa_genetic_child_samefid$mother_c, fromLast = TRUE), ]$mother_c

gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$mother_c) | 
                            duplicated(gpa_genetic_child_samefid$mother_c, fromLast = TRUE), ]

length(gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$mother_c)|
                                   duplicated(gpa_genetic_child_samefid$mother_c, fromLast = TRUE), ]$mother_c)

#### 5.5 Check for siblings on the father side ####

# F_ID_2601
father_morekid <- gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$F_ID_2601) |
                                              duplicated(gpa_genetic_child_samefid$F_ID_2601, fromLast = TRUE), ]$F_ID_2601

gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$F_ID_2601) | 
                            duplicated(gpa_genetic_child_samefid$F_ID_2601, fromLast = TRUE), ]

length(gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$F_ID_2601)|
                                   duplicated(gpa_genetic_child_samefid$F_ID_2601, fromLast = TRUE), ]$F_ID_2601)

# father_c
father_morekid_fam <- gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$father_c) |
                                                  duplicated(gpa_genetic_child_samefid$father_c, fromLast = TRUE), ]$father_c

gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$father_c) | 
                            duplicated(gpa_genetic_child_samefid$father_c, fromLast = TRUE), ]

length(gpa_genetic_child_samefid[duplicated(gpa_genetic_child_samefid$father_c)|
                                   duplicated(gpa_genetic_child_samefid$father_c, fromLast = TRUE), ]$father_c)

#### 5.6 Keep only one child record for one family ####

gpa_genetic_child_filt <- gpa_genetic_child_samefid %>%
  arrange(birthdate) %>%
  group_by(M_ID_2601) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(F_ID_2601) %>%
  slice_head(n = 1) %>%
  ungroup()

length(gpa_genetic_child_filt[duplicated(gpa_genetic_child_filt$PREG_ID_2601), ]$PREG_ID_2601)

length(gpa_genetic_child_filt[duplicated(gpa_genetic_child_filt$M_ID_2601), ]$M_ID_2601)
length(gpa_genetic_child_filt[duplicated(gpa_genetic_child_filt$mother_c), ]$mother_c)

length(gpa_genetic_child_filt[duplicated(gpa_genetic_child_filt$F_ID_2601), ]$F_ID_2601)
length(gpa_genetic_child_filt[duplicated(gpa_genetic_child_filt$father_c), ]$father_c)

# Check the number of missing data from different column
sapply(gpa_genetic_child_filt, function(x) sum(is.na(x)))

#### 6. Rearrange according to the .fam file ####

# Construct one column named famid in gpa_genetic_child_filt data frame 
# that contains a sequence of consecutive integers starting at 1 and ending at the number of rows in the data frame.
gpa_genetic_child_filt$famid = 1:nrow(gpa_genetic_child_filt)
max(gpa_genetic_child_filt$famid)  

colnames(gpa_genetic_child_filt)[colnames(gpa_genetic_child_filt) == "SENTRIX_ID"] <- "iid" 

# We want the dataframe of Trio-GCTA will be: 
# row 1:n mother
# row n+1:n+n father
# row 2n+1:2n+n child

# mother
m_arr <- m %>%
  filter(iid %in% gpa_genetic_child_filt$mother_c) %>%
  left_join(gpa_genetic_child_filt %>% select(mother_c, famid), by = c("iid" = "mother_c")) %>%
  mutate(PREG_ID_2601 = NA, BARN_NR = NA, gpa_original=NA, gpa_normalized = NA, year = NA, birthyear = NA) %>%
  select(-M_ID_2601) %>%
  arrange(famid) 

glimpse(m_arr)
max(m_arr$famid) 

# father
f_arr <- f %>%
  filter(iid %in% gpa_genetic_child_filt$father_c) %>%
  left_join(gpa_genetic_child_filt %>% select(father_c, famid), by = c("iid" = "father_c")) %>%
  mutate(PREG_ID_2601 = NA, BARN_NR = NA, gpa_original=NA, gpa_normalized = NA, year = NA, birthyear = NA) %>%
  select(-F_ID_2601) %>%
  arrange(famid) 

glimpse(f_arr)
max(f_arr$famid) 

# child
# Here we don't use the sex from fam file because there are some missing values.
c_arr <- c %>%
  filter(iid %in% gpa_genetic_child_filt$iid) %>%
  select(-sex) %>%
  left_join(gpa_genetic_child_filt %>% select(iid, famid, sex, gpa_original, gpa_normalized, year, birthyear), by = "iid") %>%
  arrange(famid)  

glimpse(c_arr)
max(c_arr$famid) 

table(c_arr$birthyear)

gpa_genetic_child_trios <-  dplyr::bind_rows(m_arr, f_arr, c_arr)

names(gpa_genetic_child_trios)
str(gpa_genetic_child_trios)
glimpse(gpa_genetic_child_trios)
head(gpa_genetic_child_trios)

summary(gpa_genetic_child_trios)

sapply(gpa_genetic_child_trios, function(x) sum(is.na(x)))

#### 7. Write file ####
size <- "gpa_28748obs"
write.csv(gpa_genetic_child_trios, paste0("data/", genefile, "_", size, ".moba"), row.names = F)