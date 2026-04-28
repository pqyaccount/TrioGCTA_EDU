# Qiyuan Peng 
# Merge Grade 8th Reading data (prepared in 1_02_data_preprocessing_nationaltest_EDUCATION_NASJONALE_PROVER_5.0)  
# with MoBa genetic data to construct the trios dataset  

rm(list = ls())

#### 1. Library and Set up (Obtain the relevant data file) ####

library(tidyverse)
library(data.table)
library(haven)

#### 1.1 Load focal individual's phenotype data ####

# Define the file name of the registered, cleaned, and normalized outcome data (e.g., reading08)
phenotype_dataset <- "1_02_data_reading08_normalized.rda"
phenotype_path <- paste0("M:/p805-qiyuanp/TrioGCTA/data/", phenotype_dataset)

load(phenotype_path)

reading08 <- reading08_data
rm(reading08_data)

#### 1.2 Get individual SNPs information from .fam file ####

# Define the path and filename of the genotype data
genedir <- "N:/durable/data/genetics/MoBaPsychGen_v1_1mil"
genefile <- "MoBaPsychGen_v1_1m"

# Load basic family and trio structure from the .fam file (PLINK format)
fam <- read.table(paste0(genedir, "/", genefile, ".fam"),
                  sep = "", na.strings = "0",
                  col.names = c("fid", "iid", "father", "mother", "sex", "phenotype"))

#### 1.3 Retrieve sample identifiers from MoBa genotyping registry ####

# Set MoBa genotyping registry directory
mobadir <- "N:/durable/data/moba/"

# Load genotyping identifiers for child, mother, and father
ce_child <- read.csv(paste0(mobadir, "/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Child_20220829.csv"), 
                     na.strings = c("", " ")) 
ce_mother <- read.csv(paste0(mobadir, "/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Mother_20220829.csv"), 
                      na.strings = c("", " ")) 
ce_father <- read.csv(paste0(mobadir, "/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Father_20220829.csv"), 
                      na.strings = c("", " "))  

# Convert IDs to character to ensure compatibility in merging
ce_child <- ce_child %>%
  mutate(PREG_ID_2601 = as.character(PREG_ID_2601),
         BARN_NR = as.character(BARN_NR))

#### 2. Merge valid genotyped data from .fam with MoBa registry to identify trio members ####

c <- inner_join(fam, ce_child, by = c("iid" = "SENTRIX_ID"))
m <- inner_join(fam, ce_mother, by = c("iid" = "SENTRIX_ID"))
f <- inner_join(fam, ce_father, by = c("iid" = "SENTRIX_ID"))

# Merge the data from c, m and f to get the family data frame
# Merge the data from c and m
cm <- inner_join(c, m, by = c("mother" = "iid"), suffix = c("", "_m"))  

# Merge the dataframe cm and f, then, we get the family dataframe
# One row, we get the information of child, his or her mother and father
cmf <- inner_join(cm, f, by = c("father" = "iid"), suffix = c("_c", "_f"))  

colnames(cmf)[colnames(cmf) == "father"] <- "father_c" 

#### 3. Link registered, cleaned and normailized reading08 data with MoBa IDs ####

linkage <- fread(paste0(mobadir, "/linkage/merged_IDs/MoBa_SSB_IDs_20250317.csv"))  

linkage <- linkage %>%
  mutate(PREG_ID_2601 = as.character(PREG_ID_2601),
         CHILD_NR = as.character(CHILD_NR))

reading08_moba <- inner_join(linkage, reading08, by = "w19_0634_lnr")  

#### 4. Link the reading08 data with the genetic data from MoBa based on the child ("PREG_ID_2601", "BARN_NR") ####

reading08_genetic_child <- inner_join(reading08_moba, cmf,
                                by =c ("SENTRIX_ID" = "iid",
                                       "PREG_ID_2601" = "PREG_ID_2601",
                                       "CHILD_NR" = "BARN_NR")
                                )

na_counts <- colSums(is.na(reading08_genetic_child))
print(na_counts)

missing_sex_c <- reading08_genetic_child %>%
  filter(is.na(sex_c))
print(missing_sex_c)

POPULATION_FASTE_OPPLYSNINGER <- read.csv("N:/durable/data/registers/SSB/01_data/data_v5.0/CORE/csv/POPULATION_FASTE_OPPLYSNINGER_reduced.csv", 
            na = c("NA", ""))
demographics <- POPULATION_FASTE_OPPLYSNINGER
rm(POPULATION_FASTE_OPPLYSNINGER)

# Match sex and birthyear information from demographics to reading08_genetic_child
demographics <- demographics %>%
  rename(
    sex = kjoenn,          
    birthdate = foedsels_aar_mnd 
  ) %>%
  mutate(
    birthyear = substr(birthdate, 1, 4)
  )
head(demographics)

reading08_genetic_child <- reading08_genetic_child %>%
  left_join(demographics, by = "w19_0634_lnr") %>%
  select(-sex_c)

na_counts <- colSums(is.na(reading08_genetic_child))
print(na_counts)

#### 5. Filter for only one child per family to fit the Trio model ####

#### 5.1 Check if there are totally duplicated individuals and delete the duplicated individuals ####

any_duplicated <- any(duplicated(reading08_genetic_child) | duplicated(reading08_genetic_child, fromLast = TRUE))
print(any_duplicated) 

#### 5.2 Only keep the families who have the fid_c=fid_m=fid_f ####

indices <- which(reading08_genetic_child$fid_c != reading08_genetic_child$fid_m | 
                   reading08_genetic_child$fid_c != reading08_genetic_child$fid_f |
                   reading08_genetic_child$fid_m != reading08_genetic_child$fid_f)

observations_with_different_values <- reading08_genetic_child[indices, ]

reading08_genetic_child_samefid <- reading08_genetic_child[-indices, ]  

#### 5.3 Check for multiple children in a single birth ####

pregnancy_morekid <- reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$PREG_ID_2601) | 
                                                    duplicated(reading08_genetic_child_samefid$PREG_ID_2601, fromLast = TRUE), ]$PREG_ID_2601

reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$PREG_ID_2601) | 
                               duplicated(reading08_genetic_child_samefid$PREG_ID_2601, fromLast = TRUE), ]

length(reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$PREG_ID_2601)|
                                      duplicated(reading08_genetic_child_samefid$PREG_ID_2601, fromLast = TRUE), ]$PREG_ID_2601)

#### 5.4 Check for siblings on the mother side ####

# M_ID_2601
mother_morekid <- reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$M_ID_2601) |
                                                 duplicated(reading08_genetic_child_samefid$M_ID_2601, fromLast = TRUE), ]$M_ID_2601

reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$M_ID_2601) | 
                               duplicated(reading08_genetic_child_samefid$M_ID_2601, fromLast = TRUE), ]

length(reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$M_ID_2601)|
                                      duplicated(reading08_genetic_child_samefid$M_ID_2601, fromLast = TRUE), ]$M_ID_2601)

# mother_c
mother_morekid_fam <- reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$mother_c) |
                                                     duplicated(reading08_genetic_child_samefid$mother_c, fromLast = TRUE), ]$mother_c

reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$mother_c) | 
                               duplicated(reading08_genetic_child_samefid$mother_c, fromLast = TRUE), ]

length(reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$mother_c)|
                                      duplicated(reading08_genetic_child_samefid$mother_c, fromLast = TRUE), ]$mother_c)

#### 5.5 Check for siblings on the father side ####

# F_ID_2601
father_morekid <- reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$F_ID_2601) |
                                                 duplicated(reading08_genetic_child_samefid$F_ID_2601, fromLast = TRUE), ]$F_ID_2601

reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$F_ID_2601) | 
                               duplicated(reading08_genetic_child_samefid$F_ID_2601, fromLast = TRUE), ]

length(reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$F_ID_2601)|
                                      duplicated(reading08_genetic_child_samefid$F_ID_2601, fromLast = TRUE), ]$F_ID_2601)

# father_c
father_morekid_fam <- reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$father_c) |
                                                     duplicated(reading08_genetic_child_samefid$father_c, fromLast = TRUE), ]$father_c

reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$father_c) | 
                               duplicated(reading08_genetic_child_samefid$father_c, fromLast = TRUE), ]

length(reading08_genetic_child_samefid[duplicated(reading08_genetic_child_samefid$father_c)|
                                      duplicated(reading08_genetic_child_samefid$father_c, fromLast = TRUE), ]$father_c)

#### 5.6 Keep only one child record for one family ####

reading08_genetic_child_filt <- reading08_genetic_child_samefid %>%
  arrange(birthdate) %>%
  group_by(M_ID_2601) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(F_ID_2601) %>%
  slice_head(n = 1) %>%
  ungroup()

length(reading08_genetic_child_filt[duplicated(reading08_genetic_child_filt$PREG_ID_2601), ]$PREG_ID_2601)

length(reading08_genetic_child_filt[duplicated(reading08_genetic_child_filt$M_ID_2601), ]$M_ID_2601)
length(reading08_genetic_child_filt[duplicated(reading08_genetic_child_filt$mother_c), ]$mother_c)

length(reading08_genetic_child_filt[duplicated(reading08_genetic_child_filt$F_ID_2601), ]$F_ID_2601)
length(reading08_genetic_child_filt[duplicated(reading08_genetic_child_filt$father_c), ]$father_c)

sapply(reading08_genetic_child_filt, function(x) sum(is.na(x)))

#### 6. Rearrange according to the .fam file ####

# Construct one column named famid in reading08_genetic_child_filt data frame 
# that contains a sequence of consecutive integers starting at 1 and ending at the number of rows in the data frame.
reading08_genetic_child_filt$famid = 1:nrow(reading08_genetic_child_filt)
max(reading08_genetic_child_filt$famid)  

colnames(reading08_genetic_child_filt)[colnames(reading08_genetic_child_filt) == "SENTRIX_ID"] <- "iid" 

# We want the dataframe of Trio-GCTA will be: 
# row 1:n mother
# row n+1:n+n father
# row 2n+1:2n+n child

# mother
m_arr <- m %>%
  filter(iid %in% reading08_genetic_child_filt$mother_c) %>%
  left_join(reading08_genetic_child_filt %>% select(mother_c, famid), by = c("iid" = "mother_c")) %>%
  mutate(PREG_ID_2601 = NA, BARN_NR = NA, reading08=NA, reading08_z = NA, testyear = NA, birthyear = NA) %>%
  select(-M_ID_2601) %>%
  arrange(famid)

glimpse(m_arr)
max(m_arr$famid) 

# father
f_arr <- f %>%
  filter(iid %in% reading08_genetic_child_filt$father_c) %>%
  left_join(reading08_genetic_child_filt %>% select(father_c, famid), by = c("iid" = "father_c")) %>%
  mutate(PREG_ID_2601 = NA, BARN_NR = NA, reading08=NA, reading08_z = NA, testyear = NA, birthyear = NA) %>%
  select(-F_ID_2601) %>%
  arrange(famid) 

glimpse(f_arr)
max(f_arr$famid) 

# child
# Here we don't use the sex from fam file because there are some missing values.
c_arr <- c %>%
  filter(iid %in% reading08_genetic_child_filt$iid) %>%
  select(-sex) %>%
  left_join(reading08_genetic_child_filt %>% select(iid, famid, sex, reading08, reading08_z, testyear, birthyear), by = "iid") %>%
  arrange(famid) 

glimpse(c_arr)
max(c_arr$famid) 

table(c_arr$birthyear)

reading08_genetic_child_trios <-  dplyr::bind_rows(m_arr, f_arr, c_arr)

names(reading08_genetic_child_trios)
str(reading08_genetic_child_trios)
glimpse(reading08_genetic_child_trios)
head(reading08_genetic_child_trios)

summary(reading08_genetic_child_trios)

sapply(reading08_genetic_child_trios, function(x) sum(is.na(x)))

#### 7. Write file ####
size <- "reading08_33943obs"
write.csv(reading08_genetic_child_trios, paste0("M:/p805-qiyuanp/TrioGCTA/data/", genefile, "_", size, ".moba"), row.names = F)