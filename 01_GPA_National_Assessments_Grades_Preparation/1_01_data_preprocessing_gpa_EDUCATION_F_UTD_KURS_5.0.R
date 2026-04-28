# Qiyuan Peng
# R scripts for preparing GPA data
  
#### 1. Libraries and Set Up ####

rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

gpa_data <- data.table::fread("N:/durable/data/registers/SSB/01_data/data_v5.0/csv/EDUCATION_F_UTD_KURS.csv", 
                        na.strings = c("NA", "")) 

gpa_data <- gpa_data %>%
  rename(
    gpa_original = GRUNNSKOLEPOENG, 
    end_date = AVGDATO, 
    register_date = REGDATO
  )
  
head(gpa_data)
str(gpa_data)
summary(gpa_data)

#### 2. Selection/Check: Missing Data, Invalid Data (Outliers), and Duplicate Data ####

# Find missing data, impossible data and infinite data
gpa_data %>% summarise_all(list(~sum(is.na(.)))) 
gpa_data %>% summarise_all(list(~sum(is.nan(.)))) 
gpa_data %>% summarise_all(list(~sum(is.infinite(.)))) 

# Keep the observations with gpa and completion date for primary school
gpa_data <- gpa_data %>%
  filter(!is.na(gpa_original) & !is.na(end_date)) 

# Find invalid/strange data: <0, ==0, 0<gpa<10, ==60, >60
gpa_data %>%
  summarise(
    gpa_lessthan0 = sum(gpa_original < 0, na.rm = TRUE),
    gpa_equalto0 = sum(gpa_original == 0, na.rm = TRUE),
    gpa_between0and10 = sum(gpa_original > 0 & gpa_original < 10, na.rm = TRUE),
    gpa_equalto60 = sum(gpa_original == 60, na.rm = TRUE),
    gpa_greaterthan60 = sum(gpa_original > 60, na.rm = TRUE)
  )  

# 0<gpa<10
gpa_between0and10 <- gpa_data %>%
  filter(gpa_original > 0 & gpa_original < 10)
gpa_between0and10 <- gpa_between0and10 %>%
  arrange(desc(end_date))

# ==60
gpa_equalto60 <- gpa_data %>%
  filter(gpa_original == 60)
gpa_equalto60 <- gpa_equalto60 %>%
  arrange(desc(end_date))

# >60
gpa_greaterthan60 <- gpa_data %>%
  filter(gpa_original > 60)
gpa_greaterthan60 <- gpa_greaterthan60 %>%
  arrange(desc(end_date))

# Only keep valid gpa that gpa >=10 & gpa <=60
gpa_clean <- gpa_data %>%
  filter(!is.na(gpa_original) & gpa_original >= 10 & gpa_original <= 60) 

head(gpa_clean)
str(gpa_clean)
summary(gpa_clean)

rm(gpa_between0and10)
rm(gpa_equalto60)
rm(gpa_greaterthan60)
rm(gpa_data)

# Check whether there are duplicated individuals
any_duplicated <- any(duplicated(gpa_clean$w19_0634_lnr))
print(any_duplicated)  

# Find all the duplicated individuals
# Identifies locations with duplicate values, either from the first or last duplicate value
duplicates_first <- duplicated(gpa_clean$w19_0634_lnr)  
duplicates_last <- duplicated(gpa_clean$w19_0634_lnr, fromLast = TRUE)
duplicates <- gpa_clean[duplicates_first | duplicates_last, ]  
table(duplicated(gpa_clean$w19_0634_lnr) | duplicated(gpa_clean$w19_0634_lnr, fromLast = TRUE)) 

duplicates <- duplicates %>%
  arrange(desc(w19_0634_lnr),desc(end_date))
glimpse(duplicates)
head(duplicates)
rm(duplicates)

# After checking the structure and details of duplicates
# Only keep the latest GPA of the duplicated individuals
gpa_latest <- gpa_clean %>%
  group_by(w19_0634_lnr) %>%
  arrange(desc(end_date)) %>% 
  slice(1) %>% 
  ungroup() 
rm(gpa_clean)

gpa_latest$year <- as.factor(substr(gpa_latest$end_date, 1, 4))

gpa_latest <- gpa_latest %>%
  select(w19_0634_lnr, gpa_original, year)

head(gpa_latest)
glimpse(gpa_latest)
summary(gpa_latest)

#### 3. Standardize the GPA ####

# Data visualization
# Box plot
ggplot(gpa_latest, aes(x = year, y = gpa_original)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Original GPA distribution plot(several plots)
ggplot(gpa_latest, aes(x = gpa_original)) +
  geom_density(fill = "blue", alpha = 0.5) +
  theme_minimal() + 
  labs(title = "GPA Distribution by Registration Year", x = "Original GPA", y = "Density") +  
  facet_wrap(~ year, scales = "free", ncol = 4)  

# Original GPA distribution plot(different years' distribution on the same plot)
ggplot(gpa_latest, aes(x = gpa_original, color = year)) +
  geom_density(alpha = 0, linewidth = 0.5) +  
  theme_minimal() + 
  labs(title = "Overlay of GPA Distribution by Registration Year", x = "Original GPA", y = "Density") 

# Descriptive stats
gpa_stats <- gpa_latest %>%
  group_by(year) %>%
  summarise(mean_gpa_original = mean(gpa_original, na.rm = TRUE),
            sd_gpa_original = sd(gpa_original, na.rm = TRUE),
            median_gpa_original = median(gpa_original, na.rm = TRUE))
gpa_stats

# The number of children for each year's tests
exam_counts <- gpa_latest %>%
  group_by(year) %>%
  summarise(count = n())
print(exam_counts)

# Standardize the GPA by year
gpa_normalizeddata <- gpa_latest %>%
  group_by(year) %>%
  mutate(gpa_normalized = as.vector(scale(gpa_original, center = TRUE, scale = TRUE)))

# Check 
results <- gpa_normalizeddata %>%
  group_by(year) %>%
  summarise(mean_gpa = mean(gpa_normalized, na.rm = TRUE),
            sd_gpa = sd(gpa_normalized, na.rm = TRUE))

gpa_normalizeddata <- gpa_normalizeddata %>%
  dplyr::select(w19_0634_lnr, gpa_original, gpa_normalized, year)

glimpse(gpa_normalizeddata)
head(gpa_normalizeddata)

# Standardized GPA distribution plot(different years' distribution on the same plot)
ggplot(gpa_normalizeddata, aes(x = gpa_normalized, color = year)) +
  geom_density(alpha = 0, linewidth = 0.5) +  
  theme_minimal() + 
  labs(title = "Overlay of Normalized GPA Distribution by Registration Year", x = "normalized GPA", y = "Density") 

#### 4. Save data ####
save(gpa_normalizeddata, file = "M:/p805-qiyuanp/TrioGCTA/data/1_01_data_gpa_normalized.rda")

#### 5. Descriptive stats #####
# Basic descriptive stats to get the min, median, mean, max, etc.
summary(gpa_normalizeddata$gpa_normalized)

# Histogram
hist(gpa_normalizeddata$gpa_normalized, 
     main = "Histogram of Normalized GPA", xlab = "Normalized GPA", breaks = 30)

# Density map
plot(density(gpa_normalizeddata$gpa_normalized), 
     main = "Density Plot of Normalized GPA", xlab = "Normalized GPA")

