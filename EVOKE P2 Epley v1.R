# ----------------------- non-Anonymous R-coder --------------------------#
# Ekstra skriftyper til plots
# ---- Libraries  ----
library(tidyverse) 
library(knitr) 
library(tableone)
library(readr) 
library(ggpubr)
library(ggsignif)
library(patchwork)
library(effsize)
library(mirt)

# ---- Load data ----
df <- read_delim("L:/LovbeskyttetMapper01/Motion Sickness - ENS/EVOKE 1 (EpleyOmniax)/REDCap/Interrim data 1.csv", 
                       delim = ";", escape_double = FALSE, col_types = cols(redcap_survey_identifier = col_skip(), 
                                                                            basisoplysninger_timestamp = col_skip()), 
                       trim_ws = TRUE)
df1 <- df
## ---- remove outliers ----

# Clean the CPR column but keep it as a character to preserve leading zeros
clean_cpr <- function(cpr) {
  gsub("[^0-9]", "", cpr)
}
df1$cpr <- sapply(df1$cpr, clean_cpr)

# Assign gender based on the last digit of the cleaned CPR
assign_gender <- function(cpr) {
  last_digit <- as.numeric(substr(cpr, nchar(cpr), nchar(cpr)))
  if (last_digit %% 2 == 0) {
    return("Female")
  } else {
    return("Male")
  }
}
df1$gender <- sapply(df1$cpr, assign_gender)

# Remove the last 4 numbers and convert to date
df1 <- df1 %>%
  mutate(
    dob = substr(cpr, 1, nchar(cpr) - 4),
    dob = dmy(dob))

# Calculate age from date of birth
df1 <- df1 %>%
  mutate(age = as.integer(floor(interval(dob, today()) / years(1))))

# Relocate the Gender column next to the CPR column using dplyr
df1 <- df1 %>%
  select(record_id, cpr, navn, gender, everything())

# Rename 'valsalva' to first or second
df1 <- df1 |> 
  mutate(valsalva = ifelse(valsalva == 1, "First", "Second"))

# Mapping tobak
mapping_tob <- c("0" = "Never smoked", "1" = "Previous smoker", "2" = "Currently smoker")
df1 <- df1 %>%
  mutate(tobak = mapping_tob[as.character(tobak)])

# Mapping alkohol
mapping_alk <- c("1" = "Never consumed alcohol", "2" = "Former Drinkers", "3" = "Current Drinkers")
df1 <- df1 %>%
  mutate(alkohol = mapping_alk[as.character(alkohol)])

#rename pr30_msaq_sum
df1 <- df1 %>%
  rename(pr30_msaq_sum_r2 = pr30_msaq_sum_r12)

# MSSQ-dem der mangler
df1 <- df1 %>%
  mutate(across(ends_with("_barn") | ends_with("_voksen"), ~{
    suppressWarnings(as.numeric(.))
  }))

df1 <- df1 %>%
  rowwise() %>%
  mutate(
    MSA_temp = sum(across(ends_with("_barn"), ~as.numeric(.)), na.rm = TRUE),
    MSB_temp = sum(across(ends_with("_voksen"), ~as.numeric(.)), na.rm = TRUE),
    mssq_raw = if_else(is.na(mssq_raw), MSA_temp + MSB_temp, mssq_raw)
  ) %>%
  ungroup() %>%
  select(-MSA_temp, -MSB_temp)

#rotationer der mangler
df1 <- df1 %>%
  mutate(rotationer = ifelse(is.na(rotationer), 80, rotationer))

# ---- Reshape data - MSAQ ----
## ---- MSAQ - Reshape the data to long format ----

df_long_main <- df1 %>%
  pivot_longer(cols = c(pre_msaq_sum_r1, pr0_msaq_sum_r1, pr10_msaq_sum_r1, pr30_msaq_sum_r1, pr60_msaq_sum_r1,
                        pre_msaq_sum_r2, pr0_msaq_sum_r2, pr10_msaq_sum_r2, pr30_msaq_sum_r2, pr60_msaq_sum_r2,), 
               names_to = "Time", values_to = "MSAQ_score")

# Create a new column for the overall category based on the presence of _r2
df_long_main <- df_long_main %>%
  mutate(Category = ifelse(grepl("_r2", Time), "Second Rotation", "First Rotation"))

# Rename the levels of the Time factor without including (r1)
df_long_main$Time <- factor(df_long_main$Time, levels = c("pre_msaq_sum_r1", "pr0_msaq_sum_r1", "pr10_msaq_sum_r1", "pr30_msaq_sum_r1", "pr60_msaq_sum_r1",
                                                "pre_msaq_sum_r2", "pr0_msaq_sum_r2", "pr10_msaq_sum_r2", "pr30_msaq_sum_r2", "pr60_msaq_sum_r2"),
                       labels = c("Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation", "60 min post rotation",
                                  "Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation", "60 min post rotation"))

# Add the valsalva column and record_id
df_long_main <- df_long_main %>%
  mutate(record_id = rep(df1$record_id, each = 10),
         valsalva_dr = rep(df1$valsalva, each = 10))

# Set the valsalva column based on the specified conditions
df_long_main <- df_long_main %>%
  mutate(valsalva_dr = case_when(
    Category == "First Rotation" & valsalva == "First" ~ "Yes",
    Category == "Second Rotation" & valsalva == "First" ~ "No",
    Category == "First Rotation" & valsalva == "Second" ~ "No",
    Category == "Second Rotation" & valsalva == "Second" ~ "Yes"
  ))

# Select and arrange the final columns
df_final_main <- df_long_main %>%
  select(record_id, MSAQ_score, Time, gender, valsalva, Category, valsalva_dr, mssq_raw, rotationer, hads_sum_angst, hads_sum_depresssion)

#delta MSAQ
df_final_main <- df_final_main %>%
  group_by(record_id, Category) %>%
  mutate(baseline_value = MSAQ_score[Time == 'Baseline'],
         delta_from_baseline = MSAQ_score - baseline_value) %>%
  ungroup() %>%
  select(-baseline_value) 

#Only men
df_final_main_men <- df_final_main %>%
  filter(gender == "Male")

#Only women
df_final_main_women <- df_final_main %>%
  filter(gender == "Female")

## ---- Long dataframe with group names ----
# Select and arrange the final columns
df_long_group <- df_long_main %>%
  select(record_id, MSAQ_score, Time, gender, valsalva, Category, valsalva, mssq_raw)
df_long_group <- df_long_group %>%
  mutate(Group = case_when(
    valsalva == "Yes" & Category == "First Rotation" ~ "A",
    valsalva == "No" & Category == "Second Rotation" ~ "B",
    valsalva == "No" & Category == "First Rotation" ~ "C",
    valsalva == "Yes" & Category == "Second Rotation" ~ "D",
    TRUE ~ NA_character_
  ))

# ---- Tests ----
# ---- Histograms ----
#Histogram for MSAQ
ggplot(df_final_main, aes(x = MSAQ_score)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of MSAQ Scores",
       x = "MSAQ Score",
       y = "Frequency") +
  theme_minimal()

#Histogram for delta MSAQ
ggplot(df_final_main, aes(x = delta_from_baseline)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of delta MSAQ Scores",
       x = "Delta MSAQ Score",
       y = "Frequency") +
  theme_minimal()

#Histogram for MSSQ
ggplot(df_final_main, aes(x = MSSQ_raw)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of MSSQ Scores",
       x = "MSSQ Score",
       y = "Frequency") +
  theme_minimal()

#Histogram for MSSQ (wide)
ggplot(df_wide, aes(x = MSSQ_raw)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of MSSQ Scores",
       x = "MSSQ Score",
       y = "Frequency") +
  theme_minimal()

# ---- Shapiro-Wilk test for normality ----
shapiro.test(df_final_main$MSAQ_score)
shapiro.test(df_final_main$delta_from_baseline)
shapiro.test(df_final_main$MSSQ_raw)
shapiro.test(df_wide$MSSQ_raw)

#If normal dist
t_test_main_msaq <- t.test(MSAQ_score ~ valsalva, data = df_final_main)
print(t_test_main_msaq)

# If non normal dist
wilcox_test_main_msaq <- wilcox.test(MSAQ_score ~ Category, data = df_final_main, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)
wilcox_test_main_msaq <- wilcox.test(MSAQ_score ~ valsalva, data = df_final_main, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)
wilcox_test_main_msaq <- wilcox.test(MSAQ_score ~ valsalva, data = df_final_main, paired = FALSE, exact = FALSE)
print(wilcox_test_main_msaq)
wilcox_test_main_msaq <- wilcox.test(MSAQ_score ~ valsalva_dr, data = df_final_main, paired = FALSE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ Category, data = df_final_main, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)
wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_final_main,paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)
wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_final_main, paired = FALSE, exact = FALSE)
print(wilcox_test_main_msaq)

### Wilcoxon MSAQ at 0 min-----------------
df_filtered_0 <- subset(df_final_main, Time == "0 min post rotation")

wilcox_test_main_msaq <- wilcox.test(MSAQ_score ~ Category, data = df_filtered_0, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(MSAQ_score ~ valsalva, data = df_filtered_0, paired = FALSE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(MSAQ_score ~ valsalva_dr, data = df_filtered_0, paired = FALSE, exact = FALSE)
print(wilcox_test_main_msaq)

### Wilcoxon delta at baseline ---------------
df_filtered_base <- subset(df_final_main, Time == "Baseline")

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ Category, data = df_filtered_base, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_filtered_base, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_filtered_base, paired = FALSE, exact = FALSE)
print(wilcox_test_main_msaq)

### Wilcoxon delta at 0 min-----------------
df_filtered_0 <- subset(df_final_main, Time == "0 min post rotation")

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ Category, data = df_filtered_0, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_filtered_0, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_filtered_0, paired = FALSE, exact = FALSE)
print(wilcox_test_main_msaq)

### Wilcoxon delta at 10 min-----------------
df_filtered_10 <- subset(df_final_main, Time == "10 min post rotation")

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ Category, data = df_filtered_10, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_filtered_10, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_filtered_10, paired = FALSE, exact = FALSE)
print(wilcox_test_main_msaq)

### Wilcoxon delta at 30 min-----------------
df_filtered_30 <- subset(df_final_main, Time == "30 min post rotation")

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ Category, data = df_filtered_30, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_filtered_30, paired = TRUE, exact = FALSE)
print(wilcox_test_main_msaq)

wilcox_test_main_msaq <- wilcox.test(delta_from_baseline ~ valsalva, data = df_filtered_30, paired = FALSE, exact = FALSE)
print(wilcox_test_main_msaq)

#MSSQ vs valsalva
wilcox_test_main_mssq <- wilcox.test(MSSQ_raw ~ valsalva, data = df_final_main, paired = FALSE, exact = FALSE)
print(wilcox_test_main_mssq)

# Effect size -----------------
### Effect size - MSAQ -----------------
cliff.delta(df_final_main$MSAQ_score[df_final_main$Time == "Baseline" & df_final_main$valsalva == "Yes"], 
            df_final_main$MSAQ_score[df_final_main$Time == "Baseline" & df_final_main$valsalva == "No"])
cliff.delta(df_final_main$MSAQ_score[df_final_main$Time == "0 min post rotation" & df_final_main$valsalva == "Yes"], 
            df_final_main$MSAQ_score[df_final_main$Time == "0 min post rotation" & df_final_main$valsalva == "No"])
cliff.delta(df_final_main$MSAQ_score[df_final_main$Time == "10 min post rotation" & df_final_main$valsalva == "Yes"], 
            df_final_main$MSAQ_score[df_final_main$Time == "10 min post rotation" & df_final_main$valsalva == "No"])
cliff.delta(df_final_main$MSAQ_score[df_final_main$Time == "30 min post rotation" & df_final_main$valsalva == "Yes"], 
            df_final_main$MSAQ_score[df_final_main$Time == "30 min post rotation" & df_final_main$valsalva == "No"])
ggplot(df_final_main, aes(x = MSAQ_score)) + 
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") + 
  labs(title = "Histogram of MSAQ Scores", x = "MSAQ Score", y = "Frequency") +
  theme_minimal()

### Effect size - delta MSAQ -----------------
cliff.delta(df_final_main$delta_from_baseline[df_final_main$Time == "Baseline" & df_final_main$valsalva == "Yes"], 
            df_final_main$delta_from_baseline[df_final_main$Time == "Baseline" & df_final_main$valsalva == "No"])
cliff.delta(df_final_main$delta_from_baseline[df_final_main$Time == "0 min post rotation" & df_final_main$valsalva == "Yes"], 
            df_final_main$delta_from_baseline[df_final_main$Time == "0 min post rotation" & df_final_main$valsalva == "No"])
cliff.delta(df_final_main$delta_from_baseline[df_final_main$Time == "10 min post rotation" & df_final_main$valsalva == "Yes"], 
            df_final_main$delta_from_baseline[df_final_main$Time == "10 min post rotation" & df_final_main$valsalva == "No"])
cliff.delta(df_final_main$delta_from_baseline[df_final_main$Time == "30 min post rotation" & df_final_main$valsalva == "Yes"], 
            df_final_main$delta_from_baseline[df_final_main$Time == "30 min post rotation" & df_final_main$valsalva == "No"])
ggplot(df_final_main, aes(x = delta_from_baseline)) + 
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") + 
  labs(title = "Histogram of MSAQ Scores", x = "MSAQ Score", y = "Frequency") +
  theme_minimal()

## Models -----------------

# Make sure categorical variables are factors
df_final_main <- df_final_main %>%
  mutate(
    valsalva = as.factor(valsalva),
    valsalva_dr = as.factor(valsalva_dr),
    gender = as.factor(gender)
  )

# Fit the model - with psyh - no effect
model <- lm(MSAQ_score ~ rotationer + valsalva + mssq_raw + gender + valsalva_dr + hads_sum_angst + hads_sum_depresssion, data = df_final_main)
summary(model)

# Simplify model
model <- lm(MSAQ_score ~ rotationer + valsalva + mssq_raw + gender, data = df_final_main)
summary(model)

### Model MSAQ -----------------
model_msaq <- lm(MSAQ_score ~ MSSQ_raw + valsalva, data = df_final_main)
summary(model_msaq)
plot(model_msaq)

### Model MSAQ delta -----------------
model_delta <- lm(delta_from_baseline ~ MSSQ_raw + valsalva, data = df_filtered_0)
summary(model_delta)
plot(model_delta)

### Model MSAQ delta -----------------
log_delta <- log(df_final_main$delta_from_baseline)
model_delta <- lm(delta_from_baseline ~ MSSQ_raw + valsalva, data = df_filtered_0)
summary(model_delta)
plot(model_delta)


### MSAQ v. MSSQ plot -----------------
ggplot(df_filtered_0, aes(x = MSAQ_score, y = mssq_raw)) +
  geom_point() + 
  geom_smooth(method = "lm", col = "blue") + 
  labs(title = "Change in MSAQ at 0 min post rotation vs Motion Sickness Susceptibilty", x = "MSAQ", y = "MSSQ") + 
  theme_minimal() + 
  stat_cor(method = "pearson", label.x = 30, label.y = 30)

### delta MSAQ v. MSSQ plot -----------------
ggplot(df_final_main, aes(x = delta_from_baseline, y = MSSQ_raw)) +
  geom_point() + 
  geom_smooth(method = "lm", col = "blue") + 
  labs(title = "MSAQ vs MSSQ", x = "Change in MSAQ", y = "MSSQ") + 
  theme_minimal() + 
  stat_cor(method = "pearson", label.x = 30, label.y = 30)

### Linear regression MSAQ at 0 min post rotation -----------------
model1 <- lm(MSAQ_score ~ valsalva + mssq_raw + Category, data = df_final_main)

# Summarize the model
summary(model1)
plot(model1)

### Linear regression MSAQ at full MSAQ -----------------
# Fit the model with MSSQ as a covariate
model2 <- lm(MSAQ_score ~ valsalva + MSSQ_raw, data = df_final_main)

# Summarize the model
summary(model2)
plot(model2)

### Linear regression MSAQ at 0 min post rotation -----------------
# Fit the model with MSSQ as a covariate
model3 <- lm(delta_from_baseline ~ valsalva + MSSQ_raw, data = df_filtered_0)

# Summarize the model
summary(model3)
plot(model3)

### Linear regression MSAQ at full MSAQ -----------------
# Fit the model with MSSQ as a covariate
model4 <- lm(delta_from_baseline ~ valsalva + MSSQ_raw, data = df_final_main)

# Summarize the model
summary(model4)
plot(model4)

# ---- Plots ----
# ---- Lineplots MSAQ ----
## ---- Line plot Valsalva during first rotation ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Timing of Valsalva",
       fill = "Timing of Valsalva") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

## ---- Line plot Valsalva during first rotation ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva_dr, group = valsalva_dr)) +
  geom_smooth(aes(fill = valsalva_dr), se = TRUE, alpha = 0.2) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during rotation",
       fill = "Valsalva during rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

## ---- Line plot Valsalva during first rotation ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva_dr, group = valsalva_dr)) +
  geom_smooth(aes(fill = valsalva_dr), se = TRUE, alpha = 0.2) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during rotation",
       fill = "Valsalva during rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

## ---- Line plot Valsalva during first rotation - men ----
ggplot(df_final_main_men, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

## ---- Line plot for poster ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(text = element_text(family = "Neo Sans"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(hjust = 0.5, size = 18),
        legend.text = element_text(size = 16),
        legend.position = c(0.5, 0.88),  # Adjust the position inside the plot
        legend.background = element_rect(fill = "white", color = "black"))

# Line plot Valsalva during first rotation
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Line plot Valsalva during first rotation only men
ggplot(df_final_main_men, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  labs(title = "MSAQ Scores by Time and Valsalva - Men",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Line plot delta from baseline Valsalva during first rotation
ggplot(df_final_main, aes(x = Time, y = delta_from_baseline, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  labs(title = "Changes in MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "Change in MSAQ Score from baseline",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Line plot Valsalva during first rotation
ggplot(df_final_main, aes(x = Time, y = delta_from_baseline, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "Changes in MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "Change in MSAQ Score from baseline",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))
colors()

## ---- Line plot (new label) MSAQ by timing and rotation ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  scale_color_manual(name = "Valsalva during first rotation", 
                     values = c("No" = "salmon", "Yes" = "turquoise3"),
                     labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  scale_fill_manual(name = "Valsalva during first rotation", 
                    values = c("No" = "salmon", "Yes" = "turquoise3"),
                    labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

df_final_main$valsalva <- factor(df_final_main$valsalva, levels = c("Yes", "No"))

## ---- Line plot (new label) MSAQ by timing ----
# Line plot Valsalva during first rotation
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  scale_color_manual(name = "Valsalva during first rotation", 
                     values = c("No" = "salmon", "Yes" = "turquoise3"),
                     labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  scale_fill_manual(name = "Valsalva during first rotation", 
                    values = c("No" = "salmon", "Yes" = "turquoise3"),
                    labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# By gender
## ---- Line plot MSAQ by gender and rotation ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = gender, group = gender)) +
  geom_smooth(aes(fill = gender), se = TRUE, alpha = 0.2) +
  facet_wrap(~ gender + Category) +
  labs(title = "MSAQ Scores by Time and Gender (Closed First)",
       x = "Time",
       y = "MSAQ Score",
       color = "Gender",
       fill = "Gender") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# ---- Line plot delta MSAQ (new label) dMSAQ by timing and rotation ----
## ---- Line plot (new label) delta MSAQ by timing and rotation ----
ggplot(df_final_main, aes(x = Time, y = delta_from_baseline, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "Change in MSAQ Score from baseline",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  scale_color_manual(name = "Valsalva during first rotation", 
                     values = c("No" = "salmon", "Yes" = "turquoise3"),
                     labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  scale_fill_manual(name = "Valsalva during first rotation", 
                    values = c("No" = "salmon", "Yes" = "turquoise3"),
                    labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

## ---- Line plot (new label) delta MSAQ by timing  ----
ggplot(df_final_main, aes(x = Time, y = delta_from_baseline, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  labs(title = "Normalized MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "Change in MSAQ Score from baseline",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  scale_color_manual(name = "Valsalva during first rotation", 
                     values = c("No" = "salmon", "Yes" = "turquoise3"),
                     labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  scale_fill_manual(name = "Valsalva during first rotation", 
                    values = c("No" = "salmon", "Yes" = "turquoise3"),
                    labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# ---- Boxplots ----
## ---- Boxplot Valsalva during first rotation ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva, group = interaction(Time, valsalva))) +
  geom_boxplot(aes(fill = valsalva), alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

## ---- Boxplot MSSQ by gender ----
ggplot(df1, aes(x = mssq_raw, y = gender, fill = gender)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
labs(title = "Motion sickness susceptibilty by gender", x = "MSSQ Raw Score", y = "Gender") +  
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  
  )

## ---- Boxplot MSSQ by group ----
ggplot(df1, aes(x = mssq_raw, y = valsalva, fill = valsalva)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Motion sickness susceptibilty by gender", x = "MSSQ Raw Score", y = "Gender") +  
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  
  )

## ---- Boxplot MSAQ A vs. B ----
ggplot(df_A_B, aes(x = Time, y = MSAQ_score, color = Group, group = interaction(Time, Group))) +
  geom_boxplot(aes(fill = Group), alpha = 0.2) +
  facet_wrap(~ Group) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Group",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

## ---- Boxplot MSAQ C vs. D ----
ggplot(df_C_D, aes(x = Time, y = MSAQ_score, color = Group, group = interaction(Time, Group))) +
  geom_boxplot(aes(fill = Group), alpha = 0.2) +
  facet_wrap(~ Group) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Group",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# ---- Violin plots ----
## ---- Violin MSAQ by timing ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, fill = valsalva)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Creates a violin plot
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +  # Adds individual points
  stat_compare_means(aes(group = valsalva), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test for non-parametric data
  theme_minimal() +
  labs(title = "Violin Plot of MSAQ Over Time",
       x = "Time Point",
       y = "MSAQ Score",
       fill = "Timing of Valsalva") +
  theme(legend.position = "top")


## ---- Violin delta MSAQ by timing ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, fill = valsalva)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +  # Adds individual points
  stat_compare_means(aes(group = valsalva), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test for non-parametric data
  theme_minimal() +
  labs(title = "Violin Plot of change in MSAQ Over Time",
       x = "Time Point",
       y = "Change in MSAQ Score from baseline",
       fill = "Valsalva during first rotation") +
  theme(legend.position = "top")

## ---- Violin delta MSAQ by timing ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, fill = Category)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +  # Adds individual points
  stat_compare_means(aes(group = Category), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test for non-parametric data
  theme_minimal() +
  labs(title = "Violin Plot of change in MSAQ Over Time",
       x = "Time Point",
       y = "Change in MSAQ Score from baseline",
       fill = "Rotation") +
  theme(legend.position = "top")

## ---- Violin delta MSAQ by timing ----
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, fill = valsalva_dr)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +  # Adds individual points
  stat_compare_means(aes(group = valsalva_dr), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test for non-parametric data
  theme_minimal() +
  labs(title = "Violin Plot of change in MSAQ Over Time",
       x = "Time Point",
       y = "Change in MSAQ Score from baseline",
       fill = "Valslva during rotation") +
  theme(legend.position = "top")

## ---- Violin delta MSAQ without baseline  ----
df_filtered <- df_final_main %>% filter(Time != 'Baseline')

## ---- Violin plot for Valsalva timing - delta MSAQ  ----
ggplot(df_filtered, aes(x = Time, y = MSAQ_score, fill = valsalva)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Creates a violin plot
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +  # Adds individual points
  stat_compare_means(aes(group = valsalva), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test for non-parametric data with significance lines
  geom_signif(comparisons = list(c("Yes", "No")), map_signif_level = TRUE, step_increase = 0.1) +  # Adds significance lines
  labs(title = "Change in MSAQ Score from baseline - by Time and Timing og Valsalva",
       x = "Time",
       y = "Change in MSAQ Score from baseline",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(legend.position = "top")

## ---- Violin plot for Valsalva - delta MSAQ  ----
ggplot(df_filtered, aes(x = Time, y = delta_from_baseline, fill = valsalva)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Creates a violin plot
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +  # Adds individual points
  stat_compare_means(aes(group = valsalva), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test for non-parametric data with significance lines
  geom_signif(comparisons = list(c("Yes", "No")), map_signif_level = TRUE, step_increase = 0.1) +  # Adds significance lines
  labs(title = "Change in MSAQ Score from baseline - by time and Valsalva",
       x = "Time",
       y = "Change in MSAQ Score from baseline",
       fill = "Performed Valsalva maneuver during rotation") +
  theme_minimal() +
  theme(legend.position = "top")

## ---- Violin plot for rotation - delta MSAQ  ----
ggplot(df_filtered, aes(x = Time, y = delta_from_baseline, fill = Category)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Creates a violin plot
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +  # Adds individual points
  stat_compare_means(aes(group = Category), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test for non-parametric data with significance lines
  geom_signif(comparisons = list(c("Yes", "No")), map_signif_level = TRUE, step_increase = 0.1) +  # Adds significance lines
  labs(title = "Change in MSAQ Score from baseline by rotation",
       x = "Time",
       y = "Change in MSAQ Score from baseline",
       fill = "Rotation") +
  theme_minimal() +
  theme(legend.position = "top") 

## ---- Combined violin plot  ----
v1 <- ggplot(df_final_main2, aes(x = Time, y = MSAQ_score, fill = valsalva)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Creates a violin plot
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +  # Adds individual points
  stat_compare_means(aes(group = valsalva), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test for non-parametric data with significance lines
  geom_signif(comparisons = list(c("Yes", "No")), map_signif_level = TRUE, step_increase = 0.1) +  # Adds significance lines
  labs(title = "MSAQ Score by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       fill = "Performed Valsalva maneuver during rotation") +
  theme_minimal() +
  theme(legend.position = "top")

# Violin plot for rotation first or second
v2 <-ggplot(df_final_main2, aes(x = Time, y = MSAQ_score, fill = Category)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Creates a violin plot
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +  # Adds individual points
  stat_compare_means(aes(group = Category), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test for non-parametric data with significance lines
  geom_signif(comparisons = list(c("Yes", "No")), map_signif_level = TRUE, step_increase = 0.1) +  # Adds significance lines
  labs(title = "MSAQ Score by Time and Rotation",
       x = "Time",
       y = "MSAQ Score",
       fill = "Rotation") +
  theme_minimal() +
  theme(legend.position = "top")

v3 <- v1 / v2
plot(v3)

# ---- Tables Main MSAQ ----
summary(df1$MSSQ_raw)
summary(df1$age)
summary(df1$hads_sum_angst)
summary(df1$hads_sum_depresssion)
summary_df <- df_tb1 %>%
  group_by(valsalva) %>%
  summarize(
    count = n(),
    mean = mean(MSSQ_raw, na.rm = TRUE),
    sd = sd(MSSQ_raw, na.rm = TRUE),
    median = median(MSSQ_raw, na.rm = TRUE),
    min = min(MSSQ_raw, na.rm = TRUE),
    max = max(MSSQ_raw, na.rm = TRUE)
  )

print(summary_df)

# ---- Table 1 ----
# Create table 1 dataframe
df_tb1 <- df1 %>%
  select(record_id, alkohol, tobak, valsalva, hads_sum_angst, hads_sum_depresssion, MSA, MSB, MSSQ_raw, age, gender) %>%
  mutate(Valsalva_first = ifelse(valsalva == 1, "Yes", "No"))

# Create Table 1
tbl1 <- CreateTableOne(vars = c("age", "gender", "MSSQ_raw", "MSA", "MSB", "hads_sum_depresssion", "hads_sum_angst"), 
                       strata = "valsalva", 
                       data = df_tb1, 
                       factorVars = c("gender"))
# Print Table 1
print(tbl1, showAllLevels = TRUE)

# Print Table 1 with test information
print(tbl1, showAllLevels = TRUE, test = TRUE)

# Convert Table 1 to a data frame
tbl1_df <- as.data.frame(print(tbl1, showAllLevels = TRUE, quote = FALSE))

# Get detailed summary with test information
summary(tbl1)

# First v. second rotation median, min, max, sd table #
# Calculate summary statistics
summary_stats <- data.frame(
  Time_Point = c("Baseline", "0 min Postrotation", "10 min Postrotation", "30 min Postrotation"),
  Median = c(median(df1$pre_msaq_sum_r1), median(df1$pr0_msaq_sum_r1), median(df1$pr10_msaq_sum_r1), median(df1$pr30_msaq_sum_r1)),
  Min = c(min(df1$pre_msaq_sum_r1), min(df1$pr0_msaq_sum_r1), min(df1$pr10_msaq_sum_r1), min(df1$pr30_msaq_sum_r1)),
  Max = c(max(df1$pre_msaq_sum_r1), max(df1$pr0_msaq_sum_r1), max(df1$pr10_msaq_sum_r1), max(df1$pr30_msaq_sum_r1)),
  SD = c(sd(df1$pre_msaq_sum_r1), sd(df1$pr0_msaq_sum_r1), sd(df1$pr10_msaq_sum_r1), sd(df1$pr30_msaq_sum_r1))
)
print(summary_stats)

# ---- Table 2 - 1st and 2nd rotation ----
# Combined first and second rotation
# Calculate summary statistics for the first rotation
summary_stats_first_rotation <- data.frame(
  Time_Point = c("Baseline", "0 min Postrotation", "10 min Postrotation", "30 min Postrotation"),
  Median_First = c(median(df1$pre_msaq_sum_r1), median(df1$pr0_msaq_sum_r1), median(df1$pr10_msaq_sum_r1), median(df1$pr30_msaq_sum_r1)),
  Min_First = c(min(df1$pre_msaq_sum_r1), min(df1$pr0_msaq_sum_r1), min(df1$pr10_msaq_sum_r1), min(df1$pr30_msaq_sum_r1)),
  Max_First = c(max(df1$pre_msaq_sum_r1), max(df1$pr0_msaq_sum_r1), max(df1$pr10_msaq_sum_r1), max(df1$pr30_msaq_sum_r1)),
  SD_First = c(sd(df1$pre_msaq_sum_r1), sd(df1$pr0_msaq_sum_r1), sd(df1$pr10_msaq_sum_r1), sd(df1$pr30_msaq_sum_r1))
)

# Calculate summary statistics for the second rotation
summary_stats_second_rotation <- data.frame(
  Time_Point = c("Baseline", "0 min Postrotation", "10 min Postrotation", "30 min Postrotation"),
  Median_Second = c(median(df1$pre_msaq_sum_r2), median(df1$pr0_msaq_sum_r2), median(df1$pr30_msaq_sum_r2), median(df1$pr30_msaq_sum_r2)),
  Min_Second = c(min(df1$pre_msaq_sum_r2), min(df1$pr0_msaq_sum_r2), min(df1$pr30_msaq_sum_r2), min(df1$pr30_msaq_sum_r2)),
  Max_Second = c(max(df1$pre_msaq_sum_r2), max(df1$pr0_msaq_sum_r2), max(df1$pr30_msaq_sum_r2), max(df1$pr30_msaq_sum_r2)),
  SD_Second = c(sd(df1$pre_msaq_sum_r2), sd(df1$pr0_msaq_sum_r2), sd(df1$pr30_msaq_sum_r2), sd(df1$pr30_msaq_sum_r2))
)

# Combine the data frames
summary_stats_combined <- cbind(summary_stats_first_rotation, summary_stats_second_rotation[, -1])

# Print the combined summary table
print(summary_stats_combined)

# ---- Table 3 - Valsalva v. no-valsalva ----
# Valslva median, min, max, sd table #
#Valsalva table
val_summary_stats <- df_final_main %>%
  group_by(Time, valsalva) %>%
  summarise(
    Median = median(MSAQ_score, na.rm = TRUE),
    Min = min(MSAQ_score, na.rm = TRUE),
    Max = max(MSAQ_score, na.rm = TRUE),
    SD = sd(MSAQ_score, na.rm = TRUE)
  )
# Calculate summary statistics for valsalva = "yes"
summary_stats_yes <- df_final_main %>%
  filter(valsalva == "Yes") %>%
  group_by(Time) %>%
  summarise(
    Median_Yes = median(MSAQ_score, na.rm = TRUE),
    Min_Yes = min(MSAQ_score, na.rm = TRUE),
    Max_Yes = max(MSAQ_score, na.rm = TRUE),
    SD_Yes = sd(MSAQ_score, na.rm = TRUE)
  )

# Calculate summary statistics for valsalva = "no"
summary_stats_no <- df_final_main %>%
  filter(valsalva == "No") %>%
  group_by(Time) %>%
  summarise(
    Median_No = median(MSAQ_score, na.rm = TRUE),
    Min_No = min(MSAQ_score, na.rm = TRUE),
    Max_No = max(MSAQ_score, na.rm = TRUE),
    SD_No = sd(MSAQ_score, na.rm = TRUE)
  )

# Combine the data frames
val_summary_stats_combined <- left_join(summary_stats_yes, summary_stats_no, by = "Time")
print(val_summary_stats_combined)

# ---- Data for article - Valsalva v. no-valsalva ----
summary(df1$MSSQ_raw)

# ---- MSSS ----
## ---- Reshape the data to long format ----
df_long_msss <- df1 %>%
  pivot_longer(cols = c(msss_pre, msss_pr0, msss_pr10, msss_pr30,
                        msss_pre_v2, msss_pr0_v2, msss_pr10_v2, msss_pr30_v2), 
               names_to = "Time", values_to = "msss")

# Create a new column for the overall category based on the presence of _v2
df_long_msss <- df_long_msss %>%
  mutate(Category = ifelse(grepl("_v2", Time), "Second Rotation", "First Rotation"))

# Rename the levels of the Time factor without including (v2)
df_long_msss$Time <- factor(df_long_msss$Time, levels = c("msss_pre", "msss_pr0", "msss_pr10", "msss_pr30",
                                                          "msss_pre_v2", "msss_pr0_v2", "msss_pr10_v2", "msss_pr30_v2"),
                            labels = c("Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation",
                                       "Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation"))

# Add the valsalva column and record_id
df_long_msss <- df_long_msss %>%
  mutate(record_id = rep(df1$record_id, each = 8),
         valsalva = rep(df1$valsalva, each = 8))

# Set the valsalva column based on the specified conditions
df_long_msss <- df_long_msss %>%
  mutate(valsalva = case_when(
    Category == "First Rotation" & valsalva == "Yes" ~ "Yes",
    Category == "Second Rotation" & valsalva == "Yes" ~ "No",
    Category == "First Rotation" & valsalva == "No" ~ "No",
    Category == "Second Rotation" & valsalva == "No" ~ "Yes"
  ))

# Select and arrange the final columns
df_long_msss <- df_long_msss %>%
  select(record_id, msss, Time, valsalva, Category, valsalva)

#Histogram for data distribution msss
ggplot(df_long_msss, aes(x = msss)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of msss Scores",
       x = "msss Score",
       y = "Frequency") +
  theme_minimal()

# Plot msss
ggplot(df_long_msss, aes(x = Time, y = msss, fill = valsalva)) +
  geom_violin(trim = FALSE, alpha = 0.5) + 
  geom_jitter(width = 0.1, alpha = 0.5) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif") + 
  theme_minimal() +
  labs(title = paste("msss Scores"),
       x = "Time",
       y = "msss")

# Line plot Valsalva during first rotation
ggplot(df_long_msss, aes(x = Time, y = msss, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "msss Scores by Time and Valsalva",
       x = "Time",
       y = "msss Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Line plot Valsalva during first rotation
ggplot(df_long_msss, aes(x = Time, y = msss, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  labs(title = "msss Scores by Time and Valsalva",
       x = "Time",
       y = "msss Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# Boxplot msss and Valsalva
ggplot(df_long_msss, aes(x = Time, y = msss, fill = valsalva)) +
  geom_boxplot() +
  facet_wrap(~ Category) +
  labs(title = "msss Scores by Time and Valsalva",
       x = "Time",
       y = "msss Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Bar plot Valsalva during first rotation
ggplot(df_long_msss, aes(x = Time, y = msss, fill = valsalva)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  facet_wrap(~ Category) +
  labs(title = "Mean MSSQ Scores by Time and Valsalva",
       x = "Time",
       y = "Mean msss Score",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# ---- MSAQ G ----

# MSAQ - Reshape the data to long format
df_long_main_g <- df1 %>%
  pivot_longer(cols = c(sum_msaq_g, sum_msaq_g_pr0, sum_msaq_g_pr10, sum_msaq_g_pr30,
                        sum_msaq_g_v2, sum_msaq_g_pr0_v2, sum_msaq_g_pr10_v2, sum_msaq_g_pr30_v2), 
               names_to = "Time", values_to = "MSAQ_score")
df_long_main_g <- df_long_main_g %>%
  mutate(Category = ifelse(grepl("_v2", Time), "Second Rotation", "First Rotation"))
df_long_main_g$Time <- factor(df_long_main_g$Time, levels = c("sum_msaq_g", "sum_msaq_g_pr0", "sum_msaq_g_pr10", "sum_msaq_g_pr30",
                                                          "sum_msaq_g_v2", "sum_msaq_g_pr0_v2", "sum_msaq_g_pr10_v2", "sum_msaq_g_pr30_v2"),
                            labels = c("Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation",
                                       "Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation"))
df_long_main_g <- df_long_main_g %>%
  mutate(record_id = rep(df1$record_id, each = 8),
         valsalva = rep(df1$valsalva, each = 8))
df_long_main_g <- df_long_main_g %>%
  mutate(valsalva = case_when(
    Category == "First Rotation" & valsalva == "Yes" ~ "Yes",
    Category == "Second Rotation" & valsalva == "Yes" ~ "No",
    Category == "First Rotation" & valsalva == "No" ~ "No",
    Category == "Second Rotation" & valsalva == "No" ~ "Yes"
  ))
df_final_main_g <- df_long_main_g %>%
  select(record_id, MSAQ_score, Time, valsalva, Category, valsalva)

#MSAQ GPlots
# Line plot Valsalva during first rotation
ggplot(df_final_main_g, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ G Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ G Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# ---- MSAQ S ----

# MSAQ - Reshape the data to long format
df_long_main_s <- df1 %>%
  pivot_longer(cols = c(sum_msaq_s, sum_msaq_s_pr0, sum_msaq_s_pr10, sum_msaq_s_pr30,
                        sum_msaq_s_v2, sum_msaq_s_pr0_v2, sum_msaq_s_pr10_v2, sum_msaq_s_pr30_v2), 
               names_to = "Time", values_to = "MSAQ_score")

# Create a new column for the overall category based on the presence of _v2
df_long_main_s <- df_long_main_s %>%
  mutate(Category = ifelse(grepl("_v2", Time), "Second Rotation", "First Rotation"))

# Rename the levels of the Time factor without including (v2)
df_long_main_s$Time <- factor(df_long_main_s$Time, levels = c("sum_msaq_s", "sum_msaq_s_pr0", "sum_msaq_s_pr10", "sum_msaq_s_pr30",
                                                              "sum_msaq_s_v2", "sum_msaq_s_pr0_v2", "sum_msaq_s_pr10_v2", "sum_msaq_s_pr30_v2"),
                              labels = c("Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation",
                                         "Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation"))

# Add the valsalva column and record_id
df_long_main_s <- df_long_main_s %>%
  mutate(record_id = rep(df1$record_id, each = 8),
         valsalva = rep(df1$valsalva, each = 8))

# Set the valsalva column based on the specified conditions
df_long_main_s <- df_long_main_s %>%
  mutate(valsalva = case_when(
    Category == "First Rotation" & valsalva == "Yes" ~ "Yes",
    Category == "Second Rotation" & valsalva == "Yes" ~ "No",
    Category == "First Rotation" & valsalva == "No" ~ "No",
    Category == "Second Rotation" & valsalva == "No" ~ "Yes"
  ))

# Select and arrange the final columns
df_final_main_s <- df_long_main_s %>%
  select(record_id, MSAQ_score, Time, valsalva, Category, valsalva)

# MSAQ S Plots
# Line plot Valsalva during first rotation
ggplot(df_final_main_s, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ S Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ S Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# ---- MSAQ C ----

# MSAQ - Reshape the data to long format
df_long_main_c <- df1 %>%
  pivot_longer(cols = c(sum_msaq_c, sum_msaq_c_pr0, sum_msaq_c_pr10, sum_msaq_c_pr30,
                        sum_msaq_c_v2, sum_msaq_c_pr0_v2, sum_msaq_c_pr10_v2, sum_msaq_c_pr30_v2), 
               names_to = "Time", values_to = "MSAQ_score")

# Create a new column for the overall category based on the presence of _v2
df_long_main_c <- df_long_main_c %>%
  mutate(Category = ifelse(grepl("_v2", Time), "Second Rotation", "First Rotation"))

# Rename the levels of the Time factor without including (v2)
df_long_main_c$Time <- factor(df_long_main_c$Time, levels = c("sum_msaq_c", "sum_msaq_c_pr0", "sum_msaq_c_pr10", "sum_msaq_c_pr30",
                                                              "sum_msaq_c_v2", "sum_msaq_c_pr0_v2", "sum_msaq_c_pr10_v2", "sum_msaq_c_pr30_v2"),
                              labels = c("Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation",
                                         "Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation"))

# Add the valsalva column and record_id
df_long_main_c <- df_long_main_c %>%
  mutate(record_id = rep(df1$record_id, each = 8),
         valsalva = rep(df1$valsalva, each = 8))

# Set the valsalva column based on the specified conditions
df_long_main_c <- df_long_main_c %>%
  mutate(valsalva = case_when(
    Category == "First Rotation" & valsalva == "Yes" ~ "Yes",
    Category == "Second Rotation" & valsalva == "Yes" ~ "No",
    Category == "First Rotation" & valsalva == "No" ~ "No",
    Category == "Second Rotation" & valsalva == "No" ~ "Yes"
  ))

# Select and arrange the final columns
df_final_main_c <- df_long_main_c %>%
  select(record_id, MSAQ_score, Time, valsalva, Category, valsalva)

# MSAQ C Plots
# Line plot Valsalva during first rotation
ggplot(df_final_main_c, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ C Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ C Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# ---- MSAQ P ----

# MSAQ - Reshape the data to long format
df_long_main_p <- df1 %>%
  pivot_longer(cols = c(sum_msaq_p, sum_msaq_p_pr0, sum_msaq_p_pr10, sum_msaq_p_pr30,
                        sum_msaq_p_v2, sum_msaq_p_pr0_v2, sum_msaq_p_pr10_v2, sum_msaq_p_pr30_v2), 
               names_to = "Time", values_to = "MSAQ_score")

# Create a new column for the overall category based on the presence of _v2
df_long_main_p <- df_long_main_p %>%
  mutate(Category = ifelse(grepl("_v2", Time), "Second Rotation", "First Rotation"))

# Rename the levels of the Time factor without including (v2)
df_long_main_p$Time <- factor(df_long_main_p$Time, levels = c("sum_msaq_p", "sum_msaq_p_pr0", "sum_msaq_p_pr10", "sum_msaq_p_pr30",
                                                              "sum_msaq_p_v2", "sum_msaq_p_pr0_v2", "sum_msaq_p_pr10_v2", "sum_msaq_p_pr30_v2"),
                              labels = c("Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation",
                                         "Baseline", "0 min post rotation", "10 min post rotation", "30 min post rotation"))

# Add the valsalva column and record_id
df_long_main_p <- df_long_main_p %>%
  mutate(record_id = rep(df1$record_id, each = 8),
         valsalva = rep(df1$valsalva, each = 8))

# Set the valsalva column based on the specified conditions
df_long_main_p <- df_long_main_p %>%
  mutate(valsalva = case_when(
    Category == "First Rotation" & valsalva == "Yes" ~ "Yes",
    Category == "Second Rotation" & valsalva == "Yes" ~ "No",
    Category == "First Rotation" & valsalva == "No" ~ "No",
    Category == "Second Rotation" & valsalva == "No" ~ "Yes"
  ))

# Select and arrange the final columns
df_final_main_p <- df_long_main_p %>%
  select(record_id, MSAQ_score, Time, valsalva, Category, valsalva)

# MSAQ P Plots
# Line plot Valsalva during first rotation
ggplot(df_final_main_p, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ P Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ P Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# ---- Dataframe subgroups ----
main_msaq <- df_final_main[, c("MSAQ_score", "Time", "record_id", "valsalva", "Category", "valsalva"), drop = FALSE]
c_msaq <- df_final_main_c[, c("MSAQ_score", "Time", "record_id", "valsalva", "Category", "valsalva"), drop = FALSE]
g_msaq <- df_final_main_g[, c("MSAQ_score", "Time", "record_id", "valsalva", "Category", "valsalva"), drop = FALSE]
p_msaq <- df_final_main_p[, c("MSAQ_score", "Time", "record_id", "valsalva", "Category", "valsalva"), drop = FALSE]
s_msaq <- df_final_main_s[, c("MSAQ_score", "Time", "record_id", "valsalva", "Category", "valsalva"), drop = FALSE]
main_msaq$group <- "Main"
c_msaq$group <- "Central"
g_msaq$group <- "Gastrointestinal"
p_msaq$group <- "Peripheral"
s_msaq$group <- "Sopite-related"
df_final_group <- rbind(main_msaq, c_msaq, p_msaq, g_msaq, s_msaq)
df_subgroup <- rbind(c_msaq, p_msaq, g_msaq, s_msaq)

# ---- Plot subgroups ----
ggplot(df_subgroup, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ group) +
  labs(title = "MSAQ Subgroup-scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  scale_color_manual(name = "Valsalva during first rotation", 
                     values = c("No" = "salmon", "Yes" = "turquoise3"),
                     labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  scale_fill_manual(name = "Valsalva during first rotation", 
                    values = c("No" = "salmon", "Yes" = "turquoise3"),
                    labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

df_subgroup$valsalva <- factor(df_subgroup$valsalva, levels = c("Yes", "No"))

# ---- Remove outliers ----
df2 <- subset(df1, !(row.names(df1) %in% c(5, 16)))
df_final_main$valsalva <- factor(df_final_main$valsalva, levels = c("Yes", "No"))

# Line plot Valsalva during first rotation
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  facet_wrap(~ Category) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  scale_color_manual(name = "Valsalva during first rotation", 
                     values = c("No" = "salmon", "Yes" = "royalblue3"),
                     labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  scale_fill_manual(name = "Valsalva during first rotation", 
                    values = c("No" = "salmon", "Yes" = "royalblue3"),
                    labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

df_final_main$valsalva <- factor(df_final_main$valsalva, levels = c("Yes", "No"))

# Line plot Valsalva during first rotation
ggplot(df_final_main, aes(x = Time, y = MSAQ_score, color = valsalva, group = valsalva)) +
  geom_smooth(aes(fill = valsalva), se = TRUE, alpha = 0.2) +
  labs(title = "MSAQ Scores by Time and Valsalva",
       x = "Time",
       y = "MSAQ Score",
       color = "Valsalva during first rotation",
       fill = "Valsalva during first rotation") +
  scale_color_manual(name = "Valsalva during first rotation", 
                     values = c("No" = "salmon", "Yes" = "royalblue3"),
                     labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  scale_fill_manual(name = "Valsalva during first rotation", 
                    values = c("No" = "salmon", "Yes" = "royalblue3"),
                    labels = c("No" = "No (C + D)", "Yes" = "Yes (A + B)")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# ---- Correllations IRT  ----
df1$sum_msaq <- rowSums(df1[, paste0("msaq", 1:16)])
item_total_correlation <- function(df) {
  total_score <- df$sum_msaq
  correlations <- numeric(16)
  for (i in 1:16) {
    item_col <- paste0("msaq", i)
    if (sd(df[[item_col]]) == 0) {
      correlations[i] <- NA
    } else {
      correlations[i] <- cor(df[[item_col]], total_score - df[[item_col]])
    }
  }
  return(correlations)
}

baseline_correlations <- item_total_correlation(df1)
print(baseline_correlations)

# Item discrimination index for baseline
discrimination_index <- function(df) {
  total_score <- df$sum_msaq
  discriminations <- numeric(16)
  for (i in 1:16) {
    item_col <- paste0("msaq", i)
    high_group <- df[[item_col]][total_score > median(total_score, na.rm = TRUE)]
    low_group <- df[[item_col]][total_score <= median(total_score, na.rm = TRUE)]
    discriminations[i] <- mean(high_group, na.rm = TRUE) - mean(low_group, na.rm = TRUE)
  }
  return(discriminations)
}

baseline_discriminations <- discrimination_index(df1)
print(baseline_discriminations)

# IRT analysis for baseline
library(mirt)

fit_irt <- function(df) {
  data <- df[, paste0("msaq", 1:16)]
  model <- mirt(data, 1, itemtype = "2PL")
  coef(model, IRTpars = TRUE)
}

baseline_irt_params <- fit_irt(df1)
print(baseline_irt_params)
