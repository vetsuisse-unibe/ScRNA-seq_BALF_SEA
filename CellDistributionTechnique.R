# Comparison of the distribution of the 5 cytologically distinguishable leukocytes in BALF 
# depending on the sample (BALF or cell suspension) and technique (cytology vs scRNA-seq) used

# Set up environment
library(knitr)
library(kableExtra)
library(webshot)
library(dplyr)
library(car)
library(lme4)

dcc_sc2022<-read.table("DCC_sc2022.txt",header=T,sep="\t")

# Create a new dataframe with the summarized data
summary_table <- dcc_sc2022 %>%
  group_by(horse, type) %>%
  summarise(
    lymphocytes_count = sum(lymphocytes),
    macrophages_count = sum(macrophages),
    neutrophils_count = sum(neutrophils),
    mastocytes_count = sum(mastocytes),
    eosinophils_count = sum(eosinophils),
    `lymphocytes/macrophages` = lymphocytes_count / macrophages_count
      )
# Keep the second repeat  of horse id and replace the first and third element
indices <- seq(2, nrow(summary_table), 3)
summary_table$horse[c(indices-1, indices+1)] <- " "

# Transpose the table
summary_table_transposed <- t(summary_table)

# Create the table using kable
table <- kable(summary_table_transposed, "html") %>%
  kable_styling(full_width = FALSE)

# Remove the row name
row.names(table) <- NULL

# Add bold formatting to the "type" row
table <- table %>%
  row_spec(2, bold = TRUE)

# Save the table as an HTML file
save_kable(table, "table.html")



# NORMALITY AND VARIANCE TESTING

# Select the numeric columns for the normality and variance tests
numeric_columns <- c("lymphocytes", "macrophages", "neutrophils", "mastocytes")
data_subset <- dcc_sc2022[numeric_columns]

# Perform the Shapiro-Wilk test for normality
normality_results <- lapply(data_subset, shapiro.test)

# Create a table of normality test results
normality_table <- data.frame(
  Variable = names(normality_results),
  p_value = sapply(normality_results, function(x) x$p.value),
  stringsAsFactors = FALSE
)

# Perform the Levene's test for equality of variances
variance_results <- lapply(data_subset, function(x) leveneTest(x, dcc_sc2022$type))

# Create a table of variance test results
variance_table <- data.frame(
  Variable = names(variance_results),
  p_value = sapply(variance_results, function(x) x$Pr[1]),
  stringsAsFactors = FALSE
)

# Print the normality and variance test results
print("Normality Test Results:")
print(normality_table)
cat("\n")
print("Variance Test Results:")
print(variance_table)


# STATISTICAL TEST: REPEATED MEASURES ANOVA

# Convert variables to numeric
dcc_sc2022$lymphocytes <- as.numeric(as.character(dcc_sc2022$lymphocytes))
dcc_sc2022$macrophages <- as.numeric(as.character(dcc_sc2022$macrophages))
dcc_sc2022$neutrophils <- as.numeric(as.character(dcc_sc2022$neutrophils))
dcc_sc2022$mastocytes <- as.numeric(as.character(dcc_sc2022$mastocytes))

# Conduct repeated measures ANOVA
anova_result <- Anova(lm(cbind(lymphocytes, macrophages, neutrophils, mastocytes) ~ type + (horse/type), data = dcc_sc2022), type = "III")
summary(anova_result)
print(anova_result)

# Type: The p-value for the type variable is 0.2029 > 0.05 => there is no significant differences between the different types when all other variables are considered.
# Horse: The p-value for the horse variable is 0.1151 > 0.05 => there is no significant differences between horses when all other variables are taken into account.
# Interaction (Type x Horse): The p-value for the interaction between type and horse is 0.6009 > 0.05 => there is no significant interaction effect between type and horse on the dependent variables.

