# Load required libraries
library(tidyverse)  # For data manipulation and visualization
library(readxl)     # For reading Excel files
library(SpATS)      # For spatial analysis of field trials
library(fields)     # For generating spatial trend plots
library(plyr)       # For merging multiple data frames

# Clear the R environment
rm(list = ls())

# Load phenotype data for 2018 and preprocess
pheno2018 <- read_excel("../../phenotype/BBL_phenotype_2020.xlsx", 
                        sheet = "2018", na = "NA") %>%
  mutate(
    col_f = factor(P),  # Convert column IDs to factor
    row_f = factor(R),  # Convert row IDs to factor
    CheckStatus = ifelse(ID %in% c("Eclipse", "Zenith", "Zorro"), "Check", "NonCheck") %>% as.factor()
  )

# Load phenotype data for 2019 and preprocess
pheno2019 <- read_excel("../../phenotype/BBL_phenotype_2020.xlsx", 
                        sheet = "2019", na = "NA") %>%
  mutate(
    col_f = factor(P),  # Convert column IDs to factor
    row_f = factor(R),  # Convert row IDs to factor
    CheckStatus = ifelse(ID %in% c("Eclipse", "Zenith", "Zorro"), "Check", "NonCheck") %>% as.factor()
  )

# Fit spatial models for various traits in 2018
fit.yd18 <- SpATS(response = "yield",  # Yield
                  spatial = ~ PSANOVA(P, R, nseg = c(20, 10)),
                  genotype = "ID",
                  fixed = ~ CheckStatus,  # Fixed effect for check genotypes
                  genotype.as.random = TRUE,
                  random = ~ row_f + col_f,  # Random effects for rows and columns
                  data = pheno2018,
                  control = list(tolerance = 1e-03))

fit.df18 <- SpATS(response = "df",  # Days to flowering
                  spatial = ~ PSANOVA(P, R, nseg = c(20, 10)),
                  genotype = "ID",
                  fixed = ~ CheckStatus,
                  genotype.as.random = TRUE,
                  random = ~ row_f + col_f,
                  data = pheno2018,
                  control = list(tolerance = 1e-03))

fit.dm18 <- SpATS(response = "dm",  # Days to maturity
                  spatial = ~ PSANOVA(P, R, nseg = c(20, 10)),
                  genotype = "ID",
                  fixed = ~ CheckStatus,
                  genotype.as.random = TRUE,
                  random = ~ row_f + col_f,
                  data = pheno2018,
                  control = list(tolerance = 1e-03))

fit.sw18 <- SpATS(response = "sw",  # Seed weight
                  spatial = ~ PSANOVA(P, R, nseg = c(20, 10)),
                  genotype = "ID",
                  fixed = ~ CheckStatus,
                  genotype.as.random = TRUE,
                  random = ~ row_f + col_f,
                  data = pheno2018,
                  control = list(tolerance = 1e-03))

fit.color18 <- SpATS(response = "color",  # Color
                     spatial = ~ PSANOVA(P, R, nseg = c(20, 10)),
                     genotype = "ID",
                     fixed = ~ CheckStatus,
                     genotype.as.random = TRUE,
                     random = ~ row_f + col_f,
                     data = pheno2018,
                     control = list(tolerance = 1e-03))


fit.app18 = SpATS(response = "app", # Appearance
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2018,
                 control = list(tolerance = 1e-03))


fit.L18 = SpATS(response = "L", # L*
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2018,
                 control = list(tolerance = 1e-03))

fit.a18 = SpATS(response = "a", # a*
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2018,
                 control = list(tolerance = 1e-03))

fit.b18 = SpATS(response = "b", # b*
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2018,
                 control = list(tolerance = 1e-03))

fit.text18 = SpATS(response = "text_kg", # texture
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2018,
                 control = list(tolerance = 1e-03))

fit.yd19 = SpATS(response = "yield", # Yield 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

fit.df19 = SpATS(response = "df", # DF 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

fit.dm19 = SpATS(response = "dm", # DM 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

fit.sw19 = SpATS(response = "sw", # sw 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

fit.color19 = SpATS(response = "color", # color 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

fit.app19 = SpATS(response = "app",  # app 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

fit.L19 = SpATS(response = "L", # L 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

fit.a19 = SpATS(response = "a", # a 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

fit.b19 = SpATS(response = "b", # b 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

fit.text19 = SpATS(response = "text_kg", # texture 2019
                 spatial = ~ PSANOVA (P, R, nseg = c(20,10)),
                 genotype = "ID",
                 fixed = ~ CheckStatus,
                 genotype.as.random = TRUE,
                 random = ~row_f + col_f,
                 data= pheno2019,
                 control = list(tolerance = 1e-03))

# Generate predicted values for each trait (BLUPs)
yield18 <- predict(fit.yd18, which = "ID")[,c(1,7)]
colnames(yield18) = c("ID", "yd18")

yield19 <- predict(fit.yd19, which = "ID")[,c(1,7)]
colnames(yield19) = c("ID", "yd19")

sw18 <- predict(fit.sw18, which = "ID")[,c(1,7)]
colnames(sw18) = c("ID", "sw18")
sw19 <- predict(fit.sw19, which = "ID")[,c(1,7)]
colnames(sw19) = c("ID", "sw19")

dm18 <- predict(fit.dm18, which = "ID")[,c(1,7)]
colnames(dm18) = c("ID", "dm18")
dm19 <- predict(fit.dm19, which = "ID")[,c(1,7)]
colnames(dm19) = c("ID", "dm19")

df18 <- predict(fit.df18, which = "ID")[,c(1,7)]
colnames(df18) = c("ID", "df18")
df19 <- predict(fit.df19, which = "ID")[,c(1,7)]
colnames(df19) = c("ID", "df19")

a18 <- predict(fit.a18, which = "ID")[,c(1,7)]
colnames(a18) = c("ID", "a18")
a19 <- predict(fit.a19, which = "ID")[,c(1,7)]
colnames(a19) = c("ID", "a19")

b18 <- predict(fit.b18, which = "ID")[,c(1,7)]
colnames(b18) = c("ID", "b18")
b19 <- predict(fit.b19, which = "ID")[,c(1,7)]
colnames(b19) = c("ID", "b19")

l18 <- predict(fit.L18, which = "ID")[,c(1,7)]
colnames(l18) = c("ID", "l18")
l19 <- predict(fit.L19, which = "ID")[,c(1,7)]
colnames(l19) = c("ID", "l19")

app18 <- predict(fit.app18, which = "ID")[,c(1,7)]
colnames(app18) = c("ID", "app18")
app19 <- predict(fit.app19, which = "ID")[,c(1,7)]
colnames(app19) = c("ID", "app19")

col18 <- predict(fit.color18, which = "ID")[,c(1,7)]
colnames(col18) = c("ID", "col18")
col19 <- predict(fit.color19, which = "ID")[,c(1,7)]
colnames(col19) = c("ID", "col19")

text18 <- predict(fit.text18, which = "ID")[,c(1,7)]
colnames(text18) = c("ID", "text18")
text19 <- predict(fit.text19, which = "ID")[,c(1,7)]
colnames(text19) = c("ID", "text19")

# Merge all predicted values into a single data frame
BBL_2020 <- join_all(list(yield18, yield19, sw18, sw19,             
                          dm18, dm19, df18, df19,
                          app18, app19, text18, text19, 
                          col18, col19, l18, l19, 
                          a18, a19, b18, b19),
                     by = "ID", type = 'full', 
                     match = "first")

head(BBL_2020)
