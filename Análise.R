#################################PRE-ANALYSIS###################################
####packages####

if(!require(pacman)) install.packages("pacman")
library(pacman)

pacman::p_load(dplyr, psych, car, MASS, DescTools, QuantPsyc, ggplot2, gridExtra,nlme, reshape, dplyr,
               tidyr)

####data####

# reading data

data <- read.csv2(file = "data1.csv", dec = ".") #reading data in csv with sep = ;

# excluding subjects with relevant missing values

data <- data[-c(1,3,10,18,42,46,60,63,65,79,85,95,97),]

# data transformation

glimpse(data) #visualization variables' types

cols_fac <- c(1, 7, 8, 10:22, 26:34, 83)
data[cols_fac] <- lapply(data[cols_fac], factor) #alteration to factor

data[[4]] <- as.character(data[[4]])

cols_date <- c(5,23,24,25,35:47,51,55,59,63,67,71,75,79)
data[cols_date] <- lapply(data[cols_date], as.Date) #alteration to date

data[[82]] <- (as.double(data[[82]])) #alteration to double

glimpse(data) #visualization variables' types

data$inhibitor_use <- factor(data$inhibitor_use, levels = c(0,1), labels = c("non-inhib", "inhib"))


####creating new variables####

# tac C/D at any time point

data$cd_before_1 <- round(data$serumlvl_before_1/data$dose_before_1,2)
data$cd_before_2 <- round(data$serumlvl_before_2/data$dose_before_2,2)
data$cd_before_3 <- round(data$serumlvl_before_3/data$dose_before_3,2)
data$cd_dur_1 <- round(data$serumlvl_dur_1/data$dose_dur_1,2)
data$cd_dur_2 <- round(data$serumlvl_dur_2/data$dose_dur_2,2)
data$cd_dur_3 <- round(data$serumlvl_dur_3/data$dose_dur_3,2)
data$cd_after_1 <- round(data$serumlvl_after_1/data$dose_after_1,2)
data$cd_after_2 <- round(data$serumlvl_after_2/data$dose_after_2,2)
data$cd_after_3 <- round(data$serumlvl_after_3/data$dose_after_3,2)

# change in Tac C/D values during inib co-therapy versus tac C/D in baseline conditions

data$mean_cd_baseline <- round(rowMeans(data[, c("cd_before_1", "cd_before_2", "cd_before_3")]),2)
data$mean_cd_dur <- round(rowMeans(data[, c("cd_dur_1", "cd_dur_2", "cd_dur_3")]),2)

data$var_tac <- data$mean_cd_dur - data$mean_cd_baseline

# fast e slow metabolizers

data$gen_vel <- factor(ifelse(data$gen_receptor == "2", 0,1), levels = c(0,1), labels = c("slow", "fast"))

# intoxication tac

data$mean_serumtac1 <- round(rowMeans(data[,c("serumlvl_before_1", "serumlvl_before_2", "serumlvl_before_3")]),2)
data$mean_serumtac2 <- round(rowMeans(data[,c("serumlvl_dur_1", "serumlvl_dur_2", "serumlvl_dur_3")]),2)
data$mean_serumtac3 <- round(rowMeans(data[,c("serumlvl_after_1", "serumlvl_after_2", "serumlvl_after_3")]),2)

data$tox_tac <- factor(
ifelse(data$inhibitor_use == "inhib" & data$mean_serumtac2 >= 10 | 
                   data$inhibitor_use == "non-inhib" & data$mean_serumtac1 >= 10, 1,0),
levels = c(0,1), labels = c("non-tox","tox"))

data$var_dose <- rowMeans(data[, c("dose_dur_1","dose_dur_2","dose_dur_3")]) -
  rowMeans(data[,c("dose_before_1","dose_before_2","dose_before_3")])

data$mean_dose2 <- rowMeans(data[, c("dose_dur_1","dose_dur_2","dose_dur_3")])

data$mean_dose1 <- rowMeans(data[, c("dose_before_1","dose_before_2","dose_before_3")])

data$var_hem <- rowMeans(data[, c("hem_dur_1","hem_dur_2","hem_dur_3")]) -
  rowMeans(data[,c("hem_before_1","hem_before_2","hem_before_3")])

data$mean_hem2 <- rowMeans(data[, c("hem_dur_1","hem_dur_2","hem_dur_3")])

data$mean_hem1 <- rowMeans(data[, c("hem_before_1","hem_before_2","hem_before_3")])

####exploratory analysis####

# normality

norm_results <- lapply(data[,sapply(data, is.numeric)], shapiro.test)
norm_pvalues <- sapply(norm_results, "[[", "p.value")
as.data.frame(norm_pvalues)

# general analysis

exp_anal <- as.data.frame(summary(data[,sapply(data, is.numeric)]))

####forming groups####

# inhibitor user group

inhib_group_yes <- data[data$inhibitor_use =="inhib",]

# non-inhibitor user group

inhib_group_no <- data[data$inhibitor_use == "non-inhib",]

################################ANALYSIS########################################
####binary logistic regression####
##verifying whether users of inhibitors with slow tac metabolism are associated with tac intoxication

df_reglog <- data[, c("inhibitor_use", "gen_vel", "tox_tac")]

##models
mod <- glm(tox_tac~gen_vel + inhibitor_use,
          family = binomial(link = 'logit') ,df_reglog)

mod2 <- glm(tox_tac~gen_vel*inhibitor_use,
            family = binomial(link = 'logit') ,df_reglog)

##presumptions
plot(mod, which = 5) #outliers
summary(stdres(mod)) #outliers
vif(mod) #multicollinearity > 10

##comparing models

PseudoR2(mod, which = "Nagelkerke") ##pseudo R2
PseudoR2(mod2, which = "Nagelkerke") ##pseudo R2

AIC(mod,mod2) ## the smaller the better
BIC(mod,mod2)

anova(mod, mod2, test = "Chisq") ## check if there is a difference between the models

ClassLog(mod2, data$tox_tac)

##results

summary(mod)
exp(cbind(OR = coef(mod), confint.default(mod))) ##putting the coefficient as an exponent and confidence interval
Anova(mod)

####fisher test####

data$inhib_gen <- ifelse(data$inhibitor_use == "inhib" & data$gen_vel == "slow", 1,
                    ifelse(data$inhibitor_use == "inhib" & data$gen_vel == "fast", 2,
                           ifelse(data$inhibitor_use == "non-inhib" & data$gen_vel == "slow", 3,
                                  ifelse(data$inhibitor_use == "non-inhib" & data$gen_vel == "fast", 4,0))))



fisher.test(table(data$inhib_gen, data$tox_tac))

####slow vs fast####
## test t

##all subjects

##inhibitor group

##non-inhibitor group


####linear mixed model####

df_lmm1 <- data[,c("record_id", "serumlvl_before_1", "serumlvl_before_2", "serumlvl_before_3", "serumlvl_dur_1",
                  "serumlvl_dur_2", "serumlvl_dur_3", "serumlvl_after_1", "serumlvl_after_2", "serumlvl_after_3",
                  "hem_before_1", "hem_before_2", "hem_before_3", "hem_dur_1", "hem_dur_2", "hem_dur_3", "hem_after_1",
                  "hem_after_2", "hem_after_3", "dose_before_1", "dose_before_2", "dose_before_3", "dose_dur_1",
                  "dose_dur_2", "dose_dur_3", "dose_after_1", "dose_after_2", "dose_after_3","age", "sex", "gen_vel")]

df_lmm <- reshape(df_lmm1, 
                  varying = c("dose_before_1","hem_before_1","serumlvl_before_1",
                              "dose_before_2","hem_before_2","serumlvl_before_2",
                              "dose_before_3","hem_before_3","serumlvl_before_3",
                              "dose_dur_1","hem_dur_1","serumlvl_dur_1",
                              "dose_dur_2","hem_dur_2","serumlvl_dur_2",
                              "dose_dur_3","hem_dur_3","serumlvl_dur_3",      
                              "dose_after_1","hem_after_1","serumlvl_after_1",
                              "dose_after_2","hem_after_2","serumlvl_after_2",
                              "dose_after_3","hem_after_3","serumlvl_after_3"), 
                  idvar = "record_id", v.names = c("serumlvl", "hem", "dose"),direction = "long" )

df_lmm$inhib <- ifelse(df_lmm$time == 4 | df_lmm$time == 5 | df_lmm$time == 6, "inhib", "non_inhib" )

#acrescentar mais uma variavel pra identificar grupo controle

df_lmm$inhib <- as.factor(df_lmm$inhib)

glimpse(df_lmm)

mod_lmm <- lme(serumlvl ~ inhib + hem + gen_vel, 
               data = df_lmm, random = ~ inhib|gen_vel, method = "ML")

summary(mod_lmm)
Anova(mod_lmm)
confint.default(mod_lmm)

##############################DATA VISUALIZATION################################
####graphs####






####tables####


