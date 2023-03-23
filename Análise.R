#################################PRE-ANALYSIS###################################
####packages####

if(!require(pacman)) install.packages("pacman")
library(pacman)

pacman::p_load(dplyr, psych, car, MASS, DescTools, QuantPsyc, ggplot2, gridExtra,nlme, reshape, dplyr,
               tidyr, lme4, influence.ME)

####data####

# reading data

data <- read.csv2(file = "data_inhib.csv", dec = ".") #reading data in csv with sep = ;


# excluding subjects with relevant missing values

data <- data[-c(1,3,10,18,42,46,60,63,65,79,85,95,97),]

names(data)[names(data) == "record_id"] <- "id"

# data transformation

glimpse(data) #visualization variables' types

cols_fac <- c(1, 7, 8, 10:22, 26:34, 83)
data[cols_fac] <- lapply(data[cols_fac], factor) #alteration to factor

data[[4]] <- as.character(data[[4]])

cols_date <- c(5,23,24,25,35:47,51,55,59,63,67,71,75,79)
data[cols_date] <- lapply(data[cols_date], as.Date) #alteration to date

data[[82]] <- (as.double(data[[82]])) #alteration to double
data[[61]] <- (as.double(data[[61]])) #alteration to double


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

####linear mixed model####

#dataframe for lmm analysis

df_lmm1 <- data[,c("id", "serumlvl_before_1", "serumlvl_before_2", "serumlvl_before_3", "serumlvl_dur_1",
                  "serumlvl_dur_2", "serumlvl_dur_3", "serumlvl_after_1", "serumlvl_after_2", "serumlvl_after_3",
                  "hem_before_1", "hem_before_2", "hem_before_3", "hem_dur_1", "hem_dur_2", "hem_dur_3", "hem_after_1",
                  "hem_after_2", "hem_after_3", "dose_before_1", "dose_before_2", "dose_before_3", "dose_dur_1",
                  "dose_dur_2", "dose_dur_3", "dose_after_1", "dose_after_2", "dose_after_3","age", "sex", "gen_vel",
                  "inhibitor_use")]

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
                  idvar = "id", v.names = c("serumlvl", "hem", "dose"),direction = "long" )

df_lmm$inhib <- as.factor(ifelse(df_lmm$time == 4 | df_lmm$time == 5 | df_lmm$time == 6, "using_inhib", "non_using" ))
df_lmm$cd <- df_lmm$serumlvl / df_lmm$dose #tac concentration/dose
df_lmm$cdlog <- log(df_lmm$cd)

sum(is.na(df_lmm))
df_lmm <- na.omit(df_lmm)
summary(df_lmm)

#building models


mod_lmm5 <- lmer(cdlog ~ inhib + gen_vel + hem + age + sex + (1+time|id), 
                 data = df_lmm, REML = F)

mod_lmm4 <- lmer(cdlog ~ inhib + gen_vel + hem + age + (1+time|id), 
                 data = df_lmm, REML = F)

mod_lmm3 <- lmer(cdlog ~ inhib + gen_vel + hem + (1+time|id), 
               data = df_lmm, REML = F)

mod_lmm2 <- lmer(cd ~ inhib + gen_vel + (1+time|id), 
                data = df_lmm, REML = F)

mod_lmm1 <- lmer(cd ~ inhib + (1+time|id), 
                data = df_lmm, REML = F)

mod_lmm0 <- lmer(cd ~ 1 + (1+time|id), 
                 data = df_lmm, REML = F)

summary(mod_lmm5)#$coefficients[,"Estimate"]
coef(mod_lmm3)
anova(mod_lmm1, mod_lmm0)
anova(mod_lmm2, mod_lmm1)
anova(mod_lmm3, mod_lmm2)
anova(mod_lmm4, mod_lmm3)
anova(mod_lmm5, mod_lmm4)

Anova(mod_lmm4)

coef_100 <- (exp(fixef(mod_lmm5)) -1) * 100

confint(mod_lmm5)

#presumptions tests

plot(fitted(mod_lmm5),residuals(mod_lmm5)) #Linearity and Homoskedasticity and outliers of residuals
plot(mod_lmm5, which = 5) #Linearity and Homoskedasticity and outliers of residuals
qqnorm(residuals(mod_lmm4))#normality of residuals
hist(residuals(mod_lmm4))#normality of residuals
vif(mod_lmm4) #multicollinearity > 10

##############################DATA VISUALIZATION################################
####graphs####
####tables####