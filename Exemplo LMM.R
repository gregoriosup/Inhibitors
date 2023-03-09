
if(!require(pacman)) install.packages("pacman")
library(pacman)

pacman::p_load(dplyr, psych, car, MASS, DescTools, QuantPsyc, ggplot2, gridExtra,nlme, reshape, dplyr,
               tidyr, lme4)

####data####
politeness=
  read.csv("http://www.bodowinter.com/tutorial/politeness_data.csv")

which(is.na(politeness$frequency))

####exploratory analysis####

boxplot(frequency ~ attitude*gender,
        col=c("white","lightgray"),politeness)


####lmm####

politeness.model = lmer(frequency ~ attitude +
                          (1|subject) + (1|scenario), data=politeness)
summary(politeness.model)

politeness.model = lmer(frequency ~ attitude +
                          gender + (1|subject) +
                          (1|scenario), data=politeness)

#comparing random intercept

politeness.null = lmer(frequency ~ gender +
                         (1|subject) + (1|scenario), data=politeness,
                       REML=FALSE)

politeness.model = lmer(frequency ~ attitude +
                          gender + (1|subject) + (1|scenario),
                        data=politeness, REML=FALSE)

anova(politeness.null,politeness.model)

#comparing random intercept and slopes

politeness.model = lmer(frequency ~ attitude +
                          gender + (1+attitude|subject) +
                          (1+attitude|scenario),
                        data=politeness,
                        REML=FALSE)

coef(politeness.model)

politeness.null = lmer(frequency ~ gender +
                         (1+attitude|subject) + (1+attitude|scenario),
                       data=politeness, REML=FALSE)

anova(politeness.null,politeness.model)

