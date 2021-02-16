#! /usr/bin/Rscript
library("lme4")
library("car")

data<-read.csv('../output/data/erd-ers_amplitude.csv')

model <- lmer(Amplitude ~ Epoch+(1|Subject), data = data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
results<-Anova(model, type = 2)
print(results)
