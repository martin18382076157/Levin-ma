library(tidyverse)
library(caret)
library(Publish)
library(survival)
library(glmnet)
library(ggpubr)
library(survminer)
library(rolr)
library(survIDINRI)
library(survAUC)
library(rms)
library(DCA)
library(timeROC)




dt.clinics <- read_csv('clinical_data.csv') %>% select(-1)
#clinical_data.csv means our clinical data in our study





set.seed(number)
#number means the number of seed
for(i.a in 1:10)
{
  idx.all <- 1:nrow(dt.clinics)
  idx.train <- createDataPartition(dt.clinics$Progress, p = 0.7, list = F) %>% c()
  idx.test <- setdiff(idx.all, idx.train)
  
  
  dt.train.clinics <- dt.clinics[idx.train, ]
  dt.test.clinics <- dt.clinics[idx.test, ]
  
  
  cox.test <- coxphSeries(Surv(PFS, Progress)~1, 
                          data = dt.temp, 
                          vars = colnames(dt.temp)[-c(1:2)])
  sel.names <- cox.test$Variable[(which(cox.test$Pvalue < 0.1))]
  rm(dt.temp)
  

  y <- Surv(dt.train.clinics$PFS, dt.train.clinics$Progress)
  
  cv.fit <- cv.glmnet(x, y, family = 'cox', nfolds = 10)
  fit <- glmnet(x, y, family = 'cox')
  
  fit.list[[i.a]] <- fit
  cvfit.list[[i.a]] <- cv.fit
  
  
  dt.os.train <- dt.train.clinics
  dt.os.test <- dt.test.clinics
  
  cat(i.a)
  
  dt.train.list[[i.a]] <- dt.os.train
  dt.test.list[[i.a]] <- dt.os.test
  

  res <- summary(res)
  p.list.test[[i.a]] <- res$coefficients[length(res$coefficients)]
  
  
  res <- summary(res)
  p.list.train[[i.a]] <- res$coefficients[length(res$coefficients)]
}

idx.sel <- which.min(p.list.test)
# idx.sel <- 7


fit <- fit.list[[idx.sel]]
cv.fit <- cvfit.list[[idx.sel]]
