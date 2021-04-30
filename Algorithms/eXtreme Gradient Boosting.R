library(xgboost)

rm(list=ls())
setwd(YourWorkDirectory)

f <- read.csv(file='Datasets/S5296/S5296.txt', header = TRUE, sep = '\t')
label <- 'DDGexp~PSSM+DCS+DOMH+P_L+P_FWY+P_RKDE+N_Hydro+N_Charg+SASA_pro+SASA_sol'
set.seed(100)
feature_list <- unlist(strsplit(unlist(strsplit(label,split = '~'))[2],split = '[+]'))
data <- f[,which(names(f) %in% feature_list)]
data <- data[,order(names(data))]
train_data <- xgb.DMatrix(data=as.matrix(data), label=f$DDGexp)
xgb.best <- list(max_depth = 10,eta = 0.01,gamma = 0,colsample_bytree = 0.8,min_child_weight = 2,
                 subsample = 0.7)
xgb_tune <- xgb.train(params=xgb.best, data=train_data, nrounds=1000)
f$PremPS <- predict(xgb_tune, train_data)
cor(f$DDGexp,f$PremPS)
