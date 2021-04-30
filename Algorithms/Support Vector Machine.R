library(e1071)

rm(list=ls())
setwd(YourWorkDirectory)

f <- read.csv(file='Datasets/S5296/S5296.txt', header = TRUE, sep = '\t')
label <- 'DDGexp~PSSM+DCS+DOMH+P_L+P_FWY+P_RKDE+N_Hydro+N_Charg+SASA_pro+SASA_sol'
set.seed(100)
model.svm <- svm(as.formula(label),data = f)
f$PremPS <- fitted(model.svm)
cor(f$DDGexp,f$PremPS)
