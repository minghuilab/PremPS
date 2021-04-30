library(forestFloor)
library(randomForest)

load("inputfiles/PremPS.RData")
test <- read.csv(file='test_outfeature',header = TRUE, sep = '\t') 
test_info <- read.csv(file='test_sunddg',header = TRUE, sep = '\t')
test_features <- test[,11:(ncol(test))]
test_features=as.data.frame(lapply(test_features,as.numeric))

ff <- forestFloor(premps.rf,data[-1],test_features,bootstrapFC = TRUE, calc_np = TRUE)
fc <- round(ff$FCmatrix,4)
fc <- as.data.frame(fc)
fc <- fc[,c('PSSM','DCS','DOMH','P_L','P_FWY','P_RKDE','N_Hydro','N_Charg','SASA_pro','SASA_sol','bootstrapFC')]
fc$PremPS <- rowSums(fc)

result <- cbind(test_info,fc[(nrow(fc)-nrow(test)+1):nrow(fc),])
write.table(result,file = 'test_outcontribution',sep = '\t',col.names = TRUE,row.names = FALSE,quote = FALSE)
