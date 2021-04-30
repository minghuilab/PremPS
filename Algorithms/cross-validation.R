library(randomForest)
library(stringr)

rm(list=ls())
setwd(YourWorkDirectory)

f <- read.csv(file='Datasets/S5296/S5296.txt', header = TRUE, sep = '\t')
rownames(f) <- paste(f$PDB.Id,f$Mutated.Chain,f$Mutation_PDB,f$Label,sep = '_')
label <- 'DDGexp~PSSM+DCS+DOMH+P_L+P_FWY+P_RKDE+N_Hydro+N_Charg+SASA_pro+SASA_sol'

set.seed(100)
model.rf <- randomForest(as.formula(label), data = f)
f$PremPS <- model.rf$predicted
cor(f$DDGexp,f$PremPS)

## CV1: randomly chose 80% of mutations from the S5296 set to train the model 
##      and used the remaining mutations for blind testing; 
##      the procedures were repeated 100 times.
f$Label_select <- paste(f$PDB.Id,f$Mutated.Chain,f$Mutation_PDB,sep = '_')
f.forward <- f[f$Label=='forward',]
f.reverse <- f[f$Label=='reverse',]
R_CV1 <- c()
for (i in 1:100) {
  index <- sample(nrow(f.forward),floor(nrow(f.forward)*0.8))
  train.forward <- f.forward[index,]
  test.forward <- f.forward[-index,]
  train.f <- f[f$Label_select %in% train.forward$Label_select,]
  test.f <- f[f$Label_select %in% test.forward$Label_select,]

  set.seed(i)
  model.rf <- randomForest(as.formula(label), data = train.f)
  test.f$PremPS <- predict(model.rf,test.f)
  R_CV1 <- c(R_CV1,cor(test.f$DDGexp, test.f$PremPS))
}
mean(R_CV1)


## CV2: randomly chose 50% of mutations from the S5296 set to train the model 
##      and used the remaining mutations for blind testing; 
##      the procedures were repeated 100 times.
f$Label_select <- paste(f$PDB.Id,f$Mutated.Chain,f$Mutation_PDB,sep = '_')
f.forward <- f[f$Label=='forward',]
f.reverse <- f[f$Label=='reverse',]
R_CV2 <- c()
for (i in 1:100) {
  index <- sample(nrow(f.forward),floor(nrow(f.forward)*0.5))
  train.forward <- f.forward[index,]
  test.forward <- f.forward[-index,]
  train.f <- f[f$Label_select %in% train.forward$Label_select,]
  test.f <- f[f$Label_select %in% test.forward$Label_select,]

  set.seed(i)
  model.rf <- randomForest(as.formula(label), data = train.f)
  test.f$PremPS <- predict(model.rf,test.f)
  R_CV2 <- c(R_CV2,cor(test.f$DDGexp, test.f$PremPS))
}
mean(R_CV2)


## CV3: a subset was created by randomly sampling up to 20 mutations for each protein from S5296; 
##      the procedure was repeated 10 times and resulted in 1,704 mutations in each subset. 
##      Then 80% mutations were randomly selected from each subset to train the model 
##      and the rest of mutations were used for testing, repeated 10 times.
f$Label_select <- paste(f$PDB.Id,f$Mutated.Chain,f$Mutation_PDB,sep = '_')
pdblist <- unique(f$PDB.Id)
R_CV3 <- c()
for (j in 1:10){
  data_subset <- c()
  for (pdb in pdblist){
    data_pdb <- f[f$PDB.Id==pdb,]
    num_muts <- dim(data_pdb)[1]
    if (num_muts <= 20){
      data_subset <- rbind(data_subset,data_pdb)
    }else{
      index <- sample(num_muts/2,10)
      data_pdb_subset_p <- data_pdb[data_pdb$Label == 'forward',][index,]
      data_pdb_subset_n <- data_pdb[data_pdb$Label_select %in% data_pdb_subset_p$Label_select & 
                                      data_pdb$Label=='reverse',]
      data_pdb_subset <- rbind(data_pdb_subset_p,data_pdb_subset_n)
      data_subset <- rbind(data_subset,data_pdb_subset)
    }
  }
  for (i in 1:10) {
    index <- sample(nrow(data_subset),floor(nrow(data_subset)*0.8))
    train.f <- data_subset[index,]
    test.f <- data_subset[-index,]
    
    set.seed(i)
    model.rf <- randomForest(as.formula(label), data = train.f)
    test.f$PremPS <- predict(model.rf,test.f)
    R_CV3 <- c(R_CV3,cor(test.f$DDGexp, test.f$PremPS))
  }
}
mean(R_CV3)


## CV4: trained on all mutations from 130 protein structures and 
##      the rest of protein/mutations were used to evaluate the performance. 
##      This procedure was repeated for each protein and its mutations.
pdblist <- unique(f$PDB.Id)
R_CV4 <- c()
for (pdb in pdblist) {
  train.f <- f[f$PDB.Id!=pdb,]
  test.f  <- f[f$PDB.Id==pdb,]

  set.seed(100)
  model.rf <- randomForest(as.formula(label), data = train.f)
  test.f$PremPS.CV4. <- predict(model.rf,test.f)
  R_CV4 <- rbind(R_CV4,test.f)
}
cor(R_CV4$DDGexp,R_CV4$PremPS.CV4.)


## CV5: Not only a protein in the validation set was removed from the training set, 
##      but also all other “similar proteins” with more than 25% sequence identity to this protein, 
##      repeated for each protein cluster.
f$Label_select <- rownames(f)
f$PDBChain <- paste(f$PDB.Id,f$Mutated.Chain,sep = '')
pdblist <- unique(f$PDB.Id)
R_CV5 <- c()
for (pdb in pdblist){
  test.f <- c()
  data_pdb <- f[f$PDB.Id==pdb,]
  test.f <- rbind(test.f,data_pdb)
  for (sametype in strsplit(as.character(unique(data_pdb$similar.proteins)),';')[[1]]){ 
    if (sametype!=as.character(unique(data_pdb$PDBChain))){
      test.f <- rbind(test.f,f[f$PDBChain==sametype,])
    }
  }
  train.f <- f[!f$Label_select %in% test.f$Label_select,]
  ########################## Random Forest ##
  set.seed(100)
  model.rf <- randomForest(as.formula(label), data = train.f)
  data_pdb$PremPS.CV5. <- predict(model.rf,data_pdb)
  R_CV5 <- rbind(R_CV5,data_pdb)
}
cor(R_CV5$DDGexp,R_CV5$PremPS.CV5.)
