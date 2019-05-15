
# test=as.numeric(ADSL$E)
# goldstandard=as.numeric(ADSL$F3)-1
# cutOff='youden'

AUROCs <- function(test,goldstandard,cutOff){
  
  cutOff_<-strsplit(x = cutOff,split = '_')[[1]]
  cutoffInd<-NULL
  
  fr_ <- data.frame(score = test,label = factor(goldstandard))
  fr <- na.omit(fr_)
  #prev <- Summarize(fr$label)['1','perc']/100
  prev <- table(fr$label)['1'] / NROW(fr)
  
  pred <- prediction(fr$score, fr$label)
  
  se <- performance(pred, "tpr")
  sp <- performance(pred, "tnr")
  
  ppv <- performance(pred, "ppv")
  npv <- performance(pred, "npv")
  
  #PPV <- (se@y.values[[1]]*prev)/(se@y.values[[1]]*prev+(1-sp@y.values[[1]])*(1-prev))
  #NPV <- (sp@y.values[[1]]*(1-prev))/((1-se@y.values[[1]])*prev+sp@y.values[[1]]*(1-prev))
  
  lrPos <- se@y.values[[1]]/(1-sp@y.values[[1]])
  lrNeg <- (1-se@y.values[[1]])/sp@y.values[[1]]
  
  ACC <- performance(pred, "acc")
  AUC <- performance(pred, "auc")
  
  #, direction = "<"
  myROC <- roc(goldstandard, test, ci=TRUE, direction = "<")
  sd <- sqrt(var(myROC))
  
  
  Result.Temp <- data.frame(SENS   = se@y.values[[1]],
                            SPEC   = sp@y.values[[1]],
                            PPV    = ppv@y.values[[1]],
                            NPV    = npv@y.values[[1]],
                            DA     = ACC@y.values[[1]],
                            CUTOFF = pred@cutoffs[[1]],
                            LRPos  = lrPos,
                            LRNeg  = lrNeg,
                            FP     = pred@fp[[1]],
                            TP     = pred@tp[[1]],
                            TN     = pred@tn[[1]],
                            FN     = pred@fn[[1]])
  Result.Temp <- Result.Temp[!is.infinite(Result.Temp$CUTOFF),]
  
  Result.Temp$F1_SCORE <- 2*(Result.Temp$PPV * Result.Temp$SENS) / (Result.Temp$PPV + Result.Temp$SENS)
  
  posTest <- (Result.Temp$TP + Result.Temp$FP) / NROW(fr)
  Result.Temp$NET_BENEFIT <- ((Result.Temp$TP - Result.Temp$FP) / NROW(fr)) * (posTest / (1 - posTest))
  
  if(cutOff_[1] == 'youden'){
    ### Maximizing Youden's Index ###
    Youden<-Result.Temp$SENS+Result.Temp$SPEC-1
    Youden.Ind<-which(Youden==max(Youden))
    
    cutoffInd <- Youden.Ind
  }
  if(cutOff_[1] == 'da'){
    ### Maximizing Diagnostic Accuracy ###
    DA.Ind<-which(ACC@y.values[[1]]==max(ACC@y.values[[1]]))
    
    cutoffInd <- DA.Ind
  }
  if(cutOff_[1] == 'se'){
    ### Minimizing 95% of Sensibility ###
    Se.Ind<-min(which(Result.Temp$SENS >= as.numeric(cutOff_[2])/100),na.rm=TRUE)
    
    cutoffInd <- Se.Ind
  }
  if(cutOff_[1] == 'sp'){
    ### Minimizing 95% of Specificity ###
    prem<-Result.Temp$SPEC
    ord<-order(prem)
    prem<-prem[ord]
    Sp.I<-min(which(prem >= as.numeric(cutOff_[2])/100),na.rm=TRUE)
    Sp.Ind<-which(order(ord)==Sp.I)
    
    cutoffInd <- Sp.Ind
  }
  if(cutOff_[1] == 'ppv'){
    ### Minimizing 95% of PPV ###
    prem<-Result.Temp$PPV
    ord<-order(prem)
    prem<-prem[ord]
    PPV.Ind<-which(prem >= as.numeric(cutOff_[2])/100)
    
    cutoffInd <- PPV.Ind
  }
  if(cutOff_[1] == 'npv'){
    ### Minimizing 95% of NPV ###
    prem<-Result.Temp$NPV
    ord<-order(prem)
    prem<-prem[ord]
    NPV.Ind<-which(prem >= as.numeric(cutOff_[2])/100)
    
    cutoffInd <- NPV.Ind
  }
  if(cutOff_[1] == 'fixed'){
    ### Fixed Cut Off ###
    prem<-pred@cutoffs[[1]][!is.infinite(pred@cutoffs[[1]])]
    ord<-order(prem)
    prem<-prem[ord]
    CutOffFixed.I<-min(which(prem >= as.numeric(cutOff_[2])),na.rm=TRUE)
    CutOffFixed.Ind<-which(order(ord) == CutOffFixed.I)
    
    cutoffInd <- CutOffFixed.Ind
  }
  
  
  M.Results<-matrix(NA,15,3,dimnames=list(c('N','Prev','N_Pos','AUC','SD','Cut-Off','DA','Sens','Spec','PPV','NPV','LR+','LR-','F1-Score','Net-Benefits'),
                                          c(cutOff,'Lower.CI.95','Upper.CI.95')))
  
  
  M.Results['N',cutOff] <- NROW(fr)
  M.Results['Prev',cutOff] <- prev
  M.Results['N_Pos',cutOff] <- table(fr$label)['1']
  
  M.Results['AUC',cutOff] <- myROC$auc
  M.Results['AUC','Lower.CI.95'] <- myROC$ci[1]
  M.Results['AUC','Upper.CI.95'] <- myROC$ci[3]
  
  M.Results['SD',cutOff] <- sd
  
  M.Results['Cut-Off',cutOff] <- mean(Result.Temp[cutoffInd,'CUTOFF'],na.rm = TRUE)
  
  M.Results['DA',cutOff] <- mean(Result.Temp[cutoffInd,'DA'],na.rm = TRUE)
  if(!is.na(M.Results['DA',cutOff])){
    M.Results['DA','Lower.CI.95'] <- binom.ci(M.Results['DA',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Lower']
    M.Results['DA','Upper.CI.95'] <- binom.ci(M.Results['DA',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Upper']
  }
  
  M.Results['Sens',cutOff] <- mean(Result.Temp[cutoffInd,'SENS'],na.rm = TRUE)
  if(!is.na(M.Results['Sens',cutOff])){
    M.Results['Sens','Lower.CI.95'] <- binom.ci(M.Results['Sens',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Lower']
    M.Results['Sens','Upper.CI.95'] <- binom.ci(M.Results['Sens',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Upper']
  }
  
  M.Results['Spec',cutOff] <- mean(Result.Temp[cutoffInd,'SPEC'],na.rm = TRUE)
  if(!is.na(M.Results['Spec',cutOff])){
    M.Results['Spec','Lower.CI.95'] <- binom.ci(M.Results['Spec',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Lower']
    M.Results['Spec','Upper.CI.95'] <- binom.ci(M.Results['Spec',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Upper']
  }
  
  M.Results['PPV',cutOff] <- mean(Result.Temp[cutoffInd,'PPV'],na.rm = TRUE)
  if(!is.na(M.Results['PPV',cutOff])){
    M.Results['PPV','Lower.CI.95'] <- binom.ci(M.Results['PPV',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Lower']
    M.Results['PPV','Upper.CI.95'] <- binom.ci(M.Results['PPV',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Upper']
  }
  
  M.Results['NPV',cutOff] <- mean(Result.Temp[cutoffInd,'NPV'],na.rm = TRUE)
  if(!is.na(M.Results['NPV',cutOff])){
    M.Results['NPV','Lower.CI.95'] <- binom.ci(M.Results['NPV',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Lower']
    M.Results['NPV','Upper.CI.95'] <- binom.ci(M.Results['NPV',cutOff]*NROW(fr),NROW(fr),method="exact")[1,'Upper']
  }
  
  M.Results['LR+',cutOff] <- mean(Result.Temp[cutoffInd,'LRPos'],na.rm = TRUE)
  M.Results['LR-',cutOff] <- mean(Result.Temp[cutoffInd,'LRNeg'],na.rm = TRUE)
  
  M.Results['F1-Score',cutOff] <- mean(Result.Temp[cutoffInd,'F1_SCORE'],na.rm = TRUE)
  
  M.Results['Net-Benefits',cutOff] <- mean(Result.Temp[cutoffInd,'NET_BENEFIT'],na.rm = TRUE)
  
  
  
  Result.List <- list()
  Result.List[['Tab1']] <- M.Results
  Result.List[['Tab2']] <- Result.Temp
  Result.List[['Tab3']] <- myROC
  
  return(Result.List)
}
