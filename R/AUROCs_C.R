
# clean.env()
# set.WorkingDir()
# 
# a <- read.csv2('0 - Stat Management/For_MD/dataset_FibroVet.csv', stringsAsFactors = FALSE)

# test1 = as.numeric(a$SCORE[a$TSFL == 'Y']); goldstandard1 = a$OUTCOME[a$TSFL == 'Y']; cutOff1 = 'youden';
# test2 = as.numeric(a$SCORE[a$FASFL == 'Y']); goldstandard2 = a$OUTCOME[a$FASFL == 'Y']; cutOff2 = 'se_80';

# a <- ADSL[ADSL$STUDYID == 'M118',]
# test1=as.numeric(a$FMVCTE); goldstandard1=as.numeric(a$F3)-1; cutOff1='youden'
# test2=as.numeric(a$E); goldstandard2=as.numeric(a$F3)-1; cutOff2='youden'


AUROCs_C <- function(test1,goldstandard1,cutOff1,test2,goldstandard2,cutOff2 = NULL){
  
  auroc1 <- AUROCs(test = test1,goldstandard = goldstandard1,cutOff = cutOff1)
  
  cutoff.Getted <- auroc1$Tab1['Cut-Off',1]
  
  if(is.null(cutOff2)){
    auroc2 <- AUROCs(test = test2,goldstandard = goldstandard2,cutOff = paste('fixed_', cutoff.Getted, sep = ''))
  }
  else{
    auroc2 <- AUROCs(test = test2,goldstandard = goldstandard2,cutOff = cutOff2)
  }
  
  dat1 <- data.frame(test1         = test1,
                     goldstandard1 = goldstandard1, stringsAsFactors = FALSE)
  dat2 <- data.frame(test2         = test2,
                     goldstandard2 = goldstandard2, stringsAsFactors = FALSE)
  
  RES <- data.frame(TEST1  = auroc1$Tab1[c("Prev","AUC","Cut-Off","DA","Sens","Spec","PPV","NPV"),1],
                    TEST2  = auroc2$Tab1[c("Prev","AUC","Cut-Off","DA","Sens","Spec","PPV","NPV"),1],
                    PVALUE = NA,
                    METHOD = '', stringsAsFactors = FALSE)
  
  ### Contingency table for TEST 1 ###
  auroc1_  <- tail(auroc1$Tab2[auroc1$Tab2$CUTOFF >= auroc1$Tab1['Cut-Off',1],], 1)
  auroc1_T <- data.frame(DISEASE   = c(auroc1_$TP, auroc1_$FN),
                         NODISEASE = c(auroc1_$FP, auroc1_$TN),stringsAsFactors = FALSE)
  rownames(auroc1_T) <- c('POSITIVE.TEST','NEGATIVE.TEST')
  
  ### Contingency table for TEST 2 ###
  auroc2_  <- tail(auroc2$Tab2[auroc2$Tab2$CUTOFF >= auroc2$Tab1['Cut-Off',1],], 1)
  auroc2_T <- data.frame(DISEASE   = c(auroc2_$TP, auroc2_$FN),
                         NODISEASE = c(auroc2_$FP, auroc2_$TN),stringsAsFactors = FALSE)
  rownames(auroc2_T) <- c('POSITIVE.TEST','NEGATIVE.TEST')
  
  ### paired test ###
  if((NROW(dat1) == NROW(dat2)) && (dat1$goldstandard1 == dat2$goldstandard2)){
    
    dat <- data.frame(GLD   = dat1$goldstandard1,
                      TEST1 = dat1$test1,
                      TEST2 = dat2$test2, stringsAsFactors = FALSE)
    dat$TEST1 <- ifelse(dat$TEST1 >= auroc1$Tab1['Cut-Off',1], 1, 0)
    dat$TEST2 <- ifelse(dat$TEST2 >= auroc2$Tab1['Cut-Off',1], 1, 0)
    
    
    f.factor <- function(x){ return(factor(x, levels = c('1','0')))}
    
    NO.DISEASE <- table(TEST2 = f.factor(dat$TEST2[dat$GLD == 0]), TEST1 = f.factor(dat$TEST1[dat$GLD == 0]))
    NO.DISEASE <- addmargins(NO.DISEASE)
    
    DISEASE    <- table(TEST2 = f.factor(dat$TEST2[dat$GLD == 1]), TEST1 = f.factor(dat$TEST1[dat$GLD == 1]))
    DISEASE    <- addmargins(DISEASE)
    
    ### Elements ... ###
    m011 <- NO.DISEASE['1','1'];   m010 <- NO.DISEASE['1','0'];   r21 <- NO.DISEASE['1','Sum']
    m001 <- NO.DISEASE['0','1'];   m000 <- NO.DISEASE['0','0'];   r20 <- NO.DISEASE['0','Sum']
    r11  <- NO.DISEASE['Sum','1']; r10  <- NO.DISEASE['Sum','0']; n0  <- NO.DISEASE['Sum','Sum']
    
    m111 <- DISEASE['1','1'];      m110 <- DISEASE['1','0'];      s21 <- DISEASE['1','Sum']
    m101 <- DISEASE['0','1'];      m100 <- DISEASE['0','0'];      s20 <- DISEASE['0','Sum']
    s11  <- DISEASE['Sum','1'];    s10  <- DISEASE['Sum','0'];    n1  <- DISEASE['Sum','Sum']
    
    
    ### For Diagnostic Accuracy
    dat3 <- data.frame(TEST1 = factor(ifelse(dat$GLD == dat$TEST1,1,0), levels = c('1','0')),
                       TEST2 = factor(ifelse(dat$GLD == dat$TEST2,1,0), levels = c('1','0')))
    dat3_ <- table(dat3$TEST1,dat3$TEST2)
    
    #########################################################################################################
                                  ### Sensitivity : Exact test (small sample) ###
    #########################################################################################################
    if((m110 + m101) < 20){
      if((m110 != 0) & (m101 != 0)){
        ### Sensitivity ###
        m <- m110 + m101
        k <- min(m110, m101)
        p_j <- function(j){ return(choose(n = m,k = j) * (1/2)^m)}
        pvalue <- 2 * sum(sapply(X = 1:k, FUN = p_j))
        
        RES$PVALUE[rownames(RES) == 'Sens'] <- pvalue
        RES$METHOD[row.names(RES) == c("Sens")] <- 'Exact'
      }
      else if(m110 == 0 & m101 == 0){
        pvalue <- 1
        
        RES$PVALUE[rownames(RES) == 'Sens'] <- pvalue
        RES$METHOD[row.names(RES) == c("Sens")] <- 'Exact'
      }
      else{
        X_2 <- (m110 - m101)^2/(m110 + m101)
        pvalue <- 1 - pchisq(q = X_2,df = (2-1)*(2-1))
        
        RES$PVALUE[rownames(RES) == 'Sens'] <- pvalue
        RES$METHOD[row.names(RES) == c("Sens")] <- 'McNemar'
      }
    }
    
    #########################################################################################################
                                  ### Specificity : Exact test (small sample) ###
    #########################################################################################################
    if(m010 + m001 < 20){
      if(m010 != 0 & m001 != 0){
        ### Specificity ###
        m <- m010 + m001
        k <- min(m010, m001)
        p_j <- function(j){ return(choose(n = m,k = j) * (1/2)^m)}
        pvalue <- 2 * sum(sapply(X = 1:k, FUN = p_j))
        
        RES$PVALUE[rownames(RES) == 'Spec'] <- pvalue
        RES$METHOD[row.names(RES) == c("Spec")] <- 'Exact'
      }
      else if(m010 == 0 & m001 == 0){
        pvalue <- 1
        
        RES$PVALUE[rownames(RES) == 'Spec'] <- pvalue
        RES$METHOD[row.names(RES) == c("Spec")] <- 'Exact'
      }
      else{
        ### Specificity ###
        X_2 <- (m010 - m001)^2/(m010 + m001)
        pvalue <- 1 - pchisq(q = X_2,df = (2-1)*(2-1))
        
        RES$PVALUE[rownames(RES) == 'Spec'] <- pvalue
        RES$METHOD[row.names(RES) == c("Spec")] <- 'McNemar'
      }
    }
    
    #########################################################################################################
                                  ### Sensitivity : Mc Nemar test ###
    #########################################################################################################
    if(m110 + m101 >= 20){
      ### Sensitivity ###
      X_2 <- (m110 - m101)^2/(m110 + m101)
      pvalue <- 1 - pchisq(q = X_2,df = (2-1)*(2-1))
      
      RES$PVALUE[rownames(RES) == 'Sens'] <- pvalue
      RES$METHOD[row.names(RES) == c("Sens")] <- 'McNemar'
    }
    
    #########################################################################################################
                                  ### Specificity : Mc Nemar test ###
    #########################################################################################################
    if(m010 + m001 >= 20){
      ### Specificity ###
      X_2 <- (m010 - m001)^2/(m010 + m001)
      pvalue <- 1 - pchisq(q = X_2,df = (2-1)*(2-1))
      
      RES$PVALUE[rownames(RES) == 'Spec'] <- pvalue
      RES$METHOD[row.names(RES) == c("Spec")] <- 'McNemar'
    }
    
    #########################################################################################################
                                  ### Diagnostic Accuracy : Exact test (small sample) ###
    #########################################################################################################
    if(dat3_['1','0'] + dat3_['0','1'] < 20){
      if(dat3_['1','0'] != 0 & dat3_['0','1'] != 0){
        m <- dat3_['1','0'] + dat3_['0','1']
        k <- min(dat3_['1','0'], dat3_['0','1'])
        p_j <- function(j){ return(choose(n = m,k = j) * (1/2)^m)}
        pvalue <- 2 * sum(sapply(X = 1:k, FUN = p_j))
        
        RES$PVALUE[rownames(RES) == 'DA'] <- pvalue
        RES$METHOD[row.names(RES) == c("DA")] <- 'Exact'
      }
      else if(dat3_['1','0'] == 0 & dat3_['0','1'] == 0){
        pvalue <- 1
        
        RES$PVALUE[rownames(RES) == 'DA'] <- pvalue
        RES$METHOD[row.names(RES) == c("DA")] <- 'Exact'
      }
      else{
        X_2 <- (dat3_['1','0'] - dat3_['0','1'])^2/(dat3_['1','0'] + dat3_['0','1'])
        pvalue <- 1 - pchisq(q = X_2,df = (2-1)*(2-1))
        
        RES$PVALUE[rownames(RES) == 'DA'] <- pvalue
        RES$METHOD[row.names(RES) == c("DA")] <- 'McNemar'
      }
    }
    
    #########################################################################################################
                                  ### Diagnostic Accuracy : Mc Nemar test ###
    #########################################################################################################
    if(dat3_['1','0'] + dat3_['0','1'] >= 20){
      X_2 <- (dat3_['1','0'] - dat3_['0','1'])^2/(dat3_['1','0'] + dat3_['0','1'])
      pvalue <- 1 - pchisq(q = X_2,df = (2-1)*(2-1))
      
      RES$PVALUE[rownames(RES) == 'DA'] <- pvalue
      RES$METHOD[row.names(RES) == c("DA")] <- 'McNemar'
    }
    
    #########################################################################################################
                                  ### Positive Predictive Test  ###
    #########################################################################################################
    
    ### Positive Predictive Value ###
    Z_ <- (m111 + m110 + m011 + m010)/(2 * m111 + m110 + m101 + 2 * m011 + m010 + m001)
    D_ <- (2 * m111 + m110 + m101)/(2 * m111 + m110 + m101 + 2 * m011 + m010 + m001)
    
    Numerator <- m111 * (1 - 2 * Z_) + m110 * (1 - Z_) + m101 * (0 - Z_)
    Numerator <- Numerator^2
    
    Denominator <- (1 - D_)^2 * (m111 * (1 - 2 * Z_)^2 + m110 * (1 - Z_)^2 + m101 * (0 - Z_)^2) + 
      (0 - D_)^2 * (m011 * (1 - 2 * Z_)^2 + m010 * (1 - Z_)^2 + m001 * (0 - Z_)^2)
    
    X_2 <- Numerator/Denominator
    pvalue <- 1 - pchisq(q = X_2,df = (2-1)*(2-1),lower.tail = FALSE)
    
    ### using package ###
    hsd <- tab.paired(d=dat$GLD, y1=dat$TEST1, y2=dat$TEST2)
    res <- pv.rpv(hsd)
    
    RES$PVALUE[rownames(RES) == 'PPV'] <- res$ppv$p.value
    RES$METHOD[row.names(RES) == c("PPV")] <- 'Chisq'
    
    #########################################################################################################
                                  ### Negative Predictive Test  ###
    #########################################################################################################
    
    RES$PVALUE[rownames(RES) == 'NPV'] <- res$npv$p.value
    RES$METHOD[row.names(RES) == c("NPV")] <- 'Chisq'
    
    ### Prevalence ###
    RES$PVALUE[rownames(RES) == 'Prev'] <- 1
    RES$METHOD[rownames(RES) == 'Prev'] <- 'McNemar'
    
  }
  ### unpaired test ###
  else{
    
    ### Prev ###
    N.Prev1 <- colSums(auroc1_T)['DISEASE']
    N.Tot1  <- sum(auroc1_T)
    N.Prev2 <- colSums(auroc2_T)['DISEASE']
    N.Tot2  <- sum(auroc2_T)
    
    pvalue <- chisq.test(x = rbind(c(N.Prev1,N.Tot1 - N.Prev1),
                                   c(N.Prev2,N.Tot2 - N.Prev2)))
    # pvalue <- prop.test(x = c(N.Prev1,N.Prev2),n = c(N.Tot1,N.Tot2))
    RES$PVALUE[rownames(RES) == 'Prev'] <- pvalue$p.value
    #RES$METHOD[rownames(RES) == 'Prev'] <- 'Chi2'
    
    ### DA ###
    N.DA1 <- sum(diag(x = as.matrix(auroc1_T)))
    N.DA2 <- sum(diag(x = as.matrix(auroc2_T)))
    
    pvalue <- chisq.test(x = rbind(c(N.DA1,N.Tot1 - N.DA1),
                                   c(N.DA2,N.Tot2 - N.DA2)))
    # pvalue <- prop.test(x = c(N.DA1,N.DA2),n = c(N.Tot1,N.Tot2))
    RES$PVALUE[rownames(RES) == 'DA'] <- pvalue$p.value
    
    ### SENSITIVITY ###
    N.TP1      <- auroc1_T['POSITIVE.TEST','DISEASE']
    N.DISEASE1 <- colSums(auroc1_T)['DISEASE']
    N.TP2      <- auroc2_T['POSITIVE.TEST','DISEASE']
    N.DISEASE2 <- colSums(auroc2_T)['DISEASE']
    
    pvalue <- chisq.test(x = rbind(c(N.TP1,N.DISEASE1 - N.TP1),
                                   c(N.TP2,N.DISEASE2 - N.TP2)))
    # pvalue <- prop.test(x = c(N.TP1,N.TP2),n = c(N.DISEASE1,N.DISEASE2))
    RES$PVALUE[rownames(RES) == 'Sens'] <- pvalue$p.value
    
    ### SPECIFICITY ###
    N.TN1        <- auroc1_T['NEGATIVE.TEST','NODISEASE']
    N.NODISEASE1 <- colSums(auroc1_T)['NODISEASE']
    N.TN2        <- auroc2_T['NEGATIVE.TEST','NODISEASE']
    N.NODISEASE2 <- colSums(auroc2_T)['NODISEASE']
    
    pvalue <- chisq.test(x = rbind(c(N.TN1,N.NODISEASE1 - N.TN1),
                                   c(N.TN2,N.NODISEASE2 - N.TN2)))
    # pvalue <- prop.test(x = c(N.TN1,N.TN2),n = c(N.NODISEASE1,N.NODISEASE2))
    RES$PVALUE[rownames(RES) == 'Spec'] <- pvalue$p.value
    
    ### PPV ###
    N.POSITIVE.TEST1 <- rowSums(auroc1_T)['POSITIVE.TEST']
    N.POSITIVE.TEST2 <- rowSums(auroc2_T)['POSITIVE.TEST']
    
    pvalue <- chisq.test(x = rbind(c(N.TP1,N.POSITIVE.TEST1 - N.TP1),
                                   c(N.TP2,N.POSITIVE.TEST2 - N.TP2)))
    # pvalue <- prop.test(x = c(N.TP1,N.TP2),n = c(N.POSITIVE.TEST1,N.POSITIVE.TEST2))
    RES$PVALUE[rownames(RES) == 'PPV'] <- pvalue$p.value
    
    ### NPV ###
    N.NEGATIVE.TEST1 <- rowSums(auroc1_T)['NEGATIVE.TEST']
    N.NEGATIVE.TEST2 <- rowSums(auroc2_T)['NEGATIVE.TEST']
    
    pvalue <- chisq.test(x = rbind(c(N.TN1,N.NEGATIVE.TEST1 - N.TN1),
                                   c(N.TN2,N.NEGATIVE.TEST2 - N.TN2)))
    # pvalue <- prop.test(x = c(N.TN1,N.TN2),n = c(N.NEGATIVE.TEST1,N.NEGATIVE.TEST2))
    RES$PVALUE[rownames(RES) == 'NPV'] <- pvalue$p.value
    
    ### Test chosen ###
    RES$METHOD[row.names(RES) %in% c("Prev","DA","Sens","Spec","PPV",'NPV')] <- 'Chi2'
    
  }
  
  ### AUC ###
  pvalue <- roc.test(auroc1$Tab3, auroc2$Tab3)
  
  RES$PVALUE[rownames(RES) == 'AUC'] <- pvalue$p.value
  RES$METHOD[rownames(RES) == 'AUC'] <- 'DeLong'
  
  return(RES)
  
}


