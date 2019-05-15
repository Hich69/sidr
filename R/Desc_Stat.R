

Desc_Stat <- function(x, digits = 2, percdigs = 1){

  fRound <- function(value, digits){
    #return(round(x = value, digits = digits))
    return(sprintf(paste('%.', digits, 'f', sep = ''), value))
  }

  fRoundPerc <- function(value, percdigs){
    return(sprintf(paste('%.', percdigs, 'f', sep = ''), value))
  }

  if(is.character(x) | is.factor(x)){
    if(is.character(x))
      x <- factor(x)

    tab <- table(x)
    tabPercentage <- prop.table(tab)
    tabPercentage <- tabPercentage * 100
    RES <- data.frame(STAT = names(tab),
                      RES  = paste(tab, ' (', fRoundPerc(tabPercentage, percdigs), ' %)' , sep = ''), stringsAsFactors = FALSE)
    RES <- rbind(RES,
                 data.frame(STAT = 'Total',
                            RES  = paste(sum(tab), ' (', fRoundPerc(sum(tabPercentage), 0), ' %)', sep = ''), stringsAsFactors = FALSE))
    RES <- data.frame(STATN = 1:NROW(RES),
                      RES, stringsAsFactors = FALSE)
  }
  else{
    res <- quantile(x = x, na.rm = TRUE)

    RES <- data.frame()
    RES <- rbind(RES,
                 data.frame(STAT = 'n',
                            RES  = length(x[!is.na(x)]), stringsAsFactors = FALSE),
                 data.frame(STAT = 'Mean (Std)',
                            RES  = paste(fRound(mean(x = x, na.rm = TRUE), digits), ' (', fRound(sd(x = x, na.rm = TRUE), digits), ')', sep = ''), stringsAsFactors = FALSE),
                 data.frame(STAT = 'Median',
                            RES  = fRound(res['50%'], digits), stringsAsFactors = FALSE),
                 data.frame(STAT = 'Q1 ; Q3',
                            RES  = paste(fRound(res['25%'], digits), fRound(res['75%'], digits), sep = ' ; '), stringsAsFactors = FALSE),
                 data.frame(STAT = 'IQR',
                            RES  = fRound(res['75%'] - res['25%'], digits), stringsAsFactors = FALSE),
                 data.frame(STAT = 'Min ; Max',
                            RES  = paste(fRound(res['0%'], digits), fRound(res['100%'], digits), sep = ' ; '), stringsAsFactors = FALSE))
    RES <- data.frame(STATN = 1:NROW(RES),
                      RES, stringsAsFactors = FALSE)
  }

  return(RES)

}

