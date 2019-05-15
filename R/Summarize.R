

### Summarize ###
Summarize <- function(x, addtotal = TRUE){
  if(is.character(x) | is.factor(x)){
    if(is.character(x))
      x <- factor(x)

    tab <- table(x)
    tabPercentage <- prop.table(tab) * 100
    tab <- rbind(tab, tabPercentage)

    if(isTRUE(addtotal)){
      tab <- cbind(tab, Total = rowSums(tab))
    }

    RES <- t(tab)
    colnames(RES) <- c('freq', 'perc')
  }
  else{
    res <- quantile(x = x, na.rm = TRUE)

    RES <- data.frame(n      = length(x),
                      nvalid = length(x[!is.na(x)]),
                      mean   = mean(x = x, na.rm = TRUE),
                      sd     = sd(x = x, na.rm = TRUE),
                      min    = res['0%'],
                      Q1     = res['25%'],
                      median = res['50%'],
                      Q3     = res['75%'],
                      max    = res['100%'], stringsAsFactors = FALSE)
  }

  return(RES)
}

