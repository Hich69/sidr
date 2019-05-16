\name{Desc_Stat}
\alias{Desc_Stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Descriptive statistics
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Desc_Stat(x, digits = 2, percdigs = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
x vector to be analyzed
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
digits number of digits
}
  \item{percdigs}{
%%     ~~Describe \code{percdigs} here~~
percdigs number of digits for a percentage
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, digits = 2, percdigs = 1)
{
    fRound <- function(value, digits) {
        return(sprintf(paste("\%.", digits, "f", sep = ""), value))
    }
    fRoundPerc <- function(value, percdigs) {
        return(sprintf(paste("\%.", percdigs, "f", sep = ""),
            value))
    }
    if (is.character(x) | is.factor(x)) {
        if (is.character(x))
            x <- factor(x)
        tab <- table(x)
        tabPercentage <- prop.table(tab)
        tabPercentage <- tabPercentage * 100
        RES <- data.frame(STAT = names(tab), RES = paste(tab,
            " (", fRoundPerc(tabPercentage, percdigs), " \%)",
            sep = ""), stringsAsFactors = FALSE)
        RES <- rbind(RES, data.frame(STAT = "Total", RES = paste(sum(tab),
            " (", fRoundPerc(sum(tabPercentage), 0), " \%)", sep = ""),
            stringsAsFactors = FALSE))
        RES <- data.frame(STATN = 1:NROW(RES), RES, stringsAsFactors = FALSE)
    }
    else {
        res <- quantile(x = x, na.rm = TRUE)
        RES <- data.frame()
        RES <- rbind(RES, data.frame(STAT = "n", RES = length(x[!is.na(x)]),
            stringsAsFactors = FALSE), data.frame(STAT = "Mean (Std)",
            RES = paste(fRound(mean(x = x, na.rm = TRUE), digits),
                " (", fRound(sd(x = x, na.rm = TRUE), digits),
                ")", sep = ""), stringsAsFactors = FALSE), data.frame(STAT = "Median",
            RES = fRound(res["50\%"], digits), stringsAsFactors = FALSE),
            data.frame(STAT = "Q1 ; Q3", RES = paste(fRound(res["25\%"],
                digits), fRound(res["75\%"], digits), sep = " ; "),
                stringsAsFactors = FALSE), data.frame(STAT = "IQR",
                RES = fRound(res["75\%"] - res["25\%"], digits),
                stringsAsFactors = FALSE), data.frame(STAT = "Min ; Max",
                RES = paste(fRound(res["0\%"], digits), fRound(res["100\%"],
                  digits), sep = " ; "), stringsAsFactors = FALSE))
        RES <- data.frame(STATN = 1:NROW(RES), RES, stringsAsFactors = FALSE)
    }
    return(RES)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
