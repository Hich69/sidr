\name{Summarize}
\alias{Summarize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summarize
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Summarize(x, addtotal = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{addtotal}{
%%     ~~Describe \code{addtotal} here~~
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
function (x, addtotal = TRUE)
{
    if (is.character(x) | is.factor(x)) {
        if (is.character(x))
            x <- factor(x)
        tab <- table(x)
        tabPercentage <- prop.table(tab) * 100
        tab <- rbind(tab, tabPercentage)
        if (isTRUE(addtotal)) {
            tab <- cbind(tab, Total = rowSums(tab))
        }
        RES <- t(tab)
        colnames(RES) <- c("freq", "perc")
    }
    else {
        res <- quantile(x = x, na.rm = TRUE)
        RES <- data.frame(n = length(x), nvalid = length(x[!is.na(x)]),
            mean = mean(x = x, na.rm = TRUE), sd = sd(x = x,
                na.rm = TRUE), min = res["0\%"], Q1 = res["25\%"],
            median = res["50\%"], Q3 = res["75\%"], max = res["100\%"],
            stringsAsFactors = FALSE)
    }
    return(RES)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
