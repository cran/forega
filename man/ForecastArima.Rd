\name{ForecastArima}
\alias{ForecastArima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Floating-Point Genetic Algorithms with Statistical Forecast Based Inheritance Operator}
\description{
This function is not called directly by the user but it is called by the forecasting_ga() for applying 
forecast operations using ARIMA. A new forecast function can be written by changing some lines of this 
function. The basic idea is to pass a chromosome's anchestors to the function, perform a forecast, return
the forecasted offspring if there is not error, else return the original chromosome as the offspring.
}
\usage{
ForecastArima(genes, Anchestors, MinimumForecastLength = 5, 
              forecastprob, robust.filter.width=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{Vector of values of chromosome}
  \item{Anchestors}{Anchestor matrix of the chromosome}
  \item{MinimumForecastLength}{If number of anchestors of the chromosome is less than the MinimumForecastLength, then statistical forecast will not be performed. The default is 5.}
  \item{forecastprob}{A gene will be forecasted if a uniform random variable is less than the forecastprob.}
  \item{robust.filter.width}{Robust filtering parameter. Default is 5.}
}

\details{
Anchestor matrix is an nxp matrix of historical gene values where n is the number of anchestors and p is the number of genes. 
Robust filtering is required to eliminate the sudden shifts which most probably happens in very early generations. Current version 
uses repeated medians as the robust filtering method. Any other robust filtering functions in package \pkg{robfilter} can be used 
in new forecasting functions rather than ForecastArima.
}

\value{
Returns a vector of values of a single offspring which is forecasted.
}

%\references{}

\author{
Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>\cr
}

%\seealso{}

\examples{

# Required package forega
require("forega")

set.seed(12345)

# This function has the global minimum at x_1 = pi and x_2 = exp(1)
f <- function (x){
	return( (x[1]-pi)^2 + (x[2]-2.71828)^2 )
}

# Performing a floating-point genetic algorithm search with forecast probability of 0.10
res <- forecasting_ga(evalFunc=f, chsize=2, minv=rep(-10.0,2),
                      maxv=rep(10.0,2), crossprob=0.80, mutationprob=0.01, 
                      popsize=100, maxiter=1000, MinimumForecastLength=20, 
                      ForecastFunction=ForecastArima, elitism=2, forecastprob=0.01)

# Show the first chromosome of the returned population matrix
print(res[1,])


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

