\name{forecasting_ga}
\alias{forecasting_ga}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Floating-Point Genetic Algorithms with Statistical Forecast Based Inheritance Operator}
\description{
This function performs a floating-point genetic algorithm search with 
a statistical forecasting operator that generates offspring which probably 
will be generated in future generations. Use of this operator makes floating-point
genetic algorithms hybrid, that is, forecasting future generations mimics a local
search. 
}
\usage{
forecasting_ga(evalFunc, chsize, minv, maxv, crossprob, mutationprob, 
               forecastprob, elitism, popsize, maxiter, MinimumForecastLength, 
               ForecastFunction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{evalFunc}{Function to be minimized}
  \item{chsize}{Chromosome size}
  \item{minv}{Vector of lower bounds of decision variables (genes)}
  \item{maxv}{Vector of upper bounds of decision variables (genes)}
  \item{crossprob}{Crossover probability}
  \item{mutationprob}{Mutation probability}
  \item{forecastprob}{Probability of forecasting operator applied on variables (genes)}
  \item{elitism}{Number of chromosomes that are directly copied into the next generations}
  \item{popsize}{Population size}
  \item{maxiter}{Maximum number of iterations (generations)}
  \item{MinimumForecastLength}{Minimum number of ancestors of a chromosome before applying statistical forecast}
  \item{ForecastFunction}{R function to be used in forecasting. ForecastArima is implemented as R code and any other forecast functions can be implemented instead}
}

\details{
EvalFunc is always encapsulates a minimization problem. Genetic operators satisfy the boundary constraints. Since the current forecast
method is ARIMA, the parameter ForecastFunction must be manually set to ForecastArima if there is no any other function at hand. A new user-defined forecast function can easly be implemented which is based on the other well-known forecast methods. 
}

\value{
  Returns the population at last generation ordered by cost (fitness) values. First element is the best.
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

