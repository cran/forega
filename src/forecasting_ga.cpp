#include <Rcpp.h>
#include "Population.h"
#include "Chromosome.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector forecasting_ga (Function evalFunc,
                       int chsize, 
                       DoubleVector minv, 
                       DoubleVector maxv, 
                       double crossprob, 
                       double mutationprob, 
                       double forecastprob, 
                       int elitism, 
                       int popsize, 
                       int maxiter, 
                       int MinimumForecastLength,
                       Function ForecastFunction) {
    double *dminv = (double*)malloc(sizeof(double) * chsize);
    double *dmaxv = (double*)malloc(sizeof(double) * chsize);
    for (int i=0;i < chsize; i++){
        dminv[i] = minv[i];
        dmaxv[i] = maxv[i];
    }
    Population *pop = new Population(popsize, chsize, dminv, dmaxv);
    pop->setCrossoverProb(crossprob);
    pop->setMutationProb(mutationprob);
    pop->setForecastProb(forecastprob);
    pop->setElitism(elitism);
    for (int citer=0; citer < maxiter; citer++){
        for (int i=0; i < pop->getPopulationSize() ; i++){
            double fitness = as<double> (evalFunc(pop->getChromosome(i)->getGenesAsDoubleRcppVector()));
            pop->setFitness(i, fitness);
        }
        pop->TournamentSelection(ForecastFunction, forecastprob, MinimumForecastLength);
    }
    
    for (int i=0; i < pop->getPopulationSize() ; i++){
        double fitness = as<double> (evalFunc(pop->getChromosome(i)->getGenesAsDoubleRcppVector()));
        pop->setFitness(i, fitness);
    }
    pop->sortPopulation();
	NumericMatrix res = pop->getPopulationAsRcppMatrix();
	delete pop;
    return(res);
}



