#ifndef POPULATION_H_
#define POPULATION_H_

#include <Rcpp.h>
#include "Chromosome.h"
#include <vector>


using namespace std;
using namespace Rcpp;

class Population {
private:
	vector<Chromosome*> chromosomes;
	vector<Chromosome*> tempPopulation;
	vector<Chromosome*> parents;
	vector<Chromosome*> offspring;
	int popsize;
	int chsize;
	double* maxv;
	double* minv;
	int elitism;
	double crossprob;
	double mutationprob;
	double forecastprob;
public:
	Population(int popsize, int chsize, double* maxv, double *minv);
	~Population();
	int getPopulationSize();
	void setElitism(int elitism);
	int getElitism();
	void setCrossoverProb(double p);
	double getCrossoverProb();
	void setMutationProb(double p);
	double getMutationProb();
	void setForecastProb(double p);
	double getForecastProb();
	vector<Chromosome*> performTournament();
	void dump();
	void setFitness(int index, double d);
	double getFitness(int index);
	Chromosome *getChromosome(int index);
	void setChromosome(int index, Chromosome *c);
	void sortPopulation();
	void TournamentSelection(Function ForecastFunction, double forecastprob, int MinimumForecastLength);
	NumericMatrix getPopulationAsRcppMatrix();
};



#endif /* POPULATION_H_ */
