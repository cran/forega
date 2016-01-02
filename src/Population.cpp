#include <Rcpp.h>
#include "Population.h"
#include <vector>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Rcpp;

/*
 * Non-class member comparator function for fast sorting of chromosomes
 */
bool chromosome_sorter(const Chromosome *c1, const Chromosome *c2) {
	return (c1->fitness < c2->fitness);
}

Population::~Population(){
    //cout << "Deleting " << Chromosome::chromosomesEverCreated.size() << " chromosomes "<<endl;
	for (unsigned int i=0; i < Chromosome::chromosomesEverCreated.size();i++){
	    Chromosome *c = Chromosome::chromosomesEverCreated[i];
	    delete c;
	}
	Chromosome::chromosomesEverCreated.clear();
	//cout << "OK" << endl;
}

Population::Population(int popsize, int chsize, double *maxv, double *minv) {
	//srand(time(NULL));
	chromosomes.clear();
	for (int i = 0; i < popsize; i++) {
		Chromosome *c = new Chromosome(chsize);
		c->setMaxValues(maxv);
		c->setMinValues(minv);
		c->setFather(NULL);
		c->setMother(NULL);
		c->randomize();
		chromosomes.push_back(c);
	}
	this->crossprob = 0.80;
	this->mutationprob = 0.01;
	this->elitism = 0;
	this->popsize = popsize;
	this->chsize = chsize;
}

int Population::getPopulationSize(){
    return(this->popsize);    
}

void Population::setElitism(int elitism) {
	this->elitism = elitism;
}

int Population::getElitism() {
	return (this->elitism);
}

void Population::setCrossoverProb(double p) {
	this->crossprob = p;
}

double Population::getCrossoverProb() {
	return (this->crossprob);
}

void Population::setMutationProb(double p) {
	this->mutationprob = p;
}

double Population::getMutationProb() {
	return (this->mutationprob);
}

void Population::setForecastProb(double d){
    this->forecastprob = d;
}

double Population::getForecastProb(){
    return(this->forecastprob);
}


vector<Chromosome*> Population::performTournament() {
	Chromosome *parent1, *parent2;
	vector<Chromosome*> result;
	int index1, index2, index3, index4;
	NumericVector indices = runif(4, 0, this->popsize);
	index1 = indices[0];  index2 = indices[1]; index3 = indices[2]; index4 = indices[3];
	if (this->chromosomes[index1]->getFitness()
			< this->chromosomes[index2]->getFitness()) {
		parent1 = this->chromosomes[index1];
	} else {
		parent1 = this->chromosomes[index2];
	}

	if (this->chromosomes[index3]->getFitness()
			< this->chromosomes[index4]->getFitness()) {
		parent2 = this->chromosomes[index3];
	} else {
		parent2 = this->chromosomes[index4];
	}
	result.push_back(parent1);
	result.push_back(parent2);
	return (result);
}

void Population::dump() {
	/*
	cout << "Population:" << endl;
	cout << "Size: " << popsize << " Chsize: " << chsize << " Crossprob: "
			<< crossprob << endl;
	for (int i = 0; i < this->popsize; i++) {
		this->chromosomes[i]->dump();
	}
	*/
}

void Population::setFitness(int index, double d) {
	this->chromosomes[index]->setFitness(d);
}

double Population::getFitness(int index) {
	return (this->chromosomes[index]->getFitness());
}

void Population::setChromosome(int index, Chromosome *c) {
	this->chromosomes[index] = c;
}

Chromosome *Population::getChromosome(int index) {
	return (this->chromosomes[index]);
}

void Population::sortPopulation() {
	std::sort(this->chromosomes.begin(), this->chromosomes.end(),
			chromosome_sorter);
}

void Population::TournamentSelection(Function ForecastFunction, double forecastprob, int MinimumForecastLength) {
	parents.clear();
	offspring.clear();
	tempPopulation.clear();

	this->sortPopulation();
	for (int i = 0; i < this->elitism; i++) {
		tempPopulation.push_back(this->chromosomes[i]);
	}

	/*
	 *
	 * Apply forecasting operator here
	 *
	 */
	NumericMatrix nm = this->getChromosome(0)->getFatherMatrix();
	NumericVector genes = this->getChromosome(0)->getGenesAsDoubleRcppVector();
	NumericVector result = ForecastFunction(genes, nm, MinimumForecastLength, forecastprob, 5);
	Chromosome *newch = this->getChromosome(0)->clone();
	for (int i=0;i<newch->getLength();i++){
	    newch->setGene(i, result[i]);
	}
	tempPopulation.push_back(newch);
	/*
	 * 
	 * End of forecasting operator
	 */
	while (true) {
		parents = this->performTournament();
		offspring = parents[0]->CrossOver(parents[1], this->crossprob);
		offspring[0]->Mutate(this->mutationprob);
		offspring[1]->Mutate(this->mutationprob);
		if (tempPopulation.size() < chromosomes.size()) {
			tempPopulation.push_back(offspring[0]);
		} else {
			break;
		}
		if (tempPopulation.size() < chromosomes.size()) {
			tempPopulation.push_back(offspring[1]);
		} else {
			break;
		}
	}
	std::copy(tempPopulation.begin(), tempPopulation.end(),
			chromosomes.begin());
}



NumericMatrix Population::getPopulationAsRcppMatrix(){
    NumericMatrix nm(this->popsize, this->chsize);
    for (int i=0;i<this->popsize;i++){
        for (int j=0;j<this->chsize;j++){
            nm(i,j) = this->chromosomes[i]->getGene(j);
        }
    }
    return(nm);
}
