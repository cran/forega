/*
 * Chromosome.h
 *
 *  Created on: May 27, 2015
 *      Author: hako
 */

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <vector>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

class Chromosome {
private:
	double* genes;
	int len;
	Chromosome *father;
	Chromosome *mother;
	double* minv;
	double* maxv;

public:
	double fitness;
	Chromosome(int len);
	~Chromosome();
	int getLength();
	void setGene(int i, double v);
	double getGene(int i);
	void setMinValues(double *minv);
	void setMaxValues(double *maxv);
	void setMinValues(double m);
	void setMaxValues(double m);
	double *getMinValues();
	double *getMaxValues();
	void setFitness(double f);
	double getFitness();
	void setMother(Chromosome *c);
	void setFather(Chromosome *c);
	Chromosome *getMother();
	Chromosome *getFather();
	void Mutate(double mutationProb);
	double* getGenes();
	void setGenes(double *d);
	vector<Chromosome*> CrossOver(Chromosome *other, double crossProb);
	void dump();
	void randomize();
	Chromosome *clone();
	void getFatherList(vector<Chromosome*> &v);
	void getMotherList(vector<Chromosome*> &v);
	NumericMatrix getFatherMatrix();
	NumericMatrix getMotherMatrix();
	DoubleVector getGenesAsDoubleRcppVector();
	
	static vector<Chromosome*> chromosomesEverCreated;
};


#endif /* CHROMOSOME_H_ */
