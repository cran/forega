#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "Chromosome.h"

using namespace std;
using namespace Rcpp;

vector<Chromosome*> Chromosome::chromosomesEverCreated;

Chromosome::Chromosome(int len) {
	this->len = len;
	this->fitness = 0;
	this->father = NULL;
	this->mother = NULL;
	this->genes = (double*) malloc(sizeof(double) * len);
	this->maxv = NULL;
	this->minv = NULL;
	chromosomesEverCreated.push_back(this);
}

Chromosome::~Chromosome() {
	free(genes);
    //cout << "Deleted Chromosome " << (unsigned long)this << endl;
}

int Chromosome::getLength(){
    return(this->len);
}

void Chromosome::setFitness(double f) {
	this->fitness = f;
}

double Chromosome::getFitness() {
	return (this->fitness);
}

void Chromosome::dump() {
	/*
	cout << "F(" << this->fitness << "): ";
	for (int i = 0; i < this->len; i++) {
		cout << this->genes[i] << ", ";
	}
	cout << endl;
    */
}

void Chromosome::setFather(Chromosome *c) {
	this->father = c;
}

void Chromosome::setMother(Chromosome *c) {
	this->mother = c;
}

Chromosome *Chromosome::getFather() {
	return (this->father);
}

Chromosome *Chromosome::getMother() {
	return (this->mother);
}

void Chromosome::setMaxValues(double m) {
	if (this->maxv == NULL) {
		this->maxv = (double*) malloc(this->len * sizeof(double));
	}

	for (int i = 0; i < this->len; i++) {
		this->maxv[i] = m;
	}
}

void Chromosome::setMinValues(double m) {
	if (this->minv == NULL) {
		this->minv = (double*) malloc(this->len * sizeof(double));
	}

	for (int i = 0; i < this->len; i++) {
		this->minv[i] = m;
	}
}

void Chromosome::setMaxValues(double* m) {
	this->maxv = m;
}

void Chromosome::setMinValues(double* m) {
	this->minv = m;
}

void Chromosome::randomize() {
	for (int i = 0; i < this->len; i++) {
		//this->genes[i] = minv[i] + ((double) rand() / (double) (RAND_MAX)) * (maxv[i] - minv[i]);
		this->genes[i] = minv[i] + R::runif(0,1) * (maxv[i] - minv[i]);
	}
}

void Chromosome::Mutate(double mutationProb) {
    /*
	for (int i = 0; i < this->len; i++) {
		if ((((double) rand()) / ((double) (RAND_MAX))) < mutationProb) {
			this->genes[i] = minv[i]
					+ (((double) rand()) / ((double) (RAND_MAX)))
							* (maxv[i] - minv[i]);
		}
	}
    */
	for (int i = 0; i < this->len; i++) {
		if (R::runif(0,1) < mutationProb) {
			this->genes[i] = minv[i] + R::runif(0,1) * (maxv[i] - minv[i]);
		}
	}
}

void Chromosome::setGenes(double *d) {
	if (this->genes == NULL) {
		this->genes = (double*) malloc(this->len * sizeof(double));
	}
	for (int i = 0; i < this->len; i++) {
		this->genes[i] = d[i];
	}
}

double *Chromosome::getGenes() {
	return (this->genes);
}

Chromosome *Chromosome::clone() {
	Chromosome *c = new Chromosome(this->len);
	//c->setFather(this->father);
	//c->setMother(this->mother);
	c->setMaxValues(this->getMaxValues());
	c->setMinValues(this->getMinValues());
	c->setFitness(this->getFitness());
	for (int i = 0; i < this->len; i++) {
		c->setGene(i, this->genes[i]);
	}
	return (c);
}

double *Chromosome::getMaxValues() {
	return (this->maxv);
}

double *Chromosome::getMinValues() {
	return (this->minv);
}

void Chromosome::setGene(int i, double v) {
	this->genes[i] = v;
}

double Chromosome::getGene(int i) {
	return (this->genes[i]);
}

vector<Chromosome*> Chromosome::CrossOver(Chromosome *other, double crossProb) {
	vector<Chromosome*> v;
	Chromosome *offspring1 = this->clone();
	Chromosome *offspring2 = other->clone();
	/*
	for (int i = 0; i < this->len; i++) {
		if (((double) rand() / (RAND_MAX)) < crossProb) {
			double beta = (((double) rand()) / ((double) (RAND_MAX)));
			offspring1->setGene(i,
					this->genes[i]
							- beta * (this->genes[i] - other->getGene(i)));
			offspring2->setGene(i,
					other->getGene(i)
							+ beta * (this->genes[i] - other->getGene(i)));
		}
	}
    */
	for (int i = 0; i < this->len; i++) {
		if (R::runif(0,1) < crossProb) {
			double beta = R::runif(0,1);
			offspring1->setGene(i,
					this->genes[i]
							- beta * (this->genes[i] - other->getGene(i)));
			offspring2->setGene(i,
					other->getGene(i)
							+ beta * (this->genes[i] - other->getGene(i)));
		}
	}
	offspring1->setFather(this);
	offspring1->setMother(other);
	offspring2->setFather(this);
	offspring2->setMother(other);
	v.push_back(offspring1);
	v.push_back(offspring2);
	return (v);
}

void Chromosome::getFatherList(vector<Chromosome*> &v) {
	Chromosome *c = this;
	while (true) {
		v.push_back(c);
		if (c->father == NULL) {
			std::reverse(v.begin(), v.end());
			return;
		}
		c = c->father;
	}
}

void Chromosome::getMotherList(vector<Chromosome*> &v) {
	Chromosome *c = this;
	while (true) {
		v.push_back(c);
		if (c->mother == NULL) {
			std::reverse(v.begin(), v.end());
			return;
		}
		c = c->mother;
	}
}

NumericMatrix Chromosome::getFatherMatrix(){
	vector<Chromosome*> fathers;
	this->getFatherList(fathers);
	NumericMatrix m(fathers.size(), this->len);
	for (unsigned int i=0;i<fathers.size();i++){
		for (int j=0;j<len;j++){
			m(i,j) = fathers[i]->getGene(j);
		}
	}
	return(m);
}


NumericMatrix Chromosome::getMotherMatrix(){
    vector<Chromosome*> mothers;
    this->getMotherList(mothers);
    NumericMatrix m(mothers.size(), this->len);
    for (unsigned int i=0;i<mothers.size();i++){
        for (int j=0;j<len;j++){
            m(i,j) = mothers[i]->getGene(j);
        }
    }
    return(m);
}



DoubleVector Chromosome::getGenesAsDoubleRcppVector(){
    DoubleVector dv(this->len);
    for (int i=0;i<this->len;i++){
        dv[i] = this->genes[i];
    }
    return(dv);
}
