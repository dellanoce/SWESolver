#include "Solution.h"

Solution::Solution() {

    this->data = nullptr;

    this->nVar   = 0;
    this->nNodes = 0;
}

Solution::Solution(Grid *meshPtr) {

    this->initSolutionOnMesh(meshPtr);
}

Solution::~Solution() {

    delete [] this->data;
}

void Solution::initSolutionOnMesh(Grid *meshPtr) {

    const int scalarVars = 1;
    const int vectorVars = 1;

    this->nNodes = meshPtr->getNNodes();
    this->nVar   = scalarVars + 2*vectorVars;

    this->data = new double[this->nVar*this->nNodes];
}

double &Solution::conservative1(const int &nodeID) {

    return this->data[nodeID];
}

double &Solution::conservative2(const int &nodeID) {

    return this->data[nodeID + this->nNodes];
}

double &Solution::conservative3(const int &nodeID) {

    return this->data[nodeID + 2*this->nNodes];
}

double Solution::primitive1(const int &nodeID) {

    return this->conservative1(nodeID);
}

double Solution::primitive2(const int &nodeID) {

    return this->conservative2(nodeID)/this->conservative1(nodeID);
}

double Solution::primitive3(const int &nodeID) {

    return this->conservative3(nodeID)/this->conservative1(nodeID);
}