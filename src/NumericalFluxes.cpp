#include <NumericalFluxes.h>

NumericalFluxes::NumericalFluxes() {

    this->normal        = nullptr;
    this->flux_l        = nullptr;
    this->flux_r        = nullptr;
    this->velocity_l    = nullptr;
    this->velocity_r    = nullptr;
    this->numericalFlux = nullptr;

    this->height_l = 0.0;
    this->height_r = 0.0;
}

NumericalFluxes::~NumericalFluxes() {

    delete [] this->normal;
    delete [] this->flux_l;
    delete [] this->flux_r;
    delete [] this->velocity_l;
    delete [] this->velocity_r;
    delete [] this->numericalFlux;
}

RoeFlux::RoeFlux(Config *config, Solution &sol) {

    this->flux_l        = new double [sol.nVar];
    this->flux_r        = new double [sol.nVar];
    this->diffW         = new double [sol.nVar];
    this->numericalFlux = new double [sol.nVar];
    this->lambda        = new double [sol.nVar];
    this->lambda_l      = new double [sol.nVar];
    this->lambda_r      = new double [sol.nVar];
    this->diffW         = new double [sol.nVar];
    this->normal        = new double [2];
    this->velocity_l    = new double [2];
    this->velocity_r    = new double [2];
    this->roeVelocity   = new double [2];

    this->matrixL = new double *[sol.nVar];
    this->matrixR = new double *[sol.nVar];
    for (int i = 0; i < sol.nVar; ++i) {
        this->matrixL[i] = new double [sol.nVar];
        this->matrixR[i] = new double [sol.nVar];
    }

    this->height_l  = 0.0;
    this->height_r  = 0.0;
    this->roeHeight = 0.0;

    this->entropyFix = config->getEntropyFix();
}

RoeFlux::~RoeFlux() {

    delete [] this->lambda_l;
    delete [] this->lambda_r;
    delete [] this->lambda;
    delete [] this->roeVelocity;
    delete [] this->diffW;

    for (int i = 0; i < 3; ++i) {
        delete [] this->matrixL[i];
        delete [] this->matrixR[i];
    }
    delete [] this->matrixL;
    delete [] this->matrixR;
}

void RoeFlux::computeFlux(Grid *mesh, Solution *sol, int &interfaceID, double &dt) {

    double wLeft[3] = {0.0, 0.0, 0.0}, wRight[3] = {0.0, 0.0, 0.0};
    double numRoe, denRoe, numEF, denEF, matrixRALElement;
    int leftID, rightID;

    mesh->getAdjacentNodesToInterface(interfaceID, leftID, rightID);
    mesh->getNormalToDualGridEdge(interfaceID, this->normal);

    wLeft[0] = sol->conservative1(leftID);      wRight[0] = sol->conservative1(rightID);
    wLeft[1] = sol->conservative2(leftID);      wRight[1] = sol->conservative2(rightID);
    wLeft[2] = sol->conservative3(leftID);      wRight[2] = sol->conservative3(rightID);

    this->height_l = wLeft[0];
    this->height_r = wRight[0];

    for (int i = 0; i < 2; ++i) {
        this->velocity_l[i] = wLeft[i + 1]/wLeft[0];
        this->velocity_r[i] = wRight[i + 1]/wRight[0];
    }

    for (int i = 0; i < sol->nVar; ++i) {
        this->diffW[i] = wRight[i] - wLeft[i];
    }

    this->roeHeight = (this->height_l + this->height_r)/2;

    for (int i = 0; i < 2; ++i) {
        numRoe = sqrt(this->height_l)*this->velocity_l[i] + sqrt(this->height_r)*this->velocity_r[i];
        denRoe = sqrt(this->height_l) + sqrt(this->height_r);
        this->roeVelocity[i] = numRoe/denRoe;
    }

    if (this->entropyFix) {

        ModelEquations::setProjEigenvalues(this->height_l, this->velocity_l, this->normal, this->lambda_l);
        ModelEquations::setProjEigenvalues(this->height_r, this->velocity_r, this->normal, this->lambda_r);
        ModelEquations::setProjEigenvalues(this->roeHeight, this->roeVelocity, this->normal, this->lambda);

        for (int i = 0; i < sol->nVar; ++i) {
            if (this->lambda_l[i] < 0 && this->lambda_r[i] > 0) {
                numEF = (this->lambda_l[i] + this->lambda_r[i])*this->lambda[i] - 2*this->lambda_l[i]*this->lambda_r[i];
                denEF = this->lambda_r[i] - this->lambda_l[i];
                this->lambda[i] = numEF/denEF;
            } else {
                this->lambda[i] = abs(this->lambda[i]);
            }
        }
    } else {

        ModelEquations::setProjEigenvalues(this->roeHeight, this->roeVelocity, this->normal, this->lambda);

        for (int i = 0; i < sol->nVar; ++i) {
            this->lambda[i] = abs(this->lambda[i]);
        }
    }

    ModelEquations::setProjFlux(this->height_l, this->velocity_l, this->normal, this->flux_l);
    ModelEquations::setProjFlux(this->height_r, this->velocity_r, this->normal, this->flux_r);

    ModelEquations::setProjLeftMatrix(this->roeHeight, this->roeVelocity, this->normal, this->matrixL);
    ModelEquations::setProjRightMatrix(this->roeHeight, this->roeVelocity, this->normal, this->matrixR);

    for (int i = 0; i < sol->nVar; ++i) {
        this->numericalFlux[i] = 0.5*(this->flux_l[i] + this->flux_r[i]);
        for (int j = 0; j < sol->nVar; ++j) {
            matrixRALElement = 0.0;
            for (int k = 0; k < sol->nVar; ++k) {
                matrixRALElement += this->matrixR[i][k]*this->lambda[k]*this->matrixL[k][j];
            }
            this->numericalFlux[i] -= 0.5*matrixRALElement*this->diffW[j];
        }
    }
}

LaxWendroffFlux::LaxWendroffFlux(Solution &sol) {

    this->flux_l         = new double [sol.nVar];
    this->flux_r         = new double [sol.nVar];
    this->flux_inter     = new double [sol.nVar];
    this->numericalFlux  = new double [sol.nVar];
    this->normal         = new double [2];
    this->velocity_l     = new double [2];
    this->velocity_r     = new double [2];
    this->velocity_inter = new double [2];

    this->projJacobian = new double *[sol.nVar];
    for (int i = 0; i < sol.nVar; ++i) {
        this->projJacobian[i] = new double [sol.nVar];
    }

    this->height_l     = 0.0;
    this->height_r     = 0.0;
    this->height_inter = 0.0;
}

LaxWendroffFlux::~LaxWendroffFlux() {

    delete [] this->flux_inter;
    delete [] this->velocity_inter;

    for (int i = 0; i < 3; ++i) {
        delete [] this->projJacobian[i];
    }
    delete [] this->projJacobian;
}

void LaxWendroffFlux::computeFlux(Grid *mesh, Solution *sol, int &interfaceID, double &dt) {

    double wLeft[3] = {0.0, 0.0, 0.0}, wRight[3] = {0.0, 0.0, 0.0};
    double size;
    int leftID, rightID;

    mesh->getAdjacentNodesToInterface(interfaceID, leftID, rightID);
    mesh->getNormalToDualGridEdge(interfaceID, this->normal);
    size = mesh->getDualEdgeMeasure(interfaceID);

    wLeft[0] = sol->conservative1(leftID);      wRight[0] = sol->conservative1(rightID);
    wLeft[1] = sol->conservative2(leftID);      wRight[1] = sol->conservative2(rightID);
    wLeft[2] = sol->conservative3(leftID);      wRight[2] = sol->conservative3(rightID);

    this->height_l     = wLeft[0];
    this->height_r     = wRight[0];
    this->height_inter = 0.5*(this->height_l + this->height_r);

    for (int i = 0; i < 2; ++i) {
        this->velocity_l[i]     = wLeft[i + 1]/wLeft[0];
        this->velocity_r[i]     = wRight[i + 1]/wRight[0];
        this->velocity_inter[i] = 0.5*(this->velocity_l[i] + this->velocity_r[i]);
    }

    ModelEquations::setProjFlux(this->height_l, this->velocity_l, this->normal, this->flux_l);
    ModelEquations::setProjFlux(this->height_r, this->velocity_r, this->normal, this->flux_r);

    ModelEquations::setProjJacobian(this->height_inter, this->velocity_inter, this->normal, this->projJacobian);

    for (int i = 0; i < sol->nVar; ++i) {
        this->numericalFlux[i] = 0.5*(this->flux_l[i] + this->flux_r[i]);
        for (int j = 0; j < sol->nVar; ++j) {
            this->numericalFlux[i] -= 0.5*(dt/size)*this->projJacobian[i][j]*(this->flux_r[j] - this->flux_l[j]);
        }
    }
}

HighResFlux::HighResFlux(Config *config, Solution &sol) : RoeFlux(config, sol) {

    this->fluxLimiter = config->getFluxLimiter();

    this->diffW_l      = new double [sol.nVar];
    this->diffW_r      = new double [sol.nVar];
    this->lambdaSigned = new double [sol.nVar];
    this->diffV        = new double [sol.nVar];
    this->diffV_l      = new double [sol.nVar];
    this->diffV_r      = new double [sol.nVar];
    this->diffV_up     = new double [sol.nVar];

    this->phi = 0.0;
}

HighResFlux::~HighResFlux() {

    delete [] this->diffW_l;
    delete [] this->diffW_r;
    delete [] this->lambdaSigned;
    delete [] this->diffV;
    delete [] this->diffV_l;
    delete [] this->diffV_r;
    delete [] this->diffV_up;
}

void HighResFlux::setFluxLimiter(double upDiff, double cenDiff) {

    double a, b;

    switch (this->fluxLimiter) {
        case NONE:
            this->phi = 0;
            break;
        case SUPERBEE:
            a = min(abs(cenDiff), 2*abs(upDiff));
            b = min(2*abs(cenDiff), abs(upDiff));
            this->phi = 0.5*(sgn(cenDiff) + sgn(upDiff))*max(a, b);
            break;
        case VAN_LEER:
            this->phi = (cenDiff*abs(upDiff) + abs(cenDiff)*upDiff)/(abs(cenDiff) + abs(upDiff) + 1.0e-8);
            break;
        case MINMOD:
            this->phi = 0.5*(sgn(cenDiff) + sgn(upDiff))*min(abs(cenDiff), abs(upDiff));
            break;
    }
}

void HighResFlux::computeFlux(Grid *mesh, Solution *sol, int &interfaceID, double &dt) {

    bool isFound = false;
    double wLeft[3] = {0.0, 0.0, 0.0}, wRight[3] = {0.0, 0.0, 0.0};
    double wLeftUp[3] = {0.0, 0.0, 0.0}, wRightUp[3] = {0.0, 0.0, 0.0};
    double size;
    int leftID, rightID;

    RoeFlux::computeFlux(mesh, sol, interfaceID, dt);

    mesh->getAdjacentNodesToInterface(interfaceID, leftID, rightID);
    size = mesh->getDualEdgeMeasure(interfaceID);

    wLeft[0] = sol->conservative1(leftID);      wRight[0] = sol->conservative1(rightID);
    wLeft[1] = sol->conservative2(leftID);      wRight[1] = sol->conservative2(rightID);
    wLeft[2] = sol->conservative3(leftID);      wRight[2] = sol->conservative3(rightID);

    for (int i = 0; i < mesh->nBoundaries && !isFound; ++i) {
        for (int j = 0; j < mesh->nodesIDOfBoundaryID[i].size() && !isFound; ++j) {
            if ((leftID == mesh->nodesIDOfBoundaryID[0][j] && this->normal[0] == 0) ||
                (leftID == mesh->nodesIDOfBoundaryID[3][j] && this->normal[0] == 1)) {
                isFound = true;
                for (int k = 0; k < sol->nVar; ++k) {
                    this->diffW_l[k] = this->diffW[k];
                }
            } else if ((rightID == mesh->nodesIDOfBoundaryID[1][j] && this->normal[0] == 1) ||
                       (rightID == mesh->nodesIDOfBoundaryID[2][j] && this->normal[0] == 0)) {
                isFound = true;
                for (int k = 0; k < sol->nVar; ++k) {
                    this->diffW_r[k] = this->diffW[k];
                }
            }
        }
    }

    if (!isFound) {
        if (this->normal[0] == 1) {
            wLeftUp[0] = sol->conservative1(leftID - 1);      wRightUp[0] = sol->conservative1(rightID + 1);
            wLeftUp[1] = sol->conservative2(leftID - 1);      wRightUp[1] = sol->conservative2(rightID + 1);
            wLeftUp[2] = sol->conservative3(leftID - 1);      wRightUp[2] = sol->conservative3(rightID + 1);
        } else {
            wLeftUp[0] = sol->conservative1(leftID - mesh->nNodesX);      wRightUp[0] = sol->conservative1(rightID + mesh->nNodesX);
            wLeftUp[1] = sol->conservative2(leftID - mesh->nNodesX);      wRightUp[1] = sol->conservative2(rightID + mesh->nNodesX);
            wLeftUp[2] = sol->conservative3(leftID - mesh->nNodesX);      wRightUp[2] = sol->conservative3(rightID + mesh->nNodesX);
        }
        for (int i = 0; i < sol->nVar; ++i) {
            this->diffW_r[i] = wRightUp[i] - wRight[i];
            this->diffW_l[i] = wLeft[i] - wLeftUp[i];
        }
    }

    ModelEquations::setProjEigenvalues(this->roeHeight, this->roeVelocity, this->normal, this->lambdaSigned);

    for (int i = 0; i < sol->nVar; ++i) {
        this->diffV_l[i] = 0.0;
        this->diffV_r[i] = 0.0;
        this->diffV[i]   = 0.0;
        for (int j = 0; j < sol->nVar; ++j) {
            this->diffV_l[i] += this->matrixL[i][j]*this->diffW_l[j];
            this->diffV_r[i] += this->matrixL[i][j]*this->diffW_r[j];
            this->diffV[i]   += this->matrixL[i][j]*this->diffW[j];
        }
    }

    for (int i = 0; i < sol->nVar; ++i) {
        for (int j = 0; j < sol->nVar; ++j) {
            this->diffV_up[j] = 0.5*(this->diffV_l[j] + this->diffV_r[j]) + 0.5*(this->diffV_l[j] - this->diffV_r[j])*sgn(this->lambdaSigned[j]);
            this->setFluxLimiter(this->diffV_up[j], this->diffV[j]);
            this->numericalFlux[i] += 0.5*this->matrixR[i][j]*(this->lambda[j] - dt/size*pow(this->lambdaSigned[j], 2))*this->phi;
        }
    }
}