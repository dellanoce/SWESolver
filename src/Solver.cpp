#include <Solver.h>

Solver::Solver() {

    this->solPtr     = nullptr;
    this->residual   = nullptr;
    this->consVector = nullptr;

    this->totalTime = 0.0;
    this->dt        = 0.0;
    this->cfl       = 0.0;

    this->eulerBoundaries.clear();
}

Solver::Solver(Config *config, Grid *mesh) {

    this->solPtr = new Solution(mesh);

    switch (config->getNumericalFluxMethod()) {
        case ROE:
            this->fluxPtr = make_shared<RoeFlux>(config, *this->solPtr);
            break;
        case LAX_WENDROFF:
            this->fluxPtr = make_shared<LaxWendroffFlux>(*this->solPtr);
            break;
        case HIGH_RES:
            this->fluxPtr = make_shared<HighResFlux>(config, *this->solPtr);
            break;
    }

    this->residual   = new double [mesh->nNodes*this->solPtr->nVar];
    this->consVector = new double [this->solPtr->nVar];

    this->totalTime = config->getTotalTime();
    this->cfl       = config->getCFLNumber();
    this->dt        = SWE_UNDEFINED;

    this->eulerBoundaries = config->getMarkerEuler();
}

Solver::~Solver() {

    delete [] this->residual;
    delete [] this->consVector;
}

void Solver::setTimeStep(Grid *mesh) {

    double hx, hy, dtMin, height;
    double *lambda, *velocity, *normal;

    dtMin    = 1000000;
    lambda   = new double [this->solPtr->nVar];
    velocity = new double [2];
    normal   = new double [2];

    for (int i = 0; i < mesh->nNodes; ++i) {

        height      = this->solPtr->primitive1(i);
        velocity[0] = this->solPtr->primitive2(i);
        velocity[1] = this->solPtr->primitive3(i);

        mesh->getDualElementMeasures(i, hx, hy);

        for (int j = 0; j < 2; ++j) {

            int offset  = (j == 0) ? 1 : -1;
            double size = (j == 0) ? hx : hy;

            normal[0] = j + offset;
            normal[1] = j;

            ModelEquations::setProjEigenvalues(height, velocity, normal, lambda);

            for (int k = 0; k < this->solPtr->nVar; ++k) {
                this->dt = this->cfl/max(abs(lambda[k]/size), 1e-9);

                if (this->dt < dtMin)
                    dtMin = this->dt;
            }
        }
    }

    this->dt = dtMin;

    delete [] lambda;
    delete [] velocity;
    delete [] normal;
}

void Solver::setInitialSolution(Grid *mesh, Solution &initialSol) {

    double a, b, radius;

    for (int i = 0; i < mesh->nNodes; ++i) {

        a = pow((mesh->nodes[i][0] - 0.2)/0.04, 2);
        b = pow((mesh->nodes[i][1] - 0.2)/0.04, 2);

        initialSol.conservative1(i) = 1 + 0.3*exp(-(a + b));
        initialSol.conservative2(i) = 0;
        initialSol.conservative3(i) = 0;
    }

//    for (int i = 0; i < mesh->nNodes; ++i) {
//
//        radius = sqrt(pow(mesh->nodes[i][0], 2) + pow(mesh->nodes[i][1], 2));
//
//        if (radius <= 100)
//            initialSol.conservative1(i) = 4.0;
//        else
//            initialSol.conservative1(i) = 1.0;
//
//        initialSol.conservative2(i) = 0.0;
//        initialSol.conservative3(i) = 0.0;
//    }

//    for (int i = 0; i < mesh->nNodes; ++i) {
//
//        radius = sqrt(pow(mesh->nodes[i][0] - 100, 2) + pow(mesh->nodes[i][1] - 100, 2));
//
//        if (radius <= 50)
//            initialSol.conservative1(i) = 10.0;
//        else
//            initialSol.conservative1(i) = 1.0;
//
//        initialSol.conservative2(i) = 0.0;
//        initialSol.conservative3(i) = 0.0;
//    }

    this->solPtr = &initialSol;
}

void Solver::runGodnuovSolver(Grid *mesh) {

    double time = 0, length;
    int leftID, rightID;

    while (time < this->totalTime) {

        this->setTimeStep(mesh);
        if (time + this->dt > this->totalTime)
            this->dt = this->totalTime - time;

        cout << " Computing solution at time: " << time + this->dt << " seconds..." << endl;
        cout << " ----------------------------------------------------------------" << endl << endl;

        this->initializeResidual(mesh);

        for (int i = 0; i < mesh->nInterfaces; ++i) {

            mesh->getAdjacentNodesToInterface(i, leftID, rightID);
            length = mesh->getDualEdgeMeasure(i);

            this->fluxPtr->computeFlux(mesh, this->solPtr, i, this->dt);

            this->addToResidual(leftID, length, this->fluxPtr->numericalFlux);
            this->subtractToResidual(rightID, length, this->fluxPtr->numericalFlux);
        }

        for (int i = 0; i < mesh->nBoundaries; ++i) {
            this->setBoundaryConditions(mesh, i);
        }

        this->updateSolution(mesh);

        OutputManager::printSolutionVTK(mesh, this->solPtr, time + this->dt);
        OutputManager::printSolutionCsv(mesh, this->solPtr, time + this->dt);

        time += this->dt;
    }
}

void Solver::updateSolution(Grid *mesh) const {

    for (int i = 0; i < mesh->nNodes; ++i) {

        this->consVector[0] = this->solPtr->conservative1(i);
        this->consVector[1] = this->solPtr->conservative2(i);
        this->consVector[2] = this->solPtr->conservative3(i);

        for (int j = 0; j < this->solPtr->nVar; ++j) {
            this->consVector[j] -= (this->dt/mesh->getDualElementArea(i))*this->residual[i*this->solPtr->nVar + j];
        }

        this->solPtr->conservative1(i) = this->consVector[0];
        this->solPtr->conservative2(i) = this->consVector[1];
        this->solPtr->conservative3(i) = this->consVector[2];
    }
}

void Solver::initializeResidual(Grid *mesh) const {

    for (int i = 0; i < mesh->nNodes; ++i) {
        for (int j = 0; j < this->solPtr->nVar; ++j) {
            this->residual[i*this->solPtr->nVar + j] = 0.0;
        }
    }
}

void Solver::addToResidual(const int &nodeID, const double length, const double *numericalFlux) const {

    for (int i = 0; i < this->solPtr->nVar; ++i) {
        this->residual[nodeID*this->solPtr->nVar + i] += numericalFlux[i]*length;
    }
}

void Solver::subtractToResidual(const int &nodeID, const double length, const double *numericalFlux) const {

    for (int i = 0; i < this->solPtr->nVar; ++i) {
        this->residual[nodeID*this->solPtr->nVar + i] -= numericalFlux[i]*length;
    }
}

void Solver::setBoundaryConditions(Grid *mesh, const int &boundaryID) const {

    bool isFoundEuler;
    int boundaryNodeID;
    double boundaryHeight, barHeight, h, hx, hy;
    double *boundaryVelocity, *boundaryNormal, *boundaryLambda, *boundaryFlux, *barVelocity;
    double **boundaryMatrixL, **boundaryMatrixR;
    double wEuler[3] = {0.0, 0.0, 0.0}, wBoundary[3] = {0.0, 0.0, 0.0};
    double deltaW[3] = {0.0, 0.0, 0.0}, deltaV[3] = {0.0, 0.0, 0.0}, deltaWBar[3] = {0.0, 0.0, 0.0};
    double boundaryState[3] = {0.0, 0.0, 0.0};
    Solution* eulerSol;

    boundaryVelocity = new double [2];
    boundaryNormal   = new double [2];
    barVelocity      = new double [2];
    boundaryLambda   = new double [this->solPtr->nVar];
    boundaryFlux     = new double [this->solPtr->nVar];
    boundaryMatrixL  = new double *[this->solPtr->nVar];
    boundaryMatrixR  = new double *[this->solPtr->nVar];
    for (int i = 0; i < this->solPtr->nVar; ++i) {
        boundaryMatrixL[i] = new double [this->solPtr->nVar];
        boundaryMatrixR[i] = new double [this->solPtr->nVar];
    }
    eulerSol = this->solPtr;

    isFoundEuler = swe::existsInVectorString(this->eulerBoundaries, to_string(boundaryID));

    for (int i = 0; i < mesh->nodesIDOfBoundaryID[boundaryID].size(); ++i) {

        boundaryNodeID = mesh->nodesIDOfBoundaryID[boundaryID][i];

        int offset = (boundaryID == 0 || boundaryID == 3) ? -1 : 1;

        if (boundaryID == 0 || boundaryID == 2) {
            boundaryNormal[0] = 0;
            boundaryNormal[1] = offset;
        } else if (boundaryID == 1 || boundaryID == 3) {
            boundaryNormal[0] = offset;
            boundaryNormal[1] = 0;
        }

        if (isFoundEuler) {
            if (boundaryNormal[0] == 0)
                eulerSol->conservative3(boundaryNodeID) = 0;
            else if (boundaryNormal[1] == 0)
                eulerSol->conservative2(boundaryNodeID) = 0;
        }

        wEuler[0] = eulerSol->conservative1(boundaryNodeID);      wBoundary[0] = this->solPtr->conservative1(boundaryNodeID);
        wEuler[1] = eulerSol->conservative2(boundaryNodeID);      wBoundary[1] = this->solPtr->conservative2(boundaryNodeID);
        wEuler[2] = eulerSol->conservative3(boundaryNodeID);      wBoundary[2] = this->solPtr->conservative3(boundaryNodeID);

        boundaryHeight = wBoundary[0];
        for (int j = 0; j < 2; ++j) {
            boundaryVelocity[j] = wBoundary[j + 1]/wBoundary[0];
        }

        ModelEquations::setProjEigenvalues(boundaryHeight, boundaryVelocity, boundaryNormal, boundaryLambda);
        ModelEquations::setProjLeftMatrix(boundaryHeight, boundaryVelocity, boundaryNormal, boundaryMatrixL);
        ModelEquations::setProjRightMatrix(boundaryHeight, boundaryVelocity, boundaryNormal, boundaryMatrixR);

        for (int j = 0; j < this->solPtr->nVar; ++j) {

            deltaV[j]    = 0.0;
            deltaWBar[j] = 0.0;
            deltaW[j]    = wEuler[j] - wBoundary[j];

            for (int k = 0; k < this->solPtr->nVar; ++k) {
                deltaV[j] += boundaryMatrixL[j][k]*deltaW[k];
            }

            if ((boundaryLambda[j] > 0 && (boundaryID == 1 || boundaryID == 2))
                || boundaryLambda[j] < 0 && (boundaryID == 0 || boundaryID == 3))
                deltaV[j] = 0;

            for (int k = 0; k < this->solPtr->nVar; ++k) {
                deltaWBar[j] += boundaryMatrixR[j][k]*deltaV[k];
            }

            boundaryState[j] = wBoundary[j] + deltaWBar[j];
        } 

        barHeight = boundaryState[0];
        for (int j = 0; j < 2; ++j) {
            barVelocity[j] = boundaryState[j + 1]/boundaryState[0];
        }

        ModelEquations::setProjFlux(barHeight, barVelocity, boundaryNormal, boundaryFlux);

        mesh->getDualElementMeasures(boundaryNodeID, hx, hy);
        h = (boundaryID == 0 || boundaryID == 2) ? hx : hy;

        this->addToResidual(boundaryNodeID, h, boundaryFlux);
    }

    delete [] boundaryVelocity;
    delete [] boundaryNormal;
    delete [] barVelocity;
    delete [] boundaryLambda;
    delete [] boundaryFlux;

    for (int i = 0; i < this->solPtr->nVar; ++i) {
        delete [] boundaryMatrixL[i];
        delete [] boundaryMatrixR[i];
    }
    delete [] boundaryMatrixL;
    delete [] boundaryMatrixR;
}