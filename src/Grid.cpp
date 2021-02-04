#include "Grid.h"

Grid::Grid() {

    this->x0 = 0.0;
    this->y0 = 0.0;
    this->xf = 0.0;
    this->yf = 0.0;

    this->nNodesX     = 0;
    this->nNodesY     = 0;
    this->nNodes      = 0;
    this->nInterfaces = 0;
    this->nBoundaries = 0;

    this->nodes = nullptr;

    this->nodesIDOfBoundaryID.clear();
}

Grid::Grid(Config *config) {

    this->x0      = config->getXInitialLimit();
    this->xf      = config->getXFinalLimit();
    this->y0      = config->getYInitialLimit();
    this->yf      = config->getYFinalLimit();
    this->nNodesX = config->getNNodesX();
    this->nNodesY = config->getNNodesY();

    this->nNodes      = this->getNNodes();
    this->nInterfaces = this->getNInterfaces();
    this->nBoundaries = 4;

    this->nodes = new double* [this->nNodes];
    for (int i = 0; i < this->nNodes; ++i) {
        this->nodes[i] = new double [2];
    }

    this->buildUniformGrid();
}

Grid::~Grid() {

    for (int i = 0; i < this->nNodes; ++i) {
        delete [] this->nodes[i];
    }
    delete [] this->nodes;
}

void Grid::buildUniformGrid() {

    double hx, hy;
    int bottomID, dxID, topID, sxID;

    if (this->nNodesY < 2)
        swe::printError("NODES_Y must be greater than 1 for a 2D quadrilateral mesh",
                        __PRETTY_FUNCTION__);

    hx = (this->xf - this->x0)/(this->nNodesX - 1);
    hy = (this->yf - this->y0)/(this->nNodesY - 1);

    for (int i = 0; i < this->nNodesY; ++i) {
        for (int j = 0; j < this->nNodesX; ++j) {
            this->nodes[j + i*this->nNodesX][0] = this->x0 + j*hx;
            this->nodes[j + i*this->nNodesX][1] = this->y0 + i*hy;
        }
    }

    this->nodesIDOfBoundaryID.resize(this->nBoundaries);

    for (int i = 0; i < this->nNodesX; ++i) {
        bottomID = i;
        topID    = this->nNodes - 1 - i;
        this->nodesIDOfBoundaryID[0].push_back(bottomID);
        this->nodesIDOfBoundaryID[2].push_back(topID);
    }

    for (int i = 0; i < this->nNodesY; ++i) {
        dxID = (i + 1)*this->nNodesX - 1;
        sxID = i*this->nNodesX;
        this->nodesIDOfBoundaryID[1].push_back(dxID);
        this->nodesIDOfBoundaryID[3].push_back(sxID);
    }
}

void Grid::getAdjacentNodesToInterface(const int &interfaceID, int &leftID, int &rightID) const {

    if (interfaceID > this->nInterfaces - 1)
        swe::printError("Interface " + to_string(interfaceID) + " does not exist.", __PRETTY_FUNCTION__);

    if (interfaceID < (this->nNodesX - 1)*this->nNodesY) {

        int *sxLimit, *dxLimit;

        sxLimit = new int[this->nNodesY];
        dxLimit = new int[this->nNodesY];

        for (int i = 0; i < this->nNodesY; ++i) {

            sxLimit[i] = i*(this->nNodesX - 1);
            dxLimit[i] = i*(this->nNodesX - 1) + (this->nNodesX - 2);

            if (interfaceID <= dxLimit[i] && interfaceID >= sxLimit[i]) {
                leftID  = interfaceID + i;
                rightID = interfaceID + i + 1;
                break;
            }
        }

        delete [] sxLimit;
        delete [] dxLimit;
    } else {
        leftID  = interfaceID - (this->nNodesX - 1)*this->nNodesY;
        rightID = interfaceID - (this->nNodesX - 1)*this->nNodesY + this->nNodesX;
    }
}

void Grid::getNormalToDualGridEdge(const int &dualEdgeID, double *normal) const {

    int mainLeftID, mainRightID;
    double xDiff, yDiff, length;

    if (dualEdgeID > this->nInterfaces - 1)
        swe::printError("Interface " + to_string(dualEdgeID) + " does not exist.", __PRETTY_FUNCTION__);

    this->getAdjacentNodesToInterface(dualEdgeID, mainLeftID, mainRightID);

    xDiff  = this->nodes[mainRightID][0] - this->nodes[mainLeftID][0];
    yDiff  = this->nodes[mainRightID][1] - this->nodes[mainLeftID][1];
    length = sqrt(pow(xDiff, 2) + pow(yDiff, 2));

    normal[0] = xDiff/length;
    normal[1] = yDiff/length;
}

void Grid::getDualElementMeasures(const int &dualElementID, double &hx, double &hy) {

    bool isOnBoundary = false, isOnCorner = false, isTopOrBottom = false;

    if (dualElementID > this->nNodes - 1)
        swe::printError("Node " + to_string(dualElementID) + " does not exist.", __PRETTY_FUNCTION__);

    for (int i = 0; i < this->nBoundaries && !isOnBoundary; ++i) {
        for (int j = 0; j < this->nodesIDOfBoundaryID[i].size() && !isOnBoundary; ++j) {
            if (this->nodesIDOfBoundaryID[i][j] == dualElementID) {
                isOnBoundary = true;
                if (i == 0 || i == 2)
                    isTopOrBottom = true;
                if (this->nodesIDOfBoundaryID[i][0] == dualElementID ||
                    this->nodesIDOfBoundaryID[i][this->nodesIDOfBoundaryID[i].size() - 1] == dualElementID)
                    isOnCorner = true;
            }
        }
    }

    hx = (this->xf - this->x0)/(this->nNodesX - 1);
    hy = (this->yf - this->y0)/(this->nNodesY - 1);

    if (isOnBoundary) {
        if (isOnCorner) {
            hx = hx/2;
            hy = hy/2;
        } else if (isTopOrBottom)
            hy = hy/2;
        else
            hx = hx/2;
    }
}

double Grid::getDualEdgeMeasure(const int &dualEdgeID) {

    double measure, hx, hy;
    double *normal;
    int leftID, rightID;

    if (dualEdgeID > this->nInterfaces - 1)
        swe::printError("Interface " + to_string(dualEdgeID) + " does not exist.", __PRETTY_FUNCTION__);

    normal = new double [2];

    this->getAdjacentNodesToInterface(dualEdgeID, leftID, rightID);
    this->getNormalToDualGridEdge(dualEdgeID, normal);
    this->getDualElementMeasures(leftID, hx, hy);

    measure = (normal[0] == 1) ? hy : hx;

    delete [] normal;

    return measure;
}

double Grid::getDualElementArea(const int &dualElementID) {

    double hx, hy, area;

    if (dualElementID > this->nNodes - 1)
        swe::printError("Node " + to_string(dualElementID) + " does not exist", __PRETTY_FUNCTION__);

    this->getDualElementMeasures(dualElementID, hx, hy);
    area = hx*hy;

    return area;
}

int Grid::getNNodes() const {

    return this->nNodesX*this->nNodesY;
}

int Grid::getNInterfaces() const {

    return (this->nNodesX - 1)*this->nNodesY + (this->nNodesY - 1)*this->nNodesX;
}