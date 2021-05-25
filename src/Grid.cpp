#include "Grid.h"

Grid::Grid() {

    this->x0Domain = 0.0;                   this->x0 = 0.0;
    this->y0Domain = 0.0;                   this->y0 = 0.0;
    this->xfDomain = 0.0;                   this->xf = 0.0;
    this->yfDomain = 0.0;                   this->yf = 0.0;

    this->nNodesXDomain = 0;                this->nNodesX = 0;
    this->nNodesYDomain = 0;                this->nNodesY = 0;

    this->nNodes      = 0;
    this->nInterfaces = 0;
    this->nBoundaries = 0;

    this->nodes = nullptr;

    this->nodesIDOfBoundaryID.clear();
}

Grid::Grid(Config *config) {

    this->x0Domain = config->getXInitialLimit();               this->x0 = 0.0;
    this->y0Domain = config->getYInitialLimit();               this->y0 = 0.0;
    this->xfDomain = config->getXFinalLimit();                 this->xf = 0.0;
    this->yfDomain = config->getYFinalLimit();                 this->yf = 0.0;

    this->nNodesXDomain = config->getNNodesX();                this->nNodesX = 0;
    this->nNodesYDomain = config->getNNodesY();                this->nNodesY = 0;

    this->nNodes      = 0;
    this->nInterfaces = 0;
    this->nBoundaries = 0;

    this->nodes = nullptr;

    this->nodesIDOfBoundaryID.clear();

    this->buildSubDomains();
}

Grid::~Grid() {

    for (int i = 0; i < this->nNodes; ++i) {
        delete [] this->nodes[i];
    }
    delete [] this->nodes;
}

void Grid::buildSubDomains() {

    MPI_Comm cartComm;
    int procForDim[2] = {0, 0}, periods[2] = {0, 0}, coords[2] = {0, 0};
    int size, rank;
    int nNodesXForProc, nNodesYForProc, remainderX, remainderY;
    double hx, hy;

    if (this->nNodesXDomain < 2 || this->nNodesYDomain < 2)
        swe::printError("NODES_X and NODES_Y must be greater than 1 for a 2D quadrilateral mesh",
                        __PRETTY_FUNCTION__);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Dims_create(size, 2, procForDim);
    MPI_Cart_create(MPI_COMM_WORLD, 2, procForDim, periods, 1, &cartComm);

    MPI_Cart_coords(cartComm, rank, 2, coords);

    nNodesXForProc = this->nNodesXDomain/procForDim[0];
    nNodesYForProc = this->nNodesYDomain/procForDim[1];
    remainderX = this->nNodesXDomain % procForDim[0];
    remainderY = this->nNodesYDomain % procForDim[1];

    this->nNodesX = nNodesXForProc;
    this->nNodesY = nNodesYForProc;

    if (coords[0] == procForDim[0] - 1)
        this->nNodesX += remainderX;
    if (coords[1] == procForDim[1] - 1)
        this->nNodesY += remainderY;

    hx = (this->xfDomain - this->x0Domain)/(this->nNodesXDomain - 1);
    hy = (this->yfDomain - this->y0Domain)/(this->nNodesYDomain - 1);

    this->x0 = this->x0Domain + nNodesXForProc*hx*coords[0];
    this->y0 = this->y0Domain + nNodesYForProc*hy*coords[1];
    this->xf = this->x0 + (this->nNodesX - 1)*hx;
    this->yf = this->y0 + (this->nNodesY - 1)*hy;

    this->nNodes      = this->getNNodes();
    this->nInterfaces = this->getNInterfaces();
    this->nBoundaries = 4;

    this->nodes = new double* [this->nNodes];
    for (int i = 0; i < this->nNodes; ++i) {
        this->nodes[i] = new double [2];
    }

    this->buildUniformGrid();

   // this->printRawMesh();
}

void Grid::buildUniformGrid() {

    double hx, hy;
    int bottomID, dxID, topID, sxID;

    hx = (this->xf - this->x0)/(this->nNodesX - 1);
    hy = (this->yf - this->y0)/(this->nNodesY - 1);

    if (this->nNodesX == 1 && this->nNodesY > 1)
        hx = 0;
    else if (this->nNodesY == 1 && this->nNodesX > 1)
        hy = 0;

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

void Grid::printRawMesh() const {

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < size; ++i) {

        if (i == rank) {
            cout << endl;
            cout << "Subdomain on rank: " << i << endl;
            cout << "--------------------" << endl << endl;

            for (int j = 0; j < this->nNodes; ++j) {
                cout << "ID: " << j << " Coords: (" << this->nodes[j][0] << ";" << this->nodes[j][1] << ")" << endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}