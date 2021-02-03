#include "Config.h"

Config::Config(const string &configFile) {

    this->initDefault();
    this->readConfigOptions(configFile);

    if (this->currentRank == MASTER_NODE) {
        this->printHeader();
        this->printSimulationSettings();
    }
}

void Config::initDefault() {

    MPI_Comm_rank(MPI_COMM_WORLD, &this->currentRank);

    this->nNodesX   = SWE_UNDEFINED;
    this->nNodesY   = SWE_UNDEFINED;
    this->x0Lim     = SWE_UNDEFINED;
    this->xFLim     = SWE_UNDEFINED;
    this->y0Lim     = SWE_UNDEFINED;
    this->yFLim     = SWE_UNDEFINED;
    this->totalTime = SWE_UNDEFINED;
    this->cfl       = SWE_UNDEFINED;

    this->markerEuler.clear();
    this->markerOpen.clear();
}

void Config::readConfigOptions(const string &configFile) {

    ifstream SETTINGS(configFile);

    if (SETTINGS.is_open()) {

        this->nNodesX = getOptionInteger(SETTINGS, configOptionList[NODES_X], true);
        this->nNodesY = getOptionInteger(SETTINGS, configOptionList[NODES_Y], true);

        this->x0Lim = getOptionDouble(SETTINGS, configOptionList[X_INITIAL], true);
        this->xFLim = getOptionDouble(SETTINGS, configOptionList[X_FINAL], true);
        this->y0Lim = getOptionDouble(SETTINGS, configOptionList[Y_INITIAL], true);
        this->yFLim = getOptionDouble(SETTINGS, configOptionList[Y_FINAL], true);

        this->totalTime = getOptionDouble(SETTINGS, configOptionList[TOTAL_TIME], true);
        this->cfl       = getOptionDouble(SETTINGS, configOptionList[CFL_NUMBER], true);

        this->numMethodFlux_s = getOptionString(SETTINGS, configOptionList[NUMERICAL_FLUX], numericalFLuxList,
                                                numericalFluxDefault);
        this->numMethodFlux   = getEnumFromMap(this->numMethodFlux_s, mapNumericalFlux);

        this->entropyFix_s = getOptionString(SETTINGS, configOptionList[ENTROPY_FIX], entropyFixList, entropyFixDefault);
        this->entropyFix   = getEnumFromMap(this->entropyFix_s, mapEntropyFix);

        this->fluxLimiter_s = getOptionString(SETTINGS, configOptionList[FLUX_LIMITER], fluxLimiterList,
                                              fluxLimiterDefault);
        this->fluxLimiter   = getEnumFromMap(this->fluxLimiter_s, mapFluxLimiter);

        this->markerEuler = getOptionVectorStrings(SETTINGS, configOptionList[MARKER_EULER]);
        this->markerOpen  = getOptionVectorStrings(SETTINGS, configOptionList[MARKER_OPEN]);
    }
    else {
        printError("No " + configFile + " found.");
    }
}

void Config::printSimulationSettings() {

    cout << "|---------------------------- SIMULATION SETTINGS -----------------------------|" << endl;
    cout << endl;
    cout << " Mesh definition:" << endl;
    cout << " ----------------" << endl;
    cout << endl;
    cout << "  Number of nodes on X axis: \t " << this->nNodesX << endl;
    cout << "  Number of nodes on Y axis: \t " << this->nNodesY << endl;
    cout << "  Boundary limits coordinates: \t x0 = " << this->x0Lim << endl;
    cout << "                               \t xf = " << this->xFLim << endl;
    cout << "                               \t y0 = " << this->y0Lim << endl;
    cout << "                               \t yf = " << this->yFLim << endl;
    cout << endl;
    cout << " Numerical method:" << endl;
    cout << " -----------------" << endl;
    cout << endl;
    cout << "  Total simulation time: \t " << this->totalTime << endl;
    cout << "  CFL number: \t\t\t " << this->cfl << endl;
    cout << "  Flux numerical method: \t " << this->numMethodFlux_s << endl;
    cout << "  Entropy fix: \t\t\t " << this->entropyFix_s << endl;
    cout << "  Flux limiter: \t\t " << this->fluxLimiter_s << endl;
    cout << endl;
    cout << " Boundary conditions:" << endl;
    cout << " --------------------" << endl;
    cout << endl;
    cout << "  Euler wall boundary marker: \t " << "(";
            for (int i = 0; i < this->markerEuler.size(); ++i) {
                if (i != this->markerEuler.size() - 1)
                    cout << this->markerEuler[i] << " ";
                else
                    cout << this->markerEuler[i];
            }
    cout << ")" << endl;
    cout << "  Open boundary marker: \t " << "(";
            for (int i = 0; i < this->markerOpen.size(); ++i) {
                if (i != this->markerOpen.size() - 1)
                    cout << this->markerOpen[i] << " ";
                else
                    cout << this->markerOpen[i];
            }
    cout << ")" << endl;
    cout << endl;
}

void Config::printHeader() {

    cout << "|------------------------------------------------------------------------------|" << endl;
    cout << "|                                                                              |" << endl;
    cout << "| 2D shallow water solver equations solver                                     |" << endl;
    cout << "|                                                                              |" << endl;
    cout << "|------------------------------------------------------------------------------|" << endl;
    cout << endl;
}