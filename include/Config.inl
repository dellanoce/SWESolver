#ifndef SWE_CONFIG_INL
#define SWE_CONFIG_INL

#include "Config.h"

inline int Config::getNNodesX() { return nNodesX; }

inline int Config::getNNodesY() { return nNodesY; }

inline int Config::getNumericalFluxMethod() { return numMethodFlux; }

inline int Config::getFluxLimiter() { return fluxLimiter; }

inline double Config::getXInitialLimit() { return x0Lim; }

inline double Config::getXFinalLimit() { return xFLim; }

inline double Config::getYInitialLimit() { return y0Lim; }

inline double Config::getYFinalLimit() { return yFLim; }

inline double Config::getCFLNumber() { return cfl; }

inline double Config::getTotalTime() { return totalTime; }

inline bool Config::getEntropyFix() { return entropyFix; }

inline vector<string> Config::getMarkerEuler() { return markerEuler; }

inline vector<string> Config::getMarkerOpen() { return markerOpen; }

#endif //SWE_CONFIG_INL