#ifndef SWE_CONFIG_H
#define SWE_CONFIG_H

#include "Enums.h"
#include "Utils.h"

using namespace std;
using namespace swe;

/*!
 * \class Config
 * \brief Class Config reads the configuration file and stores all the information for the simulation
 * \author Alberto Della Noce
 */
class Config {
private:
    int currentRank;
    int nNodesX; /*!< \brief Number of nodes on X-direction */
    int nNodesY; /*!< \brief Number of nodes on Y-direction */
    double x0Lim; /*!< \brief X-coordinate of the first vertex on X-axis */
    double xFLim; /*!< \brief X-coordinate of the second vertex on X-axis */
    double y0Lim; /*!< \brief Y-coordinate of the first vertex on Y-axis */
    double yFLim; /*!< \brief Y-coordinate of the second vertex on Y-axis */
    double totalTime; /*!< \brief Total simulation time */
    double cfl; /*!< \brief CFL number */
    numericalFlux_t numMethodFlux; /*!< \brief Enumerator for numerical flux options */
    entropyFix_t entropyFix; /*!< \brief Enumerator for entropy fix options */
    fluxLimiter_t fluxLimiter; /*!< \brief Enumerator for flux limiter options */
    vector<string> markerEuler; /*!< \brief Vector containing boundaries tag where Euler wall BC are applied */
    vector<string> markerOpen; /*!< \brief Vector containing boundaries tag where open BC are applied */
    string numMethodFlux_s; /*!< \brief Selected numerical flux method */
    string entropyFix_s; /*!< \brief Selected entropy fix */
    string fluxLimiter_s; /*!< \brief Selected flux limiter */

    /*!
     * \brief Initializes the default values
     */
    void initDefault();

    /*!
     * \brief Parses the config file and stores options and values
     * \param[in] configFile - Name of the config file
     */
    void readConfigOptions(const string& configFile);

    /*!
     * \brief Prints simulation settings on screen
     */
    void printSimulationSettings();

    /*!
     * \brief Prints header on screen
     */
    void printHeader();

public:

    /*!
     * \brief Constructor of the class
     * @param configFile
     */
    explicit Config(const string &configFile);

    /*!
     * \brief Destructor of the class
     */
    ~Config() = default;

    /*!
     * \brief Gets the number of nodes on X-direction
     * @return Number of nodes on X-direction
     */
    int getNNodesX();

    /*!
     * \brief Gets the number of nodes on Y-direction
     * @return Number of nodes on Y-direction
     */
    int getNNodesY();

    /*!
     * \brief Gets the numerical fluxes method
     * @return Enumerator value of the selected numerical flux method
     */
    int getNumericalFluxMethod();

    /*!
     * \brief Gets the flux limiter
     * @return Enumerator value of the selected flux limiter
     */
    int getFluxLimiter();

    /*!
     * \brief Gets the X-coordinate of the first vertex on X-axis
     * @return Coordinate of the first vertex on X-axis
     */
    double getXInitialLimit();

    /*!
     * \brief Gets the X-coordinate of the second vertex on X-axis
     * @return Coordinate of the second vertex on X-axis
     */
    double getXFinalLimit();

    /*!
     * \brief Gets the Y-coordinate of the first vertex on Y-axis
     * @return Coordinate of the first vertex on Y-axis
     */
    double getYInitialLimit();

    /*!
     * \brief Gets the Y-coordinate of the second vertex on Y-axis
     * @return Coordinate of the second vertex on Y-axis
     */
    double getYFinalLimit();

    /*!
     * \brief Gets the total simulation time
     * @return Total simulation time
     */
    double getTotalTime();

    /*!
     * \brief Gets the CFL number
     * @return CFL number
     */
    double getCFLNumber();

    /*!
     * \brief Gets the selected option for entropy fix
     * @return Enumerator value of the entropy fix
     */
    bool getEntropyFix();

    /*!
     * \brief Returns a vector containing boundaries tag where Euler BC are applied
     * @return Vector containing boundaries tag
     */
    vector<string> getMarkerEuler();

    /*!
     * \brief Returns a vector containing boundaries tag where open BC are applied
     * @return Vector containing boundaries tag
     */
    vector<string> getMarkerOpen();
};

#endif //SWE_CONFIG_H

#include "Config.inl"
