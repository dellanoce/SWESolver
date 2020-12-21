#ifndef SWE_ENUMS_H
#define SWE_ENUMS_H

#include <string>
#include <vector>
#include <map>

#define SWE_UNDEFINED -1

using namespace std;

/*!
 * \file Enums.h
 * \brief Header Enums.h contains enumerators defining configuration options and maps handling configuration
 * sub-options
 * \author Alberto Della Noce
 */

// --------------------------------------------- CONFIG OPTIONS ---------------------------------------------- //

/*!
 * \brief Enum listing config options
 */
enum configOption_t {

    NODES_X,
    NODES_Y,
    X_INITIAL,
    X_FINAL,
    Y_INITIAL,
    Y_FINAL,

    TOTAL_TIME,
    CFL_NUMBER,
    NUMERICAL_FLUX,
    ENTROPY_FIX,
    FLUX_LIMITER,

    MARKER_EULER,
    MARKER_OPEN,
};

/*!
 * \brief Vector containing options strings
 */
static const vector<string> configOptionList = {

        "NODES_X",
        "NODES_Y",
        "x0",
        "xF",
        "y0",
        "yF",

        "TOTAL_TIME",
        "CFL_NUMBER",
        "FLUX_NUM_METHOD",
        "ENTROPY_FIX",
        "FLUX_LIMITER",

        "MARKER_EULER",
        "MARKER_OPEN",
};

// ------------------------------------------- CONFIG SUB-OPTIONS -------------------------------------------- //

/*!
 * \brief Enumerator listing numerical methods for computing fluxes
 */
enum numericalFlux_t {
    ROE          = 0,
    LAX_WENDROFF = 1,
    HIGH_RES     = 2
};

/*!
 * \brief Map storing numerical methods sub-options
 */
static const map<numericalFlux_t, string> mapNumericalFlux = {
        {ROE, "ROE"},
        {LAX_WENDROFF, "LAX-WENDROFF"},
        {HIGH_RES, "HIGH_RES"}
};

/*!
 * \brief Vector containing options for numerical methods
 */
static const vector<string> numericalFLuxList = {"ROE", "LAX-WENDROFF", "HIGH_RES"};

/*!
 * \brief Default numerical method
 */
static const string numericalFluxDefault = "ROE";

// ----------------------------------------------------------------------------------------------------------- //

/*!
 * \brief Enumerator listing options for entropy fix
 */
enum entropyFix_t {
    NO  = false,
    YES = true
};

/*!
 * \brief Map storing entropy fix sub-options
 */
static const map<entropyFix_t, string> mapEntropyFix = {
        {NO, "NO"},
        {YES, "YES"}
};

/*!
 * \brief Vector containing options for entropy fix
 */
static const vector<string> entropyFixList = {"NO", "YES"};

/*!
 * \brief Default setting for entropy fix
 */
static const string entropyFixDefault = "NO";

// ----------------------------------------------------------------------------------------------------------- //

/*!
 * \brief Enumerator listing options for flux limiter
 */
enum fluxLimiter_t {
    NONE     = 0,
    SUPERBEE = 1,
    VAN_LEER = 2,
    MINMOD   = 3
};

/*!
 * \brief Map storing flux limiter sub-options
 */
static const map<fluxLimiter_t, string> mapFluxLimiter = {
        {NONE, "NONE"},
        {SUPERBEE, "SUPERBEE"},
        {VAN_LEER, "VAN_LEER"},
        {MINMOD, "MINMOD"}
};

/*!
 * \brief Vector containing options for flux limiter
 */
static const vector<string> fluxLimiterList = {"NONE", "SUPERBEE", "VAN_LEER", "MINMOD"};

/*!
 * \brief Default flux limiter
 */
static const string fluxLimiterDefault = "NONE";

#endif //SWE_ENUMS_H
