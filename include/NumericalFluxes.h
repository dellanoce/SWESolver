#ifndef SWE_NUMERICALFLUXES_H
#define SWE_NUMERICALFLUXES_H

#include "ModelEquations.h"
#include "Solution.h"

/*!
 * \class NumericalFluxes
 * \brief Class NumericalFluxes is an abstract class defining general attributes used to compute numerical
 * fluxes with different numerical methods
 * \author Alberto Della Noce
 */

class NumericalFluxes {
protected:
    double *normal; /*!< \brief Normal vector to an edge */
    double *flux_l; /*!< \brief Flux evaluated on left node */
    double *flux_r; /*!< \brief Flux evaluated on right node */
    double *velocity_l; /*!< \brief Velocity evaluated on left node */
    double *velocity_r; /*!< \brief Velocity evaluated on right node */
    double height_l; /*!< \brief Height evaluated on left node */
    double height_r; /*!< \brief Height evaluated on right node */

public:
    double *numericalFlux; /*!< \brief Numerical flux vector */

    /*!
     * \brief Constructor of the class
     */
    NumericalFluxes();

    /*!
     * \brief Destructor of the class
     */
    virtual ~NumericalFluxes();

    /*!
     * \brief Virtual method for computing numerical fluxes
     * @param[in] mesh - Pointer to the mesh
     * @param[in] sol - Pointer to the solution
     * @param[in] interfaceID - Interface (edge) where the numerical flux is computed
     * @param[in] dt - Time step
     */
    virtual void computeFlux(Grid *mesh, Solution *sol, int &interfaceID, double &dt) = 0;
};

/*!
 * \class RoeFlux
 * \brief Child class RoeFlux defines attributes and methods for computing the numerical flux with first
 * order Roe linearization
 * \author Alberto Della Noce
 */

class RoeFlux : public NumericalFluxes {
protected:
    bool entropyFix; /*!< \brief Enables or disables the entropy fix */
    double *lambda_l; /*!< \brief Eigenvalues evaluated on the left node */
    double *lambda_r; /*!< \brief Eigenvalues evaluated on the right node */
    double *lambda; /*!< \brief Eigenvalues evaluated on Roe averaged variables */
    double *roeVelocity; /*!< \brief Roe averaged velocity */
    double *diffW; /*!< \brief Difference between right and left conservative variables */
    double **matrixR; /*!< \brief Matrix of right eigenvectors */
    double **matrixL; /*!< \brief Matrix of left eigenvectors */
    double roeHeight; /*!< \brief Roe averaged height */

public:

    /*!
     * \brief Constructor of the class
     * @param[in] config - Pointer to the Config class
     * @param[in] sol - Reference to the solution
     */
    RoeFlux(Config *config, Solution &sol);

    /*!
     * \brief Destructor of the class
     */
    ~RoeFlux() override;

    /*!
     * \brief Overridden method for computing the numerical flux on interfaceID using the Roe linearization
     * @param[in] mesh - Pointer to the mesh
     * @param[in] sol - Pointer to the solution
     * @param[in] interfaceID - ID of the interface
     * @param[in] dt - Time step
     */
    void computeFlux(Grid *mesh, Solution *sol, int &interfaceID, double &dt) override;
};

/*!
 * \class LaxWendroffFlux
 * \brief Child class LaxWendroffFlux defines attributes and methods for computing the numerical flux with
 * Lax-Wendroff second order scheme
 * \author Alberto Della Noce
 */

class LaxWendroffFlux : public NumericalFluxes {
private:
    double *flux_inter; /*!< \brief Flux evaluated on the interface (edge) */
    double *velocity_inter; /*!< \brief Velocity evaluated on the interface (edge) */
    double **projJacobian; /*!< \brief Jacobian matrix projected in normal direction */
    double height_inter; /*!< \brief Height evaluated on the interface (edge) */

public:

    /*!
     * \brief Constructor of the class
     * @param[in] sol - Reference to the solution
     */
    explicit LaxWendroffFlux(Solution &sol);

    /*!
     * \brief Destructor of the class
     */
    ~LaxWendroffFlux() override;

    /*!
     * \brief Overridden method for computing the numerical flux on interfaceID using the Lax-Wendroff scheme
     * @param[in] mesh - Pointer to the mesh
     * @param[in] sol - Pointer to the solution
     * @param[in] interfaceID - ID of the interface
     * @param[in] dt - Time step
     */
    void computeFlux(Grid *mesh, Solution *sol, int &interfaceID, double &dt) override;
};

/*!
 * \class HighResFlux
 * \brief Child class HighResFlux defines attributes and methods for computing the numerical flux with an
 * high resolution scheme
 * \author Alberto Della Noce
 */

class HighResFlux : public RoeFlux {
private:
    double *diffW_l; /*!< \brief Difference between left and upwind left conservative variables */
    double *diffW_r; /*!< \brief Difference between right and upwind right conservative variables */
    double *diffV_l; /*!< \brief Difference between characteristics left and upwind left */
    double *diffV_r; /*!< \brief Difference between characteristics right and upwind right */
    double *diffV_up; /*!< \brief Difference between upwind characteristics */
    double *diffV; /*!< \brief Difference between right and left characteristics */
    double *lambdaSigned; /*!< \brief Signed eigenvalues */
    double phi; /*!< \brief Value of the flux limiter */
    int fluxLimiter; /*!< Selected flux limiter */

public:

    /*!
     * \brief Constructor of the class
     * @param[in] config - Pointer to the Config class
     * @param[in] sol - Reference to the solution
     */
    HighResFlux(Config *config, Solution &sol);

    /*!
     * \brief Destructor of the class
     */
    ~HighResFlux() override;

    /*!
     * \brief Sets the value of the flux limiter
     * @param[in] upDiff - Upwind difference between right and left variables
     * @param[in] cenDiff - Centered difference between right and left variables
     */
    void setFluxLimiter(double upDiff, double cenDiff);

    /*!
     * \brief Overridden method for computing the numerical flux on interfaceID using an high resolution scheme
     * @param[in] mesh - Pointer to the mesh
     * @param[in] sol - Pointer to the solution
     * @param[in] interfaceID - ID of the interface
     * @param[in] dt - Time step
     */
    void computeFlux(Grid *mesh, Solution *sol, int &interfaceID, double &dt) override;
};

// ---------------------------------------------- TEMPLATES ------------------------------------------------ //

/*!
 * \brief Returns the sign of the object T
 * @tparam[in] T - Datatype for floating point value
 * @param[in] val - Input value
 * @return Sign of val
 */
template <typename T> int sgn(T val) { return (val >= 0) ? 1 : -1; }

#endif //SWE_NUMERICALFLUXES_H
