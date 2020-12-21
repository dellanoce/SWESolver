#ifndef SWE_SOLVER_H
#define SWE_SOLVER_H

#include <memory>

#include "NumericalFluxes.h"
#include "OutputManager.h"

/*!
 * \class Solver
 * \brief Class Solver implements the finite volume method for solving the shallow waters equations
 * \author Alberto Della Noce
 */

class Solver {
public:
    shared_ptr<NumericalFluxes> fluxPtr; /*!< \brief Smart pointer to numerical fluxes */
    Solution* solPtr; /*!< \brief Pointer to the solution */
    double totalTime; /*!< \brief Total simulation time */
    double dt; /*!< \brief Time step */
    double cfl; /*!< \brief CFL number */
    double *residual; /*!< \brief Vector of residual */
    double *consVector; /*!< \brief Vector of conservative variables */
    vector<string> eulerBoundaries; /*!< \brief Vector of Euler wall boundaries tag */

    /*!
     * \brief Constructor of the class
     */
    Solver();

    /*!
     * \overload
     * @param[in] config - Pointer to the Config class
     * @param[in] mesh - Pointer to the mesh
     */
    Solver(Config *config, Grid *mesh);

    /*!
     * \brief Destructor of the class
     */
    ~Solver();

    /*!
     * \brief Sets the time step at every iteration with a fixed CFL number
     * @param[in] mesh - Pointer to the mesh
     */
    void setTimeStep(Grid *mesh);

    /*!
     * \brief Sets the initial condition
     * @param[in] mesh - Pointer to the mesh
     * @param[in] initialSol - Reference to the Solution
     */
    void setInitialSolution(Grid *mesh, Solution &initialSol);

    /*!
     * \brief Runs the solver
     * @param[in] mesh - Pointer to the mesh
     */
    void runGodnuovSolver(Grid *mesh);

    /*!
     * \brief Updates the solution with the new computed values
     * @param[in] mesh - Pointer to the mesh
     */
    void updateSolution(Grid *mesh) const;

    /*!
     * \brief Initializes the residual vector
     * @param[in] mesh - Pointer to the mesh
     */
    void initializeResidual(Grid *mesh) const;

    /*!
     * \brief Adds the numerical flux at an interface (edge) to the residual on a node
     * @param[in] nodeID - ID of the node
     * @param[in] length - Length of the edge
     * @param[in] numericalFlux - Numerical flux vector
     */
    void addToResidual(const int &nodeID, double length, const double *numericalFlux) const;

    /*!
     * \brief Subtracts the numerical flux at an interface (edge) to the residual on a node
     * @param[in] nodeID - ID of the node
     * @param[in] length - Length of the edge
     * @param[in] numericalFlux - Numerical flux vector
     */
    void subtractToResidual(const int &nodeID, double length, const double *numericalFlux) const;

    /*!
     * \brief Sets boundary conditions on boundaryID
     * @param[in] mesh - Pointer to the mesh
     * @param[in] boundaryID - ID of the boundary
     */
    void setBoundaryConditions(Grid *mesh, const int &boundaryID) const;
};

#endif //SWE_SOLVER_H
