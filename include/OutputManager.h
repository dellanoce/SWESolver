#ifndef SWE_OUTPUTMANAGER_H
#define SWE_OUTPUTMANAGER_H

#include "Grid.h"
#include "Solution.h"

/*!
 * \class OutputManager
 * \brief Static class OutputManager contains methods for printing the solution in different formats
 * \author Alberto Della Noce
 */

class OutputManager {
public:

    /*!
     * \brief Constructor of the class
     */
    OutputManager() = default;

    /*!
     * \brief Destructor of the class
     */
    ~OutputManager() = default;

    /*!
     * \brief Prints the solution in VTK format
     * @param[in] mesh - Pointer to the mesh
     * @param[in] sol - Pointer to the solution
     * @param[in] time - Time at which the solution is printed
     */
    static void printSolutionVTK(Grid *mesh, Solution *sol, double time);

    /*!
     * \brief Prints the solution in csv format
     * @param[in] mesh - Pointer to the mesh
     * @param[in] sol - Pointer to the solution
     * @param[in] time - Time at which the solution is printed
     */
    static void printSolutionCsv(Grid *mesh, Solution *sol, double time);
};

#endif //SWE_OUTPUTMANAGER_H
