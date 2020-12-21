#ifndef SWE_SOLUTION_H
#define SWE_SOLUTION_H

#include "Grid.h"

using namespace std;

/*!
 * \class Solution
 * \brief Class Solution defines the structure of the solution container
 * \author Alberto Della Noce
 */

class Solution {
private:
    double *data; /*!< \brief Vector containing the solution variables at each node of the mesh */
    int nNodes; /*!< \brief Number of nodes on the mesh */

    /*!
     * \brief Initializes the solution on a mesh
     * @param[in] meshPtr - Pointer to the mesh
     */
    void initSolutionOnMesh(Grid *meshPtr);
public:

    int nVar; /*!< \brief Number of variables of the problem */

    /*!
     * \brief Constuctor of the class
     */
    Solution();

    /*!
     * \overload
     * @param[in] meshPtr - Pointer to the mesh
     */
    explicit Solution(Grid* meshPtr);

    /*!
     * \brief Destructor of the class
     */
    ~Solution();

    /*!
     * \brief Returns the value of the conservative variable 1 (h) at nodeID
     * @param[in] nodeID - ID of the node
     * @return Value of conservative 1
     */
    double &conservative1(const int &nodeID);

    /*!
     * \brief Returns the value of the conservative variable 2 (h*u) at nodeID
     * @param[in] nodeID - ID of the node
     * @return Value of conservative 2
     */
    double &conservative2(const int &nodeID);

    /*!
     * \brief Returns the value of the conservative variable 3 (h*v) at nodeID
     * @param[in] nodeID - ID of the node
     * @return Value of conservative 3
     */
    double &conservative3(const int &nodeID);

    /*!
     * \brief Returns the value of the primitive variable 1 (h) at nodeID
     * @param[in] nodeID - ID of the node
     * @return Value of primitive 1
     */
    double primitive1(const int &nodeID);

    /*!
     * \brief Returns the value of the primitive variable 2 (u) at nodeID
     * @param[in] nodeID - ID of the node
     * @return Value of primitive 2
     */
    double primitive2(const int &nodeID);

    /*!
     * \brief Returns the value of the primitive variable 3 (v) at nodeID
     * @param[in] nodeID - ID of the node
     * @return Value of primitive 3
     */
    double primitive3(const int &nodeID);
};

#endif //SWE_SOLUTION_H
