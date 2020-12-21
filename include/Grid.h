#ifndef SWE_GRID_H
#define SWE_GRID_H

#include <cmath>

#include "Config.h"

using namespace std;

/*!
 * \class Grid
 * \brief Class Grid defines a 2D structured quadrilateral mesh on a rectangular domain
 * \author Alberto Della Noce
 */

class Grid {
public:
    double x0; /*!< \brief X-coordinate of the first vertex on X-axis */
    double xf; /*!< \brief X-coordinate of the second vertex on X-axis */
    double y0; /*!< \brief Y-coordinate of the first vertex on Y-axis */
    double yf; /*!< \brief Y-coordinate of the second vertex on Y-axis */
    int nNodesX; /*!< \brief Number of nodes on X-direction */
    int nNodesY; /*!< \brief Number of nodes on Y-direction */
    int nNodes; /*!< \brief Total number of nodes */
    int nInterfaces; /*!< \brief Number of interfaces (number of edges) */
    int nBoundaries; /*!< \brief Number of boundaries */
    double **nodes; /*!< \brief Matrix containing mesh nodes coordinates */
    vector<vector<int>> nodesIDOfBoundaryID; /*!< \brief Matrix containing nodes ID for each boundary */

    /*!
     * \brief Constructor of the class
     */
    Grid();

    /*!
     * \overload
     * @param[in] config - Pointer to Config class
     */
    explicit Grid(Config *config);

    /*!
     * \brief Destructor of the class
     */
    ~Grid();

    /*!
     * \brief Builds a 2D structured quadrilateral mesh on a rectangular domain
     */
    void buildUniformGrid();

    /*!
     * \brief Gets the adjacent nodes to an interface (edge)
     * @param[in] interfaceID - ID of the interface (edge)
     * @param[out] leftID - ID of the first adjacent node
     * @param[out] rightID - ID of the second adjacent node
     */
    void getAdjacentNodesToInterface(const int &interfaceID, int &leftID, int &rightID) const;

    /*!
     * \brief Gets the normal vector to an edge of the dual grid
     * @param[in] dualEdgeID - ID of the edge on the dual grid
     * @param[out] normal - Normal vector to the edge
     */
    void getNormalToDualGridEdge(const int &dualEdgeID, double *normal) const;

    /*!
     * \brief Gets measures of a dual grid element (length of both edges)
     * @param[in] dualElementID - ID of the dual element
     * @param[out] hx - Length of the edge on X-direction
     * @param[out] hy - Length of the edge on Y-direction
     */
    void getDualElementMeasures(const int &dualElementID, double &hx, double &hy);

    /*!
     * \brief Returns the measure of an edge of the dual grid
     * @param[in] dualEdgeID - ID of the edge on the dual face
     * @return Measure of the edge
     */
    double getDualEdgeMeasure(const int &dualEdgeID);

    /*!
     * \brief Returns the area of a dual element
     * @param[in] dualElementID - ID of the dual element
     * @return Area of the dual element
     */
    double getDualElementArea(const int &dualElementID);

    /*!
     * \brief Returns the number of nodes on the mesh
     * @return Number of nodes on the mesh
     */
    int getNNodes() const;

    /*!
     * \brief Returns the number of interfaces (edges) on the mesh
     * @return Number of interfaces on the mesh
     */
    int getNInterfaces() const;
};

#endif //SWE_GRID_H
