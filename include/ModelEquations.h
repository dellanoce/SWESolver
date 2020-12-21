#ifndef SWE_MODELEQUATIONS_H
#define SWE_MODELEQUATIONS_H

#include <cmath>

/*!
 * \class ModelEquations
 * \brief Static class ModelEquations contains methods defining the physical model of shallow waters
 * \author Alberto Della Noce
 */

const double g0 = 9.80665; /*!< \brief Gravitational acceleration at sea level */

class ModelEquations {
public:

    /*!
     * \brief Constructor of the class
     */
    ModelEquations() = default;

    /*!
     * \brief Destructor of the class
     */
    ~ModelEquations() = default;

    /*!
     * \brief Sets the flux vector projected on the direction defined by the normal
     * @param[in] height - Height value
     * @param[in] velocity - Velocity vector
     * @param[in] normal - Normal vector
     * @param[out] flux - Projected flux vector
     */
    static void setProjFlux(double height, double *velocity, const double *normal, double *flux);

    /*!
     * \brief Sets the Jacobian matrix projected on the direction defined by the normal
     * @param[in] height - Height value
     * @param[in] velocity - Velocity vector
     * @param[in] normal - Normal vector
     * @param[out] projJacobian - Projected Jacobian matrix
     */
    static void setProjJacobian(double height, double *velocity, const double *normal, double **projJacobian);

    /*!
     * \brief Sets the matrix of right eigenvectors projected on the direction defined by the normal
     * @param[in] height - Height value
     * @param[in] velocity - Velocity vector
     * @param[in] normal - Normal vector
     * @param[out] projRightMatrix - Projected right eigenvectors matrix
     */
    static void setProjRightMatrix(double height, const double *velocity, const double *normal, double **projRightMatrix);

    /*!
     * \brief Sets the matrix of left eigenvectors projected on the direction defined by the normal
     * @param[in] height - Height value
     * @param[in] velocity - Velocity vector
     * @param[in] normal - Normal vector
     * @param[out] projLeftMatrix - Projected left eigenvectors matrix
     */
    static void setProjLeftMatrix(double height, const double *velocity, const double *normal, double **projLeftMatrix);

    /*!
     * \brief Sets the eigenvalues projected on the direction defined by the normal
     * @param[in] height - Height value
     * @param[in] velocity - Velocity vector
     * @param[in] normal - Normal vector
     * @param[out] lambda - Vector of eigenvalues
     */
    static void setProjEigenvalues(double height, const double *velocity, const double *normal, double *lambda);
};

#endif //SWE_MODELEQUATIONS_H
