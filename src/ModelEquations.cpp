#include "ModelEquations.h"

void ModelEquations::setProjFlux(double height, double *velocity, const double *normal, double *flux) {

    flux[0] = height*velocity[0]*normal[0];
    flux[1] = (height*pow(velocity[0], 2) + 0.5*g0*pow(height, 2))*normal[0];
    flux[2] = height*velocity[0]*velocity[1]*normal[0];

    flux[0] += height*velocity[1]*normal[1];
    flux[1] += height*velocity[0]*velocity[1]*normal[1];
    flux[2] += (height*pow(velocity[1], 2) + 0.5*g0*pow(height, 2))*normal[1];
}

void ModelEquations::setProjJacobian(double height, double *velocity, const double *normal, double **projJacobian) {

    double waveSpeed = sqrt(g0*height);

    projJacobian[0][0] = 0;
    projJacobian[1][0] = (-pow(velocity[0], 2) + pow(waveSpeed, 2))*normal[0];
    projJacobian[2][0] = -velocity[0]*velocity[1]*normal[0];

    projJacobian[0][1] = normal[0];
    projJacobian[1][1] = 2*velocity[0]*normal[0];
    projJacobian[2][1] = velocity[1]*normal[0];

    projJacobian[0][2] = 0;
    projJacobian[1][2] = 0;
    projJacobian[2][2] = velocity[0]*normal[0];

    projJacobian[0][0] += 0;
    projJacobian[1][0] += -velocity[0]*velocity[1]*normal[1];
    projJacobian[2][0] += (-pow(velocity[1], 2) + pow(waveSpeed, 2))*normal[1];

    projJacobian[0][1] += 0;
    projJacobian[1][1] += velocity[1]*normal[1];
    projJacobian[2][1] += 0;

    projJacobian[0][2] += normal[1];
    projJacobian[1][2] += velocity[0]*normal[1];
    projJacobian[2][2] += 2*velocity[1]*normal[1];
}

void ModelEquations::setProjRightMatrix(double height, const double *velocity, const double *normal, double **projRightMatrix) {

    double waveSpeed = sqrt(g0*height);

    projRightMatrix[0][0] = 1;
    projRightMatrix[1][0] = velocity[0] - waveSpeed*normal[0];
    projRightMatrix[2][0] = velocity[1] - waveSpeed*normal[1];

    projRightMatrix[0][1] = 0;
    projRightMatrix[1][1] = normal[1];
    projRightMatrix[2][1] = normal[0];

    projRightMatrix[0][2] = 1;
    projRightMatrix[1][2] = velocity[0] + waveSpeed*normal[0];
    projRightMatrix[2][2] = velocity[1] + waveSpeed*normal[1];
}

void ModelEquations::setProjLeftMatrix(double height, const double *velocity, const double *normal, double **projLeftMatrix) {

    double waveSpeed = sqrt(g0*height);

    projLeftMatrix[0][0] = (velocity[0]*normal[0] + velocity[1]*normal[1] + waveSpeed)/(2*waveSpeed);
    projLeftMatrix[1][0] = -velocity[1]*normal[0] - velocity[0]*normal[1];
    projLeftMatrix[2][0] = (-velocity[0]*normal[0] - velocity[1]*normal[1] + waveSpeed)/(2*waveSpeed);

    projLeftMatrix[0][1] = -normal[0]/(2*waveSpeed);
    projLeftMatrix[1][1] = normal[1];
    projLeftMatrix[2][1] = normal[0]/(2*waveSpeed);

    projLeftMatrix[0][2] = -normal[1]/(2*waveSpeed);
    projLeftMatrix[1][2] = normal[0];
    projLeftMatrix[2][2] = normal[1]/(2*waveSpeed);
}

void ModelEquations::setProjEigenvalues(double height, const double *velocity, const double *normal, double *lambda) {

    double waveSpeed = sqrt(g0*height);

    lambda[0] = velocity[0]*normal[0];
    lambda[1] = velocity[0]*normal[0];
    lambda[2] = velocity[0]*normal[0];

    lambda[0] += velocity[1]*normal[1] - waveSpeed;
    lambda[1] += velocity[1]*normal[1];
    lambda[2] += velocity[1]*normal[1] + waveSpeed;
}