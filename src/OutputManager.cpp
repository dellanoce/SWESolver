#include <OutputManager.h>

void OutputManager::printSolutionVTK(Grid *mesh, Solution *sol, double time) {

    string rootName   = "time_";
    string pathFolder = "./VTK/";
    string timeStr;

    timeStr = to_string(time);
    timeStr.erase(remove(timeStr.begin(), timeStr.end(), '.'), timeStr.end());

    pathFolder.append(rootName + timeStr + ".vtk");

    ofstream SOLUTION(pathFolder.c_str());
    SOLUTION.precision(6);

    if (SOLUTION.is_open()) {

        // Header
        SOLUTION << "# vtk DataFile Version 3.0" << endl;
        SOLUTION << "vtk output" << endl;
        SOLUTION << "ASCII" << endl;
        SOLUTION << "DATASET STRUCTURED_GRID" << endl;
        SOLUTION << "DIMENSIONS " << mesh->nNodesX << " " << mesh->nNodesY << " " << 1 << endl;
        SOLUTION << "POINTS " << mesh->nNodes << " double" << endl;

        // Nodes
        for (int i = 0; i < mesh->nNodes; ++i) {
            for (int j = 0; j < 2; ++j) {
                SOLUTION << scientific << mesh->nodes[i][j] << "\t";
            }
            SOLUTION << scientific << "0.0" << "\t";
        }
        SOLUTION << endl;

        // Height
        SOLUTION << "POINT_DATA " << mesh->nNodes << endl << endl;
        SOLUTION << "SCALARS Height double 1" << endl;
        SOLUTION << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < mesh->nNodes; ++i) {
            SOLUTION << scientific << sol->conservative1(i) << "\t";
        }
        SOLUTION << endl;

        // Momentum
        SOLUTION << "VECTORS Momentum double" << endl;
        for (int i = 0; i < mesh->nNodes; ++i) {
            SOLUTION << scientific << sol->conservative2(i) << "\t" << sol->conservative3(i) << "\t";
            SOLUTION << scientific << "0.0" << "\t";
        }

        SOLUTION.close();
    } else
        swe::printError("Unable to write " + pathFolder + ".\n"
                        "->void OutputManager::printSolutionVTK(Grid *mesh, Solution *sol, double time)");
}

void OutputManager::printSolutionCsv(Grid *mesh, Solution *sol, double time) {

    string rootName   = "solution_";
    string pathFolder = "./solution_data/";
    string timeStr;

    timeStr = to_string(time);
    timeStr.erase(remove(timeStr.begin(), timeStr.end(), '.'), timeStr.end());

    pathFolder.append(rootName + timeStr + ".csv");

    ofstream SOLUTION(pathFolder.c_str());
    SOLUTION.precision(6);

    if (SOLUTION.is_open()) {
    
        SOLUTION << "Coordinate_x[m],Coordinate_y[m],Coordinate_z[m],Conservative_1[m],"
                    "Conservative_2[m^2/2],Conservative_3[m^2/s], \n";
                    
        for (int i = 0; i < mesh->nNodes; ++i) {
            SOLUTION << scientific << mesh->nodes[i][0] << ",",
            SOLUTION << scientific << mesh->nodes[i][1] << ",",
            SOLUTION << scientific << "0.0" << ",",
            SOLUTION << scientific << sol->conservative1(i) << ",",
            SOLUTION << scientific << sol->conservative2(i) << ",",
            SOLUTION << scientific << sol->conservative3(i) << endl;
        }
    } else
        swe::printError("Unable to write " + pathFolder + ".\n"
                        "->void OutputManager::printSolutionCsv(Grid *mesh, Solution *sol, double time)");
}