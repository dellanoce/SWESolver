#include "Solver.h"

int main(int argc, char *argv[]) {

    clock_t startSim = clock();
    double totalSimTime, meshTime, solTime;
    string configFileName;

    cout << "|------------------------------------------------------------------------------|" << endl;
    cout << "|                                                                              |" << endl;
    cout << "| 2D shallow water solver equations solver                                     |" << endl;
    cout << "|                                                                              |" << endl;
    cout << "|------------------------------------------------------------------------------|" << endl;
    cout << endl;

    // --- Check inputs --- //
    if (argc < 2)
        printError("No configuration file specified.");
    else
        configFileName = argv[1];

    // --- Parse configuration file --- //
    Config* config;
    config = new Config(configFileName);

    cout << "|----------------------------- RUN CFD SIMULATION -----------------------------|" << endl;
    cout << endl;

    cout << " Creating uniform structured mesh..." << endl;

    // -- Create uniform mesh --- //
    clock_t startMesh = clock();
    Grid grid(config);
    meshTime = (double) (clock() - startMesh)/CLOCKS_PER_SEC;
    cout << " Mesh created in " << meshTime << " seconds" << endl << endl;

    cout << " Creating solution structure..." << endl;

    // --- Create solution object --- //
    clock_t startSol = clock();
    Solution w0(&grid);
    solTime = (double) (clock() - startSol)/CLOCKS_PER_SEC;
    cout << " Solution structure created in " << solTime << " seconds" << endl << endl;

    cout << " Creating solver structure..." << endl << endl;

    // --- Instantiate the Godunov solver class --- //
    Solver god(config, &grid);

    cout << " Setting initial condition..." << endl;

    // --- Set initial solution --- //
    god.setInitialSolution(&grid, w0);
    swe::setupFolder("./solution_data", CLEAN_FOLDER);
    swe::setupFolder("./VTK", CLEAN_FOLDER);
    OutputManager::printSolutionVTK(&grid, &w0, 0.00000);
    OutputManager::printSolutionCsv(&grid, &w0, 0.00000);

    cout << endl;
    cout << " Start CFD computation..." << endl << endl;

    // --- Run simulation --- //
    god.runGodnuovSolver(&grid);
    totalSimTime = (double) (clock() - startSim)/CLOCKS_PER_SEC;
    cout << " Simulation ended in " << totalSimTime << " seconds" << endl;

    return EXIT_SUCCESS;
}
