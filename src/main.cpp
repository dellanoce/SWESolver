#include "Solver.h"

int main(int argc, char *argv[]) {

    double startSim = MPI_Wtime();
    double totalSimTime, meshTime, solTime;
    int worldSize, currentRank;
    string configFileName;

    // --- MPI init --- //
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);

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
    double startMesh = MPI_Wtime();
    Grid grid(config);
    meshTime = MPI_Wtime() - startMesh;
    cout << " Mesh created in " << meshTime << " seconds" << endl << endl;

    cout << " Creating solution structure..." << endl;

    // --- Create solution object --- //
    double startSol = MPI_Wtime();
    Solution w0(&grid);
    solTime = MPI_Wtime() - startSol;
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
    totalSimTime = MPI_Wtime() - startSim;
    cout << " Simulation ended in " << totalSimTime << " seconds on " << worldSize << " cores" << endl;

    MPI_Finalize();

    return EXIT_SUCCESS;
}
