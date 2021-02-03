#include "Solver.h"

int main(int argc, char *argv[]) {

    double startSim = MPI_Wtime();
    double totalSimTime, meshTime, solTime, startMesh, startSol;
    int worldSize, worldRank;
    string configFileName;
    stringstream ssMessage;

    // --- MPI init --- //
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // --- Check inputs --- //
    if (argc < 2)
        printError("No configuration file specified.");
    else
        configFileName = argv[1];

    // --- Parse configuration file --- //
    Config* config;
    config = new Config(configFileName);

    ssMessage << "|----------------------------- RUN CFD SIMULATION -----------------------------|" << endl;
    printMessage(ssMessage.str());

    // -- Create uniform mesh --- //
    printMessage(" Creating uniform structured mesh...");
    startMesh = MPI_Wtime();

    Grid grid(config);
    meshTime = MPI_Wtime() - startMesh;
    ssMessage.str("");
    ssMessage << " Mesh created in " << meshTime << " seconds" << endl;
    printMessage(ssMessage.str());

    // --- Create solution object --- //
    printMessage(" Creating solution structure...");
    startSol = MPI_Wtime();

    Solution w0(&grid);
    solTime = MPI_Wtime() - startSol;
    ssMessage.str("");
    ssMessage << " Solution structure created in " << solTime << " seconds" << endl;
    printMessage(ssMessage.str());

    // --- Instantiate the Godunov solver class --- //
    ssMessage.str("");
    ssMessage << " Creating solver structure..." << endl;
    printMessage(ssMessage.str());

    Solver god(config, &grid);

    // --- Set initial solution --- //
    printMessage(" Setting initial condition...");

    god.setInitialSolution(&grid, w0);
    if (worldRank == MASTER_NODE) {
        swe::setupFolder("./solution_data", CLEAN_FOLDER);
        swe::setupFolder("./VTK", CLEAN_FOLDER);
    }
    OutputManager::printSolutionVTK(&grid, &w0, 0.00000);
    OutputManager::printSolutionCsv(&grid, &w0, 0.00000);

    // --- Run simulation --- //
    ssMessage.str("");
    ssMessage << endl << " Start CFD computation..." << endl;
    printMessage(ssMessage.str());

    god.runGodnuovSolver(&grid);
    totalSimTime = MPI_Wtime() - startSim;
    ssMessage.str("");
    ssMessage << " Simulation ended in " << totalSimTime << " seconds on " << worldSize << " cores" << endl;
    printMessage(ssMessage.str());

    MPI_Finalize();

    return EXIT_SUCCESS;
}
