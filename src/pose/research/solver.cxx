#include <cstddef>
#include <string>
#include <iostream>
#include <stdexcept>

#include "andres/ilp/gurobi.hxx"

#include "io-hdf5.hxx"

#include "pose/solver.hxx"

typedef andres::ilp::Gurobi<> IlpSolver;
typedef pose::Solver<IlpSolver> Solver;

int main(int argc, char** argv) { 
    int timeLimit = 36000;
    if(argc < 4) {
        std::cerr << "solver <problem.h5> <solution.h5> <s/m>" << std::endl;
        return 1;
    }
    else if(argc == 5) {
        timeLimit = atoi(argv[4]);
    }
  
    const std::string problemFileName = argv[1];
    const std::string solutionFileName = argv[2];
    const std::string solverIndex = argv[3];
    
    Problem problem;
    loadProblem(problemFileName, problem);

    Solution solution;
    Solver solver(problem, solution, solverIndex, timeLimit);
    saveSolution(solutionFileName, solution);

    return 0;
}
