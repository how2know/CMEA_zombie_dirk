#include "writer.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <iostream>


#include "dirksolver.hpp"



int main(int argc, char** argv) {

    double T = 100;
    int N = 10;

    if (argc > 1)  {
        // Read N from command line
        // We use atof because we want to allow things like 1e7
        N = int(atof(argv[1]));
    }

    std::vector<std::vector<double> > u(3);
    std::vector<double> time(N + 1, 0);
    u[0].resize(N + 1, 0);
    u[1].resize(N + 1, 0);
    u[2].resize(N + 1, 0);

    u[0][0] = 500;
    u[1][0] = 0;
    u[2][0] = 0;

    ZombieOutbreak outbreak;

    DIRKSolver dirk(outbreak);
    dirk.solve(u, time, T, N);

    writeToFile("S.txt", u[0]);
    writeToFile("Z.txt", u[1]);
    writeToFile("R.txt", u[2]);
    writeToFile("time.txt", time);
}




