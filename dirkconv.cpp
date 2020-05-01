#include "writer.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <iostream>
#include <chrono>
#include "zombieoutbreak.hpp"
#include "dirksolver.hpp"


//----------------mainBegin----------------

int main(int argc, char** argv) {

    double T = 101;
    ZombieOutbreak outbreak(0, 0, 0, 0.03, 0.02);
    std::vector<double> u0(3);
    u0[0] = 500;
    u0[1] = 0;
    u0[2] = 0;
    
    // Compute the exact solution for the parameters above
    std::vector<double> exact = outbreak.computeExactNoZombies(T, u0[0]);

    // Initialize solver object for the parameters above
    DIRKSolver dirkSolver(outbreak);

    int minExp = 0;
    int maxExp = 12;
    int countExponents = maxExp - minExp +1;
    std::vector<double> numbers(countExponents);
    std::vector<double> walltimes(countExponents);
    std::vector<double> errors(countExponents);

// (write your solution here)

    /// Start of my solution ///

    for (int i = 0; i <= maxExp; i++) {

        // compute the new number of iteration for each step
        int N = 200 * std::pow(2, i);

        // declare vectors u and t and their initial value to use them dirksolver.solve
        std::vector<std::vector<double> > u(3);
        std::vector<double> t(N + 1, 0);
        u[0].resize(N + 1, 0);
        u[1].resize(N + 1, 0);
        u[2].resize(N + 1, 0);
        u[0][0] = u0[0];
        u[1][0] = u0[1];
        u[2][0] = u0[2];

        // begin to count the time
        auto begin = std::chrono::high_resolution_clock::now();

        // solve for the approximate values
        dirkSolver.solve(u, t, T, N);

        // stop to count the time and calculate the time elapsed
        auto end = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

        // compute the errors and store them
        errors[i] = std::abs( exact[0] - u[0][N] ) + std::abs( exact[1] - u[1][N] ) + std::abs( exact[2] - u[2][N] );

        // store the number of iterations
        numbers[i] = N;

        // store the elapsed time
        walltimes[i] = time;

        // print the approximate solution
        std::cout << "Approx. solution: " << u[0].back() << " " << u[1].back() << " " << u[2].back() << " for N = " << N << " in " << time/1e6 << " ms" << std::endl;
    }

    /// End of my solution ///

    writeToFile("numbers.txt", numbers);
    writeToFile("errors.txt", errors);
    writeToFile("walltimes.txt", walltimes);

}

//----------------mainEnd----------------
