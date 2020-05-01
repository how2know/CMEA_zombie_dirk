#pragma once
#include <Eigen/Core>
#include "zombieoutbreak.hpp"
#include <vector>

class ForwardEuler {
public:
    ForwardEuler(const ZombieOutbreak& zombieOutbreak_) :
        zombieOutbreak(zombieOutbreak_) {
        // empty
    }

    ///
    /// Compute N timesteps of Forward-Euler
    /// @param[in/out] should be a vector of size 3, where each
    ///                component is a vector of size N+1. u[i][0]
    ///                should have the initial value stored before
    ///                calling the funtion
    ///
    /// @param[out] time should be of length N+1
    ///
    /// @param[in] T the final time
    /// @param[in] N the number of timesteps
    /// @param[in] k1 the constant k1
    /// @param[in] k2 the constant k2
    ///
    void solve(std::vector<std::vector<double> >& u, std::vector<double>& time, double T, int N) {

        const double dt = T / N;
        Eigen::Vector3d Fu;

        for (int i = 1; i < N + 1; ++i) {

            Eigen::Vector3d uPrevious(u[0][i - 1], u[1][i - 1], u[2][i - 1]);
            zombieOutbreak.computeF(Fu, time[i - 1], uPrevious);
            Eigen::Vector3d uNext = uPrevious + dt * Fu;

            u[0][i] = uNext[0];
            u[1][i] = uNext[1];
            u[2][i] = uNext[2];

            time[i] = time[i - 1] + dt;
        }
    }

private:
    ZombieOutbreak zombieOutbreak;
}; // end ForwardEuler class
