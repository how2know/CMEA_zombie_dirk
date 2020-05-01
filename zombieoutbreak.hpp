#pragma once
#include <functional>
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <iostream>

class ZombieOutbreak {
public:

    ///
    /// Evaluates right hand side of the system of ODEs
    ///    S' = Pi*S - beta*S*Z - delta*S
    ///    Z' = beta*S*Z + zeta*R - alpha*S*Z
    ///    R' = delta*S + alpha*S*Z - zeta*R
    /// @param[in/out] F evaluation of RHS
    /// @param[in] U Vector3d with [S,Z,R]
    /// @param[in] t current time
    ///

    //----------------FStart----------------

    void computeF(Eigen::Vector3d& F, double t, Eigen::Vector3d U) {

        double S = U[0];
        double Z = U[1];
        double R = U[2];

        F[0] = Pi * S - beta(t) * S * Z - delta * S;
        F[1] = beta(t) * S * Z + zeta(t) * R - alpha(t) * S * Z;
        F[2] = delta * S + alpha(t) * S * Z - zeta(t) * R;
    }

    //----------------FEnd---------------

    ///
    /// Computes Jacobian matrix of F
    /// @param[in/out] J Jacobian matrix
    /// @param[in] t time
    /// @param[in] U Vector3d with [S,Z,R]
    ///

    //----------------JFStart----------------

    void computeJF(Eigen::Matrix3d& J, double t, Eigen::Vector3d U) {

// (write your solution here)

        /// Start of my solution ///

        double S = U[0];
        double Z = U[1];
        double R = U[2];

        J << Pi - beta(t)*Z - delta,   -beta(t)*S,             0,
             beta(t)*Z - alpha(t)*Z,   beta(t)*S-alpha(t)*S,   zeta(t),
             delta + alpha(t)*Z,       alpha(t)*S,             -zeta(t);

        /// End of my solution ///

    }

    //----------------JFEnd----------------

    double alpha(double t) {
        return a + b * (t / (1 + t));
    }

    double beta(double t) {
        return 2*a - (b / 2) * t / (1 + t);
    }

    double zeta(double t) {
        return t <= 5 ? 0 : c;
    }

    ZombieOutbreak() = default;

    ZombieOutbreak(double a_, double b_, double c_, double Pi_, double delta_) :
        a(a_), b(b_), c(c_), Pi(Pi_), delta(delta_) {
        // Empty
    }


    std::vector<double> computeExactNoZombies(double t, double S0) {
        std::vector<double> result(3);
        result[0] = S0 * std::exp( (Pi - delta) * t);
        result[1] = 0.;
        result[2] = (delta * S0 / (Pi - delta)) * ( std::exp( (Pi - delta) * t) - 1. );
        std::cout.precision(17);
        std::cout << "Exact solution: " << result[0] << " " << result[1]
            << " " << result[2] << std::endl;
        return result;
    }


private:
    //! Birth rate
    const double Pi = 3e-5;

    //! Death rate
    const double delta = 2e-5;

    // The magic Garland constants:
    const double a = 0.1;
    const double b = 0.1;
    const double c = 0.086;


};
