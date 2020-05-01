#pragma once
#include <Eigen/Core>
#include <vector>
#include "zombieoutbreak.hpp"

class DIRKSolver {
public:
    DIRKSolver(const ZombieOutbreak& zombieOutbreak_)
        : zombieOutbreak(zombieOutbreak_) {

    }

    ///
    /// Evaluates function G1(y) from task b)
    /// @param[out] G evaluation of function
    /// @param[in] y input to G
    /// @param[in] tn current time
    /// @param[in] Un computed value of U at time tn
    /// @param[in] dt timestep
    ///

    //----------------G1Start----------------

    void computeG1(Eigen::Vector3d& G, Eigen::Vector3d y, double tn, Eigen::Vector3d Un, double dt) {

// (write your solution here)

        /// Start of my solution ///

        // compute F(t+mu*dt, y1)
        Eigen::Vector3d F;
        zombieOutbreak.computeF(F, tn + mu * dt, y);

        // compute G1
        G = Un + dt * mu * F - y;

        /// End of my solution ///
    }

    //----------------G1End----------------

    ///
    /// Evaluates function G2(y1, y) from task b)
    /// @param[out] G evaluation of function
    /// @param[in] y input to G
    /// @param[in] tn current time
    /// @param[in] Un computed value of U at time tn
    /// @param[in] dt timestep
    /// @param[in] y1 computed value for first intermediate RK value
    ///

    void computeG2(Eigen::Vector3d& G, Eigen::Vector3d y, double tn, Eigen::Vector3d Un, double dt, Eigen::Vector3d y1) {

// (write your solution here)

        /// Start of my solution ///

        // compute F(t+mu*dt, y1)
        Eigen::Vector3d F1;
        zombieOutbreak.computeF(F1, tn + mu * dt, y1);

        // compute F(t+(mu-nu)*dt, y2)
        Eigen::Vector3d F2;
        zombieOutbreak.computeF(F2, tn + (mu-nu) * dt, y);

        // compute G2
        G = Un + dt * (-nu * F1 + mu * F2) - y;

        /// End of my solution ///
    }

    ///
    /// Find a solution to JG1*x = -G1 with Newton's method
    /// @param[out] x solution to the system
    /// @param[in] Un computed value of U at time tn
    /// @param[in] dt timestep
    /// @param[in] tn current time
    /// @param[in] tolerance if Newton increment smaller, successfully converged
    /// @param[in] maxIterations  max Newton iterations to try before failing
    ///

    void newtonSolveY1(Eigen::Vector3d& u, Eigen::Vector3d Un, double dt, double tn, double tolerance, int maxIterations) {

        Eigen::Vector3d RHSG1, Fu, x;
        Eigen::Matrix3d JG1, JFu;
        u = Un;
	
	   // Write your loop for the Newton solver here.
	   // You will need to use zombieOutbreak.computeF(...) 
	   // and zombieOutbreak.computeJF(...)
// (write your solution here)

        /// Start of my solution ///

        for (int i = 0; i < maxIterations; i++) {
/*
   *** I tried here to do something to use the parameter that are declared above but I prefered to do it in my way without using them ***

            // compute F(tn+mu*dt, y1)
            zombieOutbreak.computeF(Fu, tn + mu * dt, u);

            // compute -G(y1)
            RHSG1 = - (Un + dt * mu * Fu - u);

            // compute the residual and check if we are better than the tolerance
            Eigen::Vector3d res = - RHSG1;
            if (res.norm() <= tolerance){
                return;
            }
*/

            // compute G1
            Eigen::Vector3d G1;
            computeG1(G1, u, tn, Un, dt);

            // check if we are better than the tolerance
            if (G1.norm() <= tolerance){
                return;
            }

            // compute DF(tn+mu*dt, y1)
            zombieOutbreak.computeJF(JFu, tn + mu * dt, u);

            // compute DG1
            JG1 = dt * mu * JFu - Eigen::Matrix3d::Identity();

            // solve DG1 * x = -G1
            x = JG1.lu().solve(-G1);

            // compute the next u
            u += x;
        }

        /// End of my solution ///

        // If we reach this point, something wrong happened.
        throw std::runtime_error("Did not reach tolerance in Newton iteration");
    }

    ///
    /// Find a solution to JG2*x = -G2 with Newton's method
    /// @param[out] x solution to the system
    /// @param[in] Un computed value of U at time tn
    /// @param[in] y1 previous intermediate value for RK method
    /// @param[in] dt timestep
    /// @param[in] tn current time
    /// @param[in] tolerance if Newton increment smaller, successfully converged
    /// @param[in] maxIterations  max Newton iterations to try before failing
    ///

    //----------------NewtonG2Start----------------

    void newtonSolveY2(Eigen::Vector3d& v, Eigen::Vector3d Un, Eigen::Vector3d y1, double dt, double tn, double tolerance, int maxIterations) {

	   // Use newtonSolveY1 as a model for this
// (write your solution here)

        /// Start of my solution ///

        // initialize the vectors and matrices
        Eigen::Vector3d G2, x;
        Eigen::Matrix3d DF2, DG2;

        for (int i = 0; i <= maxIterations; i++) {

            // compute G2
            computeG2(G2, v, tn, Un, dt, y1);

            // check if we are better than the tolerance
            if (G2.norm() <= tolerance) {
                return;
            }

            // compute DF(tn+(mu-nu)*dt, y2)
            zombieOutbreak.computeJF(DF2, tn + (mu - nu) * dt, v);

            // compute DG1
            DG2 = dt * mu * DF2 - Eigen::Matrix3d::Identity();

            // solve DG1 * x = -G1
            x = DG2.lu().solve(-G2);

            // compute the next u
            v += x;
        }

        /// End of my solution ///
    }

    //----------------NewtonG2End----------------



    ///
    /// Compute N timesteps of DIRK(2,3)
    /// @param[in/out] u should be a vector of size 3, where each
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

    //----------------DirkStart----------------

    void solve(std::vector<std::vector<double> >& u, std::vector<double>& time, double T, int N) {

        const double dt = T / N;

	// Your main loop goes here. At iteration n,
	// 1) Find Y_1 with newtonSolveY1 (resp. Y2)
	// 2) Compute U^{n+1} with F(Y1), F(Y2)
	// 3) Write the values at u[...][n] and time[n]

// (write your solution here)

        /// Start of my solution ///

        // declare the tolerance and the maximum number of iteration for Newton
        double tol = 1e-10;
        double max = 100;

        // initialize some vectors
        Eigen::Vector3d y1, y2, vector, F1, F2;

        for (int i = 0; i < N; i++) {

            // store the initial value
            vector << u[0][i], u[1][i], u[2][i];

            // find y1 and y2 with Newton
            newtonSolveY1(y1, vector, dt, time[i], tol, max);
            newtonSolveY2(y2, vector, y1, dt, time[i], tol, max);

            // compute F(tn+mu*dt, y1) and F(tn+(mu-nu)*dt, y2)
            zombieOutbreak.computeF(F1, time[i] + mu * dt, y1);
            zombieOutbreak.computeF(F2, time[i] + (mu - nu) * dt, y2);

            // compute U_n+1
            vector = vector + dt * (0.5 * F1 + 0.5 * F2);

            // store the values in the vectors
            u[0][i+1] = vector[0];
            u[1][i+1] = vector[1];
            u[2][i+1] = vector[2];
            time[i+1] = time[i] + dt;
        }

        /// End of my solution ///

    }

    //----------------DirkEnd----------------

private:

    const double mu = 0.5 + 0.5 / sqrt(3);
    const double nu = 1. / sqrt(3);

    ZombieOutbreak zombieOutbreak;

}; // end class DIRKSolver
