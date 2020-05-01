#include <gtest/gtest.h>
#include <Eigen/Core>
#include "../dirksolver.hpp"
#include <iostream>

#define NEAR_TOLERANCE 1e-7

TEST(TestGFunctions, G1Check) {
    const double dt = 2*sqrt(3)/(1+sqrt(3));
    ZombieOutbreak outbreak(1,4,1,0,0);
    Eigen::Vector3d y(1,2,3), u(4,8,15), G;

    DIRKSolver solver(outbreak);
    solver.computeG1(G, y, 0, u, dt);

    ASSERT_NEAR(G(0), 1, NEAR_TOLERANCE);
    ASSERT_NEAR(G(1), 2, NEAR_TOLERANCE);
    ASSERT_NEAR(G(2), 18, NEAR_TOLERANCE);
}


TEST(TestGFunctions, G2Check) {
    const double dt = 2*sqrt(3)/(1+sqrt(3));
    ZombieOutbreak outbreak(1,4,0,0,0);
    Eigen::Vector3d y1(1,2,3), u(4,8,15), y(0,0,42), G;
    const double C = 1.0/(1+sqrt(3));

    DIRKSolver solver(outbreak);
    solver.computeG2(G, y, 0, u, dt, y1);

    ASSERT_NEAR(G(0), 4 + 4*C, NEAR_TOLERANCE);
    ASSERT_NEAR(G(1), 8 + 8*C, NEAR_TOLERANCE);
    ASSERT_NEAR(G(2), 15 - 12*C - 42, NEAR_TOLERANCE);
}


TEST(TestNewtonMethod, NewtonY1Check) {
    ZombieOutbreak outbreak;
    Eigen::Vector3d u(4,8,15), uOut;
    DIRKSolver solver(outbreak);
	solver.newtonSolveY1(uOut, u, 0.01, 1, 1e-10, 100);

    ASSERT_NEAR(uOut(0), 3.95630776094, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(1), 8.00617170496, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(2), 15.0375214702, NEAR_TOLERANCE);
}


TEST(TestNewtonMethod, NewtonY2Check) {
    ZombieOutbreak outbreak;
    Eigen::Vector3d u(4,8,15), y(16,23,42), uOut;
    DIRKSolver solver(outbreak);
	solver.newtonSolveY2(uOut, u, y, 0.01, 1, 1e-10, 100);

    ASSERT_NEAR(uOut(0), 4.32413965602, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(1), 7.95426989836, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(2), 14.7215886974, NEAR_TOLERANCE);
}

TEST(TestDirkSolver, DirkSolverCheckSizes) {
    std::vector<std::vector<double> > u(3);
    std::vector<double> time(201, 0);
    u[0].resize(201, 0);
    u[1].resize(201, 0);
    u[2].resize(201, 0);
    u[0][0] = 0;
    u[1][0] = 0;
    u[2][0] = 0;
    ZombieOutbreak outbreak;
    DIRKSolver dirk(outbreak);
    dirk.solve(u, time, 1, 200);

    ASSERT_EQ(u[0].size(), 201);
    ASSERT_EQ(u[1].size(), 201);
    ASSERT_EQ(u[2].size(), 201);
    ASSERT_EQ(time.size(), 201);
}

TEST(TestDirkSolver, DirkSolverWroteLast) {
    std::vector<std::vector<double> > u(3);
    std::vector<double> time(201, 0);
    u[0].resize(201, 0);
    u[1].resize(201, 0);
    u[2].resize(201, 0);
    u[0][0] = 1;
    u[1][0] = 1;
    u[2][0] = 1;
    ZombieOutbreak outbreak;
    DIRKSolver dirk(outbreak);
    dirk.solve(u, time, 1, 200);

    EXPECT_NE(u[0][200], 0);
    EXPECT_NE(u[1][200], 0);
    EXPECT_NE(u[2][200], 0);
    ASSERT_NE(time[200], 0);
}


TEST(TestDirkSolver, DirkSolverCheck) {
    const double T = 42;
    const int N = 1138;
    std::vector<std::vector<double> > u(3);
    std::vector<double> time(N + 1, 0);
    u[0].resize(N + 1, 0);
    u[1].resize(N + 1, 0);
    u[2].resize(N + 1, 0);
    u[0][0] = 1123;
    u[1][0] = 6536;
    u[2][0] = 5321;
    ZombieOutbreak outbreak;
    DIRKSolver dirk(outbreak);
    dirk.solve(u, time, T, N);

    ASSERT_NEAR(u[0][N-2], 2.88144266549722163e-18, NEAR_TOLERANCE);
    ASSERT_NEAR(u[1][N-2], 12734.2615307803462, NEAR_TOLERANCE);
    ASSERT_NEAR(u[2][N-2], 245.738493192421203, NEAR_TOLERANCE);
}

