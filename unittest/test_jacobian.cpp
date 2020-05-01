#include <gtest/gtest.h>
#include <Eigen/Core>
#include "../zombieoutbreak.hpp"

#define NEAR_TOLERANCE 1e-7
TEST(TestJacobian, FirstColCheck) {
    ZombieOutbreak outbreak;
    Eigen::Matrix3d JFu;
    Eigen::Vector3d u(0,1,0);
    outbreak.computeJF(JFu,1,u);
    
    ASSERT_NEAR(JFu(0,0), -0.17499, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,0), 0.025, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,0), 0.15002, NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(0,1), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,1), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,1), 0, NEAR_TOLERANCE);
    
    ASSERT_NEAR(JFu(0,2), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,2), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,2), 0, NEAR_TOLERANCE);
}

TEST(TestJacobian, SecondColCheck) {
    ZombieOutbreak outbreak;
    Eigen::Matrix3d JFu;
    Eigen::Vector3d u(1,0,0);
    outbreak.computeJF(JFu,1,u);
    
    ASSERT_NEAR(JFu(0,0), 1e-5, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,0), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,0), 2e-5, NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(0,1), -0.175, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,1), 0.025, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,1), 0.15, NEAR_TOLERANCE);
    
    ASSERT_NEAR(JFu(0,2), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,2), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,2), 0, NEAR_TOLERANCE);
}


TEST(TestJacobian, ThirdColCheck) {
    ZombieOutbreak outbreak;
    Eigen::Matrix3d JFu;
    Eigen::Vector3d u(0,0,0);
    outbreak.computeJF(JFu,10,u);
    
    ASSERT_NEAR(JFu(0,0), 1e-5, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,0), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,0), 2e-5, NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(0,1), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,1), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,1), 0, NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(0,2), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,2), 0.086, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,2), -0.086, NEAR_TOLERANCE);
}
