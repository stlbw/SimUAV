#include "../lib/googletest/include/gtest/gtest.h"
#include "../../declaredFun.h"
#include "vector"

TEST(InterpolationTest, CanDoLinearInterpolation) {
    vector<double> testVec = {0, 1, 2, 3, 4, 5};
    vector<double> testCoef = {1, 1.2, 1.4, 1.6, 1.8, 2.0};
    double valBetween = 2.5;
    double valQuasiExactUp = 1.001;
    double valQuasiExactDown = 0.998;
    double tol = 1e-3;
    EXPECT_NEAR(linearInterpolation(testVec, testCoef, valBetween), 1.5, tol);
    EXPECT_NEAR(linearInterpolation(testVec, testCoef, valQuasiExactDown), 1.2, tol);
    EXPECT_NEAR(linearInterpolation(testVec, testCoef, valQuasiExactUp), 1.2, tol);


}

TEST(TrimTest, CanComputeAlphaTrim) {
    AeroDB db;
    db = readData("dba_100.ini", true);
    //double V = 15;
    //double h = 100;
    double tol = 1e-3;
    ASSERT_NEAR(trim1(db).alpha_trim, 2.36, tol);
}

TEST(TrimTest, CanComputeDeltaTrim) {
    AeroDB db;
    db = readData("dba_100.ini", true);
    //double V = 15;
    //double h = 100;
    double tol = 1e-3;
    ASSERT_NEAR(trim1(db).deltae_trim, -2.16, tol);
}