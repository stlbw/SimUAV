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
    EXPECT_NEAR(linearInterpolationAlpha(testVec, testCoef, valBetween), 1.5, tol);
    EXPECT_NEAR(linearInterpolationAlpha(testVec, testCoef, valQuasiExactDown), 1.2, tol);
    EXPECT_NEAR(linearInterpolationAlpha(testVec, testCoef, valQuasiExactUp), 1.2, tol);


}

TEST(InterpolationTest, CanDoLinearInterpolationWithAltitude) {
    vector<double> testVec = {0, 1, 2, 3, 4, 5};
    vector<double> testCoef1 = {2, 4, 6, 8, 10, 12};
    vector<double> testCoef2 = {4, 8, 12, 16, 20, 24};
    double h1 = 50;
    double h2 = 1000;
    double h3 = 1400;
    double valBetween = 2.5;
    double valQuasiExactUp = 1.001;
    double valQuasiExactDown = 0.998;
    double tol = 1e-2;
    EXPECT_NEAR(linearInterpolation(testVec, testCoef1, testCoef2, valBetween, h1), 10.5, tol);
    EXPECT_NEAR(linearInterpolation(testVec, testCoef1, testCoef2, valQuasiExactUp, h1), 6, tol);
    EXPECT_NEAR(linearInterpolation(testVec, testCoef1, testCoef2, valQuasiExactDown, h1), 6, tol);

    EXPECT_NEAR(linearInterpolation(testVec, testCoef1, testCoef2, valBetween, h2), 14, tol);
    EXPECT_NEAR(linearInterpolation(testVec, testCoef1, testCoef2, valQuasiExactUp, h2), 8, tol);
    EXPECT_NEAR(linearInterpolation(testVec, testCoef1, testCoef2, valQuasiExactDown, h2), 8, tol);

    EXPECT_NEAR(linearInterpolation(testVec, testCoef1, testCoef2, valBetween, h3), 9.8, tol);
    EXPECT_NEAR(linearInterpolation(testVec, testCoef1, testCoef2, valQuasiExactUp, h3), 5.6, tol);
    EXPECT_NEAR(linearInterpolation(testVec, testCoef1, testCoef2, valQuasiExactDown, h3), 5.6, tol);



}

TEST(TrimTest, CanComputeAlphaTrim) {
    AeroDB db1, db2;
    db1 = readData("dba_100.ini", true);
    db2 = readData("dba_1000.ini", true);
    double V = 15;
    double h = 100;
    double tol = 0.05;
    ASSERT_NEAR(trimAngles(db1, db2, V, h).alpha_trim, 2.36, tol);
}

TEST(TrimTest, CanComputeDeltaTrim) {
    AeroDB db1, db2;
    db1 = readData("dba_100.ini", true);
    db2 = readData("dba_1000.ini", true);
    double V = 15;
    double h = 100;
    double tol = 0.05;
    ASSERT_NEAR(trimAngles(db1, db2, V, h).deltae_trim, -2.16, tol);
}