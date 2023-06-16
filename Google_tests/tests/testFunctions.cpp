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
    ASSERT_NEAR(trimAngles(db1, db2, V, h, 0).alpha_trim, 2.36, tol);
}

TEST(TrimTest, CanComputeDeltaTrim) {
    AeroDB db1, db2;
    db1 = readData("dba_100.ini", true);
    db2 = readData("dba_1000.ini", true);
    double V = 15;
    double h = 100;
    double tol = 0.05;
    ASSERT_NEAR(trimAngles(db1, db2, V, h, 0).deltae_trim, -2.16, tol);
}

TEST(TrimTest, CanConvertRpm) {
    double rpmMin = 3600;
    double rpmMax = 30000;
    double rpmAvg = 8300;
    ASSERT_NEAR(getRpm(0.1, rpmMin, rpmMax), 3600, 1);
    ASSERT_NEAR(getRpm(1, rpmMin, rpmMax), 30000, 1);
    ASSERT_NEAR(getRpm(0.229545, rpmMin, rpmMax), 7400, 1);

}

TEST(IntegrateMotion, CanMaintainTrim) {
    double V = 15;
    double h = 100;
    double vecRef[12] = {14.9873, 0, 0.617672, 0, 0, 0, 0, 0.0411898, 0, 100, 0, 0};
    double vecTrim[12] = {14.9873, 0, 0.617672, 0, 0, 0, 0, 0.0411898, 0, 100, 0, 0};
    double vecComTrim[4] = {0, -0.0376, 0, 0.260227};
    AeroDB db1, db2, db3, db4, DB1, DB2;
    db1 = readData("dba.ini", true);
    db2 = readData("dba_100.ini", true);
    db3 = readData("dba_1000.ini", true);
    db4 = readData("dba_2000.ini", true);
    BatteryDB bat0;
    EngineDB en0; // create dba objects from struct type
    PropDB prop0; // create dba objects from struct type
    // read not aerodynamic database
    bat0 = readBat("battery.ini"); // Open database, read it and save data to struct of type Not_Aerodb
    en0 = readEn("engine.ini"); // Open database, read it and save data to struct of type Not_Aerodb
    prop0 = readProp("propeller.ini"); // Open database, read it and save data to struct of type PropDB

    std::ofstream outputSim("testIntegrateMotion.txt");
    std::ofstream loggerRemainders("testIntegrateMotion_logRemainders.txt");
    std::ofstream loggerAcceleration("testIntegrateMotion_logAcceleration.txt");
    double stateMinusOne[12] = {0}; // i-1
    double rpm = getRpm(vecComTrim[3], en0.laps_min, en0.laps_max);//d_th * RPM_MAX
    for (int j = 0; j < 12; j++) {
        stateMinusOne[j] = vecTrim[j]; // during trim the (i-1)th state is the same as the trim. We assume the
        // simulation starts at the last trim state (i=0) and before that we assume infinite trim states take place
    }
    getAerodynamicDbWithAltitude(h, DB1, DB2, db1, db2, db3, db4);

    for (int i = 1; i < 100000; i++) {
        h = vecTrim[9];
        double* newStatesPointer = integrateEquationsOfMotion(DB1, DB2, en0, prop0, rpm, vecTrim, vecComTrim, stateMinusOne, 0.01, loggerRemainders, loggerAcceleration);
        double newStates[12] = {0};
        for (int j = 0; j < 12; j++) {
            newStates[j] = newStatesPointer[j]; // save the recently calculated state vector in newStates
            stateMinusOne[j] = vecTrim[j]; // save the OLD initial condition as the (i-1)th step
            vecTrim[j] = newStates[j]; // update the initial condition vector with the new states for the next loop iteration
        } // assign values to variable
        delete[] newStatesPointer; // delete pointer to avoid memory leak
        EXPECT_NEAR(vecTrim[0], vecRef[0], 1);
        EXPECT_NEAR(vecTrim[9], vecRef[9], 10);

    }

}