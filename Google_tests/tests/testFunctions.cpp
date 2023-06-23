#include "../lib/googletest/include/gtest/gtest.h"
#include "../../declaredFun.h"
#include "vector"

TEST(InterpolationTest, CanDoLinearInterpolation) {
    std::ofstream outputSim("../../Google_tests/tests/loggerOutput/CanDoLinearInterpolation.txt");
    if (!outputSim) {
        cout<< "Could not open file CanDoLinearInterpolation.txt"<<endl;
    }
    vector<double> testVec = {0, 1, 2, 3, 4, 5};
    vector<double> testCoef = {1, 1.2, 1.4, 1.6, 1.8, 2.0};
    double valBetween = 2.5;
    double valQuasiExactUp = 1.001;
    double valQuasiExactDown = 0.998;
    double tol = 1e-3;
    EXPECT_NEAR(linearInterpolationAlpha(testVec, testCoef, valBetween), 1.5, tol);
    EXPECT_NEAR(linearInterpolationAlpha(testVec, testCoef, valQuasiExactDown), 1.2, tol);
    EXPECT_NEAR(linearInterpolationAlpha(testVec, testCoef, valQuasiExactUp), 1.2, tol);
    outputSim << "TEST: CanDoLinearInterpolation" << endl;
    outputSim << "Test reference vector: [";
    for (int i = 1; i < testVec.size(); i++){
        outputSim << left << " "<< testVec[i]; //print to logger
    }
    outputSim << "] " << endl;
    outputSim << "Test coefficients to interpolate: [";
    for (int i = 1; i < testCoef.size(); i++){
        outputSim << left << " " << testCoef[i]; //print to logger
    }
    outputSim << "] " << endl;
    outputSim << "EXPECT_NEAR TEST (tolerance 1E-3):" << endl;
    outputSim << "  Reference value: " << valBetween << endl;
    outputSim << "  Expected value: " << 1.5 << endl;
    outputSim << "  Result: TRUE " << endl;
    outputSim << "EXPECT_NEAR TEST (tolerance 1E-3):" << endl;
    outputSim << "  Reference value: " << valQuasiExactDown << endl;
    outputSim << "  Expected value: " << 1.2 << endl;
    outputSim << "  Result: TRUE " << endl;
    outputSim << "EXPECT_NEAR TEST (tolerance 1E-3):" << endl;
    outputSim << "  Reference value: " << valQuasiExactUp << endl;
    outputSim << "  Expected value: " << 1.2 << endl;
    outputSim << "  Result: TRUE " << endl;





}

TEST(InterpolationTest, CanDoLinearInterpolationWithAltitude) {
    std::ofstream outputSim("../../Google_tests/tests/loggerOutput/CanDoLinearInterpolationWithAltitude.txt");
    if (!outputSim) {
        cout<< "Could not open file CanDoLinearInterpolationWithAltitude.txt"<<endl;
    }
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

    outputSim << "TEST: CanDoLinearInterpolation" << endl;
    outputSim << "Test reference vector: [";
    for (int i = 1; i < testVec.size(); i++){
        outputSim << left << " "<< testVec[i]; //print to logger
    }
    outputSim << "] " << endl;
    outputSim << "Test coefficients to interpolate database 1: [";
    for (int i = 1; i < testCoef1.size(); i++){
        outputSim << left << " " << testCoef1[i]; //print to logger
    }
    outputSim << "] " << endl;
    outputSim << "Test coefficients to interpolate database 2: [";
    for (int i = 1; i < testCoef2.size(); i++){
        outputSim << left << " " << testCoef2[i]; //print to logger
    }
    outputSim << "] " << endl;

    outputSim << "EXPECT_NEAR TEST (tolerance 1E-2):" << endl;
    outputSim << "  Altitude: " << h1 << endl;
    outputSim << "  Reference value: " << valBetween << endl;
    outputSim << "  Expected value: " << 10.5 << endl;
    outputSim << "  Result: TRUE " << endl;
    outputSim << "EXPECT_NEAR TEST (tolerance 1E-2):" << endl;
    outputSim << "  Altitude: " << h1 << endl;
    outputSim << "  Reference value: " << valQuasiExactUp << endl;
    outputSim << "  Expected value: " << 6 << endl;
    outputSim << "  Result: TRUE " << endl;
    outputSim << "EXPECT_NEAR TEST (tolerance 1E-2):" << endl;
    outputSim << "  Altitude: " << h1 << endl;
    outputSim << "  Reference value: " << valQuasiExactDown << endl;
    outputSim << "  Expected value: " << 6 << endl;
    outputSim << "  Result: TRUE " << endl;

    outputSim << "EXPECT_NEAR TEST (tolerance 1E-2):" << endl;
    outputSim << "  Altitude: " << h2 << endl;
    outputSim << "  Reference value: " << valBetween << endl;
    outputSim << "  Expected value: " << 14.0 << endl;
    outputSim << "  Result: TRUE " << endl;
    outputSim << "EXPECT_NEAR TEST (tolerance 1E-2):" << endl;
    outputSim << "  Altitude: " << h2 << endl;
    outputSim << "  Reference value: " << valQuasiExactUp << endl;
    outputSim << "  Expected value: " << 8 << endl;
    outputSim << "  Result: TRUE " << endl;
    outputSim << "EXPECT_NEAR TEST (tolerance 1E-2):" << endl;
    outputSim << "  Altitude: " << h2 << endl;
    outputSim << "  Reference value: " << valQuasiExactDown << endl;
    outputSim << "  Expected value: " << 8 << endl;
    outputSim << "  Result: TRUE " << endl;

    outputSim << "EXPECT_NEAR TEST (tolerance 1E-2):" << endl;
    outputSim << "  Altitude: " << h3 << endl;
    outputSim << "  Reference value: " << valBetween << endl;
    outputSim << "  Expected value: " << 9.8 << endl;
    outputSim << "  Result: TRUE " << endl;
    outputSim << "EXPECT_NEAR TEST (tolerance 1E-2):" << endl;
    outputSim << "  Altitude: " << h3 << endl;
    outputSim << "  Reference value: " << valQuasiExactUp << endl;
    outputSim << "  Expected value: " << 5.6 << endl;
    outputSim << "  Result: TRUE " << endl;
    outputSim << "EXPECT_NEAR TEST (tolerance 1E-2):" << endl;
    outputSim << "  Altitude: " << h3 << endl;
    outputSim << "  Reference value: " << valQuasiExactDown << endl;
    outputSim << "  Expected value: " << 5.6 << endl;
    outputSim << "  Result: TRUE " << endl;




}

TEST(TrimTest, CanComputeAlphaTrim) {
    std::ofstream outputSim("../../Google_tests/tests/loggerOutput/CanComputeAlphaTrim.txt");
    if (!outputSim) {
        cout<< "Could not open file CanComputeAlphaTrim.txt"<<endl;
    }
    AeroDB db1, db2;
    aero_condition atm;
    double rho = AtmosphereCalc(1, atm, 100);
    db1 = readDataTest("dba_100.ini", true);
    db2 = readDataTest("dba_1000.ini", true);
    double V = 15;
    double h = 100;
    double tol = 0.05;
    ASSERT_NEAR(trimAngles(db1, db2, V, h, 0, rho).alpha_trim, 2.36, tol);

    outputSim << "TEST: CanComputeAlphaTrim" << endl;
    outputSim << "Test reference databases dba_100.ini, dba_1000.ini" << endl;
    outputSim << "ASSERT_NEAR TEST (tolerance 0.05):" << endl;
    outputSim << "  Reference velocity [ m/s]: " << V << endl;
    outputSim << "  Reference altitude [m]:  " << h << endl;
    outputSim << "  Reference value alpha_trim [deg]:  " << 2.36 << endl;
    outputSim << "  Result: TRUE " << endl;
}

TEST(TrimTest, CanComputeDeltaTrim) {
    std::ofstream outputSim("../../Google_tests/tests/loggerOutput/CanComputeDeltaTrim.txt");
    if (!outputSim) {
        cout<< "Could not open file CanComputeDeltaTrim.txt"<<endl;
    }
    AeroDB db1, db2;
    aero_condition atm;
    double rho = AtmosphereCalc(1, atm, 100);
    db1 = readDataTest("dba_100.ini", true);
    db2 = readDataTest("dba_1000.ini", true);
    double V = 15;
    double h = 100;
    double tol = 0.05;
    ASSERT_NEAR(trimAngles(db1, db2, V, h, 0, rho).deltae_trim, -2.16, tol);
    outputSim << "TEST: CanComputeDeltaTrim" << endl;
    outputSim << "Test reference databases dba_100.ini, dba_1000.ini" << endl;
    outputSim << "ASSERT_NEAR TEST (tolerance 0.05):" << endl;
    outputSim << "  Reference velocity [ m/s]: " << V << endl;
    outputSim << "  Reference altitude [m]:  " << h << endl;
    outputSim << "  Reference value alpha_trim [deg]:  " << -2.16 << endl;
    outputSim << "  Result: TRUE " << endl;
}

TEST(TrimTest, CanConvertRpm) {
    std::ofstream outputSim("../../Google_tests/tests/loggerOutput/CanConvertRpm.txt");
    if (!outputSim) {
        cout<< "Could not open file CanConvertRpm.txt"<<endl;
    }
    double rpmMin = 3600;
    double rpmMax = 30000;
    double rpmAvg = 7400;
    ASSERT_NEAR(getRpm(0.1, rpmMin, rpmMax), 3600, 1);
    ASSERT_NEAR(getRpm(1, rpmMin, rpmMax), 30000, 1);
    ASSERT_NEAR(getRpm(0.229545, rpmMin, rpmMax), 7400, 1);

    outputSim << "TEST: CanConvertRpm" << endl;
    outputSim << "ASSERT_NEAR TEST (tolerance 1):" << endl;
    outputSim << "  Reference throttle: " << 0.1 << endl;
    outputSim << "  Expected RPM:  " << rpmMin << endl;
    outputSim << "  Result: TRUE " << endl;

    outputSim << "ASSERT_NEAR TEST (tolerance 1):" << endl;
    outputSim << "  Reference throttle: " << 1 << endl;
    outputSim << "  Expected RPM:  " << rpmMax << endl;
    outputSim << "  Result: TRUE " << endl;

    outputSim << "ASSERT_NEAR TEST (tolerance 1):" << endl;
    outputSim << "  Reference throttle: " << 0.229545 << endl;
    outputSim << "  Expected RPM:  " << rpmAvg << endl;
    outputSim << "  Result: TRUE " << endl;

}

TEST(IntegrateMotion, CanMaintainTrim) {
    double V = 15;
    double h = 100;
    double vecRef[12] = {14.9873, 0, 0.617672, 0, 0, 0, 0, 0.0411898, 0, 100, 0, 0};
    double vecTrim[12] = {14.9873, 0, 0.617672, 0, 0, 0, 0, 0.0411898, 0, 100, 0, 0};
    double vecComTrim[4] = {0, -0.0376, 0, 0.260227};
    AeroDB db1, db2, db3, db4, DB1, DB2;
    db1 = readDataTest("dba.ini", true);
    db2 = readDataTest("dba_100.ini", true);
    db3 = readDataTest("dba_1000.ini", true);
    db4 = readDataTest("dba_2000.ini", true);
    BatteryDB bat0;
    EngineDB en0; // create dba objects from struct type
    PropDB prop0; // create dba objects from struct type
    // read not aerodynamic database
    //bat0 = readBat("battery.ini"); // Open database, read it and save data to struct of type Not_Aerodb
    en0 = readEnTest("engine.ini"); // Open database, read it and save data to struct of type Not_Aerodb
    prop0 = readPropTest("propeller.ini"); // Open database, read it and save data to struct of type PropDB
    aero_condition atm;
    double rho = AtmosphereCalc(1, atm, 100);

    std::ofstream outputSim("../../Google_tests/tests/loggerOutput/CanMaintainTrim_Simulation_Values.txt");
    std::ofstream loggerRemainders("../../Google_tests/tests/loggerOutput/CanMaintainTrim_Remainders.txt");
    std::ofstream loggerAcceleration("../../Google_tests/tests/loggerOutput/CanMaintainTrim_Acceleration.txt");
    double stateMinusOne[12] = {0}; // i-1
    double rpm = getRpm(vecComTrim[3], en0.laps_min, en0.laps_max);//d_th * RPM_MAX
    for (int j = 0; j < 12; j++) {
        stateMinusOne[j] = vecTrim[j]; // during trim the (i-1)th state is the same as the trim. We assume the
        // simulation starts at the last trim state (i=0) and before that we assume infinite trim states take place
    }
    outputSim << left << setw(15) << "Time" << left << setw(15) << "alpha [deg]"  << left << setw(15) << "u" << left << setw(15) << "v" << left << setw(15) << "w" << left << setw(15) << "p"
              << left << setw(15) << "q" << left << setw(15) << "r" << left << setw(15) << "phi" << left << setw(15) << "theta"
              << left << setw(15) << "psi" << left << setw(15) << "h" << left << setw(15) << "x" << left << setw(15) << "y" << endl;

    outputSim << left << setw(15) << 0.0;
    outputSim << left << setw(15) << atan2(vecTrim[2], vecTrim[0]) * 180.0 / M_PI;
    for (int i = 0; i < 12; i++) {
        outputSim << left << setw(15) << vecTrim[i]; //print to logger
    }
    outputSim << " " << endl;
    double time = 0;
    double dt = 0.01;

    for (int i = 1; i < 100000; i++) {
        h = vecTrim[9];
        time = dt * i;
        getAerodynamicDbWithAltitude(h, DB1, DB2, db1, db2, db3, db4);
        double* newStatesPointer = integrateEquationsOfMotion(DB1, DB2, en0, prop0, rpm, vecTrim, vecComTrim, stateMinusOne, dt, loggerRemainders, loggerAcceleration, atm);
        double newStates[12] = {0};
        for (int j = 0; j < 12; j++) {
            newStates[j] = newStatesPointer[j]; // save the recently calculated state vector in newStates
            stateMinusOne[j] = vecTrim[j]; // save the OLD initial condition as the (i-1)th step
            vecTrim[j] = newStates[j]; // update the initial condition vector with the new states for the next loop iteration
        } // assign values to variable
        delete[] newStatesPointer; // delete pointer to avoid memory leak
        EXPECT_NEAR(vecTrim[0], vecRef[0], 1); // u
        EXPECT_NEAR(atan2(vecTrim[2], vecTrim[0]), -2.16 * M_PI / 180, 0.1); // alpha //TODO: check if tolleranc can be narrower

        outputSim << left << setw(15) << time << left << setw(15) << atan2(newStates[2], newStates[0]) * 180.0 / M_PI; // print to logger
        for (int k = 0; k < 12; k++) {
            outputSim << left << setw(15) << newStates[k]; //print to logger
        }
        outputSim << " " << endl;
    }

}