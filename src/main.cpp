#include <iostream>
#include "../declaredFun.h"
int main() {
    cout << "---------------------------------------------------------------" << endl;
    cout << "\t\t\t Flight Simulator - MH850-3" << endl;
    cout << "\t\t\t M.Sc. in Aerospace Engineering" << endl;
    cout << "\t\t\t Politecnico di Torino - Italy" << endl;
    cout << "\t\t\t Flight Simulation - June 2023" << endl;
    cout << "---------------------------------------------------------------" << endl;
    cout << "Created by: " << endl;
    cout << "---------------------------------------------------------------" << endl;

    cout << ">> Reading databases" << endl;
    // getDba("dba.ini"); // function implemented by Mario
    AeroDB dba0, dba100, dba1000, dba2000; // create dba objects from struct type

    try{
        // read aerodynamic database for different altitudes
        dba0 = readData("dba.ini"); // Open database, read it and save data to struct of type AeroDB
        dba100 = readData("dba_100.ini");
        dba1000 = readData("dba_1000.ini");
        dba2000 = readData("dba_2000.ini");

        BatteryDB bat0;
        EngineDB en0; // create dba objects from struct type
        PropDB prop0; // create dba objects from struct type

        // read not aerodynamic database
        bat0 = readBat("battery.ini"); // Open database, read it and save data to struct of type Not_Aerodb
        en0 = readEn("engine.ini"); // Open database, read it and save data to struct of type Not_Aerodb
        prop0 = readProp("propeller.ini"); // Open database, read it and save data to struct of type PropDB


        cout << ">> Finished reading all databases."<<endl;
        cout <<'\n'<<" Please select how should the information be displayed and press ENTER:" << endl;
        cout << "\t1 - Partial (simplified) version" << endl;
        cout << "\t2 - Print full database to screen" << endl;
        cout << "\t3 - Save database to file" << endl;
        cout << ">> ";
        char dataCheckSwitch = '3'; // used to store the user's choice
        bool flagCaseFound = false;
        bool printToFile;
        while(!flagCaseFound) {
            cin >> dataCheckSwitch; // get user input
            switch (dataCheckSwitch) {
                case '1':
                    cout << "Printing partial version to screen..." << endl;
                    flagCaseFound = true;
                    printToFile = false;
                    break;
                case '2':
                    cout << "Printing full database to screen..." << endl;
                    flagCaseFound = true;
                    printToFile = false;
                    break;
                case '3':
                    cout << "Saving databases to file..." << endl;
                    flagCaseFound = true;
                    printToFile = true;
                    break;
                default:
                    cerr << "Could not read the input. Please make sure to add a number from 1 to 3 and press ENTER." << endl;
                    break;
            }
        }
        // print database
        if(!printToFile) {
            printDba(dba0, "SEA LEVEL", dataCheckSwitch);
            printDba(dba100, "100 m", dataCheckSwitch);
            printDba(dba1000, "1000 m", dataCheckSwitch);
            printDba(dba2000, "2000 m", dataCheckSwitch);
            printBat(bat0, "battery",dataCheckSwitch);
            printEn(en0,"engine",dataCheckSwitch);
            printProp(prop0,"propeller",dataCheckSwitch);
        } else {
            string filePath = "../output/logDatabaseData.txt";
            cout << ">> Saving files to output directory. Path: " << filePath << endl;
            streambuf *coutbuf;
            ofstream out(filePath);
            coutbuf = cout.rdbuf();
            cout.rdbuf(out.rdbuf());
            printDba(dba0, "SEA LEVEL", dataCheckSwitch);
            printDba(dba100, "100 m", dataCheckSwitch);
            printDba(dba1000, "1000 m", dataCheckSwitch);
            printDba(dba2000, "2000 m", dataCheckSwitch);
            printBat(bat0, "battery",dataCheckSwitch);
            printEn(en0,"engine",dataCheckSwitch);
            printProp(prop0,"propeller",dataCheckSwitch);
            cout.rdbuf(coutbuf); //reset to standard output again
        }

        cout << "" << endl;
        cout << "---------------------------------------------------------------" << endl;
        cout << "" << endl;
        // TRIM SECTION:
        double V, h;
        cout << "Insert velocity [m/s]: ";
        cin >> V;
        cout << "" << endl;
        cout << "Insert altitude [m]: ";
        cin >> h;
        cout <<""<<endl;
        cout << "TRIM PARAMETERS: " << endl;
        cout <<""<<endl;

        // try and catch block used to deal with errors -> this way errors are always sent back to the main
        try {
            // declared db for interpolation
            AeroDB DB1;
            AeroDB DB2;

            // get correct dba with altitude
            getAerodynamicDbWithAltitude(h, DB1, DB2, dba0, dba100, dba1000, dba2000);

            // trim angles
            Trim_Angles a = trimAngles(DB1, DB2, V, h);
            cout << "Alpha trim [deg]: " << a.alpha_trim << endl;
            cout << "Elevator delta trim [deg]: " << a.deltae_trim << endl;
            cout << "Velocity component u [m/s]: " << a.u << endl;
            cout << "Velocity component w [m/s]: " << a.w << endl;

            //trim rpm, T and Throttle
            cout << "" << endl;
            Trim_Engine_Propeller y = trimEnginePropeller(DB1, DB2, en0, prop0, a, V, h);
            cout << "RPM trim: " << y.rpm  << endl;
            cout << "Thrust trim: " << y.T  << endl;
            cout << "Torque trim: " << y.Torque << endl;
            cout << "Throttle: " << y.Throttle  << endl;

            cout << "" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "" << endl;

            longitudinalStability(DB1, DB2,prop0,a,V,h); // this function computes the Routh criteria
            // for the aircraft and computes the longitudinal modes in case the criteria is respected

            cout << "" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "" << endl;

            //todo: decide whether the integration loop should be done in main or under integrateEquationsOfMotion
            // todo: get new DB1, DB2 at each step based on h
            double Tsim = 10.0;
            //double dt = 0.01;
            double dt = 0.02;
            int nStep = static_cast<int>(Tsim / dt);

            //initialize the initial conditions vector used for the integration of the aircraft's equations of motion
            // IMPORTANT: integrateEquationsOfMotion receives all values in SI units -> make sure angles are in RAD
            double vecCI[12] = {a.u, 0, a.w, 0, 0, 0, 0, (a.theta_trim * M_PI / 180.0), 0, h, 0, 0}; // [u, v, w, p, q, r, phi, theta, psi, h, x, y]
            // initialize the command vector
            double vecComm[4] = {0, (a.deltae_trim * M_PI / 180.0), 0, y.Throttle};

            double stateMinusOne[12] = {0}; // i-1
            for (int j = 0; j < 12; j++) {
                stateMinusOne[j] = vecCI[j]; // during trim the (i-1)th state is the same as the trim. We assume the
                // simulation starts at the last trim state (i=0) and before that we assume infinite trim states take place
            }

            double** fullStateMatrix = new double*[nStep];
            for (int i = 0; i < 12; i++) {fullStateMatrix[i] = new double[12];}

            // assign the first column as the trim condition
            cout << left << setw(25) << "T"  << left << setw(25) << "u" << left << setw(25) << "v" << left << setw(25) << "w" << left << setw(25) << "p"
                    << left << setw(25) << "q" << left << setw(25) << "r" << left << setw(25) << "phi" << left << setw(25) << "theta"
                    << left << setw(25) << "psi" << left << setw(25) << "h" << left << setw(25) << "x" << left << setw(25) << "y" << endl;
            for (int i = 0; i < 12; i++) {
                fullStateMatrix[0][i] = vecCI[i];
                cout << left << setw(25) << 0.0;
                cout << left << setw(25) << fullStateMatrix[0][i];
            }
            cout << " " << endl;


            for (int i = 1; i <= nStep; i++) {

                double* newStatesPointer = integrateEquationsOfMotion(DB1, DB2, en0, prop0, y.rpm, vecCI, vecComm, stateMinusOne, dt);

                double newStates[12] = {0};

                for (int j = 0; j < 12; j++) {
                    newStates[j] = newStatesPointer[j]; // save the recently calculated state vector in newStates
                    stateMinusOne[j] = vecCI[j]; // save the OLD initial condition as the (i-1)th step
                    vecCI[j] = newStates[j]; // update the initial condition vector with the new states for the next loop iteration
                } // assign values to variable
                delete[] newStatesPointer; // delete pointer to avoid memory leak

                // assign the recently calculated state to the fullStateMatrix at column i
                double time = i * dt;
                cout << left << setw(25) << time;
                for (int k = 0; k < 12; k++) {
                    fullStateMatrix[i][k] = newStates[k];
                    cout << left << setw(25) << fullStateMatrix[i][k];

                }
                cout << " " << endl;


            }

            /*for (int i = 0; i < 12; i++) {
                for (int j = 0; j <= nStep; j++) {
                    cout << fullStateMatrix[i][j] << "   ";
                }
                cout << " " << endl;
            }
             */

            delete[] fullStateMatrix; // delete matrix pointer to avoid memory overflow





        }
        catch (const range_error& e){
            cerr<<"Out of range: "<<e.what()<<endl; //print error
            return 1; // return from main
        }

    }
    catch (const runtime_error& e){
        cerr<<"Runtime error: "<<e.what()<<endl; //print error
        return 1; // return from main
    }




    return 0;
};

// Trim Control
//int main{
        //cout << '\n' << " Please select the input Velocity:" <<;
        //cin >> V;
        //cout << '\n' << " Please select the input Height:" <<;
        //cin>> h;
//};
