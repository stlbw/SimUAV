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
            // trim angles
            Trim_Angles a = trimAngles(dba100, V, h);
            cout << "Alpha trim [deg]: " << a.alpha_trim << endl;
            cout << "Elevator delta trim [deg]: " << a.deltae_trim << endl;
            cout << "Velocity component u [m/s]: " << a.u << endl;
            cout << "Velocity component w [m/s]: " << a.w << endl;

            //trim rpm, T and Throttle
            cout << "" << endl;
            Trim_Engine_Propeller y = trimEnginePropeller(dba100, en0, prop0, a, V, h);
            cout << "RPM trim: " << y.rpm  << endl;
            cout << "Thrust trim: " << y.T  << endl;
            cout << "Throttle: " << y.Throttle  << endl;

            //initialize the initial conditions vector used for the integration of the aircraft's equations of motion
            // IMPORTANT: integrateEquationsOfMotion receives all values in SI units -> make sure angles are in RAD
            double vecCI[10] = {a.u, 0, a.w, 0, 0, 0, 0, (a.theta_trim * M_PI / 180.0), 0, h}; // [u, v, w, p, q, r, phi, theta, psi, h]
            // initialize the command vector
            double vecComm[4] = {0, (a.deltae_trim * M_PI / 180.0), 0, y.Throttle};

            //todo: decide whether the integration loop should be done in main or under integrateEquationsOfMotion
            integrateEquationsOfMotion(dba100, en0, prop0, y.rpm, vecCI, vecComm);


            cout << "" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "" << endl;

            Modes md = longitudinalStability(dba100,prop0,a,V,h);
            cout << "PHUGOID: " << endl;
            cout << "Frequency [rad/s]: " << md.omega_ph <<endl;
            cout << "Damping ratio: " << md.zeta_ph <<endl;
            cout << "Time to half the amplitude [s]: " << md.t_dim_ph <<endl;
            cout << "Period [s]: " << md.T_ph <<endl;
            cout << "" << endl;
            cout << "SHORT PERIOD:" << endl;
            cout << "Frequency [rad/s]: " << md.omega_sp <<endl;
            cout << "Damping ratio: " << md.zeta_sp <<endl;
            cout << "Time to half the amplitude [s]: " << md.t_dim_sp <<endl;
            cout << "Period [s]: " << md.T_sp <<endl;

            cout << "" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "" << endl;



        }
        catch (const range_error& e){
            cerr<<"Error: "<<e.what()<<endl; //print error
            return 1; // return from main
        }

    }
    catch (const runtime_error& e){
        cerr<<"Error: "<<e.what()<<endl; //print error
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
