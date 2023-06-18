#include <iostream>
#include <fstream>
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
        double V_ref, h_ref;
        double gamma_0 = 0;
        // todo: get gamma
        double flagStartSimulation = 1;
        while(flagStartSimulation != 0) {
            cout << "Insert altitude [m]: ";
            cin >> h_ref;
            cout <<""<<endl;
            cout << "Insert velocity [m/s]: ";
            cin >> V_ref;
            cout <<""<<endl;

            // ask user gamma_0
            int flag_gamma;
            cout <<""<<endl;
            cout << "Simulation will be executed with gamma_0 = 0 [deg]. \nPress 0 to continue or 1 to enter a specific value: ";
            cin >> flag_gamma;
            if (flag_gamma == 1){
                cout << "Please enter the new value for gamma_0 in DEGREES [deg]: ";
                cin >> gamma_0;
            }
            cout << "" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "" << endl;
            cout << "The simulation will run with the following parameters: ";
            cout << "" << endl;
            cout << "V_ref = "<< V_ref << " [m/s]" << endl;
            cout << "h_ref = " << h_ref << " [m]" << endl;
            cout << "gamma_0 = " << gamma_0 << " [deg]" << endl;
            cout << "" << endl;
            cout << "Press 0 to start the simulation, 1 to change the parameters.";
            cout << "" << endl;
            cin >> flagStartSimulation;

        }

        // try and catch block used to deal with errors -> this way errors are always sent back to the main
        try {
            // declared db for interpolation
            AeroDB DB1;
            AeroDB DB2;

            // get correct dba with altitude
            getAerodynamicDbWithAltitude(h_ref, DB1, DB2, dba0, dba100, dba1000, dba2000);

            double Vmax = 25;
            double Vmin = computeVmin(DB1,DB2,h_ref);
            if (V_ref > Vmax){
                cout << "The entered velocity exceeds the maximum value " << endl;
                cout << "velocity is set to the maximum one 25 m/s" << endl;
                V_ref = Vmax;
            } else if (V_ref < Vmin){
                cout << "The entered velocity is lower than the minimum value " << endl;
                cout << "Velocity is set to the minimum one " << Vmin << endl;
                V_ref = Vmin;
            }


            cout << "" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "" << endl;

            cout << "" << endl;
            cout << "TRIM PARAMETERS: " << endl;
            cout <<""<<endl;

            // trim angles
            Trim_Angles a = trimAngles(DB1, DB2, V_ref, h_ref, gamma_0);
            cout << "Alpha trim [deg]: " << a.alpha_trim << endl;
            cout << "Elevator delta trim [deg]: " << a.deltae_trim << endl;
            cout << "Velocity component u [m/s]: " << a.u << endl;
            cout << "Velocity component w [m/s]: " << a.w << endl;

            //trim rpm, T and Throttle
            cout << "" << endl;
            Trim_Engine_Propeller y = trimEnginePropeller(DB1, DB2, en0, prop0, a, V_ref, h_ref, gamma_0);
            cout << "RPM trim: " << y.rpm  << endl;
            cout << "Thrust trim: " << y.T  << endl;
            cout << "Torque trim: " << y.Torque << endl;
            cout << "Throttle: " << y.Throttle  << endl;

            TrimCondition ss;
            ss.alphaDeg = a.alpha_trim;
            ss.V = sqrt(a.u * a.u + a.w * a.w);
            ss.h = h_ref;
            ss.rho = computeDensity(h_ref);

            cout << "" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "" << endl;

            longitudinalStability(DB1, DB2,prop0,a,V_ref,h_ref); // this function computes the Routh criteria
            // for the aircraft and computes the longitudinal modes in case the criteria is respected

            cout << "" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "" << endl;

            // initialize logger
            std::ofstream outputSim("../output/simulationData.txt");
            std::ofstream loggerRemainders("../output/logRemainders.txt");
            std::ofstream loggerAcceleration("../output/logAcceleration.txt");
            if (!outputSim) {
                string error = "Could not open file simulationData.txt";
                throw runtime_error(error);
            }
            if (!loggerRemainders) {
                string error = "Could not open file logRemainders.txt";
                throw runtime_error(error);
            }
            if (!loggerAcceleration) {
                string error = "Could not open file logAcceleration.txt";
                throw runtime_error(error);
            }
            printHeaderLogger(loggerRemainders, loggerAcceleration); //prints header for both loggers

            // SIMULATION FROM HERE
            double Tsim = 150.0;
            double dt = 0.01;
            int nStep = static_cast<int>(Tsim / dt)-2000;

            Path psi0;
            psi0 = read_psiref("SQUARE_psiref.txt");

            //initialize the initial conditions vector used for the integration of the aircraft's equations of motion
            // IMPORTANT: integrateEquationsOfMotion receives all values in SI units -> make sure angles are in RAD
            double vecCI[12] = {a.u, 0, a.w, 0, 0, 0, 0, (a.theta_trim * M_PI / 180.0), 0, h_ref, 0, 0}; // [u, v, w, p, q, r, phi, theta, psi, h, x, y]

            //double vecCI[12] = {a.u, 0, a.w, 0, 0, 0, 0, (a.theta_trim * M_PI / 180.0), 0, h_ref, -300, 0};
            //double vecCI[12] = {a.u, 0, a.w, 0, 0, 0, 0, (a.theta_trim * M_PI / 180.0), 0, h_ref, 0, -300};
            //double vecCI[12] = {a.u, 0, a.w, 0, 0, 0, 0, (a.theta_trim * M_PI / 180.0), psi0.Psi[2000-1], h_ref, 0, 0};
            //double vecCI[12] = {a.u, 0, a.w, 0, 0, 0, 0, (a.theta_trim * M_PI / 180.0), psi0.Psi[2000-1], h_ref, -300, 0};
            //double vecCI[12] = {a.u, 0, a.w, 0, 0, 0, 0, (a.theta_trim * M_PI / 180.0), psi0.Psi[2000-1], h_ref, 0, -300};
            //The choice of 2000-1 is because the simulation starts after 20 sec

            // initialize the command vector
            double vecComm[4] = {0, (a.deltae_trim * M_PI / 180.0), 0, y.Throttle}; //da,de,0,throttle

            double stateMinusOne[12] = {0}; // i-1
            for (int j = 0; j < 12; j++) {
                stateMinusOne[j] = vecCI[j]; // during trim the (i-1)th state is the same as the trim. We assume the
                // simulation starts at the last trim state (i=0) and before that we assume infinite trim states take place
            }

            // create state matrix to allocate the states at each step
            double** fullStateMatrix = new double*[nStep + 1];
            for (int i = 0; i < nStep + 1; i++) {fullStateMatrix[i] = new double[12];}

            // assign the header to simulation
            cout << left << setw(15) << "Time" << left << setw(15) << "alpha [deg]"  << left << setw(15) << "u" << left << setw(15) << "v" << left << setw(15) << "w" << left << setw(15) << "p"
                    << left << setw(15) << "q" << left << setw(15) << "r" << left << setw(15) << "phi" << left << setw(15) << "theta"
                    << left << setw(15) << "psi" << left << setw(15) << "h" << left << setw(15) << "x" << left << setw(15) << "y" << endl;
            //  assign first column as the trim condition
            cout << left << setw(15) << 0.0;
            cout << left << setw(15) << atan2(vecCI[2], vecCI[0]) * 180.0 / M_PI;

            // print to logger the trim step
            outputSim << left << setw(15) << "Time" << left << setw(15) << "alpha [deg]"  << left << setw(15) << "u" << left << setw(15) << "v" << left << setw(15) << "w" << left << setw(15) << "p"
                 << left << setw(15) << "q" << left << setw(15) << "r" << left << setw(15) << "phi" << left << setw(15) << "theta"
                 << left << setw(15) << "psi" << left << setw(15) << "h" << left << setw(15) << "x" << left << setw(15) << "y" << endl;

            outputSim << left << setw(15) << 0.0;
            outputSim << left << setw(15) << atan2(vecCI[2], vecCI[0]) * 180.0 / M_PI;

            // assign fullStateMatrix first column to trim and print it
            for (int i = 0; i < 12; i++) {
                fullStateMatrix[0][i] = vecCI[i];
                cout << left << setw(15) << fullStateMatrix[0][i]; // print to screen
                outputSim << left << setw(15) << fullStateMatrix[0][i]; //print to logger
            }
            cout << " " << endl;
            outputSim << " " << endl;


           // double flagPID = 0;
            double wantPrint=0;
            double runningflag=0;
            int wantPID = 0;
            cout <<'\n'<<" Do you want to implement the PID controller ? (Y/1 N/0):" << endl;
            cin >> wantPID;
            cout <<'\n'<<"The results are gonna be printed on the 'simulationData.txt' file"<< endl;
            cout << '\n'<< "Do you want to print the results also on this screen? (Y/1 N/0):"<<endl;
            cin >>wantPrint;
            // initialize error variables and integrative error to be used when PID is active - store the error at the previous step

            pidController PID_v(-0.0021, -0.00087, -0.0015);
            pidController PID_theta(-0.3, -3.25,-0.01);
            pidController PID_h(0.019, 0.0002,0.01);
            pidController PID_psi(1.5, 0.005, 0.01);
            pidController PID_phi(0.12, 0.0005, 0.001);

            PID_v.setDerivativeFilter(true, 0.0159);
            PID_theta.setDerivativeFilter(true, 0.0159);
            PID_h.setDerivativeFilter(true, 0.1592);

            PID_psi.setErrorCheck(true);
            //PID_phi.setErrorCheck(true);


            // LOOP INTEGRATE EQUATIONS OF MOTION
            int k = 0;
            for (int i = 1; i <= nStep; i++) {

                double time = i * dt;
                k ++;
                if (k > 15000) {
                    k = 2050;
                }

                if (wantPID == 1) {
                    double V_current = sqrt(vecCI[0]*vecCI[0]+vecCI[1]*vecCI[1]+vecCI[2]*vecCI[2]);

                    /*double theta_ref_test = PID_v.computePIDSimple(V_ref, vecCI[0], dt);
                    //double theta_ref_test = PID_v.computePIDSimple(V_ref, V_current , dt);
                    double delta_e_test = PID_theta.computePIDSimple(theta_ref_test, vecCI[7], dt);
                    double delta_th_test = PID_h.computePIDSimple(h_ref, vecCI[9], dt);
                    //double delta_a_test = 0;
                    //vecCI[8] = psi0.Psi[k - 1];
                    //double phi_test = PID_psi.computePIDSimple(psi0.Psi[k],vecCI[8],dt);
                    double phi_test = PID_psi.computePIDSimple(psi0.Psi[k+2000-1],vecCI[8],dt);
                    double delta_a_test = PID_phi.computePIDSimple(phi_test,vecCI[6],dt);*/

                    /*double theta_ref_test = PID_v.computeNewTest(V_ref, vecCI[0], dt);
                    //double theta_ref_test = PID_v.computeNewTest(V_ref, V_current, dt);
                    double delta_e_test = PID_theta.computeNewTest(theta_ref_test, vecCI[7], dt);
                    double delta_th_test = PID_h.computeNewTest(h_ref, vecCI[9], dt);
                    //double phi_test = PID_psi.computeNewTest(psi0.Psi[k],vecCI[8],dt);
                    double phi_test = PID_psi.computeNewTest(psi0.Psi[k+2000-1],vecCI[8],dt);
                    double delta_a_test = PID_phi.computeNewTest(phi_test,vecCI[6],dt);*/

                    /*double theta_ref_test = PID_v.compute(V_ref, vecCI[0], dt);
                    //double theta_ref_test = PID_v.computeNewTest(V_ref, V_current, dt);
                    double delta_e_test = PID_theta.compute(theta_ref_test, vecCI[7], dt);
                    double delta_th_test = PID_h.compute(h_ref, vecCI[9], dt);
                    //double phi_test = PID_psi.computeNewTest(psi0.Psi[k],vecCI[8],dt);
                    double phi_test = PID_psi.compute(psi0.Psi[k+2000-1],vecCI[8],dt);
                    double delta_a_test = PID_phi.compute(phi_test,vecCI[6],dt);*/


                    vecComm[0] = delta_a_test;//delta_a_test; //delta_aileron
                    vecComm[1] = delta_e_test; //delta_elevator
                    vecComm[3] = delta_th_test; //delta_throttle

                    // SATURAZIONI
                    if(vecComm[3] <= 0.1){
                        vecComm[3] = 0.1;
                        cout << "delta_throttle SATURATA -> delta_th = 0.1" << endl;
                        PID_h.resetIntegrativeError();
                    }
                    else if (vecComm[3] >= 1.0) {
                        vecComm[3] = 1.0;
                        cout << "delta_throttle SATURATA -> delta_th = 1.0" << endl;
                        PID_h.resetIntegrativeError();
                    }

                    if(vecComm[1] <= -10 * M_PI / 180){
                        vecComm[1] = -10 * M_PI / 180;
                        cout << "delta_e SATURATA -> delta_e min" << endl;
                        PID_theta.resetIntegrativeError();
                    }
                    else if (vecComm[1] >= 10 * M_PI / 180) {
                        vecComm[1] = 10 * M_PI / 180;
                        cout << "delta_e SATURATA -> delta_e max" << endl;
                        PID_theta.resetIntegrativeError();

                    }
                    if(vecComm[0] <= -10 * M_PI / 180){
                        vecComm[0] = -10 * M_PI / 180;
                        cout << "delta_a SATURATA -> delta_a min" << endl;
                        PID_phi.resetIntegrativeError();
                    }
                    else if (vecComm[0] >= 10 * M_PI / 180) {
                        vecComm[0] = 10 * M_PI / 180;
                        cout << "delta_a SATURATA -> delta_a max" << endl;
                        PID_phi.resetIntegrativeError();
                    }

                }

                double rpm = getRpm(vecComm[3], en0.laps_min, en0.laps_max);
                //double rpm =  3600+(30000-3600)*(vecComm[3]-0.1)/(1-0.1);

                // get correct dba with altitude
                double h = vecCI[9]; // update altitude
                getAerodynamicDbWithAltitude(h, DB1, DB2, dba0, dba100, dba1000, dba2000); // returns the correct DB1 and DB2 to use for the interpolation
                double* newStatesPointer = integrateEquationsOfMotion(DB1, DB2, en0, prop0, rpm, vecCI, vecComm, stateMinusOne, dt, loggerRemainders, loggerAcceleration, ss);

                double newStates[12] = {0};
                for (int j = 0; j < 12; j++) {
                    newStates[j] = newStatesPointer[j]; // save the recently calculated state vector in newStates
                    stateMinusOne[j] = vecCI[j]; // save the OLD initial condition as the (i-1)th step
                    vecCI[j] = newStates[j]; // update the initial condition vector with the new states for the next loop iteration
                } // assign values to variable
                delete[] newStatesPointer; // delete pointer to avoid memory leak

                double current_V = sqrt(vecCI[0] * vecCI[0] + vecCI[1] * vecCI[1] + vecCI[2] * vecCI[2]);


                // assign the recently calculated state to the fullStateMatrix at column i
                if (wantPrint==1) {
                    cout << left << setw(15) << time << left << setw(15) << atan2(newStates[2], newStates[0]) * 180.0 / M_PI;
                }
                else if (wantPrint==0 & runningflag==0){
                    cout << '\n'<< " Running..."<< endl;
                runningflag=1;
                }//print to screen
                outputSim << left << setw(15) << time << left << setw(15) << atan2(newStates[2], newStates[0]) * 180.0 / M_PI; // print to logger

               for (int k = 0; k < 12; k++) {
                   fullStateMatrix[i][k] = newStates[k];
                   if (wantPrint==1) {
                       cout << left << setw(15) << fullStateMatrix[i][k];
                   } //print to screen
                   outputSim << left << setw(15) << fullStateMatrix[i][k]; //print to logger
               }
               if (wantPrint == 1) {cout << " " << endl;}
               outputSim << " " << endl;

                if (current_V > Vmax || current_V < Vmin) { // da cambiare in Vmin?
                    string error = "Out of range: velocity is out of bounds. V = " + to_string(current_V) + " m/s";
                    throw range_error(error);
                }
            }
            //close loggers
            outputSim.close();
            loggerRemainders.close();
            loggerAcceleration.close();
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
