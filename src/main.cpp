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

    // read aerodynamic database for different altitudes
    dba0 = readData("dba.ini"); // Open database, read it and save data to struct of type AeroDB
    dba100 = readData("dba_100.ini");
    dba1000 = readData("dba_1000.ini");
    dba2000 = readData("dba_2000.ini");

    GeneralDB bat0, en0; // create dba objects from struct type
    PropDB prop0; // create dba objects from struct type

    // read not aerodynamic database
    bat0 = readBat("battery.ini"); // Open database, read it and save data to struct of type Not_Aerodb
    en0 = readEn("engine.ini"); // Open database, read it and save data to struct of type Not_Aerodb
    prop0 = readProp("propeller.ini"); // Open database, read it and save data to struct of type PropDB


    cout << ">> Finished reading all databases. Please select how should the information be displayed and press ENTER:" << endl;
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

    linearInterpolation(dba0.alpha, dba0.fx.cx_a, 3.58);
    return 0;
}


