#include <iostream>
#include "declaredFun.h"
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
    dba0 = readData("dba.ini"); // Open database
    dba100 = readData("dba_100.ini");
    dba1000 = readData("dba_1000.ini");
    dba2000 = readData("dba_2000.ini");


    return 0;
}
