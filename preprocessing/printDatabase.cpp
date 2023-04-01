//
// Created by Matheus Padilha on 01/04/23.
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../declaredFun.h"

using namespace std;

/**
 * Prints the simplified version of an aerodynamic database
 * For variables of type vector prints the first 5 and last 5 elements
 * @param db
 */
void printSimplifiedDba(AeroDB db, string dbName) {
    int displayData = 3; // qty of elements to be displayed for first and last elements
    cout << "" << endl;
    cout << "--------------------------------- AERODYNAMIC DATABASE @" << dbName << " ---------------------------------"
         << endl;
    cout << "DATA IN S.I. UNITS\t\tBODY AXES REFERENCE" << endl;
    cout << "" << endl;
    cout << left << setw(50) << "MASS" << left << setw(25) << db.Ad.Mass << endl;
    cout << left << setw(50) << "WING SPAN" << left << setw(25) << db.Ad.Wing_spann << endl;
    cout << left << setw(50) << "WING AREA" << left << setw(25) << db.Ad.Wing_area << endl;
    cout << left << setw(50) << "CHORD" << left << setw(25) << db.Ad.Chord << endl;
    cout << left << setw(50) << "MACH DRAG RISE" << left << setw(25) << db.Ad.Mach_drag_rise << endl;
    //skiped some variables
    if (db.Ad.rotary_deriv == 1) {
        cout << left << setw(50) << "ROTARY DERIVATIVES" << left << setw(25) << "PRESENT" << endl;
    }
    else { cout << left << setw(50) << "ROTARY DERIVATIVES" << left << setw(25) << "MISSING" << endl; }
    cout << left << setw(50) << "CENTER OF GRAVITY" << left << setw(25) << db.Ad.COG << endl;
    cout << left << setw(50) << "(REFERENCE LOCATION REF. TO CMAER)" << endl;
    cout << left << setw(50) << "JX" << left << setw(25) << db.Ad.Jx << endl;
    cout << left << setw(50) << "JY" << left << setw(25) << db.Ad.Jy << endl;
    cout << left << setw(50) << "JZ" << left << setw(25) << db.Ad.jz << endl;
    cout << left << setw(50) << "JXZ" << left << setw(25) << db.Ad.jxz << endl;
    // skipped some variables

    cout << "" << endl;
    cout << "\t\tDEFLECTION LIMITS" << endl;
    cout << left << setw(50) << "ELEVATOR (max)" << left << setw(25) << db.Dl.Elevator_max << endl;
    cout << left << setw(50) << "ELEVATOR (min)" << left << setw(25) << db.Dl.Elevator_min << endl;
    cout << left << setw(50) << "AILERONS (symmetrical)" << left << setw(25) << db.Dl.Ailerons << endl;
    cout << left << setw(50) << "RUDDER   (symmetrical)" << left << setw(25) << db.Dl.Rudder << endl;
    cout << left << setw(50) << "FLAP (up)" << left << setw(25) << db.Dl.Flap_up << endl;
    cout << left << setw(50) << "FLAP (down)" << left << setw(25) << db.Dl.Flap_down << endl;

    cout << "" << endl;
    cout << "\t\tFUEL MASS" << endl;
    if (db.Fm.Mass_switch == 1) { cout << left << setw(50) << "MASS SWITCH" << left << setw(25) << "VARIABLE" << endl; }
    else { cout << left << setw(50) << "MASS SWITCH" << left << setw(25) << "CONSTANT" << endl; }
    cout << left << setw(50) << "FUEL WEIGHT FRACTION (REF. TO TOW)" << left << setw(25) << db.Fm.Fuel_weight_fraction
         << endl;

    cout << "" << endl;
    cout << "\t\tSTEADY STATE COEFFICIENTS" << endl;
    cout << left << setw(20) << "ALPHA (deg)" << left << setw(20) << "CX" << left << setw(20) << "CY" << left
         << setw(20) << "CZ" << left << setw(20) << "CL" << left << setw(20) << "CM" << left << setw(20)
         << "CN" << endl;
    for (int i = 0; i < displayData; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.ss.cx[i] << left << setw(20) << db.ss.cy[i] << left
             << setw(20) << db.ss.cz[i] << left << setw(20) << db.ss.cl[i] << left << setw(20) << db.ss.cm[i] << left << setw(20)
             << db.ss.cn[i] << endl;
    }
    cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
         << setw(20) << "..."<< left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
         << "..." << endl;
    for (int i = db.alpha.size() - displayData ; i < db.alpha.size() ; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.ss.cx[i] << left << setw(20) << db.ss.cy[i] << left
             << setw(20) << db.ss.cz[i] << left << setw(20) << db.ss.cl[i] << left << setw(20) << db.ss.cm[i] << left << setw(20)
             << db.ss.cn[i] << endl;
    }

    cout << "" << endl;
    cout << "" << endl;
    cout << "\t\tAERODYNAMIC DERIVATIVES" << endl;
    cout << "" << endl;
    cout << "\t\tX  FORCE DERIVATIVES" << endl;
    cout << left << setw(20) << "ALPHA (deg)" << left << setw(20) << "CXA" << left << setw(20) << "CXAP" << left
         << setw(20) << "CXU" << left << setw(20) << "CXQ" << left << setw(20) << "CXB" << left << setw(20)
         << "CXP" << left << setw(20) << "CXR" << endl;
    for (int i = 0; i < displayData; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.fx.cx_a[i] << left << setw(20) << db.fx.cx_ap[i] << left
             << setw(20) << db.fx.cx_u[i] << left << setw(20) << db.fx.cx_q[i] << left << setw(20) << db.fx.cx_b[i] << left << setw(20)
             << db.fx.cx_p[i] << left << setw(20) << db.fx.cx_r[i] << endl;
    }
    cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
         << setw(20) << "..."<< left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
         << "..." << left << setw(20) << "..." << endl;
    for (int i = db.alpha.size() - displayData ; i < db.alpha.size() ; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.fx.cx_a[i] << left << setw(20) << db.fx.cx_ap[i] << left
             << setw(20) << db.fx.cx_u[i] << left << setw(20) << db.fx.cx_q[i] << left << setw(20) << db.fx.cx_b[i] << left << setw(20)
             << db.fx.cx_p[i] << left << setw(20) << db.fx.cx_r[i] << endl;
    }

    cout << "" << endl;
    cout << "\t\tY  FORCE DERIVATIVES" << endl;
    cout << left << setw(20) << "ALPHA (deg)" << left << setw(20) << "CYB" << left << setw(20) << "CYBP" << left
         << setw(20) << "CYP" << left << setw(20) << "CYR" << left << setw(20) << "CYA" << left << setw(20)
         << "CYQ" <<endl;
    for (int i = 0; i < displayData; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.fy.cy_b[i] << left << setw(20) << db.fy.cy_bp[i] << left
             << setw(20) << db.fy.cy_p[i] << left << setw(20) << db.fy.cy_r[i] << left << setw(20) << db.fy.cy_a[i] << left << setw(20)
             << db.fy.cy_q[i] << endl;
    }
    cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
         << setw(20) << "..."<< left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
         << "..." << endl;
    for (int i = db.alpha.size() - displayData ; i < db.alpha.size() ; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.fy.cy_b[i] << left << setw(20) << db.fy.cy_bp[i] << left
             << setw(20) << db.fy.cy_p[i] << left << setw(20) << db.fy.cy_r[i] << left << setw(20) << db.fy.cy_a[i] << left << setw(20)
             << db.fy.cy_q[i] << endl;
    }
}