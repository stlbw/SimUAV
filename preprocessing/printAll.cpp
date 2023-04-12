#include <iostream>
#include <string>
#include <iomanip>
#include "../declaredFun.h"

using namespace std;

/**
 * Prints the simplified version of an aerodynamic database
 * For variables of type vector prints the first 5 and last 5 elements
 * @param db
 * @param dbName
 * @param switchCase
 * @param printToFile
 */
void printDba(AeroDB db, string dbName, char switchCase) {
    int displayData; //qty of elements to be displayed for first and last elements
    if(switchCase == '1') {displayData = 3;} else {displayData = db.length;} // set size of vectors to be printed


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

    cout << left << setw(50) << "THRUST AXIS OFFSET (REF. TO XB ALONG X)" << left << setw(25) << db.Ad.Thrust_axis_offset_x << endl;
    cout << left << setw(50) << "THRUST AXIS OFFSET (REF. TO XB ALONG Y)" << left << setw(25) << db.Ad.Thrust_axis_offset_y << endl;
    cout << left << setw(50) << "THRUST AXIS OFFSET (REF. TO XB ALONG Z)" << left << setw(25) << db.Ad.Thrust_axis_offset_z << endl;
    cout << left << setw(50) << "THRUST AXIS ANGULAR OFFSET" << left << setw(25) << db.Ad.Thrust_axis_ang_off_xy << endl;
    cout << left << setw(50) << "(REF. TO XB / X-Y PLANE / POSITIVE RIGHT)" << endl;
    cout << left << setw(50) << "THRUST AXIS ANGULAR OFFSET" << left << setw(25) << db.Ad.Thrust_axis_ang_off_xz << endl;
    cout << left << setw(50) << "(REF. TO XB / X-Z PLANE / POSITIVE RIGHT)" << endl;

    cout << left << setw(50) << "NUMBER OF ANGLES OF ATTACK" << left << setw(25) << db.length << endl;

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

    if (db.Ad.option_cog_update == 1) {
        cout << left << setw(50) << "OPTION FOR C.G. UPDATE" << left << setw(25) << "YES" << endl;
    }
    else { cout << left << setw(50) << "OPTION FOR C.G. UPDATE" << left << setw(25) << "NO" << endl; }
    cout << left << setw(50) << "CENTER OF GRAVITY REFERENCE LOCATION (UPDATED)" << left << setw(25) << db.Ad.cog_updated << endl;
    cout << left << setw(50) << "PILOT POSITION (REF. TO CG ALONG XB)" << left << setw(25) << db.Ad.pilot_position_x << endl;
    cout << left << setw(50) << "PILOT POSITION (REF. TO CG ALONG YB)" << left << setw(25) << db.Ad.pilot_position_y << endl;
    cout << left << setw(50) << "PILOT POSITION (REF. TO CG ALONG ZB)" << left << setw(25) << db.Ad.pilot_position_z << endl;

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
    if(switchCase == '1') {
        cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
             << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
             << "..." << endl;

        for (int i = db.alpha.size() - displayData; i < db.alpha.size(); i++) {
            cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.ss.cx[i] << left << setw(20)
                 << db.ss.cy[i] << left
                 << setw(20) << db.ss.cz[i] << left << setw(20) << db.ss.cl[i] << left << setw(20) << db.ss.cm[i]
                 << left << setw(20)
                 << db.ss.cn[i] << endl;
        }
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
    if(switchCase == '1') {
        cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
             << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
             << "..." << endl;

        for (int i = db.alpha.size() - displayData; i < db.alpha.size(); i++) {
            cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.fx.cx_a[i] << left << setw(20)
                 << db.fx.cx_ap[i] << left
                 << setw(20) << db.fx.cx_u[i] << left << setw(20) << db.fx.cx_q[i] << left << setw(20) << db.fx.cx_b[i]
                 << left << setw(20)
                 << db.fx.cx_p[i] << left << setw(20) << db.fx.cx_r[i] << endl;
        }
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
    if(switchCase == '1') {
        cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
             << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
             << "..." << endl;

        for (int i = db.alpha.size() - displayData; i < db.alpha.size(); i++) {
            cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.fy.cy_b[i] << left << setw(20)
                 << db.fy.cy_bp[i] << left
                 << setw(20) << db.fy.cy_p[i] << left << setw(20) << db.fy.cy_r[i] << left << setw(20) << db.fy.cy_a[i]
                 << left << setw(20)
                 << db.fy.cy_q[i] << endl;
        }
    }
    cout << "" << endl;
    cout << "\t\tZ  FORCE DERIVATIVES" << endl;
    cout << left << setw(20) << "ALPHA (deg)" << left << setw(20) << "CZA" << left << setw(20) << "CZAP" << left
         << setw(20) << "CZU" << left << setw(20) << "CZQ" << left << setw(20) << "CZB" << left << setw(20)
         << "CZP" << left << setw(20) << "CZR" << endl;
    for (int i = 0; i < displayData; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.fz.cz_a[i] << left << setw(20) << db.fz.cz_ap[i] << left
             << setw(20) << db.fz.cz_u[i] << left << setw(20) << db.fz.cz_q[i] << left << setw(20) << db.fz.cz_b[i] << left << setw(20)
             << db.fz.cz_p[i] << left << setw(20) << db.fz.cz_r[i] << endl;
    }
    if(switchCase == '1') {
        cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
             << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
             << "..." << endl;

        for (int i = db.alpha.size() - displayData; i < db.alpha.size(); i++) {
            cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.fz.cz_a[i] << left << setw(20)
                 << db.fz.cz_ap[i] << left
                 << setw(20) << db.fz.cz_u[i] << left << setw(20) << db.fz.cz_q[i] << left << setw(20) << db.fz.cz_b[i]
                 << left << setw(20)
                 << db.fz.cz_p[i] << left << setw(20) << db.fz.cz_r[i] << endl;
        }
    }
    cout << "" << endl;
    cout << "\t\tROLLING MOMENT DERIVATIVES" << endl;
    cout << left << setw(20) << "ALPHA (deg)" << left << setw(20) << "CLB" << left << setw(20) << "CLBP" << left
         << setw(20) << "CLP" << left << setw(20) << "CLR" << left << setw(20) << "CLA" << left << setw(20)
         << "CLQ" <<endl;
    for (int i = 0; i < displayData; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.rm.cl_b[i] << left << setw(20) << db.rm.cl_bp[i] << left
             << setw(20) << db.rm.cl_p[i] << left << setw(20) << db.rm.cl_r[i] << left << setw(20) << db.rm.cl_a[i] << left << setw(20)
             << db.rm.cl_q[i] << endl;
    }
    if(switchCase == '1') {
        cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
             << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
             << "..." << endl;

        for (int i = db.alpha.size() - displayData; i < db.alpha.size(); i++) {
            cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.rm.cl_b[i] << left << setw(20)
                 << db.rm.cl_bp[i] << left
                 << setw(20) << db.rm.cl_p[i] << left << setw(20) << db.rm.cl_r[i] << left << setw(20) << db.rm.cl_a[i]
                 << left << setw(20)
                 << db.rm.cl_q[i] << endl;
        }
    }

    cout << "" << endl;
    cout << "\t\tPITCHING MOMENT DERIVATIVES" << endl;
    cout << left << setw(20) << "ALPHA (deg)" << left << setw(20) << "CMA" << left << setw(20) << "CMAP" << left
         << setw(20) << "CMU" << left << setw(20) << "CMQ" << left << setw(20) << "CMB" << left << setw(20)
         << "CMP" << left << setw(20) << "CMR" << endl;
    for (int i = 0; i < displayData; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.pm.cm_a[i] << left << setw(20) << db.pm.cm_ap[i] << left
             << setw(20) << db.pm.cm_u[i] << left << setw(20) << db.pm.cm_q[i] << left << setw(20) << db.pm.cm_b[i] << left << setw(20)
             << db.pm.cm_p[i] << left << setw(20) << db.pm.cm_r[i] << endl;
    }
    if(switchCase == '1') {
        cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
             << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
             << "..." << endl;

        for (int i = db.alpha.size() - displayData; i < db.alpha.size(); i++) {
            cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.pm.cm_a[i] << left << setw(20)
                 << db.pm.cm_ap[i] << left
                 << setw(20) << db.pm.cm_u[i] << left << setw(20) << db.pm.cm_q[i] << left << setw(20) << db.pm.cm_b[i]
                 << left << setw(20)
                 << db.pm.cm_p[i] << left << setw(20) << db.pm.cm_r[i] << endl;
        }
    }
    cout << "" << endl;
    cout << "\t\tYAWING MOMENT DERIVATIVES" << endl;
    cout << left << setw(20) << "ALPHA (deg)" << left << setw(20) << "CNB" << left << setw(20) << "CNBP" << left
         << setw(20) << "CNP" << left << setw(20) << "CNR" << left << setw(20) << "CNA" << left << setw(20)
         << "CNQ" <<endl;
    for (int i = 0; i < displayData; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.ym.cn_b[i] << left << setw(20) << db.ym.cn_bp[i] << left
             << setw(20) << db.ym.cn_p[i] << left << setw(20) << db.ym.cn_r[i] << left << setw(20) << db.ym.cn_a[i] << left << setw(20)
             << db.ym.cn_q[i] << endl;
    }
    if(switchCase == '1') {
        cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
             << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
             << "..." << endl;

        for (int i = db.alpha.size() - displayData; i < db.alpha.size(); i++) {
            cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.ym.cn_b[i] << left << setw(20)
                 << db.ym.cn_bp[i] << left
                 << setw(20) << db.ym.cn_p[i] << left << setw(20) << db.ym.cn_r[i] << left << setw(20) << db.ym.cn_a[i]
                 << left << setw(20)
                 << db.ym.cn_q[i] << endl;
        }
    }
    cout << "" << endl;
    cout << "\t\tCONTROL FORCE DERIVATIVES" << endl;
    cout << left << setw(20) << "ALPHA (deg)" << left << setw(20) << "CXDE" << left << setw(20) << "CXDLE" << left
         << setw(20) << "CZDE" << left << setw(20) << "CZDLE" << left << setw(20) << "CYDA" << left << setw(20)
         << "CYDR" <<endl;
    for (int i = 0; i < displayData; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.cf.cx_de[i] << left << setw(20) << db.cf.cx_dle[i] << left
             << setw(20) << db.cf.cz_de[i] << left << setw(20) << db.cf.cz_dle[i] << left << setw(20) << db.cf.cy_da[i] << left << setw(20)
             << db.cf.cy_dr[i] << endl;
    }
    if(switchCase == '1') {
        cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
             << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
             << "..." << endl;

        for (int i = db.alpha.size() - displayData; i < db.alpha.size(); i++) {
            cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.cf.cx_de[i] << left << setw(20)
                 << db.cf.cx_dle[i] << left
                 << setw(20) << db.cf.cz_de[i] << left << setw(20) << db.cf.cz_dle[i] << left << setw(20)
                 << db.cf.cy_da[i] << left << setw(20)
                 << db.cf.cy_dr[i] << endl;
        }
    }
    cout << "" << endl;
    cout << "\t\tCONTROL MOMENT DERIVATIVES" << endl;
    cout << left << setw(20) << "ALPHA (deg)" << left << setw(20) << "CLDA" << left << setw(20) << "CLDR" << left
         << setw(20) << "CMDE" << left << setw(20) << "CMDLE" << left << setw(20) << "CNDA" << left << setw(20)
         << "CNDR" <<endl;
    for (int i = 0; i < displayData; i++) {
        cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.cm.cl_da[i] << left << setw(20) << db.cm.cl_dr[i] << left
             << setw(20) << db.cm.cm_de[i] << left << setw(20) << db.cm.cm_dle[i] << left << setw(20) << db.cm.cn_da[i] << left << setw(20)
             << db.cm.cn_dr[i] << endl;
    }
    if(switchCase == '1') {
        cout << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left
             << setw(20) << "..." << left << setw(20) << "..." << left << setw(20) << "..." << left << setw(20)
             << "..." << endl;

        for (int i = db.alpha.size() - displayData; i < db.alpha.size(); i++) {
            cout << left << setw(20) << db.alpha[i] << left << setw(20) << db.cm.cl_da[i] << left << setw(20)
                 << db.cm.cl_dr[i] << left
                 << setw(20) << db.cm.cm_de[i] << left << setw(20) << db.cm.cm_dle[i] << left << setw(20)
                 << db.cm.cn_da[i] << left << setw(20)
                 << db.cm.cn_dr[i] << endl;
        }
    }
    cout << "" << endl;
    cout << "---------------------------------END OF AERODYNAMIC DATABASE @" << dbName << " ---------------------------------" << endl;

}

void printBat(GeneralDB db, string dbName, char switchCase) {
    if (switchCase == '1') {
        cout << "" << endl;
        cout << "--------------------------------- BATTERY FILE @" << dbName << " ---------------------------------"
             << endl;
        cout << "BATTERIE : LI-ION POLYMER Kok1200H" << endl;
        cout << "" << endl;
        cout << left << setw(50) << "CAPACITA' NOMINALE [mAh]" << left << setw(25) << db.Bat.C << endl;
        cout << left << setw(50) << "RENDIMENTO DI SCARICA" << left << setw(25) << db.Bat.mu << endl;
        cout << "" << endl;
        cout << "---------------------------------END OF BATTERY FILE @" << dbName
             << " ---------------------------------" << endl;
    }
}

void printEn(GeneralDB db, string dbName, char switchCase) {
    if (switchCase == '1') {
        cout << "" << endl;
        cout << "--------------------------------- ENGINE FILE @" << dbName << " ---------------------------------"
             << endl;
        cout << "MOTORE : Hacker B20-26S" << endl;
        cout << "DATI IN UNITA' S.I." << endl;
        cout << "" << endl;
        cout << left << setw(50) << "NUMERO DI GIRI MINIMO DEL MOTORE  [rpm]" << left << setw(25) << db.En.laps_min
             << endl;
        cout << left << setw(50) << "NUMERO DI GIRI MASSIMO DEL MOTORE [rpm] SENZA CARICO" << left << setw(25)
             << db.En.laps_max << endl;
        cout << left << setw(50) << "NOMINAL VOLTAGE" << left << setw(25) << db.En.nominal_voltage << endl;
        cout << left << setw(50) << "NO LOAD CURRENT" << left << setw(25) << db.En.no_load_current << endl;
        cout << left << setw(50) << "STALL CURRENT" << left << setw(25) << db.En.stall_current << endl;
        cout << left << setw(50) << "STALL TORQUE" << left << setw(25) << db.En.stall_torque << endl;
        cout << "" << endl;
        cout << "---------------------------------END OF ENGINE FILE @" << dbName
             << " ---------------------------------" << endl;
    }
}
void printProp(PropDB db, string dbName, char switchCase) {
    if (switchCase == '1') {
        cout << "" << endl;
        cout << "--------------------------------- Propeller FILE @" << dbName << " ---------------------------------"
             << endl;
        cout << "ELICA :  MICROHAWK (ELICA 5\"X 4\")" << endl;
        cout << "DATI IN UNITA' S.I." << endl;
        cout << "" << endl;
        cout << left << setw(50) << "DIAMETRO [m]" << left << setw(25) << db.Pg.diameter
             << endl;
        cout << left << setw(50) << "DIAMETRO OGIVA [m]" << left << setw(25)
             << db.Pg.diameter_ogive << endl;
        cout << left << setw(50) << "NUMERO DI PALE" << left << setw(25) << db.Pg.np << endl;
        cout << left << setw(50) << "INERZIA [kgm^2]" << left << setw(25) << db.Pg.inertia << endl;
        cout << left << setw(50) << "NUMERO DI STAZIONI" << left << setw(25) << db.Pg.nstation << endl;
        cout << "" << endl;
        cout << "PROFILO : CLARK Y" << endl;
        cout << left << setw(50) << "Clalfa [rad^-1]" << left << setw(25) << db.Pc.Clalpha << endl;
        cout << left << setw(50) << "Cl0" << left << setw(25) << db.Pc.Cl0 << endl;
        cout << left << setw(50) << "Ca0 [rad]" << left << setw(25) << db.Pc.a0 << endl;
        cout << left << setw(50) << "Cdalfa2 [rad^-2]" << left << setw(25) << db.Pc.Cdalpha2 << endl;
        cout << left << setw(50) << "Cdalfa  [rad^-1]" << left << setw(25) << db.Pc.Cdalpha << endl;
        cout << left << setw(50) << "Cd0" << left << setw(25) << db.Pc.Cd0 << endl;
        cout << "" << endl;
        cout << "\t\tCARATTERISTICHE DELL' ELICA" << endl;
        cout << left << setw(20) << "CSI" << left << setw(20) << "RD [m] " << left << setw(20) << "CH AD" << left
             << setw(20) << "BA [deg] " <<endl;
        for (int i = 0; i <20; i++) {
            cout << left << setw(20) << db.Ps.CSI[i] << left << setw(20) << db.Ps.RD[i] << left << setw(20) << db.Ps.CH_AD[i] << left
                 << setw(20) << db.Ps.BA[i] <<endl;}

        cout << "---------------------------------END OF PROPELLER FILE @" << dbName
             << " ---------------------------------" << endl;
    }

}
;

