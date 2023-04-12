#include <iostream>
#include <string>
#include <iomanip>
#include "../declaredFun.h"

using namespace std;

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

