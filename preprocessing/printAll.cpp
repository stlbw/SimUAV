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



