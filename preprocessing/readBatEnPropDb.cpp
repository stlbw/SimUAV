#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "../declaredFun.h"
using namespace std;

// Definition of the variables for reading the file: battery
struct BatteryDB{
    double C, mu;
};

// Definition of the variables for reading the file: engine
struct EngineDB{
    double laps_min, laps_max, nominal_voltage,
            no_load_current, stall_current, stall_torque;
    // reduction_rate is not considered
};

// Definition of the variables for reading the file: propeller
struct Propeller_Geometry{
    double diameter, diameter_ogive, np, inertia, nstation;
};

struct Profile_Clarky{
    double Clalpha, Cl0, a0, Cdalpha2, Cdalpha, Cd0;
};

struct Propeller_Specifics{
    vector <double> CSI;
    vector <double> RD;
    vector <double> CH_AD;
    vector <double> BA;
};

struct PropDB {
    int length;
    Propeller_Geometry Pg;
    Profile_Clarky Pc;
    Propeller_Specifics Ps;
};

// Function for saving variable from the file: battery
void save_battery(BatteryDB *db, string filePath) {
    string text;
    ifstream myfile; // ifstream is a file input stream that allows us to read any information contained in the file
    myfile.open(filePath); // to open the file
    if(myfile.is_open()) {
        cout << "Reading database " << filePath << endl;
        while (!myfile.eof()) {
            getline(myfile, text);
            if (text.find("CAPACITA' NOMINALE \t[mAh]") != string::npos) {
                istringstream A(text);
                A >> db->C;
            } else if (text.find("RENDIMENTO DI SCARICA") != string::npos) {
                istringstream A(text);
                A >> db->mu;
            }
        }
    }
    else {
        string errorMessage = "Could not open " + filePath;
        throw runtime_error(errorMessage);
    }
    myfile.close();
}

// Function for saving variable from the file: engine
void save_engine(EngineDB *db,string filePath) {
    string text;
    ifstream myfile;
    myfile.open(filePath); // to open the file
    if(myfile.is_open()) {
        cout << "Reading database " << filePath << endl; // print what the program is reading
        while (!myfile.eof()) {
            getline(myfile, text);
            if (text.find("NUMERO DI GIRI MINIMO DEL MOTORE  [rpm]") != string::npos) {
                istringstream A(text);
                A >> db->laps_min;
            } else if (text.find("NUMERO DI GIRI MASSIMO DEL MOTORE [rpm] SENZA CARICO") != string::npos) {
                istringstream A(text);
                A >> db->laps_max;
            } else if (text.find("NOMINAL VOLTAGE") != string::npos) {
                istringstream A(text);
                A >> db->nominal_voltage;
            } else if (text.find("NO LOAD CURRENT") != string::npos) {
                istringstream A(text);
                A >> db->no_load_current;
            } else if (text.find("STALL CURRENT") != string::npos) {
                istringstream A(text);
                A >> db->stall_current;
            } else if (text.find("STALL TORQUE") != string::npos) {
                istringstream A(text);
                A >> db->stall_torque;
            }
        }
    }else {
        string errorMessage = "Could not open " + filePath;
        throw runtime_error(errorMessage);
    }
    myfile.close();
}

// Function for saving variable from the file: propeller
void save_propeller(PropDB *db, string filePath) {
    int flag_ps = 0;
    string text;
    ifstream myfile;
    myfile.open(filePath); // to open the file
    if(myfile.is_open()) {
        std::vector<std::string> nameVar {""}; // declare vector to store variables
        cout << "Reading database " << filePath << endl;
        while (!myfile.eof()) {
            getline(myfile, text);
            if (text.find("DIAMETRO [m]") != string::npos) {
                istringstream A(text);
                A >> db->Pg.diameter;
            } else if (text.find("DIAMETRO OGIVA [m]") != string::npos) {
                istringstream A(text);
                A >> db->Pg.diameter_ogive;
            } else if (text.find("NUMERO DI PALE") != string::npos) {
                istringstream A(text);
                A >> db->Pg.np;
            } else if (text.find("INERZIA [kgm^2]") != string::npos) {
                istringstream A(text);
                A >> db->Pg.inertia;
            } else if (text.find("NUMERO DI STAZIONI") != string::npos) {
                istringstream A(text);
                A >> db->Pg.nstation;
            } else if (text.find("Clalfa [rad^-1]") != string::npos) {
                istringstream A(text);
                A >> db->Pc.Clalpha;
            } else if (text.find("Cl0") != string::npos) {
                istringstream A(text);
                A >> db->Pc.Cl0;
            } else if (text.find("a0      [rad]") != string::npos) {
                istringstream A(text);
                A >> db->Pc.a0;
            } else if (text.find("Cdalfa2 [rad^-2]") != string::npos) {
                istringstream A(text);
                A >> db->Pc.Cdalpha2;
            } else if (text.find("Cdalfa  [rad^-1]") != string::npos) {
                istringstream A(text);
                A >> db->Pc.Cdalpha;
            } else if (text.find("Cd0") != string::npos) {
                istringstream A(text);
                A >> db->Pc.Cd0;
            }
            if (text.find("CARATTERISTICHE DELL' ELICA") != string::npos) {
                flag_ps = 1;
            }
            if (flag_ps == 1) {
                if (text.find("********************************") != string::npos) { flag_ps = 0;}
                // we change the state of the flag, in this way we save only useful data
                if (text.find("CSI") != string::npos) {
                    nameVar = {"CSI", "RD", "CH AD", "BA"}; // clean variable to receive new values
                } else if (text.find("********************************") == string::npos && text.find("CSI") == string::npos) {
                    string token;
                    stringstream s(text);
                    int i = 0;
                    while (getline(s, token, '\t')) {
                        if (!token.empty()) {
                            if (nameVar[i].find("CSI") != string::npos) {
                                db->Ps.CSI.push_back(stof(token));
                                i++;
                            } else if (nameVar[i].find("RD") != string::npos) {
                                db->Ps.RD.push_back(stof(token));
                                i++;
                            } else if (nameVar[i].find("CH AD") != string::npos) {
                                db->Ps.CH_AD.push_back(stof(token));
                                i++;
                            } else if (nameVar[i].find("BA") != string::npos) {
                                db->Ps.BA.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }

        }
        cout << "\tFinished reading database." << endl;

    }else {
        string errorMessage = "Could not open " + filePath;
        throw runtime_error(errorMessage);
    }
    myfile.close();
}

// Function for reading the file: battery
BatteryDB readBat(string fileName) {
    string filePath = string(DATABASE_DIR) + "/" + fileName;
    BatteryDB db; // create db object
    save_battery(&db, filePath);
    return db;
}


// Function for reading the file: engine
EngineDB readEn(string fileName) {
    string filePath = string(DATABASE_DIR) + "/" + fileName;
    EngineDB db; // create db object
    save_engine(&db, filePath);
    return db;
}

// Function for reading the file: propeller
PropDB readProp(string fileName) {
    string filePath = string(DATABASE_DIR) + "/" + fileName;
    PropDB db; // create db object
    save_propeller(&db, filePath);
    return db;
}
