/#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include "../declaredFun.h"
using namespace std;

struct engine_description{
    double laps_min = -100, laps_max = -100, reduction_rate = -100, nominal_voltage = -100,
            no_load_current = -100, stall_current = -100, stall_torque = -100;
};

struct Flags {
    int flag_lapsmin, flag_lapsmax, flag_reduction, flag_volt, flag_noload, flag_scurrent, flag_storque;
    // declared flag
};

struct enginedb{
    engine_description En;
};

int get_engine(enginedb *db, string filePath){
    string text;
    ifstream myfile; // ifstream is a file input stream that allows us to read any information contained in the file
    myfile.open(filePath); // to open the file
    // if the file is correctly open...
    if(myfile.is_open()) {
        // initiate all flags to 0
        // f->flag_lapsmin = f->flag_lapsmax = 0;
        // f->flag_reduction = 0;
        // f->flag_volt = 0;
        // f->flag_noload = 0;
        // f->flag_scurrent = f->flag_storque = 0;
        // std::vector<std::string> nameVar{""}; // declare vector to store variables
        cout << "Reading database " << filePath << endl; // print what the program is reading
        // read line by line and based on the keyword save variables accordingly
        while (!myfile.eof()) {
            getline(myfile, text);
            // read and save all aircraft description data
            if (text.find("NUMERO DI GIRI MINIMO DEL MOTORE  [rpm]") != string::npos & db->En.laps_min == -100) {
                istringstream A(text);
                A >> db->En.laps_min;
            } else if (text.find("NUMERO DI GIRI MASSIMO DEL MOTORE [rpm] SENZA CARICO") != string::npos &
                       db->En.laps_max == -100) {
                istringstream A(text);
                A >> db->En.laps_max;
            } else if (text.find("RAPPORTO DI RIDUZIONE (<1)/MOLTIPLICAZIONE (>1)") != string::npos &
                       db->En.reduction_rate == -100) {
                istringstream A(text);
                A >> db->En.reduction_rate;
            } else if (text.find("NOMINAL VOLTAGE") != string::npos & db->En.nominal_voltage == -100) {
                istringstream A(text);
                A >> db->En.nominal_voltage;
            } else if (text.find("NO LOAD CURRENT") != string::npos & db->En.no_load_current == -100) {
                istringstream A(text);
                A >> db->En.no_load_current;
            } else if (text.find("STALL CURRENT") != string::npos & db->En.stall_current == -100) {
                istringstream A(text);
                A >> db->En.stall_current;
            } else if (text.find("STALL TORQUE") != string::npos & db->En.stall_torque == -100) {
                istringstream A(text);
                A >> db->En.stall_torque;
            }
        }
    }else {
        string errorMessage = "Could not open " + filePath;
        cerr << "Could not open " << filePath << endl;
        ::perror("");
    }
    myfile.close();

    return 0;
}
