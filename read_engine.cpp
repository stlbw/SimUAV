#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
using namespace std;

struct {
    struct {
        double laps_min = -100, laps_max = -100, reduction_rate = -100, nominal_voltage = -100,
                no_load_current = -100, stall_current = -100, stall_torque = -100;
    }engine_description;
}En;

inline void save_engine(string filePath) {
    string text;
    ifstream myfile; // ifstream is a file input stream that allows us to read any information contained in the file
    myfile.open(filePath); // to open the file
    // if the file is correctly open...
    if(myfile.is_open()) {
        cout << "Reading database " << filePath << endl; // print what the program is reading
        // read line by line and based on the keyword save variables accordingly
        while (!myfile.eof()) {
            getline(myfile, text);
            // read and save all aircraft description data
            if (text.find("NUMERO DI GIRI MINIMO DEL MOTORE  [rpm]") != string::npos & En.engine_description.laps_min == -100) {
                istringstream A(text);
                A >> En.engine_description.laps_min;
            } else if (text.find("NUMERO DI GIRI MASSIMO DEL MOTORE [rpm] SENZA CARICO") != string::npos & En.engine_description.laps_max == -100) {
                istringstream A(text);
                A >> En.engine_description.laps_max;
            } else if (text.find("RAPPORTO DI RIDUZIONE (<1)/MOLTIPLICAZIONE (>1)") != string::npos & En.engine_description.reduction_rate == -100) {
                istringstream A(text);
                A >> En.engine_description.reduction_rate;
            } else if (text.find("NOMINAL VOLTAGE") != string::npos & En.engine_description.nominal_voltage == -100) {
                istringstream A(text);
                A >> En.engine_description.nominal_voltage;
            } else if (text.find("NO LOAD CURRENT") != string::npos & En.engine_description.no_load_current == -100) {
                istringstream A(text);
                A >> En.engine_description.no_load_current;
            } else if (text.find("STALL CURRENT") != string::npos & En.engine_description.stall_current == -100) {
                istringstream A(text);
                A >> En.engine_description.stall_current;
            } else if (text.find("STALL TORQUE") != string::npos & En.engine_description.stall_torque == -100) {
                istringstream A(text);
                A >> En.engine_description.stall_torque;
            }
        }
    }else {
        string errorMessage = "Could not open " + filePath;
        cerr << "Could not open " << filePath << endl;
        ::perror("");
    }
    myfile.close();
}
