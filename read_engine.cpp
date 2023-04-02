#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
using namespace std;

struct {
    struct {
        double laps_min, laps_max, reduction_rate, nominal_voltage,
                no_load_current, stall_current, stall_torque;
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
            if (text.find("NUMERO DI GIRI MINIMO DEL MOTORE  [rpm]") != string::npos) {
                istringstream A(text);
                A >> En.engine_description.laps_min;
            } else if (text.find("NUMERO DI GIRI MASSIMO DEL MOTORE [rpm] SENZA CARICO") != string::npos) {
                istringstream A(text);
                A >> En.engine_description.laps_max;
            } else if (text.find("RAPPORTO DI RIDUZIONE (<1)/MOLTIPLICAZIONE (>1)") != string::npos) {
                istringstream A(text);
                A >> En.engine_description.reduction_rate;
            } else if (text.find("NOMINAL VOLTAGE") != string::npos) {
                istringstream A(text);
                A >> En.engine_description.nominal_voltage;
            } else if (text.find("NO LOAD CURRENT") != string::npos) {
                istringstream A(text);
                A >> En.engine_description.no_load_current;
            } else if (text.find("STALL CURRENT") != string::npos) {
                istringstream A(text);
                A >> En.engine_description.stall_current;
            } else if (text.find("STALL TORQUE") != string::npos) {
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
