#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

// Definizione delle variabili da leggere
struct  {
    float C, mu; // Capacità nominale e Rendimento di scarica
}Battery;

 void save_battery(string filePath) {
    string text;
    ifstream myfile; // ifstream is a file input stream that allows us to read any information contained in the file
    myfile.open(filePath); // to open the file
    // if the file is correctly open...
    if(myfile.is_open()) {
        cout << "Reading database " << filePath << endl; // print what the program is reading
        // read line by line and based on the keyword save variables accordingly
        while (getline(myfile, line)) {

            if (i == 5) {
                istringstream A(line);
                A >> Battery.C;
                cout << "" << line << "\n";
            }
            if (i == 6) {
                istringstream A(line);
                A >> Battery.mu;
            }
        }
    }
    else {
    string errorMessage = "Could not open " + filePath;
    cerr << "Could not open " << filePath << endl;
    ::perror("");
    }
    cout<< k<<'\n';
    myfile.close();

return 0;
}

