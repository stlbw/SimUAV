#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

/*struct Path{
    vector <double> Psi;
};


inline void save_psiRef(Path *db, string filePath) {
    string text;
    ifstream myfile; // ifstream is a file input stream that allows us to read any information contained in the file
    myfile.open(filePath); // to open the file
    // if the file is correctly open...
    if(myfile.is_open()) {
        // read line by line and based on the keyword save variables accordingly
        double i = 0;
        while (!myfile.eof()) {
            getline(myfile, text);
            // read and save all aircraft description data
                istringstream A(text);
                A >> db->Psi[i];
                i++;
        }
    }
    else {
        string errorMessage = "Could not open " + filePath;
        throw runtime_error(errorMessage);
    }
    myfile.close();
}

inline Path read_psiref(string fileName) {
    string filePath = string(DATABASE_DIR) + "/" + fileName;
    Path db; // create db object
    save_psiRef(&db, filePath);
    return db;
}*/
