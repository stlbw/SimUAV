#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "../declaredFun.h"
using namespace std;

struct Path{
    vector <double> Psi;
};


void save_psiRef(Path *db, string filePath) {
    ifstream myfile; // ifstream is a file input stream that allows us to read any information contained in the file
    myfile.open(filePath); // to open the file
    double length = 0;
    // if the file is correctly open...
    if(myfile.is_open()) {
        string text;
        while (!myfile.eof()) {
            getline(myfile, text);
            string token;
            stringstream s(text);
            while (getline(s, token, ' ')) {
                if (!token.empty()) {
                    db->Psi.push_back(stof(token));
                    length++;
                }
            }
        }
        myfile.close();

        for (int i=0; i <= length; i++) {
            if (db->Psi[i] < - M_PI){
                db->Psi[i] += 2 * M_PI;
            }
            else if (db->Psi[i] > M_PI){
                db->Psi[i] -= 2 * M_PI;
            }
        }

    }
    else {
        string errorMessage = "Could not open " + filePath;
        throw runtime_error(errorMessage);
    }

}

Path read_psiref(string fileName) {
    string filePath = string(DATABASE_DIR) + "/" + fileName;
    Path db; // create db object
    save_psiRef(&db, filePath);
    return db;
}


