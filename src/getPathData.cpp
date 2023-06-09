#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "../declaredFun.h"
using namespace std;

struct Path{
    vector <double> Psi;
    double sizePsi;
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
        db->sizePsi = length;

    }
    else {
        string errorMessage = "Could not open " + filePath;
        throw runtime_error(errorMessage);
    }

}

Path read_psiref(string fileName) {
    string filePath = "../path/" + fileName;
    Path db; // create db object
    save_psiRef(&db, filePath);
    return db;
}


