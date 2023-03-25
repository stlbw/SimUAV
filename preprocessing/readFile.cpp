#include <iostream>
#include <fstream>

using namespace std;

void openFile(string fileName) {
    string rootPath = "../database/";
    string filePath = rootPath + fileName;
    string text;
    ifstream myfile;
    myfile.open(filePath);
    if(myfile.is_open()) {
        while (!myfile.eof()) {
            getline(myfile, text);
            cout << "" << text << "\n";
        }
    } else {
        string errorMessage = "Could not open " + filePath;
        cerr << "Could not open " << filePath << endl;
        ::perror("");
    }


}