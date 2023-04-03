#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

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

inline void save_propeller(PropDB *db, string filePath) {
    int flag_ps = 0;
    string text;
    ifstream myfile;
    myfile.open(filePath); // to open the file
    // if the file is correctly open...
    if(myfile.is_open()) {
        std::vector<std::string> nameVar {""}; // declare vector to store variables
        cout << "Reading database " << filePath << endl; // print what the program is reading
        // read line by line and based on the keyword save variables accordingly
        while (!myfile.eof()) {
            getline(myfile, text);
            // read and save all aircraft description data
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
            } else if (text.find("Cdalfa [rad^-1]") != string::npos) {
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
                if (text.find("CSI") != string::npos) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                } else if (text.find("********************************") == string::npos &&
                           text.find("CSI") == string::npos) {
                    string token;
                    stringstream s(text);
                    int i = 0;
                    while (getline(s, token, ' ')) {
                        if (!token.empty()) {
                            /*if(token.find("\r") != string::npos) {
                                string lineSkip = "\r";
                                string::size_type i = token.find(lineSkip);
                                token.erase(i, lineSkip.length());
                            }*/
                            if (nameVar[i].find("CSI") != string::npos) {
                                db->Ps.CSI.push_back(stof(token));
                                i++;
                            } else if (nameVar[i].find("RD [m]") != string::npos) {
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
            cout << "\tFinished reading database.";
        }
    }else {
        string errorMessage = "Could not open " + filePath;
        cerr << "Could not open " << filePath << endl;
        ::perror("");
    }
    myfile.close();
}

//PropDB readData1(string fileName) {
    //string rootPath = "../database/";
    //string filePath = rootPath + fileName;
    //PropDB db; // create db object
    //Flags check; // create check object from struct Flag
    //int length = getAoALength(filePath); // get vector size to initiate struct
    //db.length = length;
    //save_propeller(&db, filePath);
    //return db;
//}