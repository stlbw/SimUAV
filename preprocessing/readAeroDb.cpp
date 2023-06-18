#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include "../declaredFun.h"
using namespace std;
//dclare structures to achieved nested structure from the most inner to the most outer for the static database
struct Aircraft_description {               //define all data description with an irreal value to use it as flag
    double Mass=-100, Wing_spann=-100, Wing_area=-100, Chord=-100, Mach_drag_rise=-100,
            Thrust_axis_offset_x=-100, Thrust_axis_offset_y=-100, Thrust_axis_offset_z=-100,
            Thrust_axis_ang_off_xy=-100, Thrust_axis_ang_off_xz=-100, numb_aoa=-100,
            rotary_deriv=-100, COG=-100, Jx=-100, Jy=-100, jz=-100, jxz=-100, option_cog_update=-100,
            cog_updated=-100, pilot_position_x=-100, pilot_position_y=-100, pilot_position_z=-100;
} ;
struct Deflection_limits {
    double Elevator_max = -100, Elevator_min = -100, Ailerons = -100, Rudder = -100, Flap_up = -100, Flap_down = -100;
} ;
struct Fuel_Mass {
    double Mass_switch = -100, Fuel_weight_fraction = -100;
} ;

struct  SSCoef {
    vector <double> cx;
    vector <double> cy;
    vector <double> cz;
    vector <double> cl;
    vector <double> cm;
    vector <double> cn;
};
struct  XForce {
    vector <double> cx_a ;
    vector <double> cx_ap;
    vector <double> cx_u ;
    vector <double> cx_q ;
    vector <double> cx_b ;
    vector <double> cx_p ;
    vector <double> cx_r ;
};
struct  YForce {
    vector <double> cy_b ;
    vector <double> cy_bp;
    vector <double> cy_p ;
    vector <double> cy_r ;
    vector <double> cy_a ;
    vector <double> cy_q ;
};
struct  ZForce {
    vector <double> cz_a ;
    vector <double> cz_ap;
    vector <double> cz_u ;
    vector <double> cz_q ;
    vector <double> cz_b ;
    vector <double> cz_p ;
    vector <double> cz_r ;

};
struct  RollMoment {
    vector <double> cl_b ;
    vector <double> cl_bp;
    vector <double> cl_p ;
    vector <double> cl_r ;
    vector <double> cl_a ;
    vector <double> cl_q ;

};
struct  PitchMoment {
    vector <double> cm_a ;
    vector <double> cm_ap;
    vector <double> cm_u ;
    vector <double> cm_q ;
    vector <double> cm_b ;
    vector <double> cm_p ;
    vector <double> cm_r ;

};
struct  YawMoment {
    vector <double> cn_b ;
    vector <double> cn_bp;
    vector <double> cn_p ;
    vector <double> cn_r ;
    vector <double> cn_a ;
    vector <double> cn_q ;
};
struct  ControlForce {
    vector <double> cx_de ;
    vector <double> cx_dle;
    vector <double> cz_de ;
    vector <double> cz_dle;
    vector <double> cy_da ;
    vector <double> cy_dr ;
};
struct  ControlMoment {
    vector <double> cl_da ;
    vector <double> cl_dr ;
    vector <double> cm_de ;
    vector <double> cm_dle;
    vector <double> cn_da ;
    vector <double> cn_dr ;
};
struct  Rotary {
    vector <double> cx_om;
    vector <double> cy_om;
    vector <double> cz_om;
    vector <double> cl_om;
    vector <double> cm_om;
    vector <double> cn_om;
};

struct AeroDB {
    int length;
    vector <double> alpha;
    Aircraft_description Ad;                    // Aircraft description
    Deflection_limits Dl;                       // Deflection limits
    Fuel_Mass Fm;                               //Fuel Mass
    SSCoef ss;                                  // steady state coefficients
    XForce fx;                                  // X Forces (Aerodynamic derivatives)
    YForce fy;                                  // Y Forces (Aerodynamic derivatives)
    ZForce fz;                                  // Y Forces (Aerodynamic derivatives)
    RollMoment rm;                              // Roll moment (Aerodynamic derivatives)
    PitchMoment pm;                             // Pitch moment (Aerodynamic derivatives)
    YawMoment ym;                               // Yaw moment (Aerodynamic derivatives)
    ControlForce cf;                            // Control force derivatives
    ControlMoment cm;                           // Control moment derivatives
    Rotary rt;                                  // Rotary derivatives
};
struct Flags {
    int flag_ss, flag_fx, flag_fy, flag_fz, flag_roll, flag_pitch, flag_yaw, flag_cf, flag_cm, flag_rotary;
    };

/** Gets the sampling size of angle of attack from file*/
int getAoALength(string filePath) {
    ifstream myfile;
    myfile.open(filePath);
    if(myfile.is_open()) {
        string text;
        while (!myfile.eof()) {
            getline(myfile, text);
            bool found = text.find("NUMBER OF ANGLES OF ATTACK") != string::npos; // find where AoA is written to allocate memory for the variables
            if(found) {
                string token;
                int vecLen;
                stringstream s(text);
                while (getline(s, token, ' ')) {
                    if (!token.empty()) {
                        vecLen = stoi(token);
                        myfile.close(); // close file
                        return vecLen;
                    }
                }
                cout << "" << text << "\n"; //prints out the database lines
            }
        }
        myfile.close();
    } else {
        string errorMessage = "Could not open " + filePath;
        throw runtime_error(errorMessage);
    }
}

/** Receives the datastructure containing the flags and initialize them according to the file.
 */
void updateFlag(Flags *f, string text) {
    if(text.find("STEADY STATE COEFFICIENTS") != string::npos) {
        f->flag_ss = 1;
        f->flag_fx = f->flag_fy = f->flag_fz = f->flag_roll = f->flag_pitch = f->flag_yaw = f->flag_cf = f->flag_cm = f->flag_rotary = 0;
    }
    else if(text.find("X  FORCE DERIVATIVES") != string::npos) {
        f->flag_fx = 1;
        f->flag_ss = f->flag_fy = f->flag_fz = f->flag_roll = f->flag_pitch = f->flag_yaw = f->flag_cf = f->flag_cm = f->flag_rotary = 0;
    }
    else if(text.find("Y  FORCE DERIVATIVES") != string::npos) {
        f->flag_fy = 1;
        f->flag_fx = f->flag_ss = f->flag_fz = f->flag_roll = f->flag_pitch = f->flag_yaw = f->flag_cf = f->flag_cm = f->flag_rotary = 0;
    }
    else if(text.find("Z  FORCE DERIVATIVES") != string::npos) {
        f->flag_fz = 1;
        f->flag_fx = f->flag_fy = f->flag_ss = f->flag_roll = f->flag_pitch = f->flag_yaw = f->flag_cf = f->flag_cm = f->flag_rotary = 0;
    }
    else if(text.find("ROLLING MOMENT DERIVATIVES") != string::npos) {
        f->flag_roll = 1;
        f->flag_fx = f->flag_fy = f->flag_fz = f->flag_ss = f->flag_pitch = f->flag_yaw = f->flag_cf = f->flag_cm = f->flag_rotary = 0;
    }
    else if(text.find("PITCHING MOMENT DERIVATIVES") != string::npos) {
        f->flag_pitch = 1;
        f->flag_fx = f->flag_fy = f->flag_fz = f->flag_roll = f->flag_ss = f->flag_yaw = f->flag_cf = f->flag_cm = f->flag_rotary = 0;
    }
    else if(text.find("YAWING MOMENT DERIVATIVES") != string::npos) {
        f->flag_yaw = 1;
        f->flag_fx = f->flag_fy = f->flag_fz = f->flag_roll = f->flag_pitch = f->flag_ss = f->flag_cf = f->flag_cm = f->flag_rotary = 0;
    }
    else if(text.find("CONTROL FORCE DERIVATIVES") != string::npos) {
        f->flag_cf = 1;
        f->flag_fx = f->flag_fy = f->flag_fz = f->flag_roll = f->flag_pitch = f->flag_yaw = f->flag_ss = f->flag_cm = f->flag_rotary = 0;
    }
    else if(text.find("CONTROL MOMENT DERIVATIVES") != string::npos) {
        f->flag_cm = 1;
        f->flag_fx = f->flag_fy = f->flag_fz = f->flag_roll = f->flag_pitch = f->flag_yaw = f->flag_cf = f->flag_ss = f->flag_rotary = 0;
    }
    else if(text.find("ROTARY DERIVATIVES") != string::npos) {
        f->flag_rotary = 1;
        f->flag_fx = f->flag_fy = f->flag_fz = f->flag_roll = f->flag_pitch = f->flag_yaw = f->flag_cf = f->flag_cm = f->flag_ss = 0;
    }
    else {};
}

/** Reads and saves aerodynamic database to a nested struct
 * @param db
 * @param f
 * @param filePath
 */
void saveData(AeroDB *db, Flags *f, string filePath) {
    string text;
    ifstream myfile;
    myfile.open(filePath);
    if(myfile.is_open()) {
        // initiate all flags to 0
        f->flag_ss = f->flag_fx = f->flag_fy = f->flag_fz = 0;
        f->flag_roll = f->flag_pitch = f->flag_yaw = f->flag_cf = 0;
        f->flag_cm = f->flag_rotary = 0;
        std::vector<std::string> nameVar {""}; // declare vector to store variables
        cout << "Reading database " << filePath << endl;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        // read line by line and based on the keyword save variables accordingly
        while (!myfile.eof()) {
            getline(myfile, text);
            // read and save all aircraft description data
            if (text.find("MASS")!= string::npos & db->Ad.Mass==-100){
                istringstream A(text);
                A >> db->Ad.Mass;
            }
            else if (text.find("WING SPAN")!= string::npos & db->Ad.Wing_spann==-100) {
                istringstream A(text);
                A >> db->Ad.Wing_spann;
            }
            else if (text.find("WING AREA")!= string::npos & db->Ad.Wing_area==-100) {
                istringstream A(text);
                A >> db->Ad.Wing_area;
            }
            else if (text.find("CHORD")!= string::npos & db->Ad.Chord==-100) {
                istringstream A(text);
                A >> db->Ad.Chord;
            }
            else if (text.find("MACH DRAG RISE")!= string::npos & db->Ad.Mach_drag_rise==-100) {
                istringstream A(text);
                A >> db->Ad.Mach_drag_rise;
            }
            else if (text.find("THRUST AXIS OFFSET (REF. TO XB ALONG X)")!= string::npos & db->Ad.Thrust_axis_offset_x==-100) {
                istringstream A(text);
                A >> db->Ad.Thrust_axis_offset_x;
            }
            else if (text.find("THRUST AXIS OFFSET (REF. TO XB ALONG Y)")!= string::npos & db->Ad.Thrust_axis_offset_y==-100) {
                istringstream A(text);
                A >> db->Ad.Thrust_axis_offset_y;
            }
            else if (text.find("THRUST AXIS OFFSET (REF. TO XB ALONG Z)")!= string::npos & db->Ad.Thrust_axis_offset_z==-100) {
                istringstream A(text);
                A >> db->Ad.Thrust_axis_offset_z;
            }
            else if (text.find("THRUST AXIS ANGULAR OFFSET (REF. TO XB / X-Y PLANE / POSITIVE RIGHT)")!= string::npos & db->Ad.Thrust_axis_ang_off_xy==-100) {
                istringstream A(text);
                A >> db->Ad.Thrust_axis_ang_off_xy;
            }
            else if (text.find("THRUST AXIS ANGULAR OFFSET (REF. TO XB / X-Z PLANE / POSITIVE DOWN)")!= string::npos & db->Ad.Thrust_axis_ang_off_xz==-100) {
                istringstream A(text);
                A >> db->Ad.Thrust_axis_ang_off_xz;
            }
            else if (text.find("NUMBER OF ANGLES OF ATTACK")!= string::npos & db->Ad.numb_aoa==-100) {
                istringstream A(text);
                A >> db->Ad.numb_aoa;
            }
            else if (text.find("ROTARY DERIVATIVES")!= string::npos & db->Ad.rotary_deriv==-100) {
                istringstream A(text);
                A >> db->Ad.rotary_deriv;
            }
            else if (text.find("CENTER OF GRAVITY REFERENCE LOCATION REF. TO CMAER")!= string::npos & db->Ad.COG==-100) {
                istringstream A(text);
                A >> db->Ad.COG;
            }
            else if (text.find("JX")!= string::npos & db->Ad.Jx==-100) {
                istringstream A(text);
                A >> db->Ad.Jx;
            }
            else if (text.find("JY")!= string::npos & db->Ad.Jy==-100) {
                istringstream A(text);
                A >> db->Ad.Jy;
            }
            else if (text.find("JZ")!= string::npos & db->Ad.jz==-100) {
                istringstream A(text);
                A >> db->Ad.jz;
            }
            else if (text.find("JXZ")!= string::npos & db->Ad.jxz==-100) {
                istringstream A(text);
                A >> db->Ad.jxz;
            }
            else if (text.find("OPTION FOR C.G. UPDATE (0/NO 1/YES)")!= string::npos & db->Ad.option_cog_update==-100) {
                istringstream A(text);
                A >> db->Ad.option_cog_update;
            }
            else if (text.find("CENTER OF GRAVITY REFERENCE LOCATION (UPDATED)")!= string::npos & db->Ad.cog_updated==-100) {
                istringstream A(text);
                A >> db->Ad.cog_updated;
            }
            else if (text.find("PILOT POSITION (REF. TO CG ALONG XB)")!= string::npos & db->Ad.pilot_position_x==-100) {
                istringstream A(text);
                A >> db->Ad.pilot_position_x;
            }
            else if (text.find("PILOT POSITION (REF. TO CG ALONG YB)")!= string::npos & db->Ad.pilot_position_y==-100) {
                istringstream A(text);
                A >> db->Ad.pilot_position_y;
            }
            else if (text.find("PILOT POSITION (REF. TO CG ALONG ZB)")!= string::npos & db->Ad.pilot_position_z==-100) {
                istringstream A(text);
                A >> db->Ad.pilot_position_z;
            }
            else if (text.find("ELEVATOR (max)")!= string::npos & db->Dl.Elevator_max==-100) {
                istringstream A(text);
                A >> db->Dl.Elevator_max;
            }
            else if (text.find("ELEVATOR (min)")!= string::npos & db->Dl.Elevator_min==-100) {
                istringstream A(text);
                A >> db->Dl.Elevator_min;
            }
            else if (text.find("AILERONS")!= string::npos & db->Dl.Ailerons==-100) {
                istringstream A(text);
                A >> db->Dl.Ailerons;
            }
            else if (text.find("RUDDER")!= string::npos & db->Dl.Rudder==-100) {
                istringstream A(text);
                A >> db->Dl.Rudder;
            }
            else if (text.find("FLAP (up)")!= string::npos & db->Dl.Flap_up==-100) {
                istringstream A(text);
                A >> db->Dl.Flap_up;
            }
            else if (text.find("FLAP (down)")!= string::npos & db->Dl.Flap_down==-100) {
                istringstream A(text);
                A >> db->Dl.Flap_down;
            }
            else if (text.find("MASS SWITCH")!= string::npos &db->Fm.Mass_switch==-100) {
                istringstream A(text);
                A >> db->Fm.Mass_switch;
            }
            else if (text.find("FUEL WEIGHT FRACTION")!= string::npos & db->Fm.Fuel_weight_fraction==-100) {
                istringstream A(text);
                A >> db->Fm.Fuel_weight_fraction;
            }
            updateFlag(f, text); // changes each flag's status to cope with database, giving information on
            // which section is being currently read and saved
            if(f->flag_ss == 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
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
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() <= db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("CX") != string::npos) {
                                db->ss.cx.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CY") != string::npos) {
                                db->ss.cy.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CZ") != string::npos) {
                                db->ss.cz.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CL") != string::npos) {
                                db->ss.cl.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CM") != string::npos) {
                                db->ss.cm.push_back(stof(token));
                                i++;
                            }
                            else if(nameVar[i].find("CN") != string::npos) {
                                db->ss.cn.push_back(stof(token));
                                i++;
                            }

                            //db->length = stoi(token);
                        }
                    }
                }
            } 
            else if(f->flag_fx == 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
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
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() < db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("CXA") != string::npos && nameVar[i].find("CXAP") == string::npos) {
                                db->fx.cx_a.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CXAP") != string::npos) {
                                db->fx.cx_ap.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CXU") != string::npos) {
                                db->fx.cx_u.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CXQ") != string::npos) {
                                db->fx.cx_q.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CXB") != string::npos) {
                                db->fx.cx_b.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CXP") != string::npos) {
                                db->fx.cx_p.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CXR") != string::npos) {
                                db->fx.cx_r.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }
            else if(f->flag_fy == 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }
                    else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
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
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() < db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("CYB") != string::npos && nameVar[i].find("CYBP") == string::npos) {
                                db->fy.cy_b.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CYBP") != string::npos) {
                                db->fy.cy_bp.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CYP") != string::npos) {
                                db->fy.cy_p.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CYR") != string::npos) {
                                db->fy.cy_r.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CYA") != string::npos) {
                                db->fy.cy_a.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CYQ") != string::npos) {
                                db->fy.cy_q.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }
            else if(f->flag_fz == 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }
                    else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
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
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() < db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("CZA") != string::npos && nameVar[i].find("CZAP") == string::npos) {
                                db->fz.cz_a.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CZAP") != string::npos) {
                                db->fz.cz_ap.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CZU") != string::npos) {
                                db->fz.cz_u.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CZQ") != string::npos) {
                                db->fz.cz_q.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CZB") != string::npos) {
                                db->fz.cz_b.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CZP") != string::npos) {
                                db->fz.cz_p.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CZR") != string::npos) {
                                db->fz.cz_r.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }
            else if(f->flag_roll== 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }
                else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
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
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() < db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("CLB") != string::npos && nameVar[i].find("CLBP") == string::npos) {
                                db->rm.cl_b.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CLBP") != string::npos) {
                                db->rm.cl_bp.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CLR") != string::npos) {
                                db->rm.cl_r.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CLP") != string::npos) {
                                db->rm.cl_p.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CLA") != string::npos) {
                                db->rm.cl_a.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CLQ") != string::npos) {
                                db->rm.cl_q.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }
            else if(f->flag_pitch == 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }
                else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
                    string token;
                    stringstream s(text);
                    int i = 0;
                    while (getline(s, token, ' ')) {
                        if (!token.empty()) {
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() < db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("CMA") != string::npos && nameVar[i].find("CMAP") == string::npos) {
                                db->pm.cm_a.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CMAP") != string::npos) {
                                db->pm.cm_ap.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CMU") != string::npos) {
                                db->pm.cm_u.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CMQ") != string::npos) {
                                db->pm.cm_q.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CMB") != string::npos) {
                                db->pm.cm_b.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CMP") != string::npos) {
                                db->pm.cm_p.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CMR") != string::npos) {
                                db->pm.cm_r.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }
            else if(f->flag_yaw == 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
                    string token;
                    stringstream s(text);
                    int i = 0;
                    while (getline(s, token, ' ')) {
                        if (!token.empty()) {
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() < db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("CNB") != string::npos && nameVar[i].find("CNBP") == string::npos) {
                                db->ym.cn_b.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CNBP") != string::npos) {
                                db->ym.cn_bp.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CNP") != string::npos) {
                                db->ym.cn_p.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CNR") != string::npos) {
                                db->ym.cn_r.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CNA") != string::npos) {
                                db->ym.cn_a.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CNQ") != string::npos) {
                                db->ym.cn_q.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }
            else if(f->flag_cf == 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
                    string token;
                    stringstream s(text);
                    int i = 0;
                    while (getline(s, token, ' ')) {
                        if (!token.empty()) {
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() < db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("CXDE") != string::npos) {
                                db->cf.cx_de.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CXDLE") != string::npos) {
                                db->cf.cx_dle.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CZDE") != string::npos) {
                                db->cf.cz_de.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CZDLE") != string::npos) {
                                db->cf.cz_dle.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CYDA") != string::npos) {
                                db->cf.cy_da.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CYDR") != string::npos) {
                                db->cf.cy_dr.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }
            else if(f->flag_cm == 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
                    string token;
                    stringstream s(text);
                    int i = 0;
                    while (getline(s, token, ' ')) {
                        if (!token.empty()) {
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() < db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("CLDA") != string::npos) {
                                db->cm.cl_da.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CLDR") != string::npos) {
                                db->cm.cl_dr.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CMDE") != string::npos) {
                                db->cm.cm_de.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CMDLE") != string::npos) {
                                db->cm.cm_dle.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CNDA") != string::npos) {
                                db->cm.cn_da.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("CNDR") != string::npos) {
                                db->cm.cn_dr.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }
            else if(f->flag_rotary == 1) {
                if(text.find("ALPHA") != string::npos ) {
                    nameVar = {}; // clean variable to receive new values
                    string textToken;
                    stringstream st(text);
                    while (getline(st, textToken, ' ')) {
                        if (!textToken.empty()) {
                            nameVar.push_back(textToken);
                        }
                    }
                }
                else if(text.find("*") == string::npos && text.find("ALPHA") == string::npos) {
                    string token;
                    stringstream s(text);
                    int i = 0;
                    while (getline(s, token, ' ')) {
                        if (!token.empty()) {
                            if(nameVar[i].find("ALPHA") != string::npos) {
                                if(db->alpha.size() < db->length) {
                                    db->alpha.push_back(stof(token));
                                }
                                i++;
                            } else if(nameVar[i].find("Cxom") != string::npos) {
                                db->rt.cx_om.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("Cyom") != string::npos) {
                                db->rt.cy_om.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("Czom") != string::npos) {
                                db->rt.cz_om.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("Clom") != string::npos) {
                                db->rt.cl_om.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("Cmom") != string::npos) {
                                db->rt.cm_om.push_back(stof(token));
                                i++;
                            } else if(nameVar[i].find("Cnom") != string::npos) {
                                db->rt.cn_om.push_back(stof(token));
                                i++;
                            }
                        }
                    }
                }
            }
            bool found = text.find("NUMBER OF ANGLES OF ATTACK") !=string::npos; // find where AoA is written to allocate memory for the variables
            //cout << "" << text << "\n"; //prints out the database lines
        }
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        cout << "\tFinished reading database. Time took: " << chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms." << endl;

    } else {
        string errorMessage = "Could not open " + filePath;
        throw runtime_error(errorMessage);
    }
}

/** Read a string containing a file name and its extension and check for its existence on the database directory.
 If the file is contained in database, print its content to the screen.
 */
AeroDB readData(string fileName, bool t = false) {
    string filePath = string(DATABASE_DIR) + "/" + fileName;
    AeroDB db; // create db object
    Flags check; // create check object from struct Flag
    int length = getAoALength(filePath); // get vector size to initiate struct
    db.length = length;
    saveData(&db, &check, filePath);
    return db;


}


