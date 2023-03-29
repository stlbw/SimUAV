#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
struct {
    struct {               //define all data description with an irreal value to use it as flag
        double Mass=-100, Wing_spann=-100, Wing_area=-100, Chord=-100, Mach_drag_rise=-100,
                Thrust_axis_offset_x=-100, Thrust_axis_offset_y=-100, Thrust_axis_offset_z=-100,
                Thrust_axis_ang_off_xy=-100, Thrust_axis_ang_off_xz=-100, numb_aoa=-100,
                rotary_deriv=-100, COG=-100, Jx=-100, Jy=-100, jz=-100, jxz=-100, option_cog_update=-100,
                cog_updated=-100, pilot_position_x=-100, pilot_position_y=-100, pilot_position_z=-100;
    } Aircraft_description;
    struct {                // definizione struttura deflessioni
        double Elevator_max, Elevator_min, Ailerons, Rudder, Flap_up, Flap_down;
    } Deflection_limits;
    struct {                // definizione struttura proprietà di massa
        double Mass_switch, Fuel_weight_fraction;
    } Fuel_Mass;
    double Alpha[106];
    struct {
        double CX[106], CY[106], CZ[106], CL[106], CM[106], CN[106];
    } Steady_state_Coefficients;
    struct {
        struct {
            double CXA[106], CXAP[106], CXU[106], CXQ[106], CXB[106], CXP[106], CXR[106];
        } Xforce;
        struct {
            double CYB[106], CYBP[106], CYP[106], CYR[106], CYA[106], CYQ[106];
        } Yforce;
        struct {
            double CZA[106], CZAP[106], CZU[106], CZQ[106], CZB[106], CZP[106], CZR[106];
        } Zforce;
        struct {
            double CLB[106], CLBP[106], CLP[106], CLR[106], CLA[106], CLQ[106];
        } RollMoment;
        struct {
            double CMA[106], CMAP[106], CMU[106], CMQ[106], CMB[106], CMP[106], CMR[106];
        } PitchMoment;
        struct {
            double CNB[106], CNBP[106], CNP[106], CNR[106], CNA[106], CNQ[106];
        } YawMoment;
        struct{
            double CXDE[106],CXDLE[106],CZDE[106],CZDLE[106],CYDA[106],CYDR[106];
        }ControlForce;
        struct{
            double CLDA[106],CLDR[106],CMDE[106],CMDLE[106],CNDA[106],CNDR[106];
        }ControlMoment;
        struct{
            double CXOM[106],CYOM[106],CZOM[106],CLOM[106],CMOM[106],CNOM[106];
        }Rotary;
    } Derivate;
}DBA_0;    //Struct defintion

int getDba(string  filename){
    int i=1;
    int k=0;
    void open(const char *filename);
    string rootPath = "../database/";
    string filePath = rootPath + filename;
    string line;
    ifstream myfile;
    myfile.open(filePath);
    if(myfile.is_open()) {
        while (!myfile.eof()) {
            (getline(myfile, line));
            if (line.find("MASS")!= string::npos & DBA_0.Aircraft_description.Mass==-100){
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Mass;
            }
            if (line.find("WING SPANN")!= string::npos & DBA_0.Aircraft_description.Wing_spann==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Wing_spann;
            }
            if (line.find("WING AREA")!= string::npos & DBA_0.Aircraft_description.Wing_area==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Wing_area;
            }
            if (line.find("CHORD")!= string::npos & DBA_0.Aircraft_description.Chord==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Chord;
            }
            if (line.find("MACH DRAG RISE")!= string::npos & DBA_0.Aircraft_description.Mach_drag_rise==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Mach_drag_rise;
            }
            if (line.find("THRUST AXIS OFFSET (REF. TO XB ALONG X)")!= string::npos & DBA_0.Aircraft_description.Thrust_axis_offset_x==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Thrust_axis_offset_x;
            }
            if (line.find("THRUST AXIS OFFSET (REF. TO XB ALONG Y)")!= string::npos & DBA_0.Aircraft_description.Thrust_axis_offset_y==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Thrust_axis_offset_y;
            }
            if (line.find("THRUST AXIS OFFSET (REF. TO XB ALONG Z)")!= string::npos & DBA_0.Aircraft_description.Thrust_axis_offset_z==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Thrust_axis_offset_z;
            }
            if (line.find("THRUST AXIS ANGULAR OFFSET (REF. TO XB / X-Y PLANE / POSITIVE RIGHT)")!= string::npos & DBA_0.Aircraft_description.Thrust_axis_ang_off_xy==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Thrust_axis_ang_off_xy;
            }
            if (line.find("THRUST AXIS ANGULAR OFFSET (REF. TO XB / X-Z PLANE / POSITIVE RIGHT)")!= string::npos & DBA_0.Aircraft_description.Thrust_axis_ang_off_xz==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Thrust_axis_ang_off_xz;
            }
            if (line.find("NUMBER OF ANGLES OF ATTACK")!= string::npos & DBA_0.Aircraft_description.numb_aoa==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.numb_aoa;
            }
            if (line.find("ROTARY DERIVATIVES")!= string::npos & DBA_0.Aircraft_description.rotary_deriv==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.rotary_deriv;
            }
            if (line.find("CENTER OF GRAVITY REFERENCE LOCATION REF. TO CMAER")!= string::npos & DBA_0.Aircraft_description.COG==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.COG;
            }
            if (line.find("JX")!= string::npos & DBA_0.Aircraft_description.Jx==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Jx;
            }
            if (line.find("JY")!= string::npos & DBA_0.Aircraft_description.Jy==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Jy;
            }
            if (line.find("JZ")!= string::npos & DBA_0.Aircraft_description.jz==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.jz;
            }
            if (line.find("JXZ")!= string::npos & DBA_0.Aircraft_description.jxz==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.jxz;
            }
            if (line.find("OPTION FOR C.G. UPDATE (0/NO 1/YES)")!= string::npos & DBA_0.Aircraft_description.option_cog_update==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.option_cog_update;
            }
            if (line.find("CENTER OF GRAVITY REFERENCE LOCATION (UPDATED)")!= string::npos & DBA_0.Aircraft_description.cog_updated==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.cog_updated;
            }
            if (line.find("PILOT POSITION (REF. TO CG ALONG XB)")!= string::npos & DBA_0.Aircraft_description.pilot_position_x==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.pilot_position_x;
            }
            if (line.find("PILOT POSITION (REF. TO CG ALONG YB)")!= string::npos & DBA_0.Aircraft_description.pilot_position_y==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.pilot_position_y;
            }
            if (line.find("PILOT POSITION (REF. TO CG ALONG ZB)")!= string::npos & DBA_0.Aircraft_description.pilot_position_z==-100) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.pilot_position_z;
            }
            if (line.find("ELEVATOR (max)")!= string::npos & DBA_0.Deflection_limits.Elevator_max==-100) {
                istringstream A(line);
                A >> DBA_0.Deflection_limits.Elevator_max;
            }
            if (line.find("ELEVATOR (min)")!= string::npos & DBA_0.Deflection_limits.Elevator_min==-100) {
                istringstream A(line);
                A >> DBA_0.Deflection_limits.Elevator_min;
            }
            if (line.find("AILERONS")!= string::npos & DBA_0.Deflection_limits.Ailerons==-100) {
                istringstream A(line);
                A >> DBA_0.Deflection_limits.Ailerons;
            }
            if (line.find("RUDDER")!= string::npos & DBA_0.Deflection_limits.Rudder==-100) {
                istringstream A(line);
                A >> DBA_0.Deflection_limits.Rudder;
            }
            if (line.find("FLAP (up)")!= string::npos & DBA_0.Deflection_limits.Flap_up==-100) {
                istringstream A(line);
                A >> DBA_0.Deflection_limits.Flap_up;
            }
            if (line.find("FLAP (down)")!= string::npos & DBA_0.Deflection_limits.Flap_down==-100) {
                istringstream A(line);
                A >> DBA_0.Deflection_limits.Flap_down;
            }
            if (line.find("MASS SWITCH")!= string::npos & DBA_0.Fuel_Mass.Mass_switch==-100) {
                istringstream A(line);
                A >> DBA_0.Fuel_Mass.Mass_switch;
            }
            if (line.find("FUEL WEIGHT FRACTION")!= string::npos & DBA_0.Fuel_Mass.Fuel_weight_fraction==-100) {
                istringstream A(line);
                A >> DBA_0.Fuel_Mass.Fuel_weight_fraction;
            }
                i++;

            if (i>44 && i<151){
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Steady_state_Coefficients.CX[k]>>
                DBA_0.Steady_state_Coefficients.CY[k]>>DBA_0.Steady_state_Coefficients.CZ[k]>>
                DBA_0.Steady_state_Coefficients.CL[k]>>DBA_0.Steady_state_Coefficients.CM[k]>>
                DBA_0.Steady_state_Coefficients.CN[k];
                k++;
            }
            if (i>156 && i<263) {
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Derivate.Xforce.CXA[k]>>DBA_0.Derivate.Xforce.CXAP[k]>>
                DBA_0.Derivate.Xforce.CXU[k]>>DBA_0.Derivate.Xforce.CXQ[k]>>DBA_0.Derivate.Xforce.CXB[k]>>
                DBA_0.Derivate.Xforce.CXP[k]>>DBA_0.Derivate.Xforce.CXR[k];
            k++;
            }
            if (i>266 && i<373) {
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Derivate.Yforce.CYB[k]>>DBA_0.Derivate.Yforce.CYBP[k]>>
                DBA_0.Derivate.Yforce.CYP[k]>>DBA_0.Derivate.Yforce.CYR[k]>>DBA_0.Derivate.Yforce.CYA[k]>>DBA_0.Derivate.Yforce.CYQ[k];
            k++;
            }
            if (i>376 && i<483) {
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Derivate.Zforce.CZA[k]>>DBA_0.Derivate.Zforce.CZAP[k]>>DBA_0.Derivate.Zforce.CZU[k]>>
                DBA_0.Derivate.Zforce.CZQ[k]>>DBA_0.Derivate.Zforce.CZB[k]>>DBA_0.Derivate.Zforce.CZP[k]>>DBA_0.Derivate.Zforce.CZR[k];
                k++;
            }
            if (i>486 && i<593) {
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Derivate.RollMoment.CLB[k]>>DBA_0.Derivate.RollMoment.CLBP[k]>>DBA_0.Derivate.RollMoment.CLP[k]>>
                DBA_0.Derivate.RollMoment.CLR[k]>>DBA_0.Derivate.RollMoment.CLA[k]>>DBA_0.Derivate.RollMoment.CLQ[k];
            k++;
            }
            if (i>596 && i<703) {
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Derivate.PitchMoment.CMA[k]>>DBA_0.Derivate.PitchMoment.CMAP[k]>>DBA_0.Derivate.PitchMoment.CMU[k]>>
                DBA_0.Derivate.PitchMoment.CMQ[k]>>DBA_0.Derivate.PitchMoment.CMB[k]>>DBA_0.Derivate.PitchMoment.CMP[k]>>DBA_0.Derivate.PitchMoment.CMR[k];
            k++;
            }
            if (i>706 && i<813) {
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Derivate.YawMoment.CNB[k]>>DBA_0.Derivate.YawMoment.CNBP[k]>>DBA_0.Derivate.YawMoment.CNP[k]>>
                DBA_0.Derivate.YawMoment.CNR[k]>>DBA_0.Derivate.YawMoment.CNA[k]>>DBA_0.Derivate.YawMoment.CNQ[k];
            k++;
            }
            if (i>816 && i<923) {
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Derivate.ControlForce.CXDE[k]>>DBA_0.Derivate.ControlForce.CXDLE[k]>>DBA_0.Derivate.ControlForce.CZDE[k]>>
                DBA_0.Derivate.ControlForce.CZDLE[k]>>DBA_0.Derivate.ControlForce.CYDA[k]>>DBA_0.Derivate.ControlForce.CYDR[k];
            k++;
            }
            if (i>926 && i<1033) {
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Derivate.ControlMoment.CLDA[k]>>DBA_0.Derivate.ControlMoment.CLDR[k]>>DBA_0.Derivate.ControlMoment.CMDE[k]>>
                DBA_0.Derivate.ControlMoment.CMDLE[k]>>DBA_0.Derivate.ControlMoment.CNDA[k]>>DBA_0.Derivate.ControlMoment.CLDR[k];
            k++;
            }
            if (i>1036 && i<1143) {
                istringstream a(line);
                a >> DBA_0.Alpha[k]>>DBA_0.Derivate.Rotary.CXOM[k]>>DBA_0.Derivate.Rotary.CYOM[k]>>DBA_0.Derivate.Rotary.CZOM[k]>>
                DBA_0.Derivate.Rotary.CLOM[k]>>DBA_0.Derivate.Rotary.CMOM[k]>>DBA_0.Derivate.Rotary.CNOM[k];
            k++;
            }
            if (k==106){
                k=0;
            }
            }
        }
    else {
        string errorMessage = "Could not open " + filePath;
        cerr << "Could not open " << filePath << endl;
        ::perror("");
    }
    myfile.close();

    return 0;
}
