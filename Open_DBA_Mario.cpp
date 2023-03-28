#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
struct {
    struct {               // definzione struttura proprietà generali
        double Mass, Wing_spann, Wing_area, Chord, Mach_drag_rise,
                Thrust_axis_offset_x, Thrust_axis_offset_y, Thrust_axis_offset_z,
                Thrust_axis_ang_off_xy, Thrust_axis_ang_off_xz, numb_aoa,
                rotary_deriv, COG, Jx, Jy, jz, jxz, option_cog_update,
                cog_updated, pilot_position_x, pilot_position_y, pilot_position_z;
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
}DBA_0;

int getDba(string filename){
    int i=1;
    int k=0;
    void open(const char *filename);
    string rootPath = "../database/";
    string filePath = rootPath + filename;
    string line;
    ifstream myfile;
    myfile.open(filePath);
    if(myfile.is_open()) {
        while (getline(myfile, line)) {

            if (i==4) {
                istringstream A(line);
                A >> DBA_0.Aircraft_description.Mass;
            }
            if (i==5){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Wing_spann;
            }
            if (i==6){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Wing_area;
            }
            if (i==7){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Chord;
            }
            if (i==8){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Mach_drag_rise;
            }
            if (i==9){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Thrust_axis_offset_x;
            }
            if (i==10){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Thrust_axis_offset_y;
            }
            if (i==11){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Thrust_axis_offset_z;
            }
            if (i==12){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Thrust_axis_ang_off_xy;
            }
            if (i==13){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Thrust_axis_ang_off_xz;
            }
            if (i==14){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.numb_aoa;
            }
            if (i==15){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.rotary_deriv;
            }
            if (i==16){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.COG;
            }
            if (i==17){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Jx;
            }
            if (i==18){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.Jy;
            }
            if (i==19){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.jz;
            }
            if (i==20){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.jxz;
            }
            if (i==21){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.option_cog_update;
            }
            if (i==22){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.cog_updated;
            }
            if (i==23){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.pilot_position_x;
            }
            if (i==24){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.pilot_position_y;
            }
            if (i==25){
                istringstream A(line);
                A >>DBA_0.Aircraft_description.pilot_position_z;
            }
// struttura deflessioni
            if (i==29){
                istringstream A(line);
                A >>DBA_0.Deflection_limits.Elevator_max;
            }
            if (i==30){
                istringstream A(line);
                A >>DBA_0.Deflection_limits.Elevator_min;
            }
            if (i==31){
                istringstream A(line);
                A >>DBA_0.Deflection_limits.Ailerons;
            }
            if (i==32){
                istringstream A(line);
                A >>DBA_0.Deflection_limits.Rudder;
            }
            if (i==33){
                istringstream A(line);
                A >>DBA_0.Deflection_limits.Flap_up;
            }
            if (i==34){
                istringstream A(line);
                A >>DBA_0.Deflection_limits.Flap_down;
            }
// struttura proprietà di massa
            if (i==38){
                istringstream A(line);
                A>> DBA_0.Fuel_Mass.Mass_switch;
            }
            if (i==39){
                istringstream A(line);
                A>> DBA_0.Fuel_Mass.Fuel_weight_fraction;
            }
            i++;
// struttura derivate
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
   // cout.precision(40);
   // cout << fixed <<DBA_0.Aircraft_description.jxz<<'\n';
    return 0;
}
