#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

struct  {               // definzione struttura proprietà generali
    double  Mass, Wing_spann, Wing_area, Chord, Mach_drag_rise,
            Thrust_axis_offset_x, Thrust_axis_offset_y, Thrust_axis_offset_z,
            Thrust_axis_ang_off_xy, Thrust_axis_ang_off_xz, numb_aoa,
            rotary_deriv, COG, Jx, Jy, jz, jxz, option_cog_update,
            cog_updated, pilot_position_x, pilot_position_y, pilot_position_z;
}Aircraft_description;
struct {                // definizione struttura deflessioni
    double Elevator_max, Elevator_min, Ailerons, Rudder, Flap_up, Flap_down;
}Deflection_limits;
struct {                // definizione struttura proprietà di massa
    double Mass_switch, Fuel_weight_fraction;
}Fuel_Mass;
struct {                // definizione struttura derivate
   vector<double> Alpha;
}Steady_state_Coefficients;

int getDba(string filename){
    int i=1;
    int k=1;
    void open(const char *filename,ios::openmode mode);
    string rootPath = "../database/";
    string filePath = rootPath + filename;
    string line;
    ifstream myfile;
    myfile.open(filePath);
    if(myfile.is_open()) {
        while (getline(myfile, line)) {

            if (i==4) {
                istringstream A(line);
                A >> Aircraft_description.Mass;
                // cout <<line<< '\n';
                cout << "" << line << "\n";
            }
            if (i==5){
                istringstream A(line);
                A >>Aircraft_description.Wing_spann;
            }
            if (i==6){
                istringstream A(line);
                A >>Aircraft_description.Wing_area;
            }
            if (i==7){
                istringstream A(line);
                A >>Aircraft_description.Chord;
            }
            if (i==8){
                istringstream A(line);
                A >>Aircraft_description.Mach_drag_rise;
            }
            if (i==9){
                istringstream A(line);
                A >>Aircraft_description.Thrust_axis_offset_x;
            }
            if (i==10){
                istringstream A(line);
                A >>Aircraft_description.Thrust_axis_offset_y;
            }
            if (i==11){
                istringstream A(line);
                A >>Aircraft_description.Thrust_axis_offset_z;
            }
            if (i==12){
                istringstream A(line);
                A >>Aircraft_description.Thrust_axis_ang_off_xy;
            }
            if (i==13){
                istringstream A(line);
                A >>Aircraft_description.Thrust_axis_ang_off_xz;
            }
            if (i==14){
                istringstream A(line);
                A >>Aircraft_description.numb_aoa;
            }
            if (i==15){
                istringstream A(line);
                A >>Aircraft_description.rotary_deriv;
            }
            if (i==16){
                istringstream A(line);
                A >>Aircraft_description.COG;
            }
            if (i==17){
                istringstream A(line);
                A >>Aircraft_description.Jx;
            }
            if (i==18){
                istringstream A(line);
                A >>Aircraft_description.Jy;
            }
            if (i==19){
                istringstream A(line);
                A >>Aircraft_description.jz;
            }
            if (i==20){
                istringstream A(line);
                A >>Aircraft_description.jxz;
            }
            if (i==21){
                istringstream A(line);
                A >>Aircraft_description.option_cog_update;
            }
            if (i==22){
                istringstream A(line);
                A >>Aircraft_description.cog_updated;
            }
            if (i==23){
                istringstream A(line);
                A >>Aircraft_description.pilot_position_x;
            }
            if (i==24){
                istringstream A(line);
                A >>Aircraft_description.pilot_position_y;
            }
            if (i==25){
                istringstream A(line);
                A >>Aircraft_description.pilot_position_z;
            }
// struttura deflessioni
            if (i==29){
                istringstream A(line);
                A >>Deflection_limits.Elevator_max;
            }
            if (i==30){
                istringstream A(line);
                A >>Deflection_limits.Elevator_min;
            }
            if (i==31){
                istringstream A(line);
                A >>Deflection_limits.Ailerons;
            }
            if (i==32){
                istringstream A(line);
                A >>Deflection_limits.Rudder;
            }
            if (i==33){
                istringstream A(line);
                A >>Deflection_limits.Flap_up;
            }
            if (i==34){
                istringstream A(line);
                A >>Deflection_limits.Flap_down;
            }
// struttura proprietà di massa
            if (i==38){
                istringstream A(line);
                A>> Fuel_Mass.Mass_switch;
            }
            if (i==39){
                istringstream A(line);
                A>> Fuel_Mass.Fuel_weight_fraction;
            }
            i++;
// struttura derivate
            if (i>44 && i<150){
                istringstream a(line);
                a >> Steady_state_Coefficients.Alpha[k];
                k++;
            }
        }

    } else {
        string errorMessage = "Could not open " + filePath;
        cerr << "Could not open " << filePath << endl;
        ::perror("");
    }
    //cout << Aircraft_description.Mass<< '\n';
    //cout<< Aircraft_description.Mach_drag_rise;
    cout<< Steady_state_Coefficients.Alpha[k]<< '\n';
    cout<< k<<'\n';
    myfile.close();


    return 0;
}
