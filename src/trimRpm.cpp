//
// Created by Matheus Padilha on 29/04/23.
//

#include "../declaredFun.h"
#include <cmath>

struct  Trim_Engine_Propeller{
    double rpm = 0;
    double T = 0;
    double Throttle = 0;
};

Trim_Engine_Propeller trimEnginePropeller(AeroDB db, EngineDB endb, PropDB pdb, Trim_Angles angles, double V, double h){
    Trim_Engine_Propeller enginePerformanceTrim;
    Propel propelResult;

    double rpm;
    double delta_rpm = 100; // [giri/min]
    double rpm_min = endb.laps_min; // [giri/min] 3600
    double rpm_max = endb.laps_max; // [giri/min] 30000

    double rho = computeDensity(h); // [kg/m^3] 1.2132

    double alpha_trim = angles.alpha_trim; // [deg] 2.36
    double deltae_trim= angles.deltae_trim; // [deg]-2.16
    double g = 9.81; // [m/s^2]
    double gamma_0 = 0;
    double res = 1; // [N]

    double S = db.Ad.Wing_area; //0.24704
    double cx_alpha = linearInterpolation(db.alpha, db.fx.cx_a, alpha_trim); //0.2115
    double cx_de = linearInterpolation(db.alpha, db.cf.cx_de, alpha_trim); // 0.04029
    double cx_ss = linearInterpolation(db.alpha, db.ss.cx, alpha_trim); // -0.010595
    double Cx_tot= cx_ss+cx_alpha*alpha_trim/180*M_PI+ cx_de*deltae_trim/180*M_PI;
    enginePerformanceTrim.T = db.Ad.Mass*g*sin(alpha_trim/180*M_PI+gamma_0/180*M_PI)-0.5*Cx_tot*rho*S*V*V; //mass 0.9726

    for(rpm = rpm_min; rpm <= rpm_max; rpm += delta_rpm){

        propelResult = getPropellerPerformance(db, endb, pdb, angles.alpha_trim, angles.deltae_trim, V, h, rpm);

        if(abs(enginePerformanceTrim.T-propelResult.T) < res){
            enginePerformanceTrim.rpm = rpm;
        }
    }


    double n = enginePerformanceTrim.rpm / 60; // [giri/s]
    double omega = n * 2 * M_PI; // [rad/s]
    double P_max = 160; // [W]
    double P_d = propelResult.Q * omega / 1000; // [W] Q = Torque

    if (P_d < P_max){
        enginePerformanceTrim.Throttle = P_d / P_max;
    }

    return enginePerformanceTrim;
}