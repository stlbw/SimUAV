#include "../declaredFun.h"
#include <cmath>

struct  Trim_Engine_Propeller{
    double rpm = 0;
    double T = 0;
    double Throttle = 0;
    double Torque = 0;
};

double getThrottle(double rpm_trim, double rpm_min, double rpm_max){
    double Throttle_min = 0.1;
    double Throttle_max = 1;
    double m = (Throttle_max - Throttle_min) / (rpm_max - rpm_min);
    double Throttle_trim = m * (rpm_trim - rpm_max) + Throttle_max;

    return Throttle_trim;
}

double getRpm(double throttle, double rpmMin, double rpmMax) {
    double throttleMin = 0.1;
    double throttleMax = 1;
    double m = (rpmMax - rpmMin) / (throttleMax - throttleMin);
    double rpm = m * (throttle - throttleMin) + rpmMin ;

    return rpm;
}

Trim_Engine_Propeller trimEnginePropeller(AeroDB db1, AeroDB db2, EngineDB endb, PropDB pdb, Trim_Angles angles, double V, double h, double gamma_0, double delta_rpm){
    Trim_Engine_Propeller enginePerformanceTrim;
    Propel propelResult;

    double rpm;
    double rpm_min = endb.laps_min; // [rpm] 3600
    double rpm_max = endb.laps_max; // [rpm] 30000

    double rho = computeDensity(h); // [kg/m^3]

    double alpha_trim = angles.alpha_trim; // [deg]
    double deltae_trim = angles.deltae_trim; // [deg]
    double g = 9.81; // [m/s^2]
    double res = 0.02; // [N]

    double S = db1.Ad.Wing_area; // [m^2]
    double cx_alpha = linearInterpolation(db1.alpha, db1.fx.cx_a, db2.fx.cx_a, alpha_trim, h);
    double cx_de = linearInterpolation(db1.alpha, db1.cf.cx_de, db2.cf.cx_de, alpha_trim, h);
    double cx_ss = linearInterpolation(db1.alpha, db1.ss.cx, db2.ss.cx, alpha_trim, h);
    double Cx_tot = cx_ss+cx_alpha*alpha_trim/180*M_PI+ cx_de * deltae_trim / 180*M_PI;
    enginePerformanceTrim.T = db1.Ad.Mass * g * sin(alpha_trim / 180 * M_PI+gamma_0 / 180 * M_PI) - 0.5 * Cx_tot * rho * S * V * V;

    for (rpm = rpm_min; rpm <= rpm_max; rpm += delta_rpm){

        propelResult = getPropellerPerformance(db1, db2, endb, pdb, angles.alpha_trim, angles.deltae_trim, V, h, rpm);

        if (abs(enginePerformanceTrim.T - propelResult.T) < res){
            enginePerformanceTrim.rpm = rpm;
        }
    }


    double n = enginePerformanceTrim.rpm / 60; // [rpm]
    double omega = n * 2 * M_PI; // [rad/s]
    double P_max = 160; // [W]
    double P_d = propelResult.Q * omega / 1000; // [W]
    enginePerformanceTrim.Torque = propelResult.Q;

    enginePerformanceTrim.Throttle = getThrottle(enginePerformanceTrim.rpm, rpm_min, rpm_max);

    return enginePerformanceTrim;
}

