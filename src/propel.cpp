//
// Created by Matheus Padilha on 29/04/23.
//
#include "../declaredFun.h"
#include <cmath>
#include<iostream>

struct Propel {
    double T = 0; //thrust
    double Q = 0; // torque
    double ct = 0; // thrust coefficient
    double cq = 0; //torque coefficient
    double J = 0; // advance ratio
    double eta = -1; // propeller's efficiency - initialized as -1 as 0 is a possible value
};
/**
 * Computes propeller performance using a simpolified Glauert's Blade Element Theory for each station of the propeller
 * defined in the propeller database and returns a struct of type Propel
 * @param db
 * @param endb
 * @param pdb
 * @param angles
 * @param V
 * @param h
 * @param rpm
 * @return
 */
Propel getPropellerPerformance(AeroDB db1, AeroDB db2, EngineDB endb, PropDB pdb, double alpha, double delta_e, double V, double h, double rpm) {
    Propel result;
    int lenVec = pdb.Pg.nstation; // [-] 100
    double diam = pdb.Pg.diameter; // [m] 0.125
    double diam_ogiva = pdb.Pg.diameter_ogive; // [m] 0.0125
    double radius = diam/2; // [m]
    vector<double> chord;
    chord.assign(lenVec+1,0);
    for (int i = 0; i < lenVec+1; ++i) {
        chord[i] = pdb.Ps.CH_AD[i]*radius; // [m]
    }
    double pitch_propeller = 10; // [deg]
    double delta_rpm = 100; // [rpm]
    double rpm_min = endb.laps_min; // [rpm] 3600
    double rpm_max = endb.laps_max; // [rpm] 30000
    double pitch_tip = pdb.Ps.BA.back(); // [deg] 9.2322
    double pitch_hub = pdb.Ps.BA.front(); // [deg] 39.1010
    double xt = radius; // [m]
    vector<double> CSI = pdb.Ps.CSI; // [-]
    double xs = CSI[0]*radius; // [m]

    double rho = computeDensity(h); // [kg/m^3] 1.2132

    double coef1 = (pitch_tip-pitch_hub)/(xt-xs); // [deg/m]
    double coef2 = pitch_hub-coef1*xs+pitch_propeller; // [deg]
    double r_step = (xt-xs)/lenVec; // [-]
    vector<double> r1;
    r1.assign(lenVec+1,0);
    for (int i = 0; i < lenVec+1; ++i) {
        r1[i] = CSI[i]*radius; // [m]
    }
    double alpha1;
    vector<double> t2;
    t2.assign(lenVec+1,0);
    vector<double> a2;
    a2.assign(lenVec+1,0);
    vector <double> b2;
    b2.assign(lenVec+1,0);
    double th, theta1, phi1, eff, DtDr, DqDr, cl, cd, CT, CQ, tem1, tem2;
    double a, anew;
    double b, bnew;
    int finished = 0;
    double rad;
    double Vlocal,V0,V2;
    double T = 0.0;
    double Torque = 0.0;

    double alpha_trim = alpha; // [deg] 2.36
    double deltae_trim = delta_e; // [deg]-2.16

    double nBlade = pdb.Pg.np; // [-] 2
    double cl0 = pdb.Pc.Cl0; // [-] 0.3832
    double cd0 = pdb.Pc.Cd0; // [-] 0.0235
    double cl_alpha = pdb.Pc.Clalpha; // [rad^-1] 3.9849
    double cd_alpha = pdb.Pc.Cdalpha; // [rad^-1] 0.154
    double cd_alpha_2 = pdb.Pc.Cdalpha2; // [rad^-2] 1.0476
    double a_0 = pdb.Pc.a0; // [rad]
    double omega;

    double S = db1.Ad.Wing_area; //0.24704
    double cx_alpha = linearInterpolation(db1.alpha, db1.fx.cx_a, db2.fx.cx_a, alpha_trim, h); //0.2115
    double cx_de = linearInterpolation(db1.alpha, db1.cf.cx_de, db2.cf.cx_de, alpha_trim, h); // 0.04029
    double cx_ss = linearInterpolation(db1.alpha, db1.ss.cx, db2.ss.cx, alpha_trim, h); // -0.010595
    double Cx_tot = cx_ss+cx_alpha*alpha_trim/180*M_PI+ cx_de*deltae_trim/180*M_PI;

    double n = rpm/60; // [giri/s]
    omega = n*2*M_PI; // [rad/s]

    for(int j=0; j<lenVec+1; j++){
        rad = r1[j]; // [m] jth-station distance from hub
        theta1 =  pdb.Ps.BA[j]+pitch_propeller; // [deg]
        t2[j] = theta1; // jth-station twist angle
        th = theta1/180.0*M_PI; // [rad] twist angle
        a = 0.1; // [-] axial inflow factor
        b = 0.01; // [-] angular inflow (swirl) factor
        finished = 0;
        int sum = 1;
        while (finished == 0){
            V0 = V*(1+a); // [m/s] flow component ~ equal to aircraft forward speed (Vinf), augmented by axial inflow factor
            V2 = omega*rad*(1-b); // [m/s] flow component ~ equal to angular velocity of the blade section (omega*rad), reduced by angular inflow factor
            phi1 = atan2(V0,V2); // [rad] angle between the two flow components V0 and V2
            alpha1 = th-phi1-a_0; // [rad] angle of attack of the blade jth section
            cl = cl0+cl_alpha*alpha1; // [-] CL lift coefficient
            cd = cd0+cd_alpha*alpha1+cd_alpha_2*alpha1*alpha1; // [-] CD drag coefficient
            Vlocal = sqrt(V0*V0+V2*V2); // [m/s] local flow velocity
            CT = cl*cos(phi1)-cd*sin(phi1); // [-] CT dimensionless thrust coefficient
            DtDr = 0.5*rho*Vlocal*Vlocal*nBlade*chord[j]*CT; // [N/m]
            CQ = cd*cos(phi1)+cl*sin(phi1); // [-] CQ dimensionless torque coefficient
            DqDr = 0.5*rho*Vlocal*Vlocal*nBlade*chord[j]*rad*CQ; // [N]
            tem1= DtDr/(4.0*M_PI*rad*rho*V*V*(1+a)); // [-] corrective factor
            tem2 = DqDr/(4.0*M_PI*rad*rad*rad*rho*V*(1+a)*omega); // [-] corrective factor
            anew = 0.5*(a+tem1);
            bnew = 0.5*(b+tem2);

            if (fabs(anew-a)<1.0/100000){
                if (fabs(bnew-b)<1.0/100000){
                    finished = 1;
                }
            }
            a = anew;
            b = bnew;
            sum = sum+1;
            if (sum > 500){
                finished = 1;
            }
        }
        a2[j] = a;
        b2[j] = b;
        T = T+DtDr*r_step; // [N] thrust
        Torque = Torque+DqDr*r_step; //[N*m] torque
    }

    double t = T/(rho*n*n*diam*diam*diam*diam); // dimensionless thrust
    double q_torque= Torque/(rho*n*n*diam*diam*diam*diam*diam); // dimensionless torque
    double J = V/(n*diam); // advance ratio
    if (t < 0){
        result.eta = 0.0; // blade efficiency
    }else{
        result.eta = t / q_torque*J/(2.0*M_PI);
    }

    result.T = T;
    result.Q = Torque;
    result.ct = t;
    result.cq = q_torque;
    result.J = J;

    return result;

}