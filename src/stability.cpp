//
// Created by Matheus Padilha on 29/04/23.
//

#include "../declaredFun.h"
#include <cmath>
#include <math.h>
#include<iostream>

// MODO FUGOIDE E CORTO PERIODO

// We have a subsonic vehicle
struct Modes{
    double omega_ph, zeta_ph, T_ph, t_dim_ph ;
    double omega_sp,zeta_sp,T_sp, t_dim_sp ;
};
Modes md; // initialize the struct of type Modes

/**
 * Computes Routh coefficients and Rputh's discriminant and evaluate both Static and Dynamic stability for the
 * longitudinal motion using Routh's criteria
 * @param Iy
 * @param rho
 * @param mu
 * @param aeroCoef
 * @return
 */
void routhCriteria(double Iy, double rho, double mu, double V, double chord, double aeroCoef[10]) {
    double Cw_e = aeroCoef[0];
    double Cl_e = aeroCoef[1];
    double Cd_e = aeroCoef[2];
    double Ct_u = aeroCoef[3];
    double Cl_a = aeroCoef[4];
    double Cd_a = aeroCoef[5];
    double Cm_a = aeroCoef[6];
    double Cm_q = aeroCoef[7];
    double Cl_ap = aeroCoef[8];
    double Cm_ap = aeroCoef[9];

    double A1 = 2 * mu * Iy * (2 * mu + Cl_ap);

    // splitted B1, C1 and D1 into smaller parts to ease the writing
    double b1 = 2 * mu * Iy * (Cl_a + Cd_e - Ct_u);
    double b2 = Iy * Ct_u * Cl_ap;
    double b3 = 2 * mu * Cm_q * Cl_ap;
    double b4 = 4 * pow(mu, 2) * (Cm_q + Cm_ap);
    double B1 =  b1 - b2 - b3 - b4;

    double c1 = 2 * mu * (Cm_q * (Ct_u - Cl_a - Cd_e) - 2 * mu * Cm_a + Cm_ap * Ct_u);
    double c2 = Iy * (2 * Cw_e * (Cw_e - Cd_a) + Ct_u * Cl_a + Cd_e * Cl_a);
    double c3 = Cm_q * Cl_ap * Ct_u;
    double C1 = c1 + c2 + c3;

    double d1 = 2 * pow(Cw_e, 2) * Cm_ap;
    double d2 = 2 * mu * Ct_u * Cm_a;
    double d3 = Ct_u * Cm_q * Cl_a;
    double d4 = 2 * Cw_e * Cm_q * (Cl_e - Cd_a);
    double d5 = 2 * Cd_e * Cm_q * Ct_u;
    double D1 = - d1 + d2 + d3 - d4 + d5;

    double E1 = - 2 * pow(Cw_e, 2) * Cm_a; // Cm_u = Cd_u = 0

    // Static stability
    bool checkStaticStability = (Cm_a < 0 && E1 > 0); // if both conditions are met, the LONGITUDINAL static stability holds

    double R = D1 * (B1 * C1 - A1 * D1) - pow(B1, 2) * E1; // delta, also called Routh's discriminant
    bool checkNegativeRealPart = (A1 > 0 && B1 > 0 && D1 > 0 && E1 > 0); // criteria to establish if the solutions have negative real part
    bool checkComplexConjugate = R > 0;
    bool checkDynamicStability = (checkComplexConjugate && checkNegativeRealPart);

    string staticCheck, dynamicCheck;
    if(checkStaticStability) { staticCheck = "verified";} else { staticCheck = "not verified";}
    if(checkDynamicStability) { dynamicCheck = "verified";} else { dynamicCheck = "not verified";}
    cout << "ROUTH CRITERIA FOR LONGITUDINAL STABILITY:" << endl;
    cout << "Longitudinal Static Stability: "<< staticCheck << endl;
    cout << "Longitudinal Dynamic Stability: "<< dynamicCheck << endl;
    cout << endl;

    double routhCoef[6] = {A1, B1, C1, D1, E1, R}; // todo: understand if it needs to be returned

    // evaluate the short and phugoid modes using the complete
    if (checkDynamicStability) {
        // iff the dynamic stability criteria is respected (complex conjugate solutions), the short period and phugoid exist and can be calculated
        double t_char = chord / (2 * V); // time constant (tempo caratteristico)

        // PHUGOID
        double ph_Re = - ((C1 * D1 - B1 * E1) / pow(C1, 2)) / 2;
        double ph_Im = sqrt(abs(pow(((C1 * D1 - B1 * E1) / pow(C1, 2)), 2) - 4 * 1 * (E1 / C1))) / 2; // the absolute value is needed to avoid having to use complex functions

        double eigenvalueModule_ph = sqrt(pow(ph_Re, 2) + pow(ph_Im, 2));
        double omega_n_ph = eigenvalueModule_ph / t_char; // natural frequency [rad/s]
        double damp_ph = - ph_Re / eigenvalueModule_ph; // damping factor
        double T_ph = (2 * M_PI) / ph_Im * t_char; // period [s]
        double timeHalfAmplitude_ph = log(0.5) / ph_Re * t_char; // time to half amplitude [s]

        // SHORT PERIOD
        double sp_Re = - (B1 / A1) / 2;
        double sp_Im = sqrt(abs(pow((B1 / A1), 2) - 4 * 1 * (C1 / A1))) / 2;

        double eigenvalueModule_sp = sqrt(pow(sp_Re, 2) + pow(sp_Im, 2));
        double omega_n_sp = eigenvalueModule_sp / t_char; // natural frequency [rad/s]
        double damp_sp = - sp_Re / eigenvalueModule_sp; // damping factor
        double T_sp = (2 * M_PI) / sp_Im * t_char; // period [s]
        double timeHalfAmplitude_sp = log(0.5) / sp_Re * t_char; // time to half amplitude [s]

        //todo: decide whether this part has to be printed and where (main?)
        cout << "PHUGOID: " << endl;
        cout << "Frequency [rad/s]: " << omega_n_ph <<endl;
        cout << "Damping ratio: " << damp_ph <<endl;
        cout << "Time to half the amplitude [s]: " << timeHalfAmplitude_ph <<endl;
        cout << "Period [s]: " << T_ph <<endl;
        cout << "" << endl;
        cout << "SHORT PERIOD:" << endl;
        cout << "Frequency [rad/s]: " << omega_n_sp <<endl;
        cout << "Damping ratio: " << damp_sp <<endl;
        cout << "Time to half the amplitude [s]: " << timeHalfAmplitude_sp <<endl;
        cout << "Period [s]: " << T_sp <<endl;
        cout << "" << endl;
        cout << "" << endl;

    }

}

Modes longitudinalStability (AeroDB db, PropDB pdb, Trim_Angles angles, double V, double h) {
    double C_Du = 0, C_mu = 0, C_Lu = 0; // these derivatives are 0 because it is a subsonic vehicle
    double g = 9.81;
    //double V = 15;
    //double h = 100;
    double rho = computeDensity(h);
    //double C_We = 0.2842;
    double C_We = (db.Ad.Mass * g) / (0.5 * rho * pow(V, 2) * db.Ad.Wing_area);
    double C_Le = C_We;
    double C_Xe = linearInterpolation(db.alpha, db.ss.cx, angles.alpha_trim); // CX coef in steady-state for the trim condition
    double C_De = -C_Xe;
    double k = 0.0382; //every flight condition
    double Cl1_alpha = 3.25; // discordanza valori tra slide 53 e 57
    double C_Tu = -0.0382;
    //double C_Tu = -3 * C_De;
    double Cz_alpha = linearInterpolation(db.alpha, db.fz.cz_a, angles.alpha_trim);
    double Cx_alpha = linearInterpolation(db.alpha, db.fx.cx_a, angles.alpha_trim);
    double Cl_alpha = Cx_alpha * sin(angles.alpha_trim * M_PI / 180.0) - Cz_alpha * cos(angles.alpha_trim * M_PI / 180.0) ;
    double Cd_alpha = 2 * k * Cl_alpha * C_Le;
    double Cz1_alpha = linearInterpolation(db.alpha, db.fz.cz_ap, angles.alpha_trim);
    double Cm_alpha = linearInterpolation(db.alpha, db.pm.cm_a, angles.alpha_trim);
    double Cm_q = linearInterpolation(db.alpha, db.pm.cm_q, angles.alpha_trim);
    double Cm1_alpha = linearInterpolation(db.alpha, db.pm.cm_ap, angles.alpha_trim) ;
    double aeroCoefVec[10] = {C_We, C_Le, C_De, C_Tu, Cl_alpha, Cd_alpha, Cm_alpha, Cm_q, Cl1_alpha, Cm1_alpha};
    //**************************************************
    double chord = db.Ad.Chord;
    double Iy = 8 * db.Ad.Jy / (rho * db.Ad.Wing_area * pow(chord, 3)); // dimensionless inertia
    double mu = 2 * db.Ad.Mass / (rho * db.Ad.Wing_area * chord); //dimensionless mass parameter
    //**************************************************

    // Routh Criteria
    routhCriteria(Iy, rho, mu, V, chord, aeroCoefVec);

    // Approximated Solution


    // PHUGOID
    double t_car = chord / (2 * V);
    double omega_ph_ad = C_We / (sqrt(2) * mu); // dimensionless omega - phugoid mode
    double omega_ph_n = omega_ph_ad / t_car; // natural frequency - phugoid mode
    md.zeta_ph = -C_Tu / (2 * sqrt(2) * C_We);
    double ph_Re = -md.zeta_ph * omega_ph_n;
    md.omega_ph = omega_ph_n * sqrt(fabs(md.zeta_ph * md.zeta_ph - 1));

    md.t_dim_ph = abs(log(0.5) / ph_Re);
    md.T_ph = 2 * M_PI / md.omega_ph;

    //**************************************************
    // SHORT PERIOD
    double omega_sp_ad = sqrt(- Cm_alpha/ Iy);
    double omega_sp_n = omega_sp_ad / t_car;
    md.zeta_sp = (Iy*Cl_alpha-2*mu*(Cm_q+Cm1_alpha))/(2*sqrt(-2*mu*Iy*(2*mu*Cm_alpha+Cm_q*Cl_alpha)));
    double sp_Re = -md.zeta_sp * omega_sp_n;
    md.omega_sp = omega_sp_n * sqrt(fabs(md.zeta_sp * md.zeta_sp - 1));
    md.t_dim_sp = log(0.5) / sp_Re;
    md.T_sp = 2 * M_PI / md.omega_sp;

    //**************************************************


    return md;
}