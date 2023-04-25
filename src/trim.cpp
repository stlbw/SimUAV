//
// Created by Matheus Padilha on 05/04/23.
//

#include "../declaredFun.h"
#include "aerodyn_std_fun.cpp" // contains function to compute density
#include <cmath>
#include <math.h>
#include<iostream>

using namespace std;

struct Trim_Angles{
    double alpha_trim, deltae_trim, u, w;
};

// Prima approssimazione
Trim_Angles trim1(AeroDB db, double V, double h) {
    // Primo tentativo gamma_0=0 --> alpha = theta

    double alpha_min;
    double alpha_max;
    alpha_min = db.alpha.front();
    alpha_max = db.alpha.back();

    double deltae_min;
    double deltae_max;
    deltae_min = db.Dl.Elevator_min;
    deltae_max = db.Dl.Elevator_max;

    //double V = 15;
    //double h = 100; // velocità e altezza che poi saranno definite da utente
    double beta = 0, delta_a = 0, gamma_0 = 0, p = 0, q = 0, r = 0;
    double g = 9.81;

    double rho = computeDensity(h);

    double alpha_int;
    double deltae_int;
    double incr_alpha = 0.02;
    double res = 0.5;
    Trim_Angles angles;
    //double deltae_trim;
    double Cz_tot;
    double Cz_ss;
    double Cz_alpha;
    double Cz_deltae;
    double Cm_ss;
    double Cm_alpha;
    double Cm_deltae;

    for (alpha_int = alpha_min; alpha_int <= alpha_max; alpha_int += incr_alpha) {
        Cm_ss = linearInterpolation(db.alpha, db.ss.cm, alpha_int);
        Cm_alpha = linearInterpolation(db.alpha, db.pm.cm_a, alpha_int);
        Cm_deltae = linearInterpolation(db.alpha, db.cm.cm_de, alpha_int);
        deltae_int = -(Cm_ss + Cm_alpha * alpha_int) / Cm_deltae;
        Cz_ss = linearInterpolation(db.alpha, db.ss.cz, alpha_int);
        Cz_alpha = linearInterpolation(db.alpha, db.fz.cz_a, alpha_int);
        Cz_deltae = linearInterpolation(db.alpha, db.cf.cz_de, alpha_int);
        Cz_tot = Cz_ss + Cz_alpha * alpha_int / 180 * M_PI + Cz_deltae * deltae_int / 180 * M_PI;
        if (abs(db.Ad.Mass * g * cos(alpha_int / 180 * M_PI) + 0.5 * Cz_tot * rho * db.Ad.Wing_area * pow(V, 2)) <
            res){
            angles.alpha_trim = alpha_int;
            Cm_ss = linearInterpolation(db.alpha, db.ss.cm, alpha_int);
            Cm_alpha = linearInterpolation(db.alpha, db.pm.cm_a, alpha_int);
            Cm_deltae = linearInterpolation(db.alpha, db.cm.cm_de, alpha_int);
            angles.deltae_trim = (-(Cm_ss + Cm_alpha * angles.alpha_trim* M_PI /180) / Cm_deltae) *180 / M_PI;
        }
    }

    angles.u = V*cos(angles.alpha_trim);
    angles.w = V*sin(angles.alpha_trim);

    return angles;
}

struct  Trim2{
    double rpm_trim, T_trim, Throttle;
};

Trim2 trim2(AeroDB db, EngineDB endb, PropDB pdb, Trim_Angles angles, double V, double h){
    Trim2 y;
    int lenVec = pdb.Pg.nstation; // [-]
    double diam = pdb.Pg.diameter; // [m]
    double radius = diam/2; // [m]
    // double chord = db.Ad.Chord; // [m]
    vector<double> chord;
    chord.assign(lenVec+1,0);
    for (int i = 0; i < lenVec+1; ++i) {
        chord[i] = pdb.Ps.CH_AD[i]*radius; // [m]
    }
    double pitch_propeller = 0;
    double rpm;
    double delta_rpm = 100; // [giri/min]
    double rpm_min = endb.laps_min; // [giri/min]
    double rpm_max = endb.laps_max; // [giri/min]
    double pitch_tip = pdb.Ps.BA.back(); // [deg]
    double pitch_hub = pdb.Ps.BA.front(); // [deg]
    double xt = radius; // [m]
    vector<double> CSI = pdb.Ps.CSI; // [-]
    double xs = CSI[0]*radius; // [m]
    //double h = 100; // [m]
    double rho = computeDensity(h); // [kg/m^3]
    double coef1 = (pitch_tip-pitch_hub)/(xt-xs); // [deg/m]
    double coef2 = pitch_hub-coef1*xs+pitch_propeller; // [deg]
    double r_step = (xt-xs)/lenVec; //calcolo step
    vector<double> r1;
    r1.assign(lenVec+1,0);
    for (int i = 0; i < lenVec+1; ++i) {
        r1[i] = CSI[i]*radius; // [m]
    }
    double alpha1;
    vector<double> t2 = pdb.Ps.BA; // [deg]
    vector<double> a2;
    a2.assign(lenVec+1,0);
    vector <double> b2;
    b2.assign(lenVec+1,0);
    double th, theta1, phi1, eff, DtDr, DqDr, cl, cd, CT, CQ, tem1, tem2;
    double a, anew;
    double b, bnew;
    int finished = 0;
    double rad;
    //double V = 15; // [m/s]
    double Vlocal,V0,V2;
    double T = 0.0; // inizializzazione vettore spinta
    double Torque = 0.0;// inizializzazione vettore coppia

    double alpha_trim = angles.alpha_trim; // [deg]
    double deltae_trim= angles.deltae_trim; // [deg]
    double g = 9.81; // [m/s^2]
    double gamma_0 = 0;
    double nBlade = pdb.Pg.np; // [-]
    double cl0 = pdb.Pc.Cl0; // [-]
    double cd0 = pdb.Pc.Cd0; // [-]
    double cl_alpha = pdb.Pc.Clalpha; // [rad^-1]
    double cd_alpha = pdb.Pc.Cdalpha; // [rad^-1]
    double cd_alpha_2 = pdb.Pc.Cdalpha2; // [rad^-2]
    double res = 1; // [N]
    double omega;

    for(rpm = rpm_min; rpm <= rpm_max; rpm += delta_rpm){
        double n = rpm/60; // [giri/s]
        omega = n*2*M_PI; // [rad/s]
        for(int j=0; j<lenVec+1; j++){
            rad = r1[j]; // [m] distanza da hub j-esima stazione
            th = t2[j]/180.0*M_PI; // [rad] angolo di svergolamento
            //theta1 = coef1*rad+coef2; //calcolo angolo di svergolamento della j-esima stazione
            //t2[j] = theta1; //angolo di svergolamento della j-esima stazione (-> BA su propeller.txt)
            //th = theta1/180.0*M_PI; //angolo di svergolamento [rad]
            a = 0.1; // [-] inizializzazione axial inflow factor (vedi pag.4 PROPEL.pdf)
            b = 0.01; // [-] inizializzazione angular inflow (swirl) factor (vedi pag.4 PROPEL.pdf)
            finished=0; // inizializzione flag
            int sum = 1; // inizializzione variabile di supporto
            while (finished == 0){
                V0 = V*(1+a); // [m/s] componente del flusso all'incirca uguale alla velocità di avanzamento del velivolo (Vinf), aumentata tramite l'axial inflow factor
                V2 = omega*rad*(1-b); // [m/s] componente del flusso all'incirca uguale alla velocità angolare della sezione della pala (omega*rad), ridotta tramite l'angular inflow factor
                phi1 = atan2(V0,V2); // [rad] angolo tra le due componenti del flusso V0 e V2
                alpha1 = th-phi1; // [rad] angolo di attacco raltivo alla j-esima sezione della pala
                cl = cl0+cl_alpha*alpha1; // [-] CL coefficiente di portanza
                cd = cd0+cd_alpha*alpha1+cd_alpha_2*alpha1*alpha1; // [-] CD coefficiente di resistenza CD = CD0+CD1*CL+CD2*CL^2 (NB nel nostro caso, CD = CD0+CD_alpha*alpha+CD_alpha2*alpha^2 -> slide lezione 2)
                Vlocal = sqrt(V0*V0+V2*V2); // [m/s] velocità locale del flusso
                CT = cl*cos(phi1)-cd*sin(phi1); // [-] CT coefficiente di spinta adimensionale
                DtDr = 0.5*rho*Vlocal*Vlocal*nBlade*chord[j]*CT; // [N/m]=[kg/s^2] contributo di spinta della j-esima sezione
                CQ = cd*cos(phi1)+cl*sin(phi1); // [-] CQ coefficiente di coppia adimensionale
                DqDr = 0.5*rho*Vlocal*Vlocal*nBlade*chord[j]*rad*CQ; // [N]=[kg*m/s^2] contributo di coppia della j-esima sezione
                tem1= DtDr/(4.0*M_PI*rad*rho*V*V*(1+a)); // [-] fattore correttivo del coefficiente "a"
                tem2 = DqDr/(4.0*M_PI*rad*rad*rad*rho*V*(1+a)*omega); // [-] però sono [rad^-1] fattore correttivo del coefficiente "b"
                anew = 0.5*(a+tem1); //nuovo valore coefficiente "a"
                bnew = 0.5*(b+tem2); //nuovo valore coefficiente "b"
                //processo iterativo per arrivare a convergenza
                if (fabs(anew-a)<1.0/100000){
                    if (fabs(bnew-b)<1.0/100000){
                        finished = 1;
                    }
                }
                a = anew; //definizione valore finale coefficiente "a"
                b = bnew; //definizione valore finale coefficiente "b"
                sum = sum+1;
                if (sum > 500){
                    finished = 1;
                }
            }
            a2[j] = a; //definizione valore finale coefficiente "a" per la j-esima stazione
            b2[j] = b; //definizione valore finale coefficiente "b" per la j-esima stazione
            T = T+DtDr*r_step; // [N] sommatoria dei contributi di spinta dalla stazione 1 alla stazione j
            Torque = Torque+DqDr*r_step; //[N*m] sommatoria dei contributi di coppia dalla stazione 1 alla stazione j
        }

        double t = T/(rho*n*n*diam*diam*diam*diam); //coefficiente di spinta adimensionale
        double q_torque= Torque/(rho*n*n*diam*diam*diam*diam*diam); //coefficiente di coppia adimensionale
        double J = V/(n*diam); //rapporto di avanzamento;
        if (t < 0){
            eff = 0.0; //efficienza elica
        }else{
            eff = t / q_torque*J/(2.0*M_PI); //efficienza elica
        }
        double S = db.Ad.Wing_area;
        double cx_alpha = linearInterpolation(db.alpha, db.fx.cx_a, alpha_trim);
        double cx_de = linearInterpolation(db.alpha, db.cf.cx_de, alpha_trim);
        double cx_ss = linearInterpolation(db.alpha, db.ss.cx, alpha_trim);
        double Cx_tot= cx_ss+cx_alpha*alpha_trim/180*M_PI+ cx_de*deltae_trim/180*M_PI;
        y.T_trim = db.Ad.Mass*g*sin(alpha_trim/180*M_PI+gamma_0/180*M_PI)-0.5*Cx_tot*rho*S*V*V;
        if(abs(y.T_trim-T) < res){
            y.rpm_trim = rpm;
        }
    }

    double P_max = 160; // [W]
    double P_d = Torque*omega/1000; // [W]

    if (P_d < P_max){
        y.Throttle = P_d/P_max;
    }

   return y;
}

// Da fare: far inserire le velocità e la quota di volo all'utente
// Slide 50 3 opzioni per il deltae

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

Modes phugoidShortPeriod (AeroDB db, PropDB pdb, Trim_Angles angles, double V, double h) {
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

