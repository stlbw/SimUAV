//
// Created by Matheus Padilha on 05/04/23.
//

#include "../declaredFun.h"
#include "aerodyn_std_fun.cpp" // contains function to compute density
#include <cmath>
#include<iostream>

using namespace std;

struct Trim_Angles{
    double alpha_trim, deltae_trim;
};

// Prima approssimazione
Trim_Angles trim1(AeroDB db) {
    // Primo tentativo gamma_0=0 --> alpha = theta

    double alpha_min;
    double alpha_max;
    alpha_min = db.alpha.front();
    alpha_max = db.alpha.back();

    double deltae_min;
    double deltae_max;
    deltae_min = db.Dl.Elevator_min;
    deltae_max = db.Dl.Elevator_max;

    double V = 15;
    double h = 100; // velocità e altezza che poi saranno definite da utente
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

    return angles;
}


double trim2(AeroDB db, EngineDB endb, PropDB pdb, Trim_Angles angles){
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
    double rpm_trim, rpm;
    double delta_rpm = 100; // [giri/min]
    double rpm_min = endb.laps_min; // [giri/min]
    double rpm_max = endb.laps_max; // [giri/min]
    double pitch_tip = pdb.Ps.BA.back(); // [deg]
    double pitch_hub = pdb.Ps.BA.front(); // [deg]
    double xt = radius; // [m]
    vector<double> CSI = pdb.Ps.CSI; // [-]
    double xs = CSI[0]*radius; // [m]
    double h = 100; // [m]
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
    double V = 15; // [m/s]
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

    for(rpm = rpm_min; rpm <= rpm_max; rpm += delta_rpm){
        double n = rpm/60; // [giri/s]
        double omega = n*2*M_PI; // [rad/s]
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
        double T_trim = db.Ad.Mass*g*sin(alpha_trim/180*M_PI+gamma_0/180*M_PI)-0.5*Cx_tot*rho*S*V*V;
        if(abs(T_trim-T) < 1){
            rpm_trim = rpm;
        }
    }
   return rpm_trim;
}

// Da fare: far inserire le velocità e la quota di volo all'utente
// Slide 50 3 opzioni per il deltae

// MODO FUGOIDE E CORTO PERIODO

// We have a subsonic vehicle
struct Modes{
    double omega_ph_n, zeta_ph, T_ph, t_dim_ph ;
    double omega_sp_n,zeta_sp,T_sp, t_dim_sp ;
};
Modes md; // initialize the struct of type Modes

Modes phugoidShortPeriod (AeroDB db, PropDB pdb, Trim_Angles angles) {
    double C_Du = 0, C_mu = 0, C_Lu = 0; // these derivatives are 0 because it is a subsonic vehicle
    double g = 9.81;
    double V = 15;
    double h = 100;
    double rho = computeDensity(h);
    //double C_We = 0.2842;
    double C_We = (db.Ad.Mass * g) / (0.5 * rho * pow(V, 2) * db.Ad.Wing_area);
    double C_Le = C_We;
    //double C_Xe = 0.0127;
    double C_Xe = linearInterpolation(db.alpha, db.ss.cx, angles.alpha_trim); // CX coef in steady-state for the trim condition
    double C_De = -C_Xe;
    double k = 0.0382; //every flight condition
    //double Cl1_alpha = 3.25; // discordanza valori tra slide 53 e 57
    double C_Tu = -3 * C_De;
    double Cz_alpha = linearInterpolation(db.alpha, db.fz.cz_a, angles.alpha_trim);
    double Cx_alpha = linearInterpolation(db.alpha, db.fx.cx_a, angles.alpha_trim);
    double Cl_alpha = Cx_alpha * sin(angles.alpha_trim) - Cz_alpha * cos(angles.alpha_trim) ;
    double Cd_alpha = 2 * k * Cl_alpha * C_Le;

    // Approximated Solution

    //**************************************************
    double rho_SL = 1.225;
    double grad = 0.0065;
    double m = 4.2561;
    double T_0 = 288.15;
    double ans;
    //double g = 9.81;
    //double h = 100;
    //double V = 15;
    ans = ((T_0 - grad * h) / T_0);
    //double rho = rho_SL * pow(ans, m);
    double chord= 0.1;
    double Iy = 8*db.Ad.Jy/rho*db.Ad.Wing_area*chord*chord*chord;
    //**************************************************
    // PHUGOID
    double mu = 2 * db.Ad.Mass/rho * db.Ad.Wing_area * chord;
    double omega_ph_ad = (sqrt(2)*g/mu);
    md.omega_ph_n = omega_ph_ad*(2*V/chord);
    md.zeta_ph = -C_Tu/(2*sqrt(2) * C_Le);
    double ph_Re = -md.zeta_ph*md.omega_ph_n;
    double ph_Im = md.omega_ph_n * sqrt(fabs(md.zeta_ph * md.zeta_ph - 1));
    double t_dim_ph = log(0.5)/ph_Re * ph_Im;
    double T_ph= 2*M_PI/ph_Im;

    //**************************************************
    // SHORT PERIOD
    double Cm_alpha = linearInterpolation(db.alpha, db.pm.cm_a, angles.alpha_trim);
    md.omega_sp_n = sqrt(- Cm_alpha/ Iy);
    double omega_sp = md.omega_sp_n * 2 * V/ chord;
    double Cm_q = linearInterpolation(db.alpha, db.pm.cm_q, angles.alpha_trim);
    double Cm1_alpha = linearInterpolation(db.alpha, db.pm.cm_ap, angles.alpha_trim) ;
    md.zeta_sp = (Iy*Cl_alpha-2*mu*(Cm_q+Cm1_alpha))/(2*sqrt(-2*mu*Iy*(2*mu*Cm_alpha+Cm_q*Cl_alpha)));
    double sp_Re = -md.zeta_sp*md.omega_sp_n;
    double sp_Im = md.omega_sp_n * sqrt(fabs(md.zeta_sp * md.zeta_sp - 1));
    double t_dim_sp = log(0.5)/sp_Re * sp_Im;
    double T_sp= 2*M_PI/sp_Im;

    //**************************************************


     return md;
}
