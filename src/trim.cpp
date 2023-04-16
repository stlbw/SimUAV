//
// Created by Matheus Padilha on 05/04/23.
//

#include "../declaredFun.h"
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
    double rho = 1.226;
    double grad = 0.0065;
    double m = 4.2561;
    double T_0 = 288.15;
    double ans;
    double g = 9.81;
    ans = ((T_0 - grad * h) / T_0);
    rho = rho * pow(ans, m);

    double alpha_int;
    double deltae_int;
    double incr_alpha = 0.2;
    Trim_Angles a;
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
            1) {
            a.alpha_trim = alpha_int;
            Cm_ss = linearInterpolation(db.alpha, db.ss.cm, alpha_int);
            Cm_alpha = linearInterpolation(db.alpha, db.pm.cm_a, alpha_int);
            Cm_deltae = linearInterpolation(db.alpha, db.cm.cm_de, alpha_int);
            a.deltae_trim = -(Cm_ss + Cm_alpha * a.alpha_trim) / Cm_deltae;
        }
    }

    return a;
}


double trim2(AeroDB db, EngineDB endb, PropDB pdb, Trim_Angles a){
    double V = 15;
    double alpha_trim = a.alpha_trim;
    double deltae_trim= a.deltae_trim;
    double h = 100; // velocità e altezza che poi saranno definite da utente
    double rho = 1.226;
    double grad = 0.0065;
    double m = 4.2561;
    double T_0 = 288.15;
    double ans;
    double g = 9.81;
    ans = ((T_0 - grad * h) / T_0);
    rho = rho * pow(ans, m);
    double gamma_0 = 0;


    // Per il calcolo degli RPM
    double rpm_min = endb.laps_min;
    double rpm_ref= rpm_min;
    double rpm_max = endb.laps_max;
    double chord= 0.1; //
    double diam = pdb.Pg.diameter;
    double xt = diam/2;
    double xs= pdb.Ps.RD.front(); // 
    double nBlade = pdb.Pg.np;
    double pitch_tip = pdb.Ps.BA.back();
    double pitch_hub = pdb.Ps.BA.front();
    double pitch_propeller = 0;

    double n = rpm_ref/60;
    double omega = n*2*M_PI;
    double coef1 = (pitch_tip-pitch_hub)/(xt-xs);
    double coef2 = pitch_hub-coef1*xs+pitch_propeller;
    int lenVec = pdb.Pg.nstation;
    vector<double> CSI = pdb.Ps.CSI;
    vector<double> t2 = pdb.Ps.BA;
    double cl0 = pdb.Pc.Cl0;
    double cd0 = pdb.Pc.Cd0;
    double cl_alpha = pdb.Pc.Clalpha;
    double cd_alpha = pdb.Pc.Cdalpha;
    double cd_alpha_2 = pdb.Pc.Cdalpha2;
    double n_step = (xt-xs)/lenVec; //calcolo step
    //double theta1;
    double alpha1;
    vector<double> a2, b2;
    double th, phi1, eff, DtDr, DqDr, cl, cd, CT, CQ, tem1, tem2;
    double a, anew;
    double b, bnew;
    int finished = 0;
    double Vlocal,V0,V2;
    double T = 0.0; //inizializzazione vettore spinta
    double Torque = 0.0;//inizializzazione vettore coppia

    for(int j=0; j<lenVec; j++){
        //theta1=coef1*CSI[j]+coef2; //calcolo angolo di svergolamento della j-esima stazione
        //t2[j]=theta1; //angolo di svergolamento della j-esima stazione (-> BA su propeller.txt)
        th = t2[j]/180.0*M_PI; //angolo di svergolamento [rad]
        a = 0.1; //inizializzazione axial inflow factor (vedi pag.4 PROPEL.pdf)
        b = 0.01; //inizializzazione angular inflow (swirl) factor (vedi pag.4 PROPEL.pdf)
        finished=0; //inizializzione flag
        int sum = 1; //inizializzione variabile di supporto
        while (finished == 0){
            V0 = V*(1+a); //componente del flusso all'incirca uguale alla velocità di avanzamento del velivolo (Vinf), aumentata tramite l'axial inflow factor
            V2 = omega*CSI[j]*(1-b); //componente del flusso all'incirca uguale alla velocità angolare della sezione della pala (omega*rad), ridotta tramite l'angular inflow factor
            phi1 = atan2(V0,V2); //angolo tra le due componenti del flusso V0 e V2
            alpha1 = th-phi1; //angolo di attacco raltivo alla j-esima sezione della pala
            cl = cl0+cl_alpha*alpha1; //L coefficiente di portanza
            cd = cd0+cd_alpha*alpha1+cd_alpha_2*pow(alpha1,2); // CD coefficiente di resistenza CD = CD0+CD1*CL+CD2*CL^2 (NB nel nostro caso, CD = CD0+CD_alpha*alpha+CD_alpha2*alpha^2 -> slide lezione 2)
            Vlocal = sqrt(V0*V0+V2*V2); // velocità locale del flusso
            CT = cl*cos(phi1)-cd*sin(phi1); //CT coefficiente di spinta adimensionale
            DtDr = 0.5*rho*Vlocal*Vlocal*nBlade*chord*CT; //contributo di spinta della j-esima sezione
            CQ = cd*cos(phi1)+cl*sin(phi1); //CQ coefficiente di coppia adimensionale
            DqDr = 0.5*rho*Vlocal*Vlocal*nBlade*chord *M_PI*CSI[j]*CQ; //contributo di coppia della j-esima sezione
            tem1=DtDr/(4.0*M_PI*CSI[j]*rho*V*V*(1+a)); //fattore correttivo del coefficiente "a"
            tem2 = DqDr/(4.0*M_PI*CSI[j]*CSI[j]*CSI[j]*rho*V*(1+a)*omega); //fattore correttivo del coefficiente "b"
            anew = 0.5*(a+tem1); //nuovo valore coefficiente "a"
            bnew = 0.5*(b+tem2); //nuovo valore coefficiente "b"
            //processo iterativo per arrivare a convergenza
            if (fabs(anew-a)<1/100000){
                if (fabs(bnew-b)<1/100000){
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
        T = T+DtDr*n_step; //sommatoria dei contributi di spinta dalla stazione 1 alla stazione j
        Torque = Torque+DqDr*n_step; //sommatoria dei contributi di coppia dalla stazione 1 alla stazione j
    }

    double t = T/(rho*n*n*diam*diam*diam*diam); //coefficiente di spinta adimensionale
    double q_torque= Torque/(rho*n*n*diam*diam*diam*diam*diam); //coefficiente di coppia adimensionale
    double J = V/(n*diam); //rapporto di avanzamento;
    if (t < 0){
        eff = 0.0; //efficienza elica
    }else{
        eff = t / q_torque*J/(2.0*M_PI); //efficienza elica
    }

    double rpm_trim;
    double delta_rpm = 100;
    double S = db.Ad.Wing_area;
    double cx_alpha = linearInterpolation(db.alpha, db.fx.cx_a, alpha_trim);
    double cx_de = linearInterpolation(db.alpha, db.cf.cx_de, alpha_trim);
    double cx_ss = linearInterpolation(db.alpha, db.ss.cx, alpha_trim);
    double Cx_tot= cx_ss+cx_alpha*alpha_trim+ cx_de*deltae_trim;
    double T_trim = m*g*sin(alpha_trim+gamma_0)-0.5*Cx_tot*rho*S*V*V;
    for(double rpm = rpm_min; rpm <= rpm_max; rpm += delta_rpm){
        if(abs(T_trim-T) < 1){
            rpm_trim = rpm;
        }
    }

   return rpm_trim;

}

