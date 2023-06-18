//
// Created by Matheus Padilha on 05/04/23.
//

#include "../declaredFun.h"
#include <math.h>
#include<iostream>

using namespace std;

struct Trim_Angles{
    double alpha_trim, deltae_trim, theta_trim, u, w;
};


// Prima approssimazione
Trim_Angles trimAngles(AeroDB db1, AeroDB db2, double V, double h, double gamma_0) {
    // Primo tentativo gamma_0=0 --> alpha = theta

    double alpha_min;
    double alpha_max;
    alpha_min = db1.alpha.front();
    alpha_max = db1.alpha.back();

    double deltae_min;
    double deltae_max;
    deltae_min = db1.Dl.Elevator_min;
    deltae_max = db1.Dl.Elevator_max;

    double beta = 0, delta_a = 0, p = 0, q = 0, r = 0;
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

    bool foundAlpha = false;
    alpha_int = alpha_min;

    while (alpha_int <= alpha_max) {
        Cm_ss = linearInterpolation(db1.alpha, db1.ss.cm, db2.ss.cm, alpha_int, h);
        Cm_alpha = linearInterpolation(db1.alpha, db1.pm.cm_a, db2.pm.cm_a, alpha_int, h);
        Cm_deltae = linearInterpolation(db1.alpha, db1.cm.cm_de, db2.cm.cm_de, alpha_int, h);
        deltae_int = -(Cm_ss + Cm_alpha * alpha_int) / Cm_deltae;
        Cz_ss = linearInterpolation(db1.alpha, db1.ss.cz, db2.ss.cz, alpha_int, h);
        Cz_alpha = linearInterpolation(db1.alpha, db1.fz.cz_a, db2.fz.cz_a, alpha_int, h);
        Cz_deltae = linearInterpolation(db1.alpha, db1.cf.cz_de,  db2.cf.cz_de, alpha_int, h);
        Cz_tot = Cz_ss + Cz_alpha * alpha_int / 180 * M_PI + Cz_deltae * deltae_int / 180 * M_PI;

        // MASS and WING AREA are also CONSTANT in the database, no need to interpolate
        if (abs(db1.Ad.Mass * g * cos(alpha_int / 180 * M_PI) + 0.5 * Cz_tot * rho * db1.Ad.Wing_area * V * V) <
            res){
            angles.alpha_trim = alpha_int;
            foundAlpha = true; // change the flag's state once a possible trim condition is found

            Cm_ss = linearInterpolation(db1.alpha, db1.ss.cm, db2.ss.cm, alpha_int, h);
            Cm_alpha = linearInterpolation(db1.alpha, db1.pm.cm_a, db2.pm.cm_a, alpha_int, h);
            Cm_deltae = linearInterpolation(db1.alpha, db1.cm.cm_de, db2.cm.cm_de, alpha_int, h);
            angles.deltae_trim = (-(Cm_ss + Cm_alpha * angles.alpha_trim* M_PI /180) / Cm_deltae) *180 / M_PI;
        }
        alpha_int += incr_alpha;
    }

    // throw error if alpha_trim cannot be found
    if(!foundAlpha) {
        string error = "Could not find alpha between alpha_min = " + to_string(alpha_min) + " [deg] and alpha_max = " +
                                                                                            to_string(alpha_max) +
                                                  " [deg] respecting the trim condition with V = "+ to_string(V) + " [m/s] and h = "+
                                                                                                                   to_string(h) +" [m]";
        throw range_error(error);
    }

    // if alpha is correct, check if delta_trim respects its boundaries, otherwise throw error.
    if(angles.deltae_trim > deltae_max || angles.deltae_trim < deltae_min) {
        string error = "Delta_trim is out of bounds [" + to_string(deltae_min) + ", " + to_string(deltae_max) + "]. Delta_trim = " +
                to_string(angles.deltae_trim) + " [deg].";
        throw range_error(error);
    }

    angles.theta_trim = angles.alpha_trim + gamma_0; // theta = alpha + gamma [DEG]
    angles.u = V*cos(angles.alpha_trim * M_PI / 180.0);
    angles.w = V*sin(angles.alpha_trim * M_PI / 180.0);

    return angles;
}

// Slide 50 3 opzioni per il deltae



