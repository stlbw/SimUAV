//
// Created by Matheus Padilha on 05/04/23.
//

#include "../declaredFun.h"
#include <cmath>
#include<iostream>

using namespace std;

// Prima approssimazione
double trim1(AeroDB db) {
    // Primo tentativo gamma_0=0 --> alpha = theta

    double alpha_min, alpha_max;
    alpha_min = db.alpha.front();
    alpha_max = db.alpha.back();

    double deltae_min, deltae_max;
    deltae_min = db.Dl.Elevator_min;
    deltae_max = db.Dl.Elevator_max;

    double V = 15,h = 100; // velocità e altezza che poi saranno definite da utente
    double beta = 0, delta_a = 0, gamma_0 = 0, p = 0, q = 0, r = 0;
    double rho = 1.226, grad = 0.0065, m = 4.2561, T_0 = 288.15, ans, g = 9.81;
    ans = ((T_0 - grad * h) / T_0);
    rho = rho * pow(ans, m);
    double alpha_int = alpha_min, deltae_int = deltae_min, incr_alpha = 0.2, alpha_trim, deltae_trim, incr_deltae = 0.2;
    double Cz_tot, Cz_ss, Cz_alpha, Cz_deltae;

    while (alpha_int <= alpha_max) {
            Cz_ss = linearInterpolation(db.alpha, db.ss.cz, alpha_int);
            Cz_alpha = linearInterpolation(db.alpha, db.fz.cz_a, alpha_int);
            Cz_deltae = linearInterpolation(db.alpha, db.cf.cz_de, alpha_int);
            Cz_tot = Cz_ss + Cz_alpha * alpha_int/180*M_PI + Cz_deltae * deltae_int/180*M_PI;
            if (abs(db.Ad.Mass * g * cos(alpha_int/180*M_PI) + 0.5 * Cz_tot * rho * db.Ad.Wing_area * pow(V, 2)) < 1) {
                alpha_trim = alpha_int;
                alpha_int = alpha_int + incr_alpha;
                deltae_int = deltae_int + incr_deltae;
            } else {
                alpha_int = alpha_int + incr_alpha;
                deltae_int = deltae_int + incr_alpha;
            }
        }

    double Cm_ss, Cm_alpha, Cm_deltae;
    Cm_ss = linearInterpolation(db.alpha, db.ss.cm, alpha_trim);
    Cm_alpha = linearInterpolation(db.alpha, db.pm.cm_a, alpha_trim);
    Cm_deltae = linearInterpolation(db.alpha, db.cm.cm_de, alpha_trim);
    deltae_trim = -(Cm_ss + Cm_alpha * alpha_trim) / Cm_deltae;

    /*double X_trim, Y_trim, Z_trim, L_trim, M_trim, N_trim;
    double cx_alfa, cx_beta, cx_p, cx_q, cx_r;
    double cy_alfa, cy_beta, cy_p, cy_q, cy_r;
    double cx_de, cy_da;
    cx_alfa = linearInterpolation(db.alpha, db.fx.cx_a, 3.58);
    cx_beta = linearInterpolation(db.alpha, db.fx.cx_b, 3.58);
    cx_p = linearInterpolation(db.alpha, db.fx.cx_p, 3.58);
    cx_q = linearInterpolation(db.alpha, db.fx.cx_q, 3.58);
    cx_r = linearInterpolation(db.alpha, db.fx.cx_r, 3.58);
    cx_de = linearInterpolation(db.alpha, db.cf.cx_de, 3.58);
    // donde sta cy_de?
    X_trim = 0.5 * rho * db.Ad.Wing_area * pow(V, 2) *
             (cx_alfa * alpha_int + cx_beta * beta + cx_p * p + cx_q  * q + cx_r * r +
             cx_de * delta_e);
    // donde sta cz_da?
    /*Y_trim = 0.5 * rho * db.Ad.Wing_area * pow(V, 2) *
             (db.fy.cy_a * alpha_int + db.fy.cy_b * beta + db.fy.cy_p * p + db.fy.cy_q * q + db.fy.cy_r * r +
              db.cf.cy_da * delta_e);
    Z_trim = 0.5 * rho * db.Ad.Wing_area * pow(V, 2) *
             (db.fz.cz_a * alpha_int + db.fz.cz_b * beta + db.fz.cz_p * p + db.fz.cz_q * q + db.fz.cz_r * r +
              db.cf.cz_de * delta_e);
    // donde sta cz_da?
    L_trim = 0.5 * rho * db.Ad.Wing_area * pow(V, 2) *
             (db.rm.cl_a * alpha_int + db.rm.cl_b * beta + db.rm.cl_p * p + db.rm.cl_q * q + db.rm.cl_r * r +
              db.cm.cl_da * delta_a);
    // donde sta delta_e
    M_trim = 0.5 * rho * db.Ad.Wing_area * pow(V, 2) *
             (db.pm.cm_a * alpha_int + db.pm.cm_b * beta + db.pm.cm_p * p + db.pm.cm_q * q + db.pm.cm_r * r +
              db.cm.cm_de * delta_e);
    N_trim = 0.5 * rho * db.Ad.Wing_area * pow(V, 2) *
             (db.ym.cn_a * alpha_int + db.ym.cn_b * beta + db.ym.cn_p * p + db.ym.cn_q * q + db.ym.cn_r * r +
              db.c_da * delta_a);
    // donde sta c_de
    }*/

    return alpha_trim;
}

