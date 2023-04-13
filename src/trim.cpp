//
// Created by Matheus Padilha on 05/04/23.
//

#include "../declaredFun.h"
#include <cmath>
#include<iostream>

using namespace std;

// Prima approssimazione
double trim1(AeroDB db) {

    double V = 15; // velocity
    double h = 100; // height
    double alpha_int = 3.58, beta = 0, delta_e = 0, delta_a = 0, gamma_0 = 0;
    double p = 0, q = 0, r = 0;

    double rho_SL = 1.226, grad = 0.0065, m = 4.2561, T_0 = 288.15;
    double ans = ((T_0 - grad * h) / T_0), rho = rho_SL * pow(ans, m);
    double X_trim, Y_trim, Z_trim, L_trim, M_trim, N_trim;
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

// Primo tentativo gamma_0=0 --> alpha = theta
// Calcolo alpha trim
    double alpha_min, alpha_max;
    alpha_min = db.alpha[0];
    alpha_max = db.alpha[0];
    alpha_max = db.alpha[0];
    int i;
    for (i = 0; i < std::size(db.alpha); i++) {
        if (db.alpha[i] < alpha_min)
            alpha_min = db.alpha[i];
        else if (db.alpha[i] > alpha_max)
            alpha_max = db.alpha[i];
    }*/
    return X_trim;
}

