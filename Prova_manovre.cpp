// Funzioni PID
#include <cmath>
#include "../declaredFun.h"
#include "integratemotion.cpp"


void longitudinalController() {

        double err_v;
        double theta_prev;
        double kp_v = -0.0021;
        double kd_v = -0.0015;
        double ki_v = -0.00087;
        double tau_v = 0.0016;
        double N_v = 1/ tau_v;
        double V_prev; //= V_trim;
        double V_current;
        int i ;
        for (i= 0; i < N_v; i++) {
          err_v = V_prev - V_current;
          theta_prev = errV * PID(kp_v, ki_v, kd_v, tau_v, s_v);
          V_prev = V_current;
        }
        double err_theta;
        double dth;
        double kp_theta = -0.3;
        double kd_theta = -0.01;
        double ki_theta = -3.5;
        double tau_theta = 0.0016;
        err_theta = theta_prev - theta_current;
        dth = err_theta * PID(kp_theta, ki_theta, kd_theta, tau_theta, s_theta);

        double err_h;
        double de;
        double kp_h = 0.019;
        double kd_h = 0.01;
        double ki_h = 0.0002;
        double tau_h = 0.1592;
        err_h = h_prev - h_current;
        de = PID(kp_h, ki_h, kd_h, tau_h, s_h);

    }
    void lateralController() {
        double err_psi;
        double kp_psi = 1.5;
        double kd_psi = 0.01;
        double ki_psi = 0.005;


        double err_phi;
        double kp_phi = 0.12;
        double kd_phi = 0.001;
        double ki_phi = 0.0005;


    }
}
double PID( double kp, double ki , double kd, double tau, double s, double err_prev, double err_current, double dt){
 double N = 1 / tau;
 double pid;
 pid = kp * err_current + (ki + dt * ( ki * err_prev) ) + (- kd * N + kd * N * err_current) ;
 return pid;
  }

