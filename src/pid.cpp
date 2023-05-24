// Funzioni PID \
#include <cmath>
#include "../declaredFun.h"


void longitudinalController(vector <double> currentState){
    double err_v;
    double theta_prev;
    double kp_v = -0.0021;
    double kd_v = -0.0015;
    double ki_v = -0.00087;
    double tau_v = 0.0016;
    double N_v = 1/ tau_v;
    double V_prev = 13.5;            //poi correggere da input
    double V_current = sqrt(pow(currentState[0], 2) + pow(currentState[1], 2) + pow(currentState[2], 2));

    err_v = V_prev - V_current;
    theta_prev = err_v * PID(kp_v, ki_v, kd_v, tau_v, s_v);

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
double PID( double kp, double ki , double kd, double tau, double s, double err_prev, double err_current, double dt, double time) {
    double N = 1 / tau;
    double pid;
    double I, D, T = 0;
    int k = 0;

    double der_err = (err_current - err_prev) / dt * time;    // è così o non è così, questo è il dilemma
if (time==0){
    I = ki;
    D = kd;
} else if {
    I += dt * (ki * err_prev);
    D = kd * N * der_err - D * N;
}
    pid = kp * der_err + (I) + (-kd * N + kd * N * der_err);

    return pid;
}


