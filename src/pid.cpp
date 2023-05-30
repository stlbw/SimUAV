#include <cmath>
#include "../declaredFun.h"

using namespace std;

double PID( double kp, double ki , double kd, double tau, double err_prev, double err_current, double dt, double time, double flagPID, double N) {
    double pid;

    double der_err = (err_current - err_prev) / dt * time;    // è così o non è così, questo è il dilemma
    double I[2];
    double D[2];
    if (flagPID == 0)
    {
        I[0] = I[1] = ki;
        D[0] = I[1] = kd;

    }
    else if (flagPID == 1)
    {   
        I[0] = I[1];
        D[0] = D[1];
        I[1] += dt * (ki * err_prev);
        D[1] = kd * N * der_err - D[0] * N;
    }
    pid = kp * der_err + I[1] + D[1];

    return pid;
}

double* longitudinalController(double currentState[12], double dt, double flagPID, double time){

    double theta_ref;
    double kp_v = -0.0021;
    double kd_v = -0.0015;
    double ki_v = -0.00087;
    double tau_v = 0.0016;
    double N_v = 1/ tau_v;
    double V_ref = 13.5;            //poi correggere da input
    double V_current = sqrt(pow(currentState[0], 2) + pow(currentState[1], 2) + pow(currentState[2], 2));
    double err_v [2];

    if (flagPID == 0)
    {
        err_v[0] = 0 ;
        err_v[1] = err_v[0];
    }
    else if (flagPID == 1){
        err_v[0] = err_v[1] ;
        err_v[1] = V_ref - V_current;
    }

    theta_ref = err_v[1] * PID(kp_v, ki_v, kd_v, tau_v, err_v[0], err_v[1], dt,time, flagPID, N_v);

    double err_theta[2];
    double dth;

    double kp_theta = -0.3;
    double kd_theta = -0.01;
    double ki_theta = -3.5;
    double tau_theta = 0.0016;
    double N_theta = 1 / tau_theta;
    double theta_current = currentState[7];

    if (flagPID == 0)
    {
        err_theta[0] = 0 ;
        err_theta[1] = err_theta[0];
    }
    else if (flagPID == 1){
        err_theta[0] = err_theta[1] ;
        err_theta[1] = theta_ref - theta_current;
    }
    dth = err_theta[1] * PID(kp_theta, ki_theta, kd_theta, tau_theta, err_theta[0], err_theta[1], dt,time, flagPID, N_theta);

    double err_h[2];
    double de;
    double kp_h = 0.019;
    double kd_h = 0.01;
    double ki_h = 0.0002;
    double tau_h = 0.1592;
    double N_h = 1 / tau_h;
    double h_ref = 100;       // da inserire utente
    double h_current = currentState[9];

    if (flagPID == 0)
    {
        err_h[0] = 0;
        err_h[1] = err_h[0];
    }
    else if (flagPID == 1){
        err_h[0] = err_h[1] ;
        err_h[1] = h_ref - h_current;
    }
    de = err_h[1] * PID(kp_h, ki_h, kd_h, tau_h, err_h[0], err_h[1], dt,time, flagPID, N_h);

    double long_commands[2];

    long_commands[0] = dth;
    long_commands[1] = de;

    return long_commands;

}

double lateralController(double state[12], double dt, double flagPID, double time){
    double err_psi[2];
    double kp_psi = 1.5;
    double kd_psi = 0.01;
    double ki_psi = 0.005;
    double psi_ref = 0;  // dovrebbe essere un vettore preso da Matlab
    double phi_ref;
    double psi_current = state[8];


    if (flagPID == 0)
    {
        err_psi[0] = 0 ;
        err_psi[1] = err_psi[0];
    }
    else if (flagPID == 1){
        err_psi[0] = err_psi[1] ;
        err_psi[1] = psi_ref - psi_current;
    }
    phi_ref = err_psi[1] * (kp_psi*err_psi[1]+kd_psi*(err_psi[1]-err_psi[0])/dt +ki_psi*err_psi[1]*dt);


    double err_phi[2];
    double kp_phi = 0.12;
    double kd_phi = 0.001;
    double ki_phi = 0.0005;
    double phi_current =  state[6];
    double da;

    if (flagPID == 0)
    {
        err_phi[0] = 0 ;
        err_phi[1] = err_phi[0];
    }
    else if (flagPID == 1){
        err_phi[0] = err_phi[1] ;
        err_phi[1] = phi_ref - phi_current;
    }

    da = err_phi[1] * (kp_phi*err_phi[1]+kd_phi*(err_phi[1]-err_phi[0])/dt +ki_phi*err_phi[1]*dt);



double lat_commands = da;

return lat_commands;
}




