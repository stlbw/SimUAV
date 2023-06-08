#include <cmath>
#include "../declaredFun.h"

using namespace std;

double* PID( double kp, double ki , double kd, double tau, double err_prev, double err_current, double dt, double time, double flagPID, double N, double I_previous) {
    double P=0;
    double I=0;
    double D=0;
    double counter = 0;

    I = I_previous + ki*(err_current)*dt;
    D = D - D * N * dt + kd * N * (err_current - err_prev);
    P = kp * err_current;

    double* pid = new double[2];
    pid[0] = P + I + D;
    pid[1] = I;
    return pid;
}

double* longitudinalController(double V_ref, double h_ref, double currentState[12], double dt, double flagPID,
                               double time, double err_v_previous, double err_theta_previous, double err_h_previous,
                               double I_v_previous, double I_theta_previous, double I_h_previous){

    double theta_ref;
    double kp_v = -0.0021;
    double kd_v = -0.0015; //D_V
    double ki_v = -0.00087;
    double tau_v = 0.0159;
    double N_v = 1/ tau_v;
    //double V_ref = 13.5;            //poi correggere da input
    double V_current = sqrt(pow(currentState[0], 2) + pow(currentState[1], 2) + pow(currentState[2], 2));
    double err_v [2];

    err_v[0] = err_v_previous ;
    err_v[1] = V_ref - V_current;
    /*if (flagPID == 0)
    {
        err_v[0] = 0 ;
        err_v[1] = V_ref - V_current;
    }
    else if (flagPID == 1){
        err_v[0] = err_v[1] ;
        err_v[1] = V_ref - V_current;
    }*/


    double* theta_pid_pointer = PID(kp_v, ki_v, kd_v, tau_v, err_v[0], err_v[1], dt,time, flagPID, N_v, I_v_previous);
    theta_ref = err_v[1] * theta_pid_pointer[0];
    double I_v = theta_pid_pointer[1];
    delete[] theta_pid_pointer;

    double err_theta[2];
    double de;

    double kp_theta = -0.3;
    double kd_theta = -0.01;
    double ki_theta = -3.25;
    double tau_theta = 0.0159;
    double N_theta = 1 / tau_theta;
    double theta_current = currentState[7];

    err_theta[0] = err_theta_previous ; // at the first step, the previous error is 0 (trim condition, assume it maintains V and h)
    err_theta[1] = theta_ref - theta_current; // the current error is not 0!

    /*if (flagPID == 0)
    {
        err_theta[0] = 0 ; // at the first step, the previous error is 0 (trim condition, assume it maintains V and h)
        err_theta[1] = theta_ref - theta_current; // the current error is not 0!
    }
    else if (flagPID == 1){
        err_theta[0] = err_theta[1] ;
        err_theta[1] = theta_ref - theta_current;
    }*/
    double* de_pid_pointer = PID(kp_theta, ki_theta, kd_theta, tau_theta, err_theta[0], err_theta[1], dt,time, flagPID, N_theta, I_theta_previous);
    de = err_theta[1] * de_pid_pointer[0];
    double I_theta = de_pid_pointer[1];
    delete[] de_pid_pointer;
    // de Saturation

    if(de<= -15*(M_PI/180)){
        de=-15*(M_PI/180);
    }
    else if (de>= 15*(M_PI/180)){
        de=15*(M_PI/180);
    };

    //dth_sat 0.55 -0.45
    double err_h[2];
    double dth;
    double kp_h = 0.019;
    double kd_h = 0.01;
    double ki_h = 0.0002;
    double tau_h = 0.1592;
    double N_h = 1 / tau_h;
    //double h_ref = 100;       // da inserire utente
    double h_current = currentState[9];

    err_h[0] = err_h_previous;
    err_h[1] = h_ref - h_current;
/*
    if (flagPID == 0)
    {
        err_h[0] = 0;
        err_h[1] = h_ref - h_current;
    }
    else if (flagPID == 1){
        err_h[0] = err_h[1] ;
        err_h[1] = h_ref - h_current;
    }*/

    double* dth_pid_pointer =PID(kp_h, ki_h, kd_h, tau_h, err_h[0], err_h[1], dt,time, flagPID, N_h, I_h_previous);
    dth = err_h[1] * dth_pid_pointer[0];
    double I_h = dth_pid_pointer[1];
    delete[] dth_pid_pointer;

    //dth_sat 0.55 -0.45
    if(dth<= -0.45){
        dth=-0.45;
    }
    else if (dth>= 0.55){
        dth=0.55;
    }

    double* long_commands = new double[8];

    long_commands[0] = dth;
    long_commands[1] = de;
    // store current errors to send back to external loop
    long_commands[2] = err_v[1];
    long_commands[3] = err_theta[1];
    long_commands[4] = err_h[1];
    long_commands[5] = I_v;
    long_commands[6] = I_theta;
    long_commands[7] = I_h;

    return long_commands;

}


double* lateralController(double state[12], double dt, double flagPID, double time, double err_psi_previous, double err_phi_previous, double I_psi_previous, double I_phi_previous){
    double err_psi[2];
    double kp_psi = 1.5;
    double kd_psi = 0.01;
    double ki_psi = 0.005;
    double psi_ref = 0;  // dovrebbe essere un vettore preso da Matlab
    double phi_ref;
    double psi_current = state[8];

    err_psi[0] = err_psi_previous ;
    err_psi[1] = psi_ref - psi_current;


    /*if (flagPID == 0)
    {
        err_psi[0] = 0 ;
        err_psi[1] = psi_ref - psi_current;
    }
    else if (flagPID == 1){
        err_psi[0] = err_psi[1] ;
        err_psi[1] = psi_ref - psi_current;
    }
     */
    double I_psi = I_psi_previous + ki_psi*err_psi[1]*dt; // takes into account the previous integral increments
    double D_psi = kd_psi*(err_psi[1]-err_psi[0])/dt;
    double P_psi = kp_psi*err_psi[1];
    phi_ref = err_psi[1] * (P_psi+ D_psi + I_psi);



    double err_phi[2];
    double kp_phi = 0.12;
    double kd_phi = 0.001;
    double ki_phi = 0.0005;
    double phi_current =  state[6];

    err_phi[0] = err_phi_previous ;
    err_phi[1] = phi_ref - phi_current;
    /*if (flagPID == 0)
    {
        err_phi[0] = 0 ;
        err_phi[1] = phi_ref - phi_current;
    }
    else if (flagPID == 1){
        err_phi[0] = err_phi[1] ;
        err_phi[1] = phi_ref - phi_current;
    }*/
    double I_phi = I_phi_previous + ki_phi*err_phi[1]*dt;
    double D_phi = kd_phi*(err_phi[1]-err_phi[0])/dt;
    double P_phi = kp_phi*err_phi[1];
    double da = err_phi[1] * (P_phi+ D_phi + I_phi);

    double* lat_commands = new double[5];
    lat_commands[0] = da;
    lat_commands[1] = err_psi[1];
    lat_commands[2] = err_phi[1];
    lat_commands[3] = I_psi;
    lat_commands[4] = I_phi;


    return lat_commands;
}




