#include <cmath>
#include "../declaredFun.h"

/**
 * Evaluate the aerodynamic forces given a initial condition and the command vector. Uses linear interpolation and the
 * single dependency of alpha (AoA) to compute aerodynamic coefficients.
 * @param db
 * @param initialConditions
 * @param command
 * @param result
 * @returns a vector containing the 6 aerodynamic forces (3) and moments(3)
 */

struct TrimCondition {
    double alphaDeg = 0;
    double V = 0;
    double h = 0;
    double rho = 0;
};

double* getAerodynamicForces(AeroDB db1, AeroDB db2, const double initialConditions[12], const double command[4], aero_condition hb) {
    // unpack initial conditions vector
    double u     = initialConditions[0];
    double v     = initialConditions[1];
    double w     = initialConditions[2];
    double p     = initialConditions[3];
    double q     = initialConditions[4];
    double r     = initialConditions[5];
    double phi   = initialConditions[6];
    double theta = initialConditions[7];
    double psi   = initialConditions[8];
    double h     = initialConditions[9];
    double x     = initialConditions[10];
    double y     = initialConditions[11];
    
    //unpack command vector - ANGLES are IN RAD
    double delta_a = command[0];
    double delta_e = command[1];
    double delta_r = command[2];

    double V = sqrt(u * u + v * v + w * w);

    double alpha     = atan2(w, u);    // w(i)/u(i) [rad]
    double alpha_deg = alpha * 180.0 / M_PI;
    double beta      = asin(v / V);     // v(i)/V(i) [rad]

    double rho  = AtmosphereCalc(1,hb,h);
    double S    = db1.Ad.Wing_area;  //[m2] wing area
    double b    = db1.Ad.Wing_spann; // [m] wing span
    double c    = db1.Ad.Chord;      // [m] mean chord


    // make dimensionless the angular velocities p, q and r
    double pHat = p * b / (2 * V);
    double qHat = q * c / (2 * V);
    double rHat = r * b / (2 * V);

    // XForces

    double Cxss = linearInterpolation(db1.alpha, db1.ss.cx, db2.ss.cx, alpha_deg, h);

    double Cxa = linearInterpolation(db1.alpha, db1.fx.cx_a, db2.fx.cx_a, alpha_deg, h);
    double Cxb = linearInterpolation(db1.alpha, db1.fx.cx_b, db2.fx.cx_b, alpha_deg, h);
    double Cxp = linearInterpolation(db1.alpha, db1.fx.cx_p, db2.fx.cx_p, alpha_deg, h);
    double Cxq = linearInterpolation(db1.alpha, db1.fx.cx_q, db2.fx.cx_q, alpha_deg, h);
    double Cxr = linearInterpolation(db1.alpha, db1.fx.cx_r, db2.fx.cx_r, alpha_deg, h);

    double Cxdelta_a = 0;
    double Cxdelta_e = linearInterpolation(db1.alpha, db1.cf.cx_de, db2.cf.cx_de, alpha_deg, h);
    double Cxdelta_r = 0; // same as Cxdelta_a

    double CxTot  = Cxss + Cxa * alpha + Cxb * beta + Cxp * pHat + Cxq * qHat + Cxr * rHat + Cxdelta_a * delta_a + Cxdelta_e * delta_e + Cxdelta_r * delta_r;
    double XForce = 0.5 * rho * V * V * S * CxTot;
    
    
    //YForces

    double Cyss = linearInterpolation(db1.alpha, db1.ss.cy, db2.ss.cy, alpha_deg, h);

    double Cya  = linearInterpolation(db1.alpha, db1.fy.cy_a, db2.fy.cy_a, alpha_deg, h);
    double Cyb  = linearInterpolation(db1.alpha, db1.fy.cy_b, db2.fy.cy_b, alpha_deg, h);
    double Cyp  = linearInterpolation(db1.alpha, db1.fy.cy_p, db2.fy.cy_p, alpha_deg, h);
    double Cyq  = linearInterpolation(db1.alpha, db1.fy.cy_q, db2.fy.cy_q, alpha_deg, h);
    double Cyr  = linearInterpolation(db1.alpha, db1.fy.cy_r, db2.fy.cy_r, alpha_deg, h);

    double Cydelta_a = linearInterpolation(db1.alpha, db1.cf.cy_da, db2.cf.cy_da, alpha_deg, h);
    double Cydelta_e = 0;
    double Cydelta_r = linearInterpolation(db1.alpha, db1.cf.cy_dr, db2.cf.cy_dr, alpha_deg, h);

    double CyTot  = Cyss + Cya * alpha + Cyb * beta + Cyp * pHat + Cyq * qHat + Cyr * rHat + Cydelta_a * delta_a + Cydelta_e * delta_e + Cydelta_r * delta_r;
    double YForce = 0.5 * rho * V * V * S * CyTot;
    
    //ZForces

    double Czss = linearInterpolation(db1.alpha, db1.ss.cz, db2.ss.cz, alpha_deg, h);

    double Cza  = linearInterpolation(db1.alpha, db1.fz.cz_a, db2.fz.cz_a, alpha_deg, h);
    double Czb  = linearInterpolation(db1.alpha, db1.fz.cz_b, db2.fz.cz_b, alpha_deg, h);
    double Czp  = linearInterpolation(db1.alpha, db1.fz.cz_p, db2.fz.cz_p, alpha_deg, h);
    double Czq  = linearInterpolation(db1.alpha, db1.fz.cz_q, db2.fz.cz_q, alpha_deg, h);
    double Czr  = linearInterpolation(db1.alpha, db1.fz.cz_r, db2.fz.cz_r, alpha_deg, h);

    double Czdelta_a = 0;
    double Czdelta_e = linearInterpolation(db1.alpha, db1.cf.cz_de, db2.cf.cz_de, alpha_deg, h);
    double Czdelta_r = 0;

    double CzTot = Czss + Cza * alpha + Czb * beta + Czp * pHat + Czq * qHat + Czr * rHat + Czdelta_a * delta_a + Czdelta_e * delta_e + Czdelta_r * delta_r;
    double ZForce = 0.5 * rho * V * V * S * CzTot;
    
    //LMoment

    double Clss = linearInterpolation(db1.alpha, db1.ss.cl, db2.ss.cl, alpha_deg, h);

    double Cla = linearInterpolation(db1.alpha, db1.rm.cl_a, db2.rm.cl_a, alpha_deg, h);
    double Clb = linearInterpolation(db1.alpha, db1.rm.cl_b, db2.rm.cl_b, alpha_deg, h);
    double Clp = linearInterpolation(db1.alpha, db1.rm.cl_p, db2.rm.cl_p, alpha_deg, h);
    double Clq = linearInterpolation(db1.alpha, db1.rm.cl_q, db2.rm.cl_q, alpha_deg, h);
    double Clr = linearInterpolation(db1.alpha, db1.rm.cl_r, db2.rm.cl_r, alpha_deg, h);

    double Cldelta_a = linearInterpolation(db1.alpha, db1.cm.cl_da, db2.cm.cl_da, alpha_deg, h);
    double Cldelta_e = 0;
    double Cldelta_r = linearInterpolation(db1.alpha, db1.cm.cl_dr, db2.cm.cl_dr, alpha_deg, h);

    double ClTot = Clss + Cla * alpha + Clb * beta + Clp * pHat + Clq * qHat + Clr * rHat + Cldelta_a * delta_a + Cldelta_e * delta_e + Cldelta_r * delta_r;
    double LMoment = 0.5 * rho * V * V * S * ClTot * b;

    //MMoment

    double Cmss = linearInterpolation(db1.alpha, db1.ss.cm, db2.ss.cm, alpha_deg, h);

    double Cma = linearInterpolation(db1.alpha, db1.pm.cm_a, db2.pm.cm_a, alpha_deg, h);
    double Cmb = linearInterpolation(db1.alpha, db1.pm.cm_b, db2.pm.cm_b, alpha_deg, h);
    double Cmp = linearInterpolation(db1.alpha, db1.pm.cm_p, db2.pm.cm_p, alpha_deg, h);
    double Cmq = linearInterpolation(db1.alpha, db1.pm.cm_q, db2.pm.cm_q, alpha_deg, h);
    double Cmr = linearInterpolation(db1.alpha, db1.pm.cm_r, db2.pm.cm_r, alpha_deg, h);

    double Cmdelta_a = 0;
    double Cmdelta_e = linearInterpolation(db1.alpha, db1.cm.cm_de, db2.cm.cm_de, alpha_deg, h);
    double Cmdelta_r = 0;

    double CmTot = Cmss + Cma * alpha + Cmb * beta + Cmp * pHat + Cmq * qHat + Cmr * rHat + Cmdelta_a * delta_a + Cmdelta_e * delta_e + Cmdelta_r * delta_r;
    double MMoment = 0.5 * rho * V * V * S * CmTot * c;

    //NMoment

    double Cnss = linearInterpolation(db1.alpha, db1.ss.cn, db2.ss.cn, alpha_deg, h);

    double Cna = linearInterpolation(db1.alpha, db1.ym.cn_a, db2.ym.cn_a, alpha_deg, h);
    double Cnb = linearInterpolation(db1.alpha, db1.ym.cn_b, db2.ym.cn_b, alpha_deg, h);
    double Cnp = linearInterpolation(db1.alpha, db1.ym.cn_p, db2.ym.cn_p, alpha_deg, h);
    double Cnq = linearInterpolation(db1.alpha, db1.ym.cn_q, db2.ym.cn_q, alpha_deg, h);
    double Cnr = linearInterpolation(db1.alpha, db1.ym.cn_r, db2.ym.cn_r, alpha_deg, h);

    double Cndelta_a = linearInterpolation(db1.alpha, db1.cm.cn_da, db2.cm.cn_da, alpha_deg, h);
    double Cndelta_e = 0;
    double Cndelta_r = linearInterpolation(db1.alpha, db1.cm.cn_dr, db2.cm.cn_dr, alpha_deg, h);

    double CnTot   = Cnss + Cna * alpha + Cnb * beta + Cnp * pHat + Cnq * qHat + Cnr * rHat + Cndelta_a * delta_a + Cndelta_e * delta_e + Cndelta_r * delta_r;
    double NMoment = 0.5 * rho * V * V * S * CnTot * b;

    double* vecForceMoment = new double[6];

    vecForceMoment[0] = XForce;
    vecForceMoment[1] = YForce;
    vecForceMoment[2] = ZForce;
    vecForceMoment[3] = LMoment;
    vecForceMoment[4] = MMoment;
    vecForceMoment[5] = NMoment;

    return vecForceMoment;

}

double* getRemainders(const double state[10], const double command[4], const double inertiaParameters[5], const double forces[6], const double thrust, const double acceleration[6]) {
    // unpack state vector -> i-th step
    double u     = state[0];
    double v     = state[1];
    double w     = state[2];
    double p     = state[3];
    double q     = state[4];
    double r     = state[5];
    double phi   = state[6];
    double theta = state[7];
    double psi   = state[8];
    double h     = state[9];
    double x     = state[10];
    double y     = state[11];

    double g = 9.81;

    //unpack inertia parameters
    double m   = inertiaParameters[0];
    double Jx  = inertiaParameters[1];
    double Jy  = inertiaParameters[2];
    double Jz  = inertiaParameters[3];
    double Jxz = 0;

    //unpack forces and moments
    double X = forces[0];
    double Y = forces[1];
    double Z = forces[2];
    double L = forces[3];
    double M = forces[4];
    double N = forces[5];

    // acceleration = [u_dot, v_dot, w_dot, p_dot, q_dot, r_dot]
    double p_dot = acceleration[3];
    double r_dot = acceleration[5];

    double du     = (r * v - q * w) - g * sin(theta) + X / m + thrust / m;
    double dv     = (p * w - r * u) + g * sin(phi) * cos(theta) + Y / m;
    double dw     = (q * u - p * v) + g * cos(phi) * cos(theta) + Z / m;
    double dp     = - ((Jz - Jy) * q * r) / Jx + (p * q + r_dot) * Jxz / Jx + L / Jx;
    double dq     = - ((Jx - Jz) * p * r) / Jy - (p * p -  r * r) * Jxz / Jy + M / Jy;
    double dr     = - ((Jy - Jx) * p * q) / Jz - (q * r - p_dot) * Jxz / Jz + N / Jz;
    double dphi   = p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta);
    double dtheta = q * cos(phi) - r * sin(phi);
    double dpsi   = q * sin(phi) / cos(theta) + r * cos(phi) / cos(theta);
    double dh     = - u * sin(theta) + v * cos(theta) * sin(phi) + w * cos(theta) * cos(phi);
    double dx     = u * cos(theta) * cos(psi) + v * (sin(phi) * sin(theta) * cos(psi) - cos(phi) * sin(psi)) + w * (cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi));
    double dy     = u * cos(theta) * sin(psi) + v * (sin(phi) * sin(theta) * sin(psi) + cos(phi) * cos(psi)) + w * (cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi));

    double* remainder = new double[12];

    remainder[0]  = du;
    remainder[1]  = dv;
    remainder[2]  = dw;
    remainder[3]  = dp;
    remainder[4]  = dq;
    remainder[5]  = dr;
    remainder[6]  = dphi;
    remainder[7]  = dtheta;
    remainder[8]  = dpsi;
    remainder[9]  = dh;
    remainder[10] = dx;
    remainder[11] = dy;

    return remainder;

}

double* getInertialForces(AeroDB db, double acceleration[6]) {

    double inertia[6] = {0};

    inertia[0] = db.Ad.Mass;
    inertia[1] = db.Ad.Mass;
    inertia[2] = db.Ad.Mass;

    inertia[3] = db.Ad.Jx;
    inertia[4] = db.Ad.Jy;
    inertia[5] = db.Ad.jz;

    double *forces = new double[6];

    for (int i = 0; i < 6; i++) {
        forces[i] = inertia[i] * acceleration[i];
    }

    return forces;
}


double* getAcceleration(double currentState[6], double previousState[6], double dt){
    // { u v w p q r }
    double* acceleration = new double[6];
    for (int i = 0; i < 6; i++){
        acceleration[i] = (currentState[i] - previousState[i]) / dt; // as a vector field
    }
    return acceleration;

}

/**
 * Integrates the i-th step of the equations of motion using Euler's Explicit method. Requires all parameters to be given
 * using SI units. Angles in initialConditions and command have to be given in RAD.
 * @param db
 * @param endb
 * @param pdb
 * @param rpm
 * @param initialConditions
 * @param command
 * @returns the (i+1)-th integration step
 */

double* integrateEquationsOfMotion(AeroDB db1, AeroDB db2, EngineDB endb, PropDB pdb, double rpm, double initialConditions[12], double command[4], double previousState[6], double dt, std::ofstream& loggerReminder, std::ofstream& loggerAcceleration, aero_condition hb) { //double dt = 0.02
    // compute forces -> initialize vector to 0 + external function to compute it
    // get initial conditions [i-th step]
    // compute remaining -> initialize vector to 0

    Propel propellerData;
    double aeroForces[6] = {0};

    double u = initialConditions[0];
    double v = initialConditions[1];
    double w = initialConditions[2];
    double p = initialConditions[3];
    double q = initialConditions[4];
    double r = initialConditions[5];
    double h = initialConditions[9];

    double alpha       = atan2(w, u); //[rad]
    double alpha_deg   = alpha * 180.0 / M_PI; //[DEG]
    double delta_e     = command[1]; //[rad]
    double delta_e_deg = delta_e * 180 / M_PI; //[DEG]

    double V = sqrt(u * u + v * v + w * w);

    double mass = db1.Ad.Mass;

    //Aerodynamic forces

    double *aeroPointer = getAerodynamicForces(db1, db2, initialConditions, command, hb); // get aerodynamic forces and return it to aeroForces
    for (int i = 0; i < 6; i++) { aeroForces[i] = aeroPointer[i]; }
    delete[] aeroPointer; // delete pointer to avoid memory leak

    //propeller forces
    double rho = AtmosphereCalc(1,hb,h);
    propellerData = getPropellerPerformance(db1, db2, endb, pdb, alpha_deg, delta_e_deg, V, h, rpm, rho);

    // initialize velocity vectors
    double previousVelocity[6] = {0};
    double currentVelocity[6] = {0};

    // assign velocities to vectors
    for (int i = 0; i <6; i++) {
        previousVelocity[i] = previousState[i]; // (i-1)-th step
        currentVelocity[i] = initialConditions[i]; // i-th step
    }

    // initialize acceleration vectors
    double acceleration[6] = {0};
    double *accelerationPointer = getAcceleration(previousVelocity, currentVelocity, dt);
    for (int i = 0; i < 6; i++) {
        acceleration[i] = accelerationPointer[i];
        loggerAcceleration << left << setw(15) << acceleration[i]; //logger
    } // assign values to variable

    loggerAcceleration << " " << endl;
    delete[] accelerationPointer; // delete pointer to avoid memory leak

    double completeForce[6] = {0};
    for (int i = 0; i < 6; i++) {
        completeForce[i] = aeroForces[i];
    }

    double remainder[12] = {0};
    //initialize remainders vector for the i-th step

    // From database
    double inertiaParameters[5] = {0.9726, 1.088831471831585E-002, 1.192188547073001E-002 , 2.232461745527658E-002, 4.590493189389271E-005 };

    double *remainderPointer = getRemainders(initialConditions, command, inertiaParameters, completeForce,
                                             propellerData.T, acceleration);
    for (int i = 0; i < 12; i++){
        remainder[i] = remainderPointer[i];
        loggerReminder << left << setw(15) << remainder[i];
    } // assign values to variable
    loggerReminder << " " << endl;
    delete[] remainderPointer; // delete pointer to avoid memory leak
    // print remainders to logger

    // Initialize current states vector with the initial conditions (safety to avoid mistakes)

    double* currentState = new double[12];
    for (int i = 0; i < 12; i++) {
        currentState[i] = initialConditions[i];
    }

    // Euler's explicit method implementation
    for (int i = 0; i < 12; i++) {
        currentState[i] = initialConditions[i] + dt * remainder[i];
    }
    return currentState;

}
/** Here is th function to comput the minimum Velocity
 * */

double computeVmin (AeroDB db1, AeroDB db2, double h, double rho) {

    double S        = db1.Ad.Wing_area;
    double m        = db1.Ad.Mass;
    double AlphaMax = db1.alpha.back();
    double g        = 9.81;

    double CX       = linearInterpolation(db1.alpha, db1.ss.cx, db2.ss.cx, AlphaMax, h);
    double CZ       = linearInterpolation(db1.alpha, db1.ss.cz, db2.ss.cz, AlphaMax, h);
    double CLmax    = -(CX * sin(M_PI / 180 * AlphaMax) + CZ * cos(M_PI / 180 * AlphaMax));

    double Vmin = sqrt(2 * m * g / (S * rho * CLmax));

    return Vmin;
}







