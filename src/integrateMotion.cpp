//
// Created by Matheus Padilha on 30/04/23.
//
#include <cmath>

/**
 * Evaluate the aerodynamic forces given a initial condition and the command vector. Uses linear interpolation and the
 * single dependency of alpha (AoA) to compute aerodynamic coefficients.
 * @param db
 * @param initialConditions
 * @param command
 * @param result
 * @returns a vector containing the 6 aerodynamic forces (3) and moments(3)
 */
double* getAerodynamicForces(AeroDB db, const double initialConditions[10], const double command[4]) {
    // unpack initial conditions vector
    double u = initialConditions[0];
    double v = initialConditions[1];
    double w = initialConditions[2];
    double p = initialConditions[3];
    double q = initialConditions[4];
    double r = initialConditions[5];
    //double phi = initialConditions[6];
    //double theta = initialConditions[7];
    //double psi = initialConditions[8];
    double h = initialConditions[9];
    
    //unpack command vector (throttle is not relevant here) - ANGLES IN RAD
    double delta_a = command[0];
    double delta_e = command[1];
    double delta_r = command[2];

    double V = sqrt(pow(u, 2) + pow(v, 2) + pow(w, 2));
    double alpha = atan2(w, u); // w(i)/u(i) [rad]
    double alpha_deg = alpha * 180.0 / M_PI;
    double beta = asin(v / V); // v(i)/V(i) [rad]
    double rho = computeDensity(h);
    double S = db.Ad.Wing_area; //[m2] wing area
    
    // XForces
    double Cxa = linearInterpolation(db.alpha, db.fx.cx_a, alpha_deg);
    double Cxb = linearInterpolation(db.alpha, db.fx.cx_b, alpha_deg);
    double Cxp = linearInterpolation(db.alpha, db.fx.cx_p, alpha_deg);
    double Cxq = linearInterpolation(db.alpha, db.fx.cx_q, alpha_deg);
    double Cxr = linearInterpolation(db.alpha, db.fx.cx_r, alpha_deg);
    double Cxdelta_a = 0; //todo: cannot find CX_da on database
    double Cxdelta_e = linearInterpolation(db.alpha, db.cf.cx_de, alpha_deg );
    double Cxdelta_r = 0; // same as Cxdelta_a

    double CxTot = Cxa * alpha + Cxb * beta + Cxp * p + Cxq * q + Cxr * r + Cxdelta_a * delta_a + Cxdelta_e * delta_e + Cxdelta_r * delta_r;
    double XForce = 0.5 * rho * pow(V, 2) * S * CxTot;
    
    
    //YForces
    double Cya = linearInterpolation(db.alpha, db.fy.cy_a, alpha_deg);
    double Cyb = linearInterpolation(db.alpha, db.fy.cy_b, alpha_deg);
    double Cyp = linearInterpolation(db.alpha, db.fy.cy_p, alpha_deg);
    double Cyq = linearInterpolation(db.alpha, db.fy.cy_q, alpha_deg);
    double Cyr = linearInterpolation(db.alpha, db.fy.cy_r, alpha_deg);
    double Cydelta_a = linearInterpolation(db.alpha, db.cf.cy_da, alpha_deg);
    double Cydelta_e = 0; //cannot finc Cy_deltae
    double Cydelta_r = linearInterpolation(db.alpha, db.cf.cy_dr, alpha_deg);

    double CyTot = Cya * alpha + Cyb * beta + Cyp * p + Cyq * q + Cyr * r + Cydelta_a * delta_a + Cydelta_e * delta_e + Cydelta_r * delta_r;
    double YForce = 0.5 * rho * pow(V, 2) * S * CyTot;
    
    //ZForces
    double Cza = linearInterpolation(db.alpha, db.fz.cz_a, alpha_deg);
    double Czb = linearInterpolation(db.alpha, db.fz.cz_b, alpha_deg);
    double Czp = linearInterpolation(db.alpha, db.fz.cz_p, alpha_deg);
    double Czq = linearInterpolation(db.alpha, db.fz.cz_q, alpha_deg);
    double Czr = linearInterpolation(db.alpha, db.fz.cz_r, alpha_deg);
    double Czdelta_a = 0; //cannot find Cz_deltaa
    double Czdelta_e = linearInterpolation(db.alpha, db.cf.cz_de, alpha_deg);
    double Czdelta_r = 0; //cannot find Cz_deltar

    double CzTot = Cza * alpha + Czb * beta + Czp * p + Czq * q + Czr * r + Czdelta_a * delta_a + Czdelta_e * delta_e + Czdelta_r * delta_r;
    double ZForce = 0.5 * rho * pow(V, 2) * S * CzTot;
    
    //LMoment
    double Cla = linearInterpolation(db.alpha, db.rm.cl_a, alpha_deg);
    double Clb = linearInterpolation(db.alpha, db.rm.cl_b, alpha_deg);
    double Clp = linearInterpolation(db.alpha, db.rm.cl_p, alpha_deg);
    double Clq = linearInterpolation(db.alpha, db.rm.cl_q, alpha_deg);
    double Clr = linearInterpolation(db.alpha, db.rm.cl_r, alpha_deg);
    double Cldelta_a = linearInterpolation(db.alpha, db.cm.cl_da, alpha_deg);
    double Cldelta_e = 0; // cannot find Cl_deltae
    double Cldelta_r = linearInterpolation(db.alpha, db.cm.cl_dr, alpha_deg);

    double ClTot = Cla * alpha + Clb * beta + Clp * p + Clq * q + Clr * r + Cldelta_a * delta_a + Cldelta_e * delta_e + Cldelta_r * delta_r;
    double LMoment = 0.5 * rho * pow(V, 2) * S * ClTot;

    //MMoment
    double Cma = linearInterpolation(db.alpha, db.pm.cm_a, alpha_deg);
    double Cmb = linearInterpolation(db.alpha, db.pm.cm_b, alpha_deg);
    double Cmp = linearInterpolation(db.alpha, db.pm.cm_p, alpha_deg);
    double Cmq = linearInterpolation(db.alpha, db.pm.cm_q, alpha_deg);
    double Cmr = linearInterpolation(db.alpha, db.pm.cm_r, alpha_deg);
    double Cmdelta_a = 0;
    double Cmdelta_e = linearInterpolation(db.alpha, db.cm.cm_de, alpha_deg);
    double Cmdelta_r = 0;

    double CmTot = Cma * alpha + Cmb * beta + Cmp * p + Cmq * q + Cmr * r + Cmdelta_a * delta_a + Cmdelta_e * delta_e + Cmdelta_r * delta_r;
    double MMoment = 0.5 * rho * pow(V, 2) * S * CmTot;

    //NMoment
    double Cna = linearInterpolation(db.alpha, db.ym.cn_a, alpha_deg);
    double Cnb = linearInterpolation(db.alpha, db.ym.cn_b, alpha_deg);
    double Cnp = linearInterpolation(db.alpha, db.ym.cn_p, alpha_deg);
    double Cnq = linearInterpolation(db.alpha, db.ym.cn_q, alpha_deg);
    double Cnr = linearInterpolation(db.alpha, db.ym.cn_r, alpha_deg);
    double Cndelta_a = linearInterpolation(db.alpha, db.cm.cn_da, alpha_deg);
    double Cndelta_e = 0; // cannot find Cn_deltae
    double Cndelta_r = linearInterpolation(db.alpha, db.cm.cn_dr, alpha_deg);

    double CnTot = Cna * alpha + Cnb * beta + Cnp * p + Cnq * q + Cnr * r + Cndelta_a * delta_a + Cndelta_e * delta_e + Cndelta_r * delta_r;
    double NMoment = 0.5 * rho * pow(V, 2) * S * CnTot;

    double* vecForceMoment = new double[6];
    vecForceMoment[0] = XForce;
    vecForceMoment[1] = YForce;
    vecForceMoment[2] = ZForce;
    vecForceMoment[3] = LMoment;
    vecForceMoment[4] = MMoment;
    vecForceMoment[5] = NMoment;

    return vecForceMoment;






}

double* getRemainders(const double state[10], const double command[4], const double inertiaParameters[5], const double forces[6], const double thrust) {
    // unpack state vector -> i-th step
    double u = state[0];
    double v = state[1];
    double w = state[2];
    double p = state[3];
    double q = state[4];
    double r = state[5];
    double phi = state[6];
    double theta = state[7];
    double psi = state[8];
    double h = state[9];

    double g = 9.81;

    //unpack inertia parameters
    double m = inertiaParameters[0];
    double Jx = inertiaParameters[1];
    double Jy = inertiaParameters[2];
    double Jz = inertiaParameters[3];
    double Jxz = inertiaParameters[4];

    //unpack forces and moments
    double X = forces[0];
    double Y = forces[1];
    double Z = forces[2];
    double L = forces[3];
    double M = forces[4];
    double N = forces[5];

    //double du = (r * v - q * w) - g * sin(theta)

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
void integrateEquationsOfMotion(AeroDB db, EngineDB endb, PropDB pdb, double rpm, double initialConditions[10], double command[4]) {
    // compute forces -> initialize vector to 0 + external function to compute it
    // get initial conditions [i-th step]
    // compute remaining -> initialize vector to 0
    Propel propellerData;
    double aeroForces[6] = {0};
    double propForces[6] = {0};
    double inertialForces[6] = {0};
    double gravForces[6] = {0};

    double u = initialConditions[0];
    double v = initialConditions[1];
    double w = initialConditions[2];
    double h = initialConditions[9];
    double alpha = atan2(w, u); //[rad]
    double delta_e = command[1]; //[rad]
    double V = sqrt(pow(u, 2) + pow(v, 2) + pow(w, 2));

    //aerodynamic forces
    double* aeroPointer = getAerodynamicForces(db, initialConditions, command); // get aerodynamic forces and return it to aeroForces
    for(int i = 0; i < 6; i++) {aeroForces[i] = aeroPointer[i];}
    delete[] aeroPointer; // delete pointer to avoid memory leak

    //propeller forces
    propellerData = getPropellerPerformance(db, endb, pdb, alpha, delta_e, V, h, rpm);
    propForces[0] = - propellerData.T; // all other entries are set to 0
    //todo: how to compute L, M, N propeller?

    //todo: how to compute inertial forces?
    //todo: add gravitational forces -> add rotation matrix L_vb

    double completeForce[6] = {0};
    for (int i = 0; i < size(completeForce); i++) {
        completeForce[i] = aeroForces[i] + propForces[i] + inertialForces[i] + gravForces[i];
    }

    double remainder[10] = {0}; //initialize remainders vector for the i-th step
    double inertiaParameters[5] = {db.Ad.Mass, db.Ad.Jx, db.Ad.Jy, db.Ad.jz, db.Ad.jxz};

    double* remainderPointer = getRemainders(initialConditions, command, inertiaParameters, completeForce, propellerData.T);
    for(int i = 0; i < 6; i++) {remainder[i] = remainderPointer[i];} // assign values to variable
    delete[] remainderPointer; // delete pointer to avoid memory leak

    // integrate equations of motion





}