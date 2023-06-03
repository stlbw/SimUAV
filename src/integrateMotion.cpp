//
// Created by Matheus Padilha on 30/04/23.
//
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
double* getAerodynamicForces(AeroDB db1, AeroDB db2, const double initialConditions[12], const double command[4]) {
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
    double x = initialConditions[10];
    double y = initialConditions[11];
    
    //unpack command vector (throttle is not relevant here) - ANGLES IN RAD
    double delta_a = command[0];
    double delta_e = command[1];
    double delta_r = command[2];

    double V = sqrt(pow(u, 2) + pow(v, 2) + pow(w, 2));
    double alpha = atan2(w, u); // w(i)/u(i) [rad]
    double alpha_deg = alpha * 180.0 / M_PI;
    double beta = asin(v / V); // v(i)/V(i) [rad]
    double rho = computeDensity(h);
    double S = db1.Ad.Wing_area; //[m2] wing area
    double b = db1.Ad.Wing_spann; // [m] wing span
    double c = db1.Ad.Chord; // [m] mean chord


    // make dimensionless the angular velocities p, q and r
    double pHat = p * b / (2 * V);
    double qHat = q * c / (2 * V);
    double rHat = r * b / (2 * V);

    // XForces
    double Cxa = linearInterpolation(db1.alpha, db1.fx.cx_a, db2.fx.cx_a, alpha_deg, h);
    double Cxb = linearInterpolation(db1.alpha, db1.fx.cx_b, db2.fx.cx_b, alpha_deg, h);
    double Cxp = linearInterpolation(db1.alpha, db1.fx.cx_p, db2.fx.cx_p, alpha_deg, h);
    double Cxq = linearInterpolation(db1.alpha, db1.fx.cx_q, db2.fx.cx_q, alpha_deg, h);
    double Cxr = linearInterpolation(db1.alpha, db1.fx.cx_r, db2.fx.cx_r, alpha_deg, h);

    double Cxdelta_a = 0;
    double Cxdelta_e = linearInterpolation(db1.alpha, db1.cf.cx_de, db2.cf.cx_de, alpha_deg, h);
    double Cxdelta_r = 0; // same as Cxdelta_a

    double CxTot = Cxa * alpha + Cxb * beta + Cxp * pHat + Cxq * qHat + Cxr * rHat + Cxdelta_a * delta_a + Cxdelta_e * delta_e + Cxdelta_r * delta_r;
    double XForce = 0.5 * rho * pow(V, 2) * S * CxTot;
    
    
    //YForces
    double Cya = linearInterpolation(db1.alpha, db1.fy.cy_a, db2.fy.cy_a, alpha_deg, h);
    double Cyb = linearInterpolation(db1.alpha, db1.fy.cy_b, db2.fy.cy_b, alpha_deg, h);
    double Cyp = linearInterpolation(db1.alpha, db1.fy.cy_p, db2.fy.cy_p, alpha_deg, h);
    double Cyq = linearInterpolation(db1.alpha, db1.fy.cy_q, db2.fy.cy_q, alpha_deg, h);
    double Cyr = linearInterpolation(db1.alpha, db1.fy.cy_r, db2.fy.cy_r, alpha_deg, h);
    double Cydelta_a = linearInterpolation(db1.alpha, db1.cf.cy_da, db2.cf.cy_da, alpha_deg, h);

    double Cydelta_e = 0; //cannot finc Cy_deltae
    double Cydelta_r = linearInterpolation(db1.alpha, db1.cf.cy_dr, db2.cf.cy_dr, alpha_deg, h);

    double CyTot = Cya * alpha + Cyb * beta + Cyp * pHat + Cyq * qHat + Cyr * rHat + Cydelta_a * delta_a + Cydelta_e * delta_e + Cydelta_r * delta_r;
    double YForce = 0.5 * rho * pow(V, 2) * S * CyTot;
    
    //ZForces
    double Cza = linearInterpolation(db1.alpha, db1.fz.cz_a, db2.fz.cz_a, alpha_deg, h);
    double Czb = linearInterpolation(db1.alpha, db1.fz.cz_b, db2.fz.cz_b, alpha_deg, h);
    double Czp = linearInterpolation(db1.alpha, db1.fz.cz_p, db2.fz.cz_p, alpha_deg, h);
    double Czq = linearInterpolation(db1.alpha, db1.fz.cz_q, db2.fz.cz_q, alpha_deg, h);
    double Czr = linearInterpolation(db1.alpha, db1.fz.cz_r, db2.fz.cz_r, alpha_deg, h);

    double Czdelta_a = 0; //cannot find Cz_deltaa
    double Czdelta_e = linearInterpolation(db1.alpha, db1.cf.cz_de, db2.cf.cz_de, alpha_deg, h);
    double Czdelta_r = 0; //cannot find Cz_deltar

    double CzTot = Cza * alpha + Czb * beta + Czp * pHat + Czq * qHat + Czr * rHat + Czdelta_a * delta_a + Czdelta_e * delta_e + Czdelta_r * delta_r;
    double ZForce = 0.5 * rho * pow(V, 2) * S * CzTot;
    
    //LMoment
    double Cla = linearInterpolation(db1.alpha, db1.rm.cl_a, db2.rm.cl_a, alpha_deg, h);
    double Clb = linearInterpolation(db1.alpha, db1.rm.cl_b, db2.rm.cl_b, alpha_deg, h);
    double Clp = linearInterpolation(db1.alpha, db1.rm.cl_p, db2.rm.cl_p, alpha_deg, h);
    double Clq = linearInterpolation(db1.alpha, db1.rm.cl_q, db2.rm.cl_q, alpha_deg, h);
    double Clr = linearInterpolation(db1.alpha, db1.rm.cl_r, db2.rm.cl_r, alpha_deg, h);
    double Cldelta_a = linearInterpolation(db1.alpha, db1.cm.cl_da, db2.cm.cl_da, alpha_deg, h);
    double Cldelta_e = 0; // cannot find Cl_deltae
    double Cldelta_r = linearInterpolation(db1.alpha, db1.cm.cl_dr, db2.cm.cl_dr, alpha_deg, h);

    double ClTot = Cla * alpha + Clb * beta + Clp * pHat + Clq * qHat + Clr * rHat + Cldelta_a * delta_a + Cldelta_e * delta_e + Cldelta_r * delta_r;
    double LMoment = 0.5 * rho * pow(V, 2) * S * ClTot * b;

    //MMoment
    double Cma = linearInterpolation(db1.alpha, db1.pm.cm_a, db2.pm.cm_a, alpha_deg, h);
    double Cmb = linearInterpolation(db1.alpha, db1.pm.cm_b, db2.pm.cm_b, alpha_deg, h);
    double Cmp = linearInterpolation(db1.alpha, db1.pm.cm_p, db2.pm.cm_p, alpha_deg, h);
    double Cmq = linearInterpolation(db1.alpha, db1.pm.cm_q, db2.pm.cm_q, alpha_deg, h);
    double Cmr = linearInterpolation(db1.alpha, db1.pm.cm_r, db2.pm.cm_r, alpha_deg, h);

    double Cmdelta_a = 0;
    double Cmdelta_e = linearInterpolation(db1.alpha, db1.cm.cm_de, db2.cm.cm_de, alpha_deg, h);
    double Cmdelta_r = 0;

    double CmTot = Cma * alpha + Cmb * beta + Cmp * pHat + Cmq * qHat + Cmr * rHat + Cmdelta_a * delta_a + Cmdelta_e * delta_e + Cmdelta_r * delta_r;
    double MMoment = 0.5 * rho * pow(V, 2) * S * CmTot * c;

    //NMoment
    double Cna = linearInterpolation(db1.alpha, db1.ym.cn_a, db2.ym.cn_a, alpha_deg, h);
    double Cnb = linearInterpolation(db1.alpha, db1.ym.cn_b, db2.ym.cn_b, alpha_deg, h);
    double Cnp = linearInterpolation(db1.alpha, db1.ym.cn_p, db2.ym.cn_p, alpha_deg, h);
    double Cnq = linearInterpolation(db1.alpha, db1.ym.cn_q, db2.ym.cn_q, alpha_deg, h);
    double Cnr = linearInterpolation(db1.alpha, db1.ym.cn_r, db2.ym.cn_r, alpha_deg, h);
    double Cndelta_a = linearInterpolation(db1.alpha, db1.cm.cn_da, db2.cm.cn_da, alpha_deg, h);
    double Cndelta_e = 0; // cannot find Cn_deltae
    double Cndelta_r = linearInterpolation(db1.alpha, db1.cm.cn_dr, db2.cm.cn_dr, alpha_deg, h);

    double CnTot = Cna * alpha + Cnb * beta + Cnp * pHat + Cnq * qHat + Cnr * rHat + Cndelta_a * delta_a + Cndelta_e * delta_e + Cndelta_r * delta_r;
    double NMoment = 0.5 * rho * pow(V, 2) * S * CnTot * b;

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
    double x = state[10];
    double y = state[11];

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
    // acceleration = [u_dot, v_dot, w_dot, p_dot, q_dot, r_dot]
    double p_dot = acceleration[3];
    double r_dot = acceleration[5];

    double du = (r * v - q * w) - g * sin(theta) + X / m + thrust / m;
    double dv = (p * w - r * u) + g * sin(phi) * cos(theta) + Y / m;
    double dw = (q * u - p * v) + g * cos(phi) * cos(theta) + Z / m;
    double dp = -((Jz - Jy) * q * r) / Jx + (p * q + r_dot) * Jxz / Jx + L / Jx;
    double dq = -((Jx - Jz) * p * r) / Jy - (pow(p, 2) - pow(r, 2)) * Jxz / Jy + M / Jy;
    double dr = -((Jy - Jx) * p * q) / Jz - (q * r - p_dot) * Jxz / Jz + N / Jz;
    double dphi = p + q * sin(phi) * tan(theta) + r * cos(phi) * tan(theta);
    double dtheta = q * cos(phi) - r * sin(phi);
    double dpsi = q * sin(phi) / cos(theta) + r * cos(phi) / cos(theta);
    double dh = -u * sin(theta) + v * cos(theta) * sin(phi) + w * cos(theta) * cos(phi);
    double dx = u * cos(theta) * cos(psi) + v * (sin(phi) * sin(theta) * cos(psi) - cos(theta) * sin(psi)) + w * (cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi));
    double dy = u * cos(theta) * sin(psi) - v * (sin(phi) * sin(theta) * sin(psi) - cos(phi) * cos(psi)) + w * (cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi));



    double* remainder = new double[10];
    remainder[0] = du;
    remainder[1] = dv;
    remainder[2] = dw;
    remainder[3] = dp;
    remainder[4] = dq;
    remainder[5] = dr;
    remainder[6] = dphi;
    remainder[7] = dtheta;
    remainder[8] = dpsi;
    remainder[9] = dh;
    remainder[10] = dx;
    remainder[11] = dy;

    return remainder;

}

double* getGravitationalForces (const double state[10], const double mass) {

    double gravity[3] = {0};
    double g = 9.81;
    gravity[2] = g * mass ;
    double det = 0;
    double rotationMatrix[3][3] = {0};
    double inverse_matrix[3][3] = {0};
    double roll = state[6];
    double pitch = state[7];
    double yaw = state[8];
    double* gravityForces = new double[3];

    /*
    double cr = cos(roll);
    double sr = sin(roll);
    double cp = cos(pitch);
    double sp = sin(pitch);
    double cy = cos(yaw);
    double sy = sin(yaw);

    double m00, m01, m02, m10, m11, m12, m20, m21, m22;


    m00 = cp * cy;
    m01 = cp * sy;
    m02 = -sp;

    m10 = -cr * sy + sr * sp * cy;
    m11 = cr * cy + sr * sp * sy;
    m12 = sr * cp;

    m20 = sr * sy + cr * sp * cy;
    m21 = -sr * cy + cr * sp * sy;
    m22 = cr * cp;




    rotationMatrix[0][0] = m00;
    rotationMatrix[0][1] = m01;
    rotationMatrix[0][2] = m02;
    rotationMatrix[1][0] = m10;
    rotationMatrix[1][1] = m11;
    rotationMatrix[1][2] = m12;
    rotationMatrix[2][0] = m20;
    rotationMatrix[2][1] = m21;
    rotationMatrix[2][2] = m22;
     */

    double cosYaw = cos(yaw);
    double sinYaw = sin(yaw);
    double cosPitch = cos(pitch);
    double sinPitch = sin(pitch);
    double cosRoll = cos(roll);
    double sinRoll = sin(roll);

    rotationMatrix[0][0] = cosYaw * cosPitch;
    rotationMatrix[0][1] = -sinYaw * cosRoll + cosYaw * sinPitch * sinRoll;
    rotationMatrix[0][2] = sinYaw * sinRoll + cosYaw * sinPitch * cosRoll;

    rotationMatrix[1][0] = sinYaw * cosPitch;
    rotationMatrix[1][1] = cosYaw * cosRoll + sinYaw * sinPitch * sinRoll;
    rotationMatrix[1][2] = -cosYaw * sinRoll + sinYaw * sinPitch * cosRoll;

    rotationMatrix[2][0] = -sinPitch;
    rotationMatrix[2][1] = cosPitch * sinRoll;
    rotationMatrix[2][2] = cosPitch * cosRoll;

    //finding determinant of the matrix
    for(int i = 0; i < 3; i++)
        det = det + (rotationMatrix[0][i] * (rotationMatrix[1][(i+1)%3] * rotationMatrix[2][(i+2)%3] - rotationMatrix[1][(i+2)%3] * rotationMatrix[2][(i+1)%3]));
    if (det > 0)
    {
        // cout<<"\n Inverse of the matrix is: \n";
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++)
                inverse_matrix[i][j] = ((rotationMatrix[(j+1)%3][(i+1)%3] * rotationMatrix[(j+2)%3][(i+2)%3]) - (rotationMatrix[(j+1)%3][(i+2)%3] *rotationMatrix[(j+2)%3][(i+1)%3])); //finding adjoint and dividing it by determinant
        }
    }
    else std::cout<<"Inverse doesn't exist for this matrix";

    // g_FB = L_VB ^-1 * {0 0 g}:
    for(int m = 0; m < 3; m++){
        for (int n = 0; n < 3; n++){
           gravityForces[m] += (inverse_matrix[m][n] * gravity[n]);
        }
    }

    return gravityForces;

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
        forces[i] = inertia[i] * acceleration[i]; //kg * m/s^2  = kg*m/s^2 = N and kg*m^2 * 1/s^2 = kg * m * m * 1/s^2 = N * m
    }

    return forces;
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


double* getAcceleration(double currentState[6], double previousState[6], double dt){
    // { u v w p q r }
    double* acceleration = new double[6];
    for (int i = 0; i < 6; i++){
        acceleration[i] = (currentState[i] - previousState[i]) / dt; // as a vector field
    }
    return acceleration;

}


double* integrateEquationsOfMotion(AeroDB db1, AeroDB db2, EngineDB endb, PropDB pdb, double rpm, double initialConditions[12], double command[4], double previousState[6], double dt) { //double dt = 0.02
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
    double p = initialConditions[3];
    double q = initialConditions[4];
    double r = initialConditions[5];
    double h = initialConditions[9];

    double alpha = atan2(w, u); //[rad]
    double delta_e = command[1]; //[rad]
    double V = sqrt(pow(u, 2) + pow(v, 2) + pow(w, 2));

    double mass = db1.Ad.Mass;

    //aerodynamic forces
    double *aeroPointer = getAerodynamicForces(db1, db2, initialConditions,
                                               command); // get aerodynamic forces and return it to aeroForces
    for (int i = 0; i < 6; i++) { aeroForces[i] = aeroPointer[i]; }
    delete[] aeroPointer; // delete pointer to avoid memory leak

    // gravity forces
    double *gravityPointer = getGravitationalForces(initialConditions, mass); // get gravity forces and return it to aeroForces
    for (int i = 0; i < 3; i++) { gravForces[i] = gravityPointer[i]; }
    delete[] gravityPointer; // delete pointer to avoid memory leak

    //propeller forces
    propellerData = getPropellerPerformance(db1, db2, endb, pdb, alpha, delta_e, V, h, rpm);
    propForces[0] = -propellerData.T; // all other entries are set to 0
    //L, M, N propeller? -> all zero

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
    for (int i = 0; i < 6; i++) { acceleration[i] = accelerationPointer[i]; } // assign values to variable
    delete[] accelerationPointer; // delete pointer to avoid memory leak

    //compute inertial forces
    double *inertialPointer = getInertialForces(db1, acceleration);
    for (int i = 0; i < 6; i++) { inertialForces[i] = inertialPointer[i]; } // assign values to variable
    delete[] inertialPointer; // delete pointer to avoid memory leak

    double completeForce[6] = {0};
    for (int i = 0; i < size(completeForce); i++) {
        completeForce[i] = aeroForces[i] + propForces[i] + inertialForces[i] + gravForces[i];
    }


    // FORZE INERZIALI
    // MOMENTI PROPELLER

    double remainder[10] = {0};
    //initialize remainders vector for the i-th step

    // INSERITI NUOVI
    double inertiaParameters[5] = {db1.Ad.Mass, db1.Ad.Jx, db1.Ad.Jy, db1.Ad.jz, db1.Ad.jxz};
    //inertiaParameters[1] = db1.Ad.Jx * pow((cos(alpha)),2)+db1.Ad.jz*pow((sin(alpha)),2)-2*db1.Ad.jxz*cos(alpha)*sin(alpha);
    //inertiaParameters[3] =db1.Ad.jz*pow((cos(alpha)),2)+db1.Ad.Jx*pow((sin(alpha)),2)+2*db1.Ad.jxz*cos(alpha)*sin(alpha);
    //inertiaParameters[4] = db1.Ad.jxz*pow((cos(alpha)),2)-pow((sin(alpha)),2)+(db1.Ad.Jx-db1.Ad.jz)*cos(alpha)*sin(alpha);
    //

    double *remainderPointer = getRemainders(initialConditions, command, inertiaParameters, completeForce,
                                             propellerData.T);
    for (int i = 0; i < 6; i++) { remainder[i] = remainderPointer[i]; } // assign values to variable
    delete[] remainderPointer; // delete pointer to avoid memory leak

//prova

    // Initialize current states vector with the initial conditions (safety to avoid mistakes)
    double* currentState = new double[12];
    for (int i = 0; i < 12; i++) {
        currentState[i] = initialConditions[i];
    }

    // Euler's explicit method implementation
    for (int i = 0; i < 11; i++) {
        currentState[i] = initialConditions[i] + dt * remainder[i];
    }

    //todo: check V, alpha, ecc
    double deltae_min;
    double deltae_max;
    deltae_min = db1.Dl.Elevator_min/180*M_PI;
    deltae_max = db1.Dl.Elevator_max/180*M_PI;
    if(command[1]> deltae_max || command[1]< deltae_min) {
        string error = "Delta_trim is out of bounds [" + to_string(deltae_min) + ", " + to_string(deltae_max) + "]. Delta_trim = " +
                       to_string(command[1]) + " [rad].";
        throw range_error(error);
    }
    return currentState;

}

double computeVmin (AeroDB db1, AeroDB db2, double h) {
    double S = db1.Ad.Wing_area;
    double m = db1.Ad.Mass;
    double rho = computeDensity(h);
    double AlphaMax = db1.alpha.back();
    double g = 9.81;
    double CX = linearInterpolation(db1.alpha, db1.ss.cx, db2.ss.cx, AlphaMax, h);
    double CZ = linearInterpolation(db1.alpha, db1.ss.cz, db2.ss.cz, AlphaMax, h);
    double CLmax = -CX * sin(M_PI / 180 * AlphaMax) + CZ * cos(M_PI / 180 * AlphaMax);
    double Vmin = sqrt(2 * m * g / (S * rho * CLmax));

    return Vmin;
}







