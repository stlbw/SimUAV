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

    double p_dot = 0; //change
    double r_dot = 0; //change

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
    double dx = u * cos(theta) * cos(psi) + v * (sin(phi) * sin(theta) * cos(psi) - cos(phi) * sin(psi)) + w * (cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi));
    double dy = u * cos(theta) * sin(psi) + v * (sin(phi) * sin(theta) * sin(psi) + cos(phi) * cos(psi)) + w * (cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi));



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
// Gravitational forces
struct G_forces{
    double gravity_forcex, gravity_forcey, gravity_forcez;
};

G_forces Gravity_forces(const double state[10]){

    G_forces gravForces;
    gravForces.gravity_forcex = 0;
    gravForces.gravity_forcey = 0;
    gravForces.gravity_forcez = 0;
    double g = 9.81;
    double det = 0;
    double rotationMatrix[3][3] = {0};
    double inverse_matrix[3][3] = {0};
    double roll = state[6];
    double pitch = state[7];
    double yaw = state[8];
    vector<double> Gravity{0,0,g}; // vector to be multiply per rotational matrix
    vector<double> gravity(3,0); //final vector with the results



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
           gravity[m] += (inverse_matrix[m][n]*Gravity[n]);

        }

    }
    gravForces.gravity_forcex = gravity[0];
    gravForces.gravity_forcey = gravity[1];
    gravForces.gravity_forcez = gravity[2];

    return gravForces;

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
    double *aeroPointer = getAerodynamicForces(db, initialConditions,
                                               command); // get aerodynamic forces and return it to aeroForces
    for (int i = 0; i < 6; i++) { aeroForces[i] = aeroPointer[i]; }
    delete[] aeroPointer; // delete pointer to avoid memory leak

    //propeller forces
    propellerData = getPropellerPerformance(db, endb, pdb, alpha, delta_e, V, h, rpm);
    propForces[0] = -propellerData.T; // all other entries are set to 0
    //todo: how to compute L, M, N propeller?

    //todo: how to compute inertial forces?
}
    // Inertial Forces
    struct Inertia {
        double X_inert, Y_inert, Z_inert, L_inert, M_inert, N_inert;
    };

  /*  Inertia InertialForces(AeroDB db){

    Inertia Inertial_forces;
    double mass = db.Ad.Mass;
    double Ix = db.Ad.Jx;
    double Iy = db.Ad.Jy;
    double Iz = db.Ad.jz;

    Inertial_forces.X_inert = mass;
    Inertial_forces.Y_inert = mass;
    Inertial_forces.Z_inert = mass;

    Inertial_forces.L_inert = Ix;
    Inertial_forces.M_inert = Iy;
    Inertial_forces.N_inert = Iz;

    return Inertial_forces;

    }
    */
