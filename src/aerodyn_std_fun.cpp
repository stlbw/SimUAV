#include <cmath>
/**
 * Computes the air density at a certain altitude using a simplified mathematical model
 * @param h
 * @return
 */
double computeDensity(double h) {
    double rho_SL = 1.225;
    double grad = 0.0065;
    double m = 4.2561;
    double T_0 = 288.15;
    double ans = ((T_0 - grad * h) / T_0);
    double rho = rho_SL * pow(ans, m);

    return rho;
}