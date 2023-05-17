//
// Created by Matheus Padilha on 05/04/23.
//
#include <vector>
#include "../declaredFun.h"
using namespace std;

     /**
 *  Returns the variable value as a result of a linear interpolation of the angle of attack (alpha)
 * @param alpha
 * @param c
 * @param aoa
 * @return
 */
 double linearInterpolationAlpha(vector <double> alpha, vector<double> c, double aoa) {
     double interpCoef; // returning variable
     double alphaUp, alphaLow;
     int pos;

     if (aoa < alpha.front() || aoa > alpha.back()) {
         string error = "The desired angle of attack is outside the database range ALPHA = [" + to_string(alpha.front()) +
                 ", " + to_string(alpha.back()) + "]";
         throw range_error(error);
     }
     vector<double>::iterator it;
     it = find_if(alpha.begin(), alpha.end(), [aoa](double value) {return std::abs(aoa - value) < 0.0001; });
     // find_if returns the element for which the statement is first found valid. If there is no such condition, returns the alpha.end() value

     if (it != alpha.end()) {
         pos = it - alpha.begin(); // returns the position of the element for which it was found the value within the tollerance
         interpCoef = c[pos];
     } else {
         // proceeds with interpolation
         // gets the first angle of attack bigger than the aoa desired and establishes the 2 extreme interpolation points
         for(int i = 0; i < alpha.size(); i++) {
             if (alpha[i] > aoa) {
                 alphaUp= alpha[i]; // upper bound
                 alphaLow = alpha[i  - 1]; // lower bound
                 pos = i;
                 break; // once this condition is met once -> exit the loop
             }
         }
         // the interpolation is done whenever to account for rounding error due to the type variable when saved (e.g. 3.6 != 3.5999999999998)
         if (pos != 0){
             double delta = (c[pos] - c[pos - 1])/(alphaUp - alphaLow); // takes the line coefficient as [y(i) - y(i-1)]/[x(i) - x(i-1)]
             double dAlpha = aoa - alphaLow;
             interpCoef = c[pos - 1] + delta * dAlpha;
         }
     }

     return interpCoef;
 }

 /**
  * Returns the variable value as a result of a linear interpolation of the angle of attack (alpha) for
 * a predetermined altitude.
  * @param alpha
  * @param var1
  * @param var2
  * @param aoa
  * @param h
  * @return
  */
double linearInterpolation(vector <double> alpha, vector<double> var1, vector<double> var2, double aoa, double h) {
     if(var1.size() != alpha.size() || var2.size() != alpha.size()){
         string error = "Cannot interpolate database as the variables do not have the same size.";
         throw range_error(error);
     }

    double interp1 = linearInterpolationAlpha(alpha, var1, aoa); // h1
    double interp2 = linearInterpolationAlpha(alpha, var2, aoa); // h2 > h1
    double delta = 0;
    double dH = 0;
     if (0 <= h && h <= 100) {
         delta = (interp2 - interp1) / 100.0;
         dH = h;
     }
     else if (100 < h && h <= 1000) {
         delta = (interp2 - interp1) / (1000.0 - 100.0);
         dH = (h - 100.0);
    }
     else if (1000 < h && h <= 2000) {
          delta = (interp2 - interp1) / 1000.0;
          dH = (h - 1000.0);
     }
     else {
         string error = "Could not interpolate database as the altitude exceeds the maximum altitude. Declared h = " +
                 to_string(h) + " m.";
         throw range_error(error);
     }

    double interpCoef = interp1 + delta * dH;
     return interpCoef;
 }

 /**
  * Selects the correct databases to interpolate from based on the altitude. This function receives the databses for
  * sea-level, 100m, 1000m and 2000m and the altitude. Returns the range of databases for the desired altitude as db1 and db2
  * @param h
  * @param db1 lower-bound database
  * @param db2 upper-bound database
  * @param dba0 sea-level
  * @param dba100 100m
  * @param dba1000 1000m
  * @param dba2000 2000m
  */
 void getAerodynamicDbWithAltitude(double h, AeroDB& db1, AeroDB& db2, AeroDB dba0, AeroDB dba100, AeroDB dba1000, AeroDB dba2000) {
     if (0 <= h && h <= 100) {
         db1 = dba100; //dba0 is not considered
         db2 = dba100;
     }
     else if (100 < h && h <= 1000) {
         db1 = dba100;
         db2 = dba1000;
     }
     else if (1000 < h && h <= 2000) {
         db1 = dba1000;
         db2 = dba2000;
     }
     else {
         string error = "Could not get database as altitude exceeds the maximum altitude. Declared h = " +
                        to_string(h) + " m.";
         throw range_error(error);
     }
}