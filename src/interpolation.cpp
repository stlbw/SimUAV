//
// Created by Matheus Padilha on 05/04/23.
//
#include <vector>
#include "../declaredFun.h"
using namespace std;

     /**
 *  Returns the variable value as a result of a linear interpolation of the angle of attack (alpha) for
 * a predetermined altitude.
 * @param alpha
 * @param c
 * @param aoa
 * @return
 */
 double linearInterpolation(vector <double> alpha, vector<double> c, double aoa) {
     double interpCoef; // returning variable
     double alphaUp, alphaLow;
     int pos;

     if (aoa < alpha.front() || aoa > alpha.back()) {
         cerr << "ERROR: The desired angle of attack is outside the database range ALPHA = [" << alpha.front() << ", " << alpha.back() << "]" << endl;
         ::exit(1); // exits program if out of bound
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


