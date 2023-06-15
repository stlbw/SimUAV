//
// Created by giaca on 11/06/2023.
//

#include "pidController.h"
class pidController {
public:
    pidController(double kp, double ki, double kd)
        : kp_(kp), ki_(ki), kd_(kd), filterConstant_(0.0), integralErrorSum_(0.0), lastError_(0.0),
        lastErrorDerivative_(0.0), flagFilter_(false), flagErrorCheck_(false) {}

    double compute(double referencePoint, double currentValue, double dt) {
        double error = referencePoint - currentValue;
        // add check on error for PSI -> must set the flag using setErrorCheck method
        if (flagErrorCheck_) {
            if (error < - M_PI) {
                error +=  2*M_PI;
            }
            else if (error > M_PI) {
                error -= 2*M_PI;
            }
        }

        //CAMBIATO DA QUA
        //double derivative = (error - lastError_) / dt;
        integralErrorSum_ += ((error+lastError_)/2) * dt;
        double output ;

        if (flagFilter_ && filterConstant_ != 0.0) {
            double N = 1 / filterConstant_;
            double derivative = lastErrorDerivative_ - kd_ * N * (error - lastError_) - lastErrorDerivative_ * N * dt;
            integralErrorSum_ += ki_ * error * dt;

            output = kp_ * error + integralErrorSum_ + derivative;

            lastErrorDerivative_ = derivative;
        }

        lastError_ = error;
        return output;
        //END CHANGES





    }

    double computePIDSimple(double referencePoint, double currentValue, double dt) {
        double error = referencePoint - currentValue;
        double derivative = (error - lastError_) / dt;

        integralErrorSum_ += error * dt;

        double output = kp_ * error + ki_ * integralErrorSum_ + kd_ * derivative;

        lastErrorDerivative_ = derivative;

        lastError_ = error;
        return output;

    }

    double computeNewTest(double referencePoint, double currentValue, double dt) {
        double error = referencePoint - currentValue;
        if (flagErrorCheck_) {
            if (error < - M_PI) {
                error += 2 * M_PI;
            }
            else if (error > M_PI) {
                error -= 2  * M_PI;
            }
        }
        double N = 1 / filterConstant_;
        if (filterConstant_ == 0) {N = 0;}

        double P = kp_ * error;
        double D = (lastErrorDerivative_ - kd_ * N * lastError_) - N * dt * lastErrorDerivative_ + kd_ * N * error;

        lastError_ = error;
        lastErrorDerivative_ = D;
        integralErrorSum_ += ki_ * lastError_ * dt; // I

        double output;
        output = P + integralErrorSum_ + D;
        return output;
    }

    void setDerivativeFilter(bool flagFilter, double filterConstant) {
        if (flagFilter) {
            flagFilter_ = true;
            filterConstant_ = filterConstant;
        }
    }

    void setErrorCheck(bool errorCheck) {
        if (errorCheck) {
            flagErrorCheck_ = true;
        }
    }
private:
    double kp_;
    double ki_;
    double kd_;
    double filterConstant_;

    double integralErrorSum_;
    double lastError_;
    double lastErrorDerivative_;
    bool flagFilter_;
    bool flagErrorCheck_;

};