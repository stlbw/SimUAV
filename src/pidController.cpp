//
// Created by giaca on 11/06/2023.
//

#include "pidController.h"
class pidController {
public:
    pidController(double kp, double ki, double kd)
        : kp_(kp), ki_(ki), kd_(kd), filterConstant_(0.0), integralErrorSum_(0.0), lastError_(0.0),
        lastErrorDerivative_(0.0), flagFilter_(false) {}

    double compute(double referencePoint, double currentValue, double dt) {
        double error = referencePoint - currentValue;
        //CAMBIATO DA QUA
        double derivative = (error - lastError_) / dt;
        double output;
        if (flagFilter_ && filterConstant_ != 0.0) {
            double filteredDerivative = (filterConstant_ * lastErrorDerivative_ + dt * derivative) / (filterConstant_ + dt);
            output = kp_ * error + ki_ * integralErrorSum_ + kd_ * filteredDerivative;
            lastErrorDerivative_ = filteredDerivative;
        }
        else {
            output = kp_ * error + ki_ * integralErrorSum_ + kd_ * derivative;
            lastErrorDerivative_ = derivative;
        }

        integralErrorSum_ += error * dt;
        lastError_ = error;
        return output;
        //END CHANGES

        //double N = 1 / filterConstant_;
        //double derivative = lastErrorDerivative_ - kd_ * N * (error - lastError_) - lastErrorDerivative_ * N * dt;
        //integralErrorSum_ += ki_ * error * dt;

        //double output = kp_ * error + integralErrorSum_ + derivative;

        //lastErrorDerivative_ = derivative;



    }

    void setDerivativeFilter(bool flagFilter, double filterConstant) {
        if (flagFilter) {
            flagFilter_ = true;
            filterConstant_ = filterConstant;
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

};