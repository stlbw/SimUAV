//
// Created by giaca on 11/06/2023.
//

#include "pidController.h"
class pidController {
public:
    pidController(double kp, double ki, double kd, double filterConstant)
        : kp_(kp), ki_(ki), kd_(kd), filterConstant_(filterConstant), integralErrorSum_(0.0), lastError_(0.0), lastErrorDerivative_(0.0) {}

    double compute(double referencePoint, double currentValue, double dt) {
        double error = referencePoint - currentValue;
        double N = 1 / filterConstant_;
        double derivative = lastErrorDerivative_ - kd_ * N * (error - lastError_) - lastErrorDerivative_ * N * dt;
        integralErrorSum_ += ki_ * error * dt;

        double output = kp_ * error + integralErrorSum_ + derivative;

        lastErrorDerivative_ = derivative;
        lastError_ = error;

        return output;


    }
private:
    double kp_;
    double ki_;
    double kd_;
    double filterConstant_;

    double integralErrorSum_;
    double lastError_;
    double lastErrorDerivative_;

};