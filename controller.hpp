#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <Eigen/Dense>

class controller
{
private:
    /* data */
public:
    // Constructor
    controller();
    double dt_integral = 0.1;

    Eigen::Array3d Kp;
    Eigen::Array3d Ki;
    Eigen::Vector3d integral_error;
    Eigen::Array3d integral_limit;

    Eigen::Vector3d pi_controller(Eigen::VectorXd x, Eigen::Vector3d target, Eigen::Matrix3d I_SC);

};


#endif // CONTROLLER_H