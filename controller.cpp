#include "controller.hpp"

#include <Eigen/Dense>

// Constructor
controller::controller(){
    
    Kp = Eigen::Array3d(1500.0, 1800.0, 1800.0); // Discrete P gains
    Ki = Eigen::Array3d(800.0, 800.0, 800.0); // Discrete I gains
    integral_error = Eigen::Vector3d::Zero();
    integral_limit = Eigen::Array3d(100.0, 50.0, 50.0);

}

 Eigen::Vector3d controller::pi_controller(Eigen::VectorXd x, Eigen::Vector3d target, Eigen::Matrix3d I_SC){
    // Get Rates from state vector
    Eigen::Vector3d rate_rad = Eigen::Vector3d(x[11], x[12], x[13]);
    Eigen::Vector3d rate_error = target - rate_rad;

    // Update integral error
    integral_error = integral_error + rate_error*dt_integral;
    // Limit Inetgral error
    integral_error = integral_error.array().min(integral_limit).max(-integral_limit);
    // Cross coupling term
    Eigen::Vector3d cross_coupling = rate_rad.cross(I_SC * rate_rad);

    // Return u as Vector3d Torque 
    return ((Kp * rate_error.array()) + (Ki * integral_error.array()) + cross_coupling.array()).matrix();

 }