#ifndef SPACECRAFT_H // check whether the given token has been #defined earlier in included file
#define SPACECRAFT_H

#include <Eigen/Dense>
#include <tuple>
#include <cmath>

//using namespace Eigen;



class spacecraft
{
private:
    /* data */
public:
    // Constructor (Bare at the moment, keep all variables to default values)
    spacecraft();

    // Default Values (Will some to move to public class later)
    // S/C Momen tof Inertia Matrix [kg-m^2]
    //Eigen::Matrix3d I_SC;
    Eigen::Matrix3d I_SC;
    // Inverse
    Eigen::Matrix3d I_SC_inv;
    // LLO ORbit Params
    double inclination;
    double raan;
    double arg_peri;
    double nu;
    double orbit_period;
    double moon_radius;
    double moon_mu;
    double peri_alt;
    double apo_alt;
    double ra;
    double rp;
    double semi_a;
    double ecc;

    // Initial Conditoins
    Eigen::Vector3d r0;
    Eigen::Vector3d v0;
    Eigen::Quaterniond q0;
    Eigen::Vector3d w0;

    Eigen::VectorXd x_IC;
    Eigen::Vector3d u_IC;

    double dt_sim;



    // Get Inertial Position and Velocity from Orbital Elements
    std::tuple<Eigen::Vector3d, Eigen::Vector3d> get_inertial_pos_vel_from_orbit();
    // Get Time Window
    double get_time_window(double& num_orbits);
    // Get Gravity Gradient
    Eigen::Vector3d gravity_grad_torque(Eigen::Vector3d r_I, Eigen::Quaterniond q_i2b);    
    // x_dot = f(x, u) spacecraft dynamics
    Eigen::VectorXd spacecraft_dynamics(Eigen::VectorXd x, Eigen::Vector3d u);
    // Runge-Kutta 4 Step
    Eigen::VectorXd rk4_step(Eigen::VectorXd x, Eigen::Vector3d u);

};

#endif // SPACECRAFT_H