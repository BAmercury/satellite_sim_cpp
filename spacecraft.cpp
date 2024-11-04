#include "spacecraft.hpp"

#include <Eigen/Dense>
#include <tuple>
#include <cmath>


// Constructor
spacecraft::spacecraft() {

    I_SC << 2402.52, 36.33, 78.17,
            36.33, 2072.58, -67.71,
            78.17, -67.71, 2022.10;

    I_SC_inv = I_SC.inverse();

    inclination = 102.6 * (M_PI/180.0); // [rad] Inclination
    raan = -132.0 * (M_PI/180.0); // [rad] Long. of Ascending Node
    arg_peri = 146.6 * (M_PI/180.0); // [rad] Argument of Periapsis
    nu = 0.0; // [rad] True anomaly

    orbit_period = 8.2; // [hr] Orbital Period

    moon_radius = 1737.4; // [km] Moon Radius
    moon_mu = 4902.8001; // [km^3/s^2] Lunar Grav. Constant

    peri_alt = 59.9; // [km] Periapsis Alt
    apo_alt = 6016.1; // [km] Apoapsis Alt
    ra = apo_alt + moon_radius; // [km] Periapsis Radius
    rp = peri_alt + moon_radius; // [km] Apoapsis Radius
    semi_a = (ra + rp) / 2.0; // [km] Semi-Major Axis 
    ecc = (ra - rp) / (ra + rp); // Eccentricity 

    // Initial Conditions
    // State Vector:  [pos, vel, attitude quat, angular velocity]
    auto [vector1, vector2] = get_inertial_pos_vel_from_orbit();
    r0 = vector1; // [km] Poosition Inertial Frame
    v0 = vector2; // [km/s] Velocity inertial Frame
    // Initialize Quaternion
    q0 = Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0); // Eigen uses Scalar part first
    // Initial Angular Rates [rad/s]
    w0 = Eigen::Vector3d(0.01, 0.01, 0.01);

    // Assemble State Vector IC, Quaternion scalar part goes last in the state
    x_IC.resize(13);
    x_IC << r0[0], r0[1], r0[2],
            v0[0], v0[1], v0[2],
            q0.x(), q0.y(), q0.z(), q0.w(),
            w0[0], w0[1], w0[2];
    // Initial Body Torques if needed
    u_IC = Eigen::Vector3d(0.0, 0.0, 0.0);
    // Simulation Time Step
    dt_sim = 0.01; 
}

std::tuple<Eigen::Vector3d, Eigen::Vector3d> spacecraft::get_inertial_pos_vel_from_orbit(){
    // Calculate semi-latus rectum
    double l = semi_a * (1.0 - std::pow(ecc, 2.0));
    // Get Specific Angular Momentum
    double h = std::sqrt(l*moon_mu);
    // Position and Velocity Vector in Perifocal Frame
    Eigen::Vector3d temp = Eigen::Vector3d(std::cos(nu), std::sin(nu), 0.0);
    Eigen::Vector3d r_w= std::pow(h, 2) / moon_mu / (1.0 + ecc * std::cos(nu)) * temp;
    Eigen::Vector3d temp2 = Eigen::Vector3d(-std::sin(nu), ecc + std::cos(nu), 0.0);
    Eigen::Vector3d v_w = moon_mu / h * temp2;
    // DCM Perifocal Frame to MCI frame
    // Create individual rotations
    Eigen::AngleAxisd Z1(raan, Eigen::Vector3d::UnitZ()); // Z rotation
    Eigen::AngleAxisd X(inclination, Eigen::Vector3d::UnitX()); //X rotation
    Eigen::AngleAxisd Z2(arg_peri, Eigen::Vector3d::UnitZ()); // Z rotation

    // Combine rotations and convert to matrix (extrinsic rotation)
    Eigen::Matrix3d C_pf2mci = (Z2 * X * Z1).toRotationMatrix();

    Eigen::Vector3d r_I = C_pf2mci * r_w;
    Eigen::Vector3d v_I = C_pf2mci * v_w;

    return std::make_tuple(r_I, v_I);
}

double spacecraft::get_time_window(double& num_orbits){
    // Function to hlep get a time window based on number of orbits
    double period = 2.0*M_PI * std::sqrt((std::pow(semi_a, 3)) / moon_mu);
    return period*num_orbits; // In Seconds
}

Eigen::Vector3d spacecraft::gravity_grad_torque(Eigen::Vector3d r_I, Eigen::Quaterniond q_i2b){
    //  From Wertz https://link.springer.com/book/10.1007/978-94-009-9907-7
    // Assume Geometric center and center of mass of S/C is the same
    // Get Rotation inertial frame to body frame
    Eigen::Matrix3d C_i2b = q_i2b.toRotationMatrix();
    Eigen::Vector3d r_B = C_i2b * r_I;
    // Get unit vec of r_B
    Eigen::Vector3d r_B_unit = r_B.normalized();
    // Get magnitude as well
    double r_B_mag = r_B.norm();
    // Divide by zero protection
    if (r_B_mag < (1e-10))
    {
        return Eigen::Vector3d::Zero();
    }
    // Gravity Gradient Toque
    double temp = 3.0 * moon_mu / (std::pow(r_B_mag, 3));
    return temp * r_B_unit.cross(I_SC * r_B_unit);
}


Eigen::VectorXd spacecraft::spacecraft_dynamics(Eigen::VectorXd x, Eigen::Vector3d u){
    // Spacecraft Rigid Body Equations of motion
    // Extract states  x = r0, v0, q0.as_quat(), w0
    Eigen::Vector3d r(x[0], x[1], x[2]);
    Eigen::Vector3d v(x[3], x[4], x[5]);
    // Get quat, ensure it's normalized
    // q_i2b
    Eigen::Quaterniond q(x[9], x[6], x[7], x[8]); // Eigen expects the scalar first 
    q.normalize();
    Eigen::Vector3d w(x[10], x[11], x[12]);
    // Create skew symmetric matrix
    Eigen::Matrix3d w_skew;
    w_skew << 0.0, -w[2], w[1],
            w[2], 0.0, -w[0],
           -w[1], w[0], 0.0;

    // Quaternion Kinematics
    Eigen::MatrixXd q_matrix(4,3);
    q_matrix << q.w(), -q.z(), q.y(),
                q.z(), q.w(), -q.x(),
               -q.y(), q.x(), q.w(),
               -q.x(), -q.y(), -q.z();

    // Q_dot (With Scalar part last)
    Eigen::VectorXd q_dot = (0.5) * (q_matrix * w);

    // Gravity Gradient Torque
    Eigen::Vector3d L_gg_b = gravity_grad_torque(r, q); 
    // Sum moments
    Eigen::Vector3d L = u + L_gg_b;
    // Angular Acceleration (Body Frame)
    Eigen::Vector3d w_dot = I_SC_inv * (-w_skew * (I_SC * w) + L);

    // Translatoinal Dynamics
    // Graviity modeel From Wertz https://link.springer.com/book/10.1007/978-94-009-9907-7
    Eigen::Vector3d r_dot = v;
    // Model with Keplerian Motion (2-body)
    // Assume M_moon >> M_sc
    // # https://tfaws.nasa.gov/wp-content/uploads/Rickman-Presentation.pdf
    double r_I_mag = r.norm();
    Eigen::Vector3d v_dot = -moon_mu * r / (std::pow(r_I_mag, 3)); // Translational Acceleration in inertial frame

    // Concatenate
    Eigen::VectorXd x_dot(13);
    x_dot << r_dot[0], r_dot[1], r_dot[2],
             v_dot[0], v_dot[1], v_dot[2],
             q_dot[0], q_dot[1], q_dot[2], q_dot[3],
             w_dot[0], w_dot[1], w_dot[2];
    
    return x_dot;
}

Eigen::VectorXd spacecraft::rk4_step(Eigen::VectorXd x, Eigen::Vector3d u){
    // Runge-Kutta 4 Integration Step
    Eigen::VectorXd f1 = spacecraft_dynamics(x, u);
    Eigen::VectorXd f2 = spacecraft_dynamics(x + 0.5*dt_sim*f1, u);
    Eigen::VectorXd f3 = spacecraft_dynamics(x + 0.5*dt_sim*f2, u);
    Eigen::VectorXd f4 = spacecraft_dynamics(x + dt_sim*f3, u);

    return (x + (dt_sim/6.0) * (f1 + 2.0*f2 + 2.0*f3 + f4));
}

