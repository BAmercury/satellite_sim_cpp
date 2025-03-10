#include "thruster.hpp"

# include <Eigen/Dense>
# include <tuple>
#include <cmath>

// Constructor
thruster::thruster(){
    dt_control = 0.1; // FSW control timestep
    max_torque = Eigen::Array3d(320.0, 220.0, 220.0); // [X, Y, Z] Nm
    // Deadband threshold (Percentage of Max Torque per axis)
    deadband_percentage =  Eigen::Array3d(0.005, 0.005, 0.005);
    deadband = max_torque * deadband_percentage;
}

Eigen::Vector3d thruster::apply_deadband(Eigen::Vector3d torque){
    // Apply deadband to prevent thruster chattering near zero
    return (torque.array().abs() < deadband).select(0, torque);
}

std::tuple<Eigen::Vector3d, Eigen::Vector3d> thruster::compute_pulse_command(Eigen::Vector3d desired_torque){
    // Compute thruster pulse width commands based on desired torque and dt of control
    // pulse_widths: duration thruster should be on (seconds)
    // directions: direction for each axis (-1 or 1)
    Eigen::Vector3d torque = apply_deadband(desired_torque);
    Eigen::Vector3d directions = torque.array().sign();

    // Calculate required pulse width based on torque ratio
    Eigen::Vector3d magnitude_ratio = torque.array().abs() / max_torque;
    Eigen::Vector3d pulse_widths = magnitude_ratio * dt_control;

    // Clip pulse widths to dt control
    pulse_widths = pulse_widths.array().min(dt_control).max(0);

    return std::make_tuple(pulse_widths, directions);
}

Command thruster::torque_to_command(Eigen::Vector3d u){
    // Convert desired torque to thruster commands
    // Returns a Command struct 
    auto [pulse_widths, directions] = compute_pulse_command(u);

    return {pulse_widths, directions, dt_control};

}

Eigen::Vector3d thruster::thruster_model(Command commands){
    // If thruster pulse-width command is equal to or greater than the DT Control timestep
    // Then thruster should stay on continously 
    Eigen::Vector3d pulse_widths = commands.pulse_widths;
    Eigen::Vector3d directions = commands.directions;
    Eigen::Vector3d u_torque = Eigen::Vector3d::Zero();
    for (int i = 0; i < 3; ++i){
        if (pulse_widths(i) >= dt_control){
            u_torque(i) = directions(i) * max_torque(i);
        }
        else{
            // When pulsewidth command is less than the FSW timestep
            // thruster force is full on for a fraction of dt control, off for remainder of that time
            u_torque(i) = directions(i) * max_torque(i) * (pulse_widths(i) / dt_control);
        }
    }

    return u_torque;

}
