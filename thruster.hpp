#ifndef THRUSTER_H
#define THRUSTER_H

#include <Eigen/Dense>
#include <tuple>

struct Command {
    Eigen::Vector3d pulse_widths;
    Eigen::Vector3d directions;
    double timestep;
};

class thruster
{
private:
    /* data */
public:

    // using default values (later may add ability to override default values via constructor)
    thruster();

    double dt_control;
    Eigen::Array3d max_torque;
    Eigen::Array3d deadband;
    Eigen::Array3d deadband_percentage;

    Eigen::Vector3d apply_deadband(Eigen::Vector3d torque);
    std::tuple<Eigen::Vector3d, Eigen::Vector3d> compute_pulse_command(Eigen::Vector3d desired_torque);
    Command torque_to_command(Eigen::Vector3d u);
    Eigen::Vector3d thruster_model(Command commands);
    

};


#endif // THRUSTER_H