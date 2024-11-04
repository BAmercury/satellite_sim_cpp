#include "spacecraft.hpp"
#include "thruster.hpp"
#include "controller.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

 
int main()
{
  // Initialize Class Objects
  spacecraft sc_obj;
  thruster thruster_obj;
  controller controller_obj;

  // final sim time
  double num_orbits = 1.0;
  double tFinal = sc_obj.get_time_window(num_orbits);
  //double tFinal = 100.0;
  double t = 0.0;

  // Time to run next control
  double t_next_control = 0.0;

  // Store Results
  // Output file
  std::ofstream output_file("spacecraft_data.csv");
  output_file << std::setprecision(15);

  // Write Header
  output_file << "time,";
  output_file << "pos_x,pos_y,pos_z,";
  output_file << "vel_x,vel_y,vel_z,";
  output_file << "quat_x,quat_y,quat_z,quat_w,";
  output_file << "omega_x,omega_y,omega_z,";
  output_file << "torque_x,torque_y,torque_z\n";

  std::vector<double> out_times;
  std::vector<Eigen::VectorXd> out_states;
  std::vector<Eigen::Vector3d> out_controls;

  // Initialize Running Variables
  Eigen::VectorXd state = sc_obj.x_IC;
  Eigen::Vector3d u_torque = sc_obj.u_IC;
  Eigen::Vector3d u_thruster = Eigen::Vector3d::Zero();

  auto commands = thruster_obj.torque_to_command(u_torque);
  // Begin sim
  while (t < tFinal){
    // Save current state
    out_times.push_back(t);
    out_states.push_back(state);
    out_controls.push_back(u_thruster);

    std::cout << "Time: " << t << std::endl;

    // Check if it's time for next control
    if (t >= t_next_control){
      // Get desired torque command from PI controller
      u_torque = controller_obj.pi_controller(state, sc_obj.I_SC);
      // Convert torque (Nm) to thruster command (PWM)
      commands = thruster_obj.torque_to_command(u_torque);
      // Increment control timer
      t_next_control += controller_obj.dt_control;
    }

    // Run Thruster model
    u_thruster = thruster_obj.thruster_model(commands);

    // Run Rigid Body Dynamics
    state = sc_obj.rk4_step(state, u_thruster);

    // Increment time
    t += sc_obj.dt_sim;
  }

  // Output data to CSV
  for (size_t i = 0; i < out_times.size(); ++i){
    output_file << out_times[i] << ",";
    // Vehicle Position in inertial frame
    output_file << out_states[i](0) << "," << out_states[i](1) << "," << out_states[i](2) << ",";
    // Vehicle Velocity in inertial frame
    output_file << out_states[i](3) << "," << out_states[i](4) << "," << out_states[i](5) << ",";
    // Quaternion [x y z w]
    output_file << out_states[i](6) << "," << out_states[i](7) << "," << out_states[i](8) << "," << out_states[i](9) << ",";
    // Angular Velocity
    output_file << out_states[i](10) << "," << out_states[i](11) << "," << out_states[i](12) << ",";
    // Control Torques
    output_file << out_controls[i](0) << "," << out_controls[i](1) << "," << out_controls[i](2) << "\n";

  }

  output_file.close();

  std::cout << "Sim Complete. Data saved to file" << std::endl;


  return 0;
}