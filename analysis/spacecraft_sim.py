import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial.transform import Rotation as R


# Define SpaceCraft Class
# Contains parameters for Spacecraft as well as helper functions to be used for simulation
class Spacecraft:
    def __init__(self):
        # S/C Moment of Inertia Matrix [kg-m^2]
        self.I_SC = np.array([
            [2402.52, 36.33, 78.17],
            [36.33, 2072.58, -67.71],
            [78.17, -67.71, 2022.10]
        ])
        # Get inverse
        self.I_SC_inv = np.linalg.inv(self.I_SC)

        # Max Control Torques of RCS Thrusters [Nm]
        self.max_torque = np.array([320.0, 220.0, 220.0])

        # LLO Orbit Parameters
        self.inclination = np.deg2rad(102.6) # [rad] inclination
        self.raan = np.deg2rad(-132.0) # [rad] Long of Ascending Node 
        self.arg_peri = np.deg2rad(146.6) # [rad] Argument of Periapsis
        self.nu = 0 # [rad] True anomaly

        self.orbit_period =  8.2 # [hr] Orbital Period

        self.moon_radius = 1737.4 # [km] moon radius
        self.moon_mu = 4902.8001 # [km^3/s^2] Lunar gravitational constant 

        self.peri_alt = 59.9 # [km] Periapsis Alt
        self.apo_alt = 6016.1 # [km] Apoapsis Alt
        self.ra = self.apo_alt + self.moon_radius # [km] Periapsis Radius
        self.rp = self.peri_alt + self.moon_radius # [km] Apoapsis Radius
        self.semi_a = (self.ra + self.rp) / 2 # [km] Semi-Major Axis
        self.ecc = (self.ra - self.rp) / (self.ra + self.rp)


        # Set Initial Conditions State Vector: [pos, vel, attitude quat, angular velocity]

        # Convert orbital elements into position and velocity vector (In MCI Frame) as IC
        # Position and Velocity Inertial frame [m], [m/s]
        r0, v0 = self.get_inertial_pos_vel_from_orbit() # Returns as km and km/s
        self.r0 = r0  # [km] Position Inertial Frame 
        self.v0 = v0  # [km/s] Velocity Inertial Frame
        self.q0 = R.from_quat([0, 0, 0, 1]) # Scalar last is assumed by default, MCI to Body Frame
        self.w0 =  np.array([0.01, 0.01, 0.01])  # Initial angular rates [rad/s]

        self.x_IC = np.concatenate([self.r0, self.v0, self.q0.as_quat(), self.w0]) # S/C Initial State Vector
        self.u_IC = np.array([0, 0, 0]) # Initial controls if needed


        # Simulation Parameters (Static)
        self.dt_sim = 0.01 # 100 Hz

    def get_inertial_pos_vel_from_orbit(self):
        # Calculate semi-latus rectum
        l = self.semi_a * (1-self.ecc**2)
        # Get Specific Angular Momentum
        h = np.sqrt(l*self.moon_mu)
        # Position and Velocity vectory in Perifocal Frame
        r_w = h**2 / self.moon_mu / (1 + self.ecc * np.cos(self.nu)) * np.array((np.cos(self.nu), np.sin(self.nu), 0))
        v_w = self.moon_mu / h * np.array((-np.sin(self.nu), self.ecc + np.cos(self.nu), 0))
        # DCM Perifocal frame to MCI frame
        C_pf2mci = R.from_euler("ZXZ", [self.raan, self.inclination, self.arg_peri])
        r_I = C_pf2mci.as_matrix() @ r_w
        v_I = C_pf2mci.as_matrix() @ v_w

        return r_I, v_I
    
    def get_time_window(self, num_orbits):
        # Function to help get a time window based on number of orbits
        period = 2*np.pi * np.sqrt((self.semi_a**3)/self.moon_mu) # seconds
        tFinal = period*num_orbits
        return tFinal
    
    def gravity_gradient_torque(self, r_I, R_i2b):
        # From Wertz https://link.springer.com/book/10.1007/978-94-009-9907-7
        # Assume geometric center and center of mass of vehicle are the same
        # Rotation from inertial frame to body frame
        C_i2b = R_i2b.as_matrix()
        r_B = C_i2b @ r_I
        # Get unit vec of r_B
        r_mag_B = np.linalg.norm(r_B)
        if r_mag_B < 1e-10:  # Add numerical stability check
            return np.zeros(3)
        r_B_unit = r_B / r_mag_B
         # Verify unit vector (optional)
        r_B_unit = r_B_unit / np.linalg.norm(r_B_unit)

        # Gravity Gradient torque
        temp = 3 * self.moon_mu / (r_mag_B**3)
        return temp * np.cross(r_B_unit, self.I_SC @ r_B_unit)

    
    def spacecraft_dynamics(self, x, u):
        # Spacecraft Rigid Body Equations of Motion
        # Extract states x = r0, v0, q0.as_quat(), w0

        r = x[0:3] # R_I
        v = x[3:6] # v_I
        # Scipy automatically normalized quaterion
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.from_quat.html
        R_i2b = R.from_quat(x[6:10]) # qi2b
        w = x[10:13] # w_b_n_Brf
        # Skew Symmetric Matrix
        w_skew = np.array([
            [0, -w[2], w[1]],
            [w[2], 0, -w[0]],
            [-w[1], w[0], 0]
        ])

        # Quaternion Kinematics
        q = R_i2b.as_quat()
        q_matrix = np.array([
            [q[3], -q[2], q[1]],
            [q[2], q[3], -q[0]],
            [-q[1], q[0], q[3]],
            [-q[0], -q[1], -q[2]]
        ])

        q_dot = (0.5) * (q_matrix @ w)

        # Gravity Gradient Torque (In Body Frame)
        L_gg_b = self.gravity_gradient_torque(r, R_i2b)
        # Sum Moments
        L = u + L_gg_b
        # Angular Acceleration (Body Frame)
        w_dot_b_n_Brf = self.I_SC_inv @ ( (-w_skew @ (self.I_SC @ w)) + L)

        # Translational Dynamics
        # Gravity Model
        # From Wertz https://link.springer.com/book/10.1007/978-94-009-9907-7
        r_dot = v # Inertial Frame
        # Model with Keplerian Motion (2-body)
        # Assume M_moon >> M_sc
        # https://tfaws.nasa.gov/wp-content/uploads/Rickman-Presentation.pdf
        r_I_mag = np.linalg.norm(r)
        v_dot = -self.moon_mu * r / (r_I_mag**3)# Translational Acceleration in Inertial Frame

        # Concatenate
        #r_dot = r_dot.flatten()
        #v_dot = v_dot.flatten()
        #q_dot = q_dot.flatten()
        #w_dot_b_n_Brf = w_dot_b_n_Brf.flatten()
        x_dot = np.concatenate([r_dot, v_dot, q_dot, w_dot_b_n_Brf]) # S/C Initial State Vector
        return x_dot
    
    def rk4_step(self, t, x, u):
        # Runge-Kutta 4 integration step
        f1 = self.spacecraft_dynamics(x, u)
        f2 = self.spacecraft_dynamics(x + 0.5*self.dt_sim*f1, u)
        f3 = self.spacecraft_dynamics(x + 0.5*self.dt_sim*f2, u)
        f4= self.spacecraft_dynamics(x + self.dt_sim*f3, u)

        return (x + (self.dt_sim/6.0) * (f1 + 2*f2 + 2*f3 + f4))

def pi_controller(x, target, Kp, Ki, integral_error, I_SC):
    rate_rad = x[10:13]
    rate_error = target - rate_rad
    integral_limit = np.array([100.0, 100.0, 100.0])
    integral_error += rate_error * 0.1
    integral_error = np.clip(integral_error, -integral_limit, integral_limit)
    # Cross coupling term
    cross_coupling = np.cross(rate_rad, I_SC @ rate_rad)
    u = (Kp * rate_error) + (Ki * integral_error) + cross_coupling
    return u, integral_error

class ThrusterController:
    def __init__(self):
        # System parameters
        self.fsw_timestep = 0.1  # Flight software timestep (seconds)
        self.max_torque = np.array([320.0, 220.0, 220.0])  # [X, Y, Z] in N.m
        
        # Deadband threshold (5% of max torque)
        self.deadband_percentage = 0.05
        self.deadband = self.max_torque * self.deadband_percentage
        
    def apply_deadband(self, torque):
        """Apply deadband to prevent thruster chattering near zero"""
        # For x-axis command sbelow +/-16 Nm is ignored
        # For y/z axis commands below +/- 11Nm is ignored
        return np.where(np.abs(torque) < self.deadband, 0, torque)
    
    def compute_pulse_command(self, desired_torque):
        """
        Compute thruster pulse commands based on desired torque and FSW timestep
        
        Parameters:
        desired_torque: numpy.ndarray - Desired torque for each axis [X, Y, Z] in N.m
        
        Returns:
        tuple: (pulse_widths, directions)
            - pulse_widths: Duration thruster should be on at full force (in seconds)
            - directions: Direction for each axis (1 or -1)
        """
        # Apply deadband
        torque = self.apply_deadband(desired_torque)
        
        # Determine direction (1 or -1)
        directions = np.sign(torque)
        directions = np.where(directions == 0, 1, directions)  # Replace 0 with 1
        
        # Calculate required pulse width based on torque ratio
        magnitude_ratio = np.abs(torque) / self.max_torque
        pulse_widths = magnitude_ratio * self.fsw_timestep
        
        # Clip pulse widths to FSW timestep
        pulse_widths = np.clip(pulse_widths, 0, self.fsw_timestep)
        
        return pulse_widths, directions
    def torque_to_command(self, u):
        """
        Convert desired torque to thruster commands
        
        Parameters:
        u: numpy.ndarray - Desired torque for each axis [X, Y, Z] in N.m
        
        Returns:
        dict: Command parameters for each axis including timing and scaling
        """
        if not isinstance(u, np.ndarray):
            u = np.array(u)
            
        if u.shape != (3,):
            raise ValueError("Input torque must be a 3-element array [X, Y, Z]")
        
        pulse_widths, directions = self.compute_pulse_command(u)
        
        # Create command structure
        commands = {
            'pulse_widths': pulse_widths,    # Duration of pulse in seconds
            'directions': directions,         # Direction of thrust (+1 or -1)
            'timestep': self.fsw_timestep    # FSW timestep for reference
        }
        
        return commands

# Main simulation loop
if __name__ == "__main__":

    spacecraft_obj = Spacecraft()
    #spacecraft_obj.dt_sim = 0.1
    # Set time vector
    num_orbits = 1
    #tFinal = spacecraft_obj.get_time_window(num_orbits) / 4
    spacecraft_obj.dt_sim = 0.05
    tFinal = 3600.0 / 4
    t = 0.0
    dt_control = 0.1
    dt_pwm = 0.1
    t_next_control = 0.0

    # Storage for results
    times = []
    states = []
    controls = []

    # Control Params
    target = np.array([np.deg2rad(0.5), 0., 0.])  # Target angular rates

    #ts = 10
    Kp = np.array([20000.0, 900.0, 900.0])
    Ki = np.array([1000.0, 1000.0, 1000.0])

    # Derivative gain - for rate damping
    #Kd = np.diag([200.0, 800.0, 800.0])   
    #Kp =spacecraft_obj.I_SC @ Kp
    #Ki = spacecraft_obj.I_SC @ Ki
    integral_error = np.zeros(3)
    # ICs
    state = spacecraft_obj.x_IC
    u = spacecraft_obj.u_IC
    u_torque = spacecraft_obj.u_IC
    duty_cycle = np.array([0, 0, 0])
    pwm_signal = np.array([0, 0, 0])
    signs = np.array([1, 1, 1])

    converter = ThrusterController()

    commands = converter.torque_to_command(u)
    duty_cycle = commands["pulse_widths"]
    signs = commands["directions"]
    while t < tFinal:
        # Store current state
        times.append(t)
        states.append(state.copy())
        controls.append(u_torque.copy())

        print(t)

        # Controls
        if t >= t_next_control:
            # Get control torques (3 axis)
            u, integral_error = pi_controller(state, target, Kp, Ki, integral_error, spacecraft_obj.I_SC)
            commands = converter.torque_to_command(u)
            duty_cycle = commands["pulse_widths"]
            signs = commands["directions"]
            t_next_control = t_next_control + dt_control


        # Thruster model
        i = 0
        for signal in duty_cycle:
            # If thruster pulsewidth command is equal to or greater than FSW timestep
            # Thruster stays on continuously 
            if signal >= dt_control:
                u_torque[i] = signs[i]*spacecraft_obj.max_torque[i] # T max
            else:
                u_torque[i] = signs[i]*spacecraft_obj.max_torque[i]*(duty_cycle[i]/dt_control)
            i = i + 1
        # Step dynamics
        state= spacecraft_obj.rk4_step(t, state, u_torque)

        t = t + spacecraft_obj.dt_sim
        #print(t)

    # Get Angular Rates [pos, vel, attitude quat, angular velocity]
    states = np.asarray(states)
    controls = np.asarray(controls)
    times = np.array(times)
    #rates_deg = np.rad2deg(states[:, 10:13])
    rates_rad = states[:, 10:13]


    plt.figure(figsize=(12, 8))
    plt.subplot(211)
    plt.plot(times, np.rad2deg(rates_rad))
    plt.ylabel('Angular Rates (deg/s)')
    plt.legend(['X', 'Y', 'Z'])
    plt.grid(True)
    plt.subplot(212)
    plt.plot(times, controls)
    plt.xlabel('Time (s)')
    plt.ylabel('Control Torque [Nm]')
    plt.legend(['X', 'Y', 'Z'])
    plt.grid(True)

    r = states[:, 0:3]
    plt.figure(figsize=(12, 8))
    plt.plot(times, r)
    plt.ylabel('S/C Position (m)')
    plt.legend(['X', 'Y', 'Z'])
    plt.grid(True)
    plt.xlabel('Time (s)')

    v = states[:, 3:6]
    plt.figure(figsize=(12, 8))
    plt.plot(times / 3600, v)
    plt.ylabel('S/C Velocity (km/s)')
    plt.legend(['X', 'Y', 'Z'])
    plt.grid(True)
    plt.xlabel('Time (hr)')
    
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(r[:, 0], r[:, 1], label='Orbit path')
    plt.xlabel('x [km]')
    plt.ylabel('y [km]')
    plt.title('Orbital Path in MCI Frame')
    plt.axis('equal')
    plt.grid()

    # Plot velocity magnitude over time
    vel_magnitude = np.linalg.norm(v, axis=1)
    plt.subplot(1, 2, 2)
    plt.plot(times / 3600, vel_magnitude)  # Convert time to hours
    plt.xlabel('Time [hours]')
    plt.ylabel('Velocity Magnitude [km/s]')
    plt.title('Velocity Magnitude over Time')
    plt.grid()


    plt.tight_layout()

    plt.show()


            
        

    





