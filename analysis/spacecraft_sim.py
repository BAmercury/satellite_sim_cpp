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
        # Max Control Torques of RCS Thrusters [Nm]
        self.max_torque = np.array([320.0, 220.0, 220.0])

        # LLO Orbit Parameters
        self.inclination = np.deg2rad(102.6) # [rad] inclination
        self.raan = np.deg2rad(-132.0) # [rad] Long of Ascending Node 
        self.arg_peri = np.deg2rad(146.6) # [rad] Argument of Periapsis
        self.nu = 0 # [rad] True anomaly

        self.orbit_period =  8.2 # [hr] Orbital Period

        self.moon_radius = 1737.5 # [km] moon radius
        self.moon_mu = 4902.8001 # [km^3/s^2] Lunar gravitational constant 

        self.peri_alt = 59.9 # [km] Periapsis Alt
        self.apo_alt = 6016.1 # [km] Apoapsis Alt
        self.ra = self.apo_alt + self.moon_radius # [km] Periapsis Radius
        self.rp = self.peri_alt + self.moon_radius # [km] Apoapsis Radius
        self.semi_a = (self.ra + self.rp) / 2 # [km] Semi-Major Axis
        self.ecc = (self.ra - self.rp) / (self.ra - self.rp)


        # Set Initial Conditions State Vector: [pos, vel, attitude quat, angular velocity]

        # Convert orbital elements into position and velocity vector (In MCI Frame) as IC
        # Position and Velocity Inertial frame [m], [m/s]
        r0, v0 = self.get_inertial_pos_vel_from_orbit() # Returns as km and km/s
        self.r0 = r0 * 1000 # [m] Position Inertial Frame 
        self.v0 = v0 * 1000 # [m/s] Velocity Inertial Frame
        self.q0 = R.from_quat([0, 0, 0, 1]) # Scalar last is assumed by default, MCI to Body Frame
        self.w0 =  np.array([0.01, 0.01, 0.01])  # Initial angular rates [rad/s]

        self.x_IC = np.concatenate([r0, v0, q0.as_quat(), w0]) # S/C Initial State Vector
        self.u_IC = np.array([0, 0, 0]) # Initial controls if needed


        # Simulation Parameters (Static)
        self.dt_sim = 0.01 # 100 Hz

    def get_inertial_pos_vel_from_orbit(self):
        # Calculate semi-latus rectum
        l = self.semi_a * (1-self.ecc**2)
        # Get Specific Angular Momentum
        h = np.sqrt(l*self.moon_mu)
        # Position and Velocity vectory in Perifocal Frame
        r_w = h**2 / self.moon_mu / (1 + ecc * np.cos(self.nu)) * np.array((np.cos(self.nu), np.sin(self.nu), 0))
        v_w = self.moon_mu / h * np.array((-np.sin(self.nu), self.ecc + np.cos(self.nu), 0))
        # DCM Perifocal frame to MCI frame
        C_pf2mci = Rotation.from_euler("ZXZ", [self.raan, self.inclination, self.arg_peri])
        r_I = C_pf2mci.as_matrix() @ r_w
        v_I = C_pf2mci.as_matrix() @ v_w

        return r_I, v_I
    
    def gravity_gradient_torque(self, r_I, q_i2b)




