import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial.transform import Rotation as R

# Single-Body Spacecraft Dynamics 
# Embedded wheels, thrusters, flexible modes
# Modeling full 6-DOF dynamics 

# Define SpaceCraft Class
# Contains parameters for Spacecraft as well as helper functions to be used for simulation
class Spacecraft:
    def __init__(self):
        # S/C Moment of Inertia Matrix [kg-m^2]
        # Big Black Monolith (BBM)
        self.I_SC = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 4.0]
        ])
        # Get inverse
        self.I_SC_inv = np.linalg.inv(self.I_SC)

        # Set Initial Conditions State Vector: [attitude quat, angular velocity, wheel ang vel, pos, vel]

        # Convert orbital elements into position and velocity vector (In MCI Frame) as IC
        self.q0 = R.from_quat([0, 0, 0, 1]) # Scalar last is assumed by default, MCI to Body Frame
        # Major Spin Axis
        #self.w0 =  np.array([0.1, 0.0, 5.0])  # Initial angular rates [rad/s]
        # Minor Spin Axis
        #self.w0 =  np.array([5.0, 0.1, 0.0])  # Initial angular rates [rad/s]
        # Intermediate Axis
        #self.w0 =  np.array([0.1, 5.0, 0.0])*np.pi/180.0  # Initial angular rates [rad/s]
        self.w0 = np.array([0.0, 0.0, 0.0])*np.pi/180.0

        # S/C position and velocity
        self.pos0 = np.array([10, 0, 10])
        self.vel0 = np.array([0, -0.01, 0])
        # Reference Orbit, 700-km circular
        self.mu = 3.98600e14 # m^3/sec^2
        self.OrbRad = 6378.145e3+700.0e3
        self.muR3 = self.mu/(self.OrbRad**3)
        self.OrbRate = np.sqrt(self.muR3)

        # Wheel class
        self.wheel_model = ReactionWheels()



        self.x_IC = np.concatenate([self.q0.as_quat(), self.w0, self.wheel_model.x_IC, self.pos0, self.vel0]) # S/C Initial State Vector
        self.u_IC = np.array([0, 0, 0, 0, 0, 0, 0]) # Initial controls if needed [rw1 rw2 rw3 thx thy thz]



        # Simulation Parameters (Static)
        self.dt_sim = 0.01 # 100 Hz


    
    def spacecraft_dynamics(self, t, x, u):
        # Spacecraft Rigid Body Equations of Motion
        # Extract states x = r0, v0, q0.as_quat(), w0
        # Scipy automatically normalized quaterion
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.from_quat.html
        R_i2b = R.from_quat(x[0:4]) # qi2b
        w = x[4:7] # w_b_n_Brf
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

        # Process control inputs
        Tw = u[0:4] # Torque commands for reaction wheel
        T_thr = u[5:7] # Thruster torque commands

        # Embedded Reaction wheels
        # Wheel angular acceleration
        wheel_x_dot = Tw / Wheel_model.J_w
        # Wheel torques from wheel frame to body frame
        Lw = self.wheel_model.A @ Tw # wheel torques
        # angular momentum in body frame
        w_w = x[7:11] # Get wheel angular velocity
        H = (self.I_SC @ w)  + self.wheel_model.A @ (self.wheel_model.J_w*w_w)
        # Sum Moments
        L = -(Lw + T_thr)# Reaction torque from wheel
        # Total Angular Acceleration (Body Frame)
        w_dot_b_n_Brf = self.I_SC_inv @ ( (-w_skew @ (H)) + L)

        # Translational Dynamics
        PosR = x[11:14]
        VelR = x[14:17]

        OrbPosN = self.OrbRad * np.array([np.cos(self.OrbRate*t), np.sin(self.OrbRate*t), 0]) # OrbPosN R
        r = OrbPosN + PosR

        fq = self.EnckeFQ(r, PosR)

        acc = np.array([0.0, 0.0, 0.0]) # F/m

        rdot = VelR
        vdot = acc - self.muR3*(PosR+fq*r)
        



        # Concatenate
        #q_dot = q_dot.flatten()
        #w_dot_b_n_Brf = w_dot_b_n_Brf.flatten()
        x_dot = np.concatenate([q_dot, w_dot_b_n_Brf, wheel_x_dot, rdot, vdot]) # S/C Initial State Vector
        return x_dot
    
    def EnckeFQ(self, r, PosR):
        delta = PosR
        q = np.dot(delta, (delta - 2.0*r)) / (np.dot(r, r))
        q_dot = (q*(3.0+q*(3.0+q))/(1.0+np.sqrt((1.0+q)**3)))
        return q_dot
    
    def rk4_step(self, t, x, u):
        # Runge-Kutta 4 integration step
        f1 = self.spacecraft_dynamics(t, x, u)
        f2 = self.spacecraft_dynamics(t + self.dt_sim/2, x + 0.5*self.dt_sim*f1, u)
        f3 = self.spacecraft_dynamics(t + self.dt_sim/2, x + 0.5*self.dt_sim*f2, u)
        f4 = self.spacecraft_dynamics(t + self.dt_sim, x + self.dt_sim*f3, u)

        return (x + (self.dt_sim/6.0) * (f1 + 2*f2 + 2*f3 + f4))
    
class FlexModel:
    """
    Flexible dynamics model with 2 modes
    """
    def __init__(self):
        # Flexible modes
        # 2 modes
        self.w_f = [1.2, 5.0]*2*np.pi # Natural Frequency
        self.d_f = [0.005, 0.0025] # Damping Ratio
        # one input (force node), 3 RotDOF x 2 modes
        self.flex_node_in = [[
            [0.1, -0.6],
            [-0.2, 0.5],
            [0.3, 0.4]
        ]]
        # One output node (measurement), 3 RotDOF x 2 modes
        self.flex_node_out = [[
            [0.25, 0.15],
            [0.35, 0.5],
            [-0.15, -0.25]
        ]]
        # Initial conditions for flexible dynamics
        eta0 = np.zeros(2)
        xi0 = np.zeros(2)
        self.x = np.concatenate([eta0, xi0])
    def flex_eom(self, u):
        # Get states:
        eta = self.x[0:2]
        xi = self.x[2:4]
        # Flex Deflection at Output node
        thetaf = self.flex_node_out @ eta  
        # Flex deflection at input node
        Tf = np.transpose(self.flex_node_in) @ u
        # Flex
        xi_dot = Tf - 2*self.d_f *self.w_f*xi - self.w_f*self.w_f*eta
        eta_dot = xi

        x_dot = np.concatenate([eta_dot, xi_dot])


class ReactionWheels:
    """
    4 wheels fixed in the body, unit axis vector a_hat and a rotary momen tof inertia J_w
    States:
        angular momentum of wheel h
    4-wheel pyramid configuration
    """
    def __init__(self):
        alpha = 45*np.pi/180
        eta = 30*np.pi/180
        x = np.cos(alpha)*np.cos(eta)
        y = np.sin(alpha)*np.cos(eta)
        z = np.sin(eta)
        self.A = np.array([
            [x, -x, -x, x],
            [y, y, -y, -y],
            [z, -z, z, -z]
        ])
        self.Ap = np.linalg.pinv(self.A)
        self.J_w = np.array([1.0, 1.0, 1.0, 1.0])
        self.x_IC = np.array([0.0, 0.0, 0.0, 0.0]) # Wheel angular velocity states

    
class Thrusters:
    """
    Thruster exerts torque and force on body B
    Thruster is fixed in B, at position p wrt the origin of convenience O_B
    Thrust axis defined by unit vector a_hat
    Assume thruster can exert a force magnitude f in range 0 <= f <= f_max
    Pulsed Thrusters: 0 or f_max can be modeled by averaging the force they exert over the simulation timestep to yield correct impulse delivered

    T = r x F = (p - p_B*) x (f*a_hat) [moment arm cross Force]

    Simulating 12 thrusters 
    """
    def __init__(self):
        self.CmPosB = np.array([0.0, 0.0, 0.0])  # Spacecraft Mass center (0 by default)
        self.ThrPosB = np.array([ # Thruster position in Brf
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, -1.0, 1.0]
        ])
        self.ThrAxis = np.array([ # Thruster force axis
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0]
        ])   
        self.rxA = np.zeros((12,3))# thruster moment arm
        self.f = 2.0e-4 # thruster max force
        for i in range(0,12):
            self.rxA[i,:] = np.cross(self.ThrPosB[i, :] - self.CmPosB, self.ThrAxis[i, :])
    def step(self, Tcmd):
        # Step thruster
        # Allocate torque command to appropiate thrusters
        T = np.zeros(3)
        for i in range(0, 12):
            TorxA = np.dot(Tcmd, self.rxA[i, :])
            # if thruster is activated, pulse on
            if (TorxA > 0.0):
                f = 2.0e-4
                T = T + self.f*self.rxA[i,:]
        return T


# Main simulation loop
if __name__ == "__main__":

    spacecraft_obj = Spacecraft()
    #spacecraft_obj.dt_sim = 0.1
    # Set time vector
    #num_orbits = 1
    #tFinal = spacecraft_obj.get_time_window(num_orbits) / 4
    #spacecraft_obj.dt_sim = 0.05
    #tFinal = 3600.0 / 4
    tFinal = 10.0
    t = 0.0
    # Storage for results
    times = []
    states = []
    # ICs
    state = spacecraft_obj.x_IC
    u = spacecraft_obj.u_IC

    spacecraft_obj.dt_sim = 0.1

    Thruster_Model = Thrusters()
    Wheel_model = ReactionWheels()

    while t < tFinal:
        # Store current state
        times.append(t)
        states.append(state.copy())

        #print(t)
        # Controls
        q = state[0:4]
        w = state[4:7]
        # Control law
        #Tcmd = -3.0e-2*w - 1.0e-3*q[0:3]*np.sign(q[3]) # 3 axis torque command reaction wheels
        # Psuedo-inverse control allocation (would need to include thrusters here if using that for attitude control)
        #Tw = -Wheel_model.Ap @ Tcmd # 4 wheels
        # Thruster model
        #u = Thruster_Model.step(Tcmd)
        #print(Tw)
        #u = np.concatenate([Tw, np.array([0, 0, 0])])

        # Impulse torque at t = 1
        if (t == 1.0):
            T = np.array([3,2,1]) # Nm
        else:
            T = np.array([0, 0, 0])
        Tw = np.array([0,0,0,0])
        u = np.concatenate([Tw, T])
        # Step dynamics
        state= spacecraft_obj.rk4_step(t, state, u)

        t = t + spacecraft_obj.dt_sim
    # Get Angular Rates [pos, vel, attitude quat, angular velocity]
    states = np.asarray(states)
    times = np.array(times)
    #rates_deg = np.rad2deg(states[:, 10:13])
    quats = states[:, 0:4]
    rates_rad = states[:, 4:7]

    plt.figure(figsize=(12, 8))
    plt.plot(times, quats)
    plt.ylabel('Quaternion')
    plt.legend(['q_1', 'q_2', 'q_3', 'q_4'])
    plt.grid(True)


    plt.figure(figsize=(12, 8))
    plt.plot(times, np.rad2deg(rates_rad))
    plt.ylabel('Angular Rates (deg/s)')
    plt.legend(['X', 'Y', 'Z'])
    plt.grid(True)

    pos = states[:, 11:14]
    vel = states[:, 14:17]

    plt.figure(figsize=(12, 8))
    plt.plot(times, pos)
    plt.ylabel('Pos (m)')
    plt.legend(['X', 'Y', 'Z'])

    plt.figure(figsize=(12, 8))
    plt.plot(times, vel)
    plt.ylabel('Vel (m/s)')
    plt.legend(['X', 'Y', 'Z'])


    plt.tight_layout()

    plt.show()


            
        

    





