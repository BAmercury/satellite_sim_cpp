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
        self.M_SC = 1

        # Set Initial Conditions State Vector: [attitude quat, angular velocity, wheel ang vel, pos, vel]

        # Convert orbital elements into position and velocity vector (In MCI Frame) as IC
        self.q0 = R.from_quat([0, 0, 0, 1]) # Scalar last is assumed by default, MCI to Body Frame
        # Major Spin Axis
        #self.w0 =  np.array([0.1, 0.0, 5.0])  # Initial angular rates [rad/s]
        # Minor Spin Axis
        #self.w0 =  np.array([5.0, 0.1, 0.0])  # Initial angular rates [rad/s]
        # Intermediate Axis
        #self.w0 =  np.array([0.1, 5.0, 0.0])*np.pi/180.0  # Initial angular rates [rad/s]
        self.w0 = np.array([5.0, 5.0, 5.0])*np.pi/180.0

        # S/C position and velocity
        self.pos0 = np.array([1, 0, 0]) # 10m offset in X and Z
        self.vel0 = np.array([0, -0.001, 0]) # -0.01 m/s in Y direction
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

        # Environmental torques
        T_env = 1.0e-4*np.array([1,2,3]) + 0.1e-4*np.sin(0.001*t) * np.array([0,1,0])

        # Process control inputs
        Tw = u[0:4] # Torque commands for reaction wheel, this is torque Tw exerted between wheel and the body
        # Convention: +Tw applied to the wheel, -Tw applied to the body
        T_thr = u[4:7] # Thruster torque commands
        F_thr_b = u[7:10] # Thruster force commands body frame

        # Embedded Reaction wheels
        # Wheel angular acceleration
        wheel_x_dot = Tw / Wheel_model.J_w
        # Wheel torques from wheel frame to body frame
        Lw = self.wheel_model.A @ Tw # wheel torques
        # angular momentum in body frame
        w_w = x[7:11] # Get wheel angular velocity
        H = (self.I_SC @ w)  + self.wheel_model.A @ (self.wheel_model.J_w*w_w)
        # Sum Moments
        L = -Lw + T_env + T_thr # Reaction torque from wheel, thruster torque, and environmental torque
        # Total Angular Acceleration (Body Frame)
        # I*w_dot = T - Tw - wxH
        w_dot_b_n_Brf = self.I_SC_inv @ ( (-w_skew @ (H)) + L)
        #w_dot_b_n_Brf = self.I_SC_inv @ ( (-w_skew @ (self.I_SC @ w)) + L)

        # Translational Dynamics
        PosR = x[11:14]
        VelR = x[14:17]

        OrbPosN = self.OrbRad * np.array([np.cos(self.OrbRate*t), np.sin(self.OrbRate*t), 0]) # OrbPosN R
        r = OrbPosN + PosR

        fq = self.EnckeFQ(r, PosR)

        #acc = np.array([0.0, 0.0, 0.0]) # F/m
        F_thr_N = np.transpose(R_i2b.as_matrix()) @ F_thr_b
        acc = F_thr_N / self.M_SC

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
    6 nodes:
        1. Wheels (input torque node)
        2-5. Thrusters (input force node)
        6. Sensor (output deflection node) 
    """
    def __init__(self):
        # Flexible modes
        # 2 modes
        self.w_f = np.array([1.2, 5.0])*2*np.pi # Natural Frequency
        self.d_f = np.array([0.005, 0.0025]) # Damping Ratio
        # one input (force node), 3 RotDOF x 2 modes
        self.theta = np.zeros((6,3,2)) # 6 nodes, 3 DOF, 2 modes
        # Node 1: Wheels, input torque node
        self.theta[0,:,:] = np.array([
            [0.1, -0.6],
            [-0.2, 0.5],
            [0.3, 0.4]
            ])
        # Node 6, Gyro, rotational deflection output node
        self.theta[5,:,:] = np.array([
            [0.05, 0.06],
            [0.07, 0.08],
            [0.09, 0.10]
        ])
        # Thruster nodes exert forces (input node)
        # Sensor is not sensitive to translation
        # Wheels are not exerting forces (not considering jitter in this sim)
        self.PSI = np.zeros((6,3,2))
        self.PSI[1,:,:] = np.array([
            [0.06, 0.05],
            [0.04, 0.03],
            [0.02, 0.01]
        ])
        self.PSI[2,:,:] = np.array([
            [0.06, 0.05],
            [0.04, 0.03],
            [0.02, 0.01]
        ])
        self.PSI[3,:,:] = np.array([
            [0.06, 0.05],
            [0.04, 0.03],
            [0.02, 0.01]
        ])
        self.PSI[4,:,:] = np.array([
            [0.06, 0.05],
            [0.04, 0.03],
            [0.02, 0.01]
        ])
        # self.flex_node_in = np.array([
        #     [0.1, -0.6],
        #     [-0.2, 0.5],
        #     [0.3, 0.4]
        # ])
        # # One output node (measurement), 3 RotDOF x 2 modes
        # self.flex_node_out = np.array([
        #     [0.25, 0.15],
        #     [0.35, 0.5],
        #     [-0.15, -0.25]
        # ])
        # Initial conditions for flexible dynamics
        eta0 = np.zeros(2)
        xi0 = np.zeros(2)
        self.x0 = np.concatenate([eta0, xi0])

        self.wheel_model = ReactionWheels()
        self.thrust_model = Thrusters()

        # dt
        self.dt_sim = 0.01
    def flex_dynamics(self, t, x, u):
        # Get states:
        eta = x[0:2]
        xi = x[2:4]
        # Get inputs
        Tw = u[0:4] # Torque commands for reaction wheel
        WhlTrqB = -self.wheel_model.A @ Tw
        F_thr = u[7:10] # Thruster Force commands
        # Flex deflection at input node
        # Wheel input node
        Flex_Trq = np.zeros(2)
        Flex_Trq[0] = np.transpose(self.theta[0,:,0]) @ WhlTrqB
        Flex_Trq[1] = np.transpose(self.theta[0,:,1]) @ WhlTrqB
        # Thruster input node
        Flex_Frc = np.zeros(2)
        ThrNode = np.array([1,1,1,2,2,2,3,3,3,4,4,4]);
        for i in range(0,12):
            for j in range(0, 1):
                Flex_Frc[j] = Flex_Frc[j] + np.transpose(self.PSI[ThrNode[i],:,j]) @ (self.thrust_model.f*self.thrust_model.ThrAxis[i,:])
        # Flex
        xi_dot = Flex_Frc + Flex_Trq - 2*self.d_f *self.w_f*xi - self.w_f*self.w_f*eta
        eta_dot = xi

        x_dot = np.concatenate([eta_dot, xi_dot])
        return x_dot
    def get_thetaf(self, eta):
        # Returns flex angle at output node
        return self.theta[5,:,:] @ eta
    def get_thetadotf(self, eta):
        # Returns  Flex rate at output node
        return self.theta[5,:,:] @ eta

    def rk4_step(self, t, x, u):
        # Runge-Kutta 4 integration step
        f1 = self.flex_dynamics(t, x, u)
        f2 = self.flex_dynamics(t + self.dt_sim/2, x + 0.5*self.dt_sim*f1, u)
        f3 = self.flex_dynamics(t + self.dt_sim/2, x + 0.5*self.dt_sim*f2, u)
        f4 = self.flex_dynamics(t + self.dt_sim, x + self.dt_sim*f3, u)

        return (x + (self.dt_sim/6.0) * (f1 + 2*f2 + 2*f3 + f4))


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
    # Clusters of 3 on 4 corners of a cube
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
                T = T + self.f*self.rxA[i,:]
        return T
    def momentum_unload(self, H):
        # Momentum unload with thrusters
        Tunl = -1.0e-4 * H
        TrqB = np.zeros(3)
        FrcB = np.zeros(3)
        for i in range(0,12):
            TorxA = np.dot(Tunl, self.rxA[i, :])
            # if thruster is activated, pulse on
            if (TorxA > 0.0):
                FrcB = FrcB + self.f*self.ThrAxis[i,:]
                TrqB = TrqB + self.f*self.rxA[i,:]
        return np.concatenate([TrqB, FrcB]) 


# Main simulation loop
if __name__ == "__main__":

    spacecraft_obj = Spacecraft()
    #spacecraft_obj.dt_sim = 0.1
    # Set time vector
    #num_orbits = 1
    #tFinal = spacecraft_obj.get_time_window(num_orbits) / 4
    #spacecraft_obj.dt_sim = 0.05
    #tFinal = 3600.0 / 4
    tFinal = 1000.0
    t = 0.0
    # Storage for results
    times = []
    states = []
    # ICs
    state = spacecraft_obj.x_IC
    u = spacecraft_obj.u_IC

    states_flex = []
    Flex_Model = FlexModel()
    state_flex = Flex_Model.x0

    spacecraft_obj.dt_sim = 0.01

    Thruster_Model = Thrusters()
    Wheel_model = ReactionWheels()

    event_triggered = False

    while t < tFinal:
        # Store current state
        times.append(t)
        states.append(state.copy())
        states_flex.append(state_flex.copy())

        # Get States
        q = state[0:4]
        w = state[4:7]


        # Gyro Model
        eta = state_flex[0:2]
        FlexRate = Flex_Model.get_thetadotf(eta)
        GyroRate = w + FlexRate

        state[4:7] = GyroRate # Update state vector with GyroRate

        # Control law
        Tcmd = -3.0e-2*GyroRate - 1.0e-3*q[0:3]*np.sign(q[3]) # 3 axis torque command reaction wheels
        # Psuedo-inverse control allocation (would need to include thrusters here if using that for attitude control)
        Tw = -Wheel_model.Ap @ Tcmd # 4 wheels
        # Balance Wheels
        w_w = state[7:11] # Wheel angular velocity
        h_w = Wheel_model.J_w * w_w
        Tw = Tw - 1.0e-5*h_w
        H_sys = spacecraft_obj.I_SC @ w + Wheel_model.A @ h_w
        # Thruster model
        # T_thr = [TorqueX, TorqueY, TorqueZ, ForceX, ForceY, ForceZ]
        T_thr = Thruster_Model.momentum_unload(H_sys)
        #u = Thruster_Model.step(Tcmd)
        #print(Tw)
        #u = np.concatenate([Tw, np.array([0, 0, 0])])

        # Impulse torque at t = 1
        # if (t >= 1.0 and not event_triggered):
        #     T = np.array([3,2,1]) # Nm
        #     event_triggered = True
        # else:
        #     T = np.array([0, 0, 0])
        # Tw = np.array([0,0,0,0])
        #T_thr = np.array([0,0,0,0,0,0])
        u = np.concatenate([Tw, T_thr])
        #print(u)
        # Step dynamics
        state= spacecraft_obj.rk4_step(t, state, u)
        state_flex = Flex_Model.rk4_step(t, state_flex, u)

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
    plt.title('Quaternion qi2b')
    plt.grid(True)


    plt.figure(figsize=(12, 8))
    plt.plot(times, np.rad2deg(rates_rad))
    plt.ylabel('Angular Rates (deg/s)')
    plt.legend(['X', 'Y', 'Z'])
    plt.title('Angular Body Rates')
    plt.grid(True)

    pos = states[:, 11:14]
    vel = states[:, 14:17]

    plt.figure(figsize=(12, 8))
    plt.plot(times, pos)
    plt.ylabel('Pos (m)')
    plt.legend(['X', 'Y', 'Z'])
    plt.title('Spacecraft Translational Position')
    plt.grid(True)

    plt.figure(figsize=(12, 8))
    plt.plot(times, vel)
    plt.ylabel('Vel (m/s)')
    plt.legend(['X', 'Y', 'Z'])
    plt.title('Spacecraft Translational Velocity')
    plt.grid(True)

    # Plot flex eta (Dispalcement)
    states_flex = np.asarray(states_flex)
    etas = states_flex[:, 0:2]
    xis = states_flex[:, 2:4]

    # Flex states eta can be used to verify the modal frequency and damping
    plt.figure(figsize=(12, 8))
    plt.plot(times, etas)
    plt.legend(['eta_1', 'eta_1'])
    plt.title('Eta')
    plt.grid(True)

    # Deflections theta at the output node are a linear combination of the flex modes
    thetafs = []
    # Get Thetaf
    for row in etas:
        theta_f = Flex_Model.get_thetaf(row)
        thetafs.append(theta_f.copy())

    thetafs = np.asarray(thetafs)
    plt.figure(figsize=(12, 8))
    plt.plot(times, thetafs)
    plt.legend(['theta_1', 'theta_2', 'theta_3'])
    plt.title('Thetaf')
    plt.grid(True)   

    plt.tight_layout()

    plt.show()


            
        

    





