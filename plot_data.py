import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def load_and_plot_spacecraft_data(filename='spacecraft_data.csv'):
    # Read the data
    data = pd.read_csv(filename)
    plt.style.use('seaborn')
    
    # 1. Position and Velocity Plot (as subplots)
    fig1, (ax1_pos, ax1_vel) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Position subplot
    ax1_pos.plot(data['time'], data['pos_x'], label='X')
    ax1_pos.plot(data['time'], data['pos_y'], label='Y')
    ax1_pos.plot(data['time'], data['pos_z'], label='Z')
    ax1_pos.set_xlabel('Time (s)')
    ax1_pos.set_ylabel('Position (km)')
    ax1_pos.grid(True)
    ax1_pos.legend()
    ax1_pos.set_title('Position vs Time')
    
    # Velocity subplot
    ax1_vel.plot(data['time'], data['vel_x'], label='X')
    ax1_vel.plot(data['time'], data['vel_y'], label='Y')
    ax1_vel.plot(data['time'], data['vel_z'], label='Z')
    ax1_vel.set_xlabel('Time (s)')
    ax1_vel.set_ylabel('Velocity (km/s)')
    ax1_vel.grid(True)
    ax1_vel.legend()
    ax1_vel.set_title('Velocity vs Time')
    
    plt.tight_layout()

    # 2. Angular Velocity and Control Torques (as subplots)
    fig2, (ax2_omega, ax2_torque) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Angular Velocity subplot
    ax2_omega.plot(data['time'], np.rad2deg(data['omega_x']), label='X')
    ax2_omega.plot(data['time'], np.rad2deg(data['omega_y']), label='Y')
    ax2_omega.plot(data['time'], np.rad2deg(data['omega_z']), label='Z')
    ax2_omega.set_xlabel('Time (s)')
    ax2_omega.set_ylabel('Angular Velocity (deg/s)')
    ax2_omega.grid(True)
    ax2_omega.legend()
    ax2_omega.set_title('Angular Velocity vs Time')
    
    # Control Torques subplot
    ax2_torque.plot(data['time'], data['torque_x'], label='X')
    ax2_torque.plot(data['time'], data['torque_y'], label='Y')
    ax2_torque.plot(data['time'], data['torque_z'], label='Z')
    ax2_torque.set_xlabel('Time (s)')
    ax2_torque.set_ylabel('Torque (N⋅m)')
    ax2_torque.grid(True)
    ax2_torque.legend()
    ax2_torque.set_title('Control Torques vs Time')
    
    plt.tight_layout()

    # 3. 3D Trajectory Plot
    fig3 = plt.figure(figsize=(10, 10))
    ax3 = fig3.add_subplot(111, projection='3d')
    ax3.plot3D(data['pos_x'], data['pos_y'], data['pos_z'])
    ax3.set_xlabel('X (km)')
    ax3.set_ylabel('Y (km)')
    ax3.set_zlabel('Z (km)')
    ax3.set_title('Spacecraft Trajectory')
    
    # Make the 3D plot aspect ratio equal
    max_range = np.array([data['pos_x'].max()-data['pos_x'].min(),
                         data['pos_y'].max()-data['pos_y'].min(),
                         data['pos_z'].max()-data['pos_z'].min()]).max() / 2.0
    mid_x = (data['pos_x'].max()+data['pos_x'].min()) * 0.5
    mid_y = (data['pos_y'].max()+data['pos_y'].min()) * 0.5
    mid_z = (data['pos_z'].max()+data['pos_z'].min()) * 0.5
    ax3.set_xlim(mid_x - max_range, mid_x + max_range)
    ax3.set_ylim(mid_y - max_range, mid_y + max_range)
    ax3.set_zlim(mid_z - max_range, mid_z + max_range)

    # 4. Quaternion Components
    fig4 = plt.figure(figsize=(12, 8))
    ax4 = fig4.add_subplot(111)
    ax4.plot(data['time'], data['quat_w'], label='w')
    ax4.plot(data['time'], data['quat_x'], label='x')
    ax4.plot(data['time'], data['quat_y'], label='y')
    ax4.plot(data['time'], data['quat_z'], label='z')
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('Quaternion Components')
    ax4.set_title('Quaternion Evolution vs Time')
    ax4.legend()
    ax4.grid(True)

    # 5. Quaternion Norm
    fig5 = plt.figure(figsize=(12, 8))
    ax5 = fig5.add_subplot(111)
    quat_norm = np.sqrt(data['quat_w']**2 + data['quat_x']**2 + 
                       data['quat_y']**2 + data['quat_z']**2)
    ax5.plot(data['time'], quat_norm)
    ax5.set_xlabel('Time (s)')
    ax5.set_ylabel('Quaternion Norm')
    ax5.set_title('Quaternion Norm vs Time')
    ax5.grid(True)

    # Print statistics
    print("\nSimulation Statistics:")
    print("-" * 50)
    print(f"Simulation duration: {data['time'].max():.2f} seconds")
    print(f"Maximum angular velocity (deg/s): {np.rad2deg(data[['omega_x', 'omega_y', 'omega_z']].abs().max()).max():.2f}")
    print(f"Maximum control torque (N⋅m): {data[['torque_x', 'torque_y', 'torque_z']].abs().max().max():.2f}")
    print(f"Maximum position (km): {data[['pos_x', 'pos_y', 'pos_z']].abs().max().max():.2f}")
    print(f"Maximum velocity (km/s): {data[['vel_x', 'vel_y', 'vel_z']].abs().max().max():.2f}")
    
    # Calculate orbital period if possible
    if data['time'].max() > 3600:  # If simulation runs for more than an hour
        pos_x = data['pos_x']
        zero_crossings = np.where(np.diff(np.signbit(pos_x)))[0]
        if len(zero_crossings) >= 2:
            period = 2 * np.mean(np.diff(data['time'][zero_crossings]))
            print(f"Estimated orbital period: {period/3600:.2f} hours")

    # Show all plots at once
    plt.show()

if __name__ == "__main__":
    load_and_plot_spacecraft_data()