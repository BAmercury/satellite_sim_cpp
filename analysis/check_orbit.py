import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
mu_moon = 4902.8001  # Gravitational parameter of the Moon [km^3/s^2]

# Initial conditions
# Position and velocity vectors in MCI frame (assuming given periapsis and apoapsis)
r0 = np.array([843.62154762, 1259.48406671, 965.55158334])  # [km]
v0 = np.array([1.06001641, 0.60447933, -1.71465181])        # [km/s]

# Dynamics function: computes [dr/dt, dv/dt]
def orbital_dynamics(t, y):
    r = y[:3]
    v = y[3:]
    r_norm = np.linalg.norm(r)
    drdt = v
    dvdt = -mu_moon * r / r_norm**3
    return np.concatenate((drdt, dvdt))

# Time span for the integration (one orbital period)
orbital_period = 8.2 * 3600  # [s] (8.2 hours)
t_span = (0, orbital_period)
t_eval = np.linspace(0, orbital_period, 1000)  # 1000 points for smooth orbit

# Initial state vector
y0 = np.concatenate((r0, v0))

# Integrate the equations of motion
sol = solve_ivp(orbital_dynamics, t_span, y0, method='RK45', t_eval=t_eval)

# Extract position and velocity
positions = sol.y[:3].T
velocities = sol.y[3:].T
times = sol.t

# Plot the orbit
plt.figure(figsize=(10, 5))

# Plot position
plt.subplot(1, 2, 1)
plt.plot(positions[:, 0], positions[:, 1], label='Orbit path')
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('Orbital Path in MCI Frame')
plt.axis('equal')
plt.grid()

# Plot velocity magnitude over time
vel_magnitude = np.linalg.norm(velocities, axis=1)
plt.subplot(1, 2, 2)
plt.plot(times / 3600, vel_magnitude)  # Convert time to hours
plt.xlabel('Time [hours]')
plt.ylabel('Velocity Magnitude [km/s]')
plt.title('Velocity Magnitude over Time')
plt.grid()

plt.figure(figsize=(10, 5))
plt.plot(times / 3600, velocities)
plt.ylabel('S/C Velocity (m/s)')
plt.legend(['X', 'Y', 'Z'])
plt.grid(True)
plt.xlabel('Time (hours)')


# Set up 3D plot
# Constants
mu_moon = 4902.8001  # Gravitational parameter of the Moon [km^3/s^2]
moon_radius = 1737.4  # Radius of the Moon [km]



fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([1, 1, 1])

# Plot the Moon as a sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = moon_radius * np.outer(np.cos(u), np.sin(v))
y = moon_radius * np.outer(np.sin(u), np.sin(v))
z = moon_radius * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='gray', alpha=0.5, rstride=5, cstride=5, edgecolor='k', linewidth=0.1)

# Plot the spacecraft's orbit
ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], 'b', label='Spacecraft Orbit')

# Labels and plot settings
ax.set_xlabel('X [km]')
ax.set_ylabel('Y [km]')
ax.set_zlabel('Z [km]')
ax.set_title('3D Orbit around the Moon')
ax.legend()

plt.tight_layout()
plt.show()