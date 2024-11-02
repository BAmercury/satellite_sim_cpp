import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
mu_moon = 4902.8001  # Gravitational parameter of the Moon [km^3/s^2]
moon_radius = 1737.4  # Radius of the Moon [km]

# Initial conditions (assuming periapsis and apoapsis positions)
r0 = np.array([843.62154762, 1259.48406671, 965.55158334])  # [km]
v0 = np.array([1.06001641, 0.60447933, -1.71465181])        # [km/s]

# Dynamics function for the ODE
def orbital_dynamics(t, y):
    r = y[:3]
    v = y[3:]
    r_norm = np.linalg.norm(r)
    drdt = v
    dvdt = -mu_moon * r / r_norm**3
    return np.concatenate((drdt, dvdt))

# Time span for integration (one orbital period)
orbital_period = 8.2 * 3600  # [s]
t_span = (0, orbital_period)
t_eval = np.linspace(0, orbital_period, 1000)

# Initial state vector
y0 = np.concatenate((r0, v0))

# Solve the equations of motion
sol = solve_ivp(orbital_dynamics, t_span, y0, method='RK45', t_eval=t_eval)
positions = sol.y[:3].T

# Set up 3D plot
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
plt.show()