from astropy import units as u
import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit

inclination = 27.7116 << u.deg # [deg] inclination
peri_alt = 185 << u.km # [km] Periapsis Alt
apo_alt = 40000 << u.km # [km] Apoapsis Alt
raan = 351.7394 << u.deg #[deg] Long of Ascending Node 
arg_peri = 0.0 << u.deg # [deg] Argument of Periapsis
nu = 240.0 << u.deg # [deg] True anomaly
earth_mu = 3.986e5 << (u.km**3)/(u.s**2) # [km^3/s^2]
#moon_mu = 1.62 << (u.m / u.s**2) # [m/s^2]
earth_radius = 6371 << u.km # [km] 

Ra = apo_alt + earth_radius # Radius at Apoapsis
Rp = peri_alt + earth_radius # Radius at Periapsis
a = (Ra + Rp) / 2 << u.km # Semi-major axis [km]

ecc = (Ra - Rp) / (Ra + Rp) << u.one

orb = Orbit.from_classical(Earth, a, ecc, inclination, raan, arg_peri, nu)


print(orb.period.to(u.second))

print(386303.56111470016 + 386303.56111470016) # Two orbits

