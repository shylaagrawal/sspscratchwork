# Shyla Agrawal
# 7/7 - 
# Orbit Determination Library

import numpy as np

# constants
AU = 149597870.7    # how many kilometers in one AU
DAY = 86400         # how many seconds in one day
GAUSSIAN_K = 0.0172020989484

def km_per_sec_to_AU_per_Gday(r_vec, v_vec):
    # conversions
    r_vec = r_vec / AU
    v_vec = v_vec * DAY / AU
    v_vec /= GAUSSIAN_K
    return r_vec, v_vec

def hms_to_deg(hours, minutes, seconds):
    # Converts an angle in (hours, minutes, seconds) to decimal degrees
    return (hours + minutes/60 + seconds/3600) * (360/24)

def dms_to_deg(is_pos, degrees, minutes, seconds):
    # Converts an angle in (degrees, arcminutes, arcseconds) to decimal degrees
    res = degrees + minutes/60 + seconds/3600

    if (is_pos):
        return res
    else:
        return -res

def angular_momentum(r_vec, v_vec):
    r, v = km_per_sec_to_AU_per_Gday(r_vec, v_vec)
    return np.cross(r, v)

def orbital_elements(r_vec, v_vec):
    h_vec = angular_momentum(r_vec, v_vec)
    r_vec, v_vec = km_per_sec_to_AU_per_Gday(r_vec, v_vec)
    mu = 1.0

    # magnitudes
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    h = np.linalg.norm(h_vec)

    # easy orbital elements
    a = 1 / ( (2/r) - v**2)
    e = np.sqrt(1 - h**2/(mu*a))
    i = np.arccos(h_vec[2]/h)
    Omega = np.arctan2(h_vec[0], -h_vec[1])
    if Omega < 0:
        Omega += 2*np.pi

    # U (angular distance from the ascending node to the asteroid)
    cos_U = (r_vec[0]*np.cos(Omega) + r_vec[1]*np.sin(Omega))/r
    sin_U = r_vec[2]/(r*np.sin(i))
    U = np.arctan2(sin_U, cos_U)
    if U < 0:
        U += 2*np.pi

    # true anomaly nu
    cos_nu = (a*(1-e**2)/r - 1)/e
    r_dot_v = np.dot(r_vec, v_vec)
    sin_nu = (a*(1-e**2)/(e*h)) * (r_dot_v/r)
    nu = np.arctan2(sin_nu, cos_nu)
    if nu < 0:
        nu += 2*np.pi

    # argument of perihelion
    omega = U - nu
    if omega < 0:
        omega += 2*np.pi

    # eccentric anomaly
    cos_E = (1/e)*(1-r/a)
    E = np.arccos(cos_E)
    if np.pi < nu < 2*np.pi:
        E = 2*np.pi - E

    # mean anomaly
    M = E - e*np.sin(E)
    if M < 0:
        M += 2*np.pi

    return np.array([a*AU, e, np.degrees(i), np.degrees(Omega), np.degrees(omega), np.degrees(M)])