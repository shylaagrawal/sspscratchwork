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

def ephemeris(a, e, i, Om, om, M, jd_ref, jd_target, sun_vec_eq):
    a_AU = a / AU

    # mean anomaly to target date
    n = GAUSSIAN_K / (a_AU**1.5)
    M_target = M + np.degrees(n * (jd_target - jd_ref))

    r_ecliptic = elements_to_position(a, e, i, Om, om, M_target)
    # ecliptic coordinates to equatorial coordinates
    epsilon = np.radians(23.43928)

    Rx = np.array([[1, 0, 0],
                   [0, np.cos(epsilon), -np.sin(epsilon)],
                   [0, np.sin(epsilon), np.cos(epsilon)]])

    r_equatorial = Rx @ r_ecliptic
    sun_vec_eq = sun_vec_eq / AU
    rho = sun_vec_eq + r_equatorial
    rho_mag = np.linalg.norm(rho)

    # RA
    RA = np.arctan2(rho[1], rho[0])
    if RA < 0:
        RA += 2*np.pi

    # Declination
    Dec = np.arcsin(rho[2] / rho_mag)

    return np.degrees(RA), np.degrees(Dec)

def solve_kepler(M_deg, e, tol=1e-12):
    # initial guess
    M = np.radians(M_deg)
    E = M

    # newton-Raphson iteration
    while True:
        E_new = E - (E - e*np.sin(E) - M)/(1 - e*np.cos(E))

        if abs(E_new - E) < tol:
            break

        E = E_new
    
    return E_new

def elements_to_position(a_km, e, i_deg, Om_deg, om_deg, M):
    
    E = solve_kepler(M, e)

    # conversions
    a = a_km / AU
    i = np.radians(i_deg)
    Om = np.radians(Om_deg)
    om = np.radians(om_deg)

    # position in orbital plane
    x_orbit = a * (np.cos(E) - e)
    y_orbit = a * np.sqrt(1 - e**2) * np.sin(E)
    r_orbit = np.array([x_orbit, y_orbit, 0])

    # rotations
    Rz_Om = np.array([[np.cos(Om), -np.sin(Om), 0],
                      [np.sin(Om), np.cos(Om), 0],
                      [0, 0, 1]])

    Rx_i = np.array([[1, 0, 0],
                    [0, np.cos(i), -np.sin(i)],
                    [0, np.sin(i), np.cos(i)]])

    Rz_om = np.array([[np.cos(om), -np.sin(om), 0],
                      [np.sin(om), np.cos(om), 0],
                      [0, 0, 1]])

    return Rz_Om @ Rx_i @ Rz_om @ r_orbit