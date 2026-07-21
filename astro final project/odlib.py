# Shyla Agrawal
# 7/7 - 
# Orbit Determination Library

import numpy as np

# constants
AU = 149597870.7    # how many kilometers in one AU
DAY = 86400         # how many seconds in one day
GAUSSIAN_K = 0.0172020989484

def AUday_to_AUGday(r_vec, v_vec):
    return np.array(r_vec), np.array(v_vec) / GAUSSIAN_K

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
    r_vec, v_vec = AUday_to_AUGday(r_vec, v_vec)
    return np.cross(r_vec, v_vec)

def orbital_elements(r_vec, v_vec):
    h_vec = angular_momentum(r_vec, v_vec)
    r_vec, v_vec = AUday_to_AUGday(r_vec, v_vec)
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

    return np.array([a, e, np.degrees(i), np.degrees(Omega), np.degrees(omega), np.degrees(M)])

def ephemeris(a, e, i, Om, om, M, jd_ref, jd_target, sun_vec_eq):

    # mean anomaly to target date
    n = GAUSSIAN_K / (a**1.5)
    M_target = M + np.degrees(n * (jd_target - jd_ref))

    r_ecliptic = elements_to_position(a, e, i, Om, om, M_target)
    # ecliptic coordinates to equatorial coordinates
    epsilon = np.radians(23.43928)

    Rx = np.array([[1, 0, 0],
                   [0, np.cos(epsilon), -np.sin(epsilon)],
                   [0, np.sin(epsilon), np.cos(epsilon)]])

    r_equatorial = Rx @ r_ecliptic
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

def elements_to_position(a, e, i_deg, Om_deg, om_deg, M):
    
    E = solve_kepler(M, e)

    # conversions
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

def f_and_g(r_vec, v_vec, tau, order_four):
    r = np.linalg.norm(r_vec)
    u = 1 / r**3
    z = np.dot(r_vec, v_vec) / r**2
    q = np.dot(v_vec, v_vec) / r**2 - u

    if order_four:
        f = (1 - 0.5*u*tau**2 + 0.5*u*z*tau**3 + (1/24)*(3*u*q - 15*u*z**2 + u**2)*tau**4)
        g = (tau - (1/6)*u*tau**3 + 0.25*u*z*tau**4)

    else:
        f = (1 - 0.5*u*tau**2 + 0.5*u*z*tau**3)
        g = (tau - (1/6)*u*tau**3)

    return f, g

def radec_to_decimal_to_rhohat(ra_str, dec_str):
    # RA
    h, m, s = map(float, ra_str.split())
    ra = np.radians(hms_to_deg(h, m, s))

    # Dec
    parts = dec_str.split()

    if parts[0][0] == "-":
        sign = -1
        d = abs(float(parts[0]))
    else:
        sign = 1
        d = float(parts[0])

    m = float(parts[1])
    s = float(parts[2])

    dec = np.radians(sign * (d + m/60 + s/3600))

    return np.array([
        np.cos(ra)*np.cos(dec),
        np.sin(ra)*np.cos(dec),
        np.sin(dec)
    ])

def determinent(a, b, c):
    return np.dot(np.cross(a, b), c)

def equatorial_to_ecliptic(vector):
    eps = np.radians(23.4374)

    rotation = np.array([[1, 0, 0],
                        [0, np.cos(eps), np.sin(eps)],
                        [0, -np.sin(eps), np.cos(eps)]])

    return rotation @ vector

def gauss_method1(t1, t2, t3, ra1, dec1, ra2, dec2, ra3, dec3, R1, R2, R3, order_four=True, tolerance=1e-10, max_iter=100):
    rhohat1 = radec_to_decimal_to_rhohat(ra1,dec1)
    rhohat2 = radec_to_decimal_to_rhohat(ra2,dec2)
    rhohat3 = radec_to_decimal_to_rhohat(ra3,dec3)

    tau1 = GAUSSIAN_K*(t1-t2)
    tau0 = GAUSSIAN_K*(t3-t1)
    tau3 = GAUSSIAN_K*(t3-t2)
    a1 = tau3/tau0
    a3 = -tau1/tau0

    rho1 = (a1 * determinent(R1, rhohat2, rhohat3) - determinent(R2, rhohat2, rhohat3) + a3 * determinent(R3, rhohat2, rhohat3)) / (a1 * determinent(rhohat1, rhohat2, rhohat3))
    rho2 = (a1 * determinent(rhohat1, R1, rhohat3) - determinent(rhohat1, R2, rhohat3) + a3 * determinent(rhohat1, R3, rhohat3)) / (-1 * determinent(rhohat1, rhohat2, rhohat3))
    rho3 = (a1 * determinent(rhohat2, R1, rhohat1) - determinent(rhohat2, R2, rhohat1) + a3 * determinent(rhohat2, R3, rhohat1)) / (a3 * determinent(rhohat2, rhohat3, rhohat1))

    r1_init = rho1 * rhohat1 - R1
    r2_init = rho2 * rhohat2 - R2
    r3_init = rho3 * rhohat3 - R3

    v12 = (r2_init - r1_init) / (t2 - t1)
    v23 = (r3_init - r2_init) / (t3 - t2)

    v2_init = ((t3-t2) * v12 + (t2-t1) * v23) / (t3-t1) / GAUSSIAN_K

    f1, g1 = f_and_g(r2_init, v2_init, tau1, order_four)
    f3, g3 = f_and_g(r2_init, v2_init, tau3, order_four)

    a1 = g3 / (f1 * g3 - f3 * g1)
    a3 = -g1 / (f1 * g3 - f3 * g1)

    r2_prev = r2_init

    for iteration in range(max_iter):
        d1 = determinent(rhohat1, rhohat2, rhohat3)
        d2 = determinent(rhohat2, rhohat3, rhohat1)

        rho1 = (a1 * determinent(R1, rhohat2, rhohat3) - determinent(R2, rhohat2, rhohat3) + a3 * determinent(R3, rhohat2, rhohat3)) / (a1 * d1)
        rho2 = (a1 * determinent(rhohat1, R1, rhohat3) - determinent(rhohat1, R2, rhohat3) + a3 * determinent(rhohat1, R3, rhohat3)) / (-1 * d1)
        rho3 = (a1 * determinent(rhohat2, R1, rhohat1) - determinent(rhohat2, R2, rhohat1) + a3 * determinent(rhohat2, R3, rhohat1)) / (a3 * d2)

        r1 = rho1 * rhohat1 - R1
        r2 = rho2 * rhohat2 - R2
        r3 = rho3 * rhohat3 - R3

        denom = f1 * g3 - f3 * g1
        r2_refined = (g3 * r1 - g1 * r3) / denom
        v2_refined = (f3 * r1 - f1 * r3) / (f3 * g1 - f1 * g3)

        f1, g1 = f_and_g(r2_refined, v2_refined, tau1, order_four)
        f3, g3 = f_and_g(r2_refined, v2_refined, tau3, order_four)

        a1 = g3 / (f1 * g3 - f3 * g1)
        a3 = -g1 / (f1 * g3 - f3 * g1)

        mag_r2 = np.linalg.norm(r2_refined)
        diff = np.linalg.norm(r2_refined - r2_prev) / mag_r2

        if diff < tolerance:
            break

        r2_prev = r2_refined

    return r2_refined, v2_refined



def gauss_method(t1, t2, t3, ra1, dec1, ra2, dec2, ra3, dec3, R1, R2, R3, order_four=True, tolerance=1e-10, max_iter=100):
    c_au = 173.144643267

    rhohat1 = radec_to_decimal_to_rhohat(ra1, dec1)
    rhohat2 = radec_to_decimal_to_rhohat(ra2, dec2)
    rhohat3 = radec_to_decimal_to_rhohat(ra3, dec3)

    tau1 = GAUSSIAN_K * (t1 - t2)
    tau0 = GAUSSIAN_K * (t3 - t1)
    tau3 = GAUSSIAN_K * (t3 - t2)
    a1 = tau3 / tau0
    a3 = -tau1 / tau0

    rho1 = (a1 * determinent(R1, rhohat2, rhohat3) - determinent(R2, rhohat2, rhohat3) + a3 * determinent(R3, rhohat2, rhohat3)) / (a1 * determinent(rhohat1, rhohat2, rhohat3))
    rho2 = (a1 * determinent(rhohat1, R1, rhohat3) - determinent(rhohat1, R2, rhohat3) + a3 * determinent(rhohat1, R3, rhohat3)) / (-1 * determinent(rhohat1, rhohat2, rhohat3))
    rho3 = (a1 * determinent(rhohat2, R1, rhohat1) - determinent(rhohat2, R2, rhohat1) + a3 * determinent(rhohat2, R3, rhohat1)) / (a3 * determinent(rhohat2, rhohat3, rhohat1))

    r1_init = rho1 * rhohat1 - R1
    r2_init = rho2 * rhohat2 - R2
    r3_init = rho3 * rhohat3 - R3

    v12 = (r2_init - r1_init) / (t2 - t1)
    v23 = (r3_init - r2_init) / (t3 - t2)

    v2_init = ((t3 - t2) * v12 + (t2 - t1) * v23) / (t3 - t1) / GAUSSIAN_K

    f1, g1 = f_and_g(r2_init, v2_init, tau1, order_four)
    f3, g3 = f_and_g(r2_init, v2_init, tau3, order_four)

    a1 = g3 / (f1 * g3 - f3 * g1)
    a3 = -g1 / (f1 * g3 - f3 * g1)

    r2_prev = r2_init

    for _ in range(max_iter):
        t1_corr = t1 - (rho1 / c_au)
        t2_corr = t2 - (rho2 / c_au)
        t3_corr = t3 - (rho3 / c_au)

        tau1 = GAUSSIAN_K * (t1_corr - t2_corr)
        tau3 = GAUSSIAN_K * (t3_corr - t2_corr)

        d1 = determinent(rhohat1, rhohat2, rhohat3)
        d2 = determinent(rhohat2, rhohat3, rhohat1)

        rho1 = (a1 * determinent(R1, rhohat2, rhohat3) - determinent(R2, rhohat2, rhohat3) + a3 * determinent(R3, rhohat2, rhohat3)) / (a1 * d1)
        rho2 = (a1 * determinent(rhohat1, R1, rhohat3) - determinent(rhohat1, R2, rhohat3) + a3 * determinent(rhohat1, R3, rhohat3)) / (-1 * d1)
        rho3 = (a1 * determinent(rhohat2, R1, rhohat1) - determinent(rhohat2, R2, rhohat1) + a3 * determinent(rhohat2, R3, rhohat1)) / (a3 * d2)

        r1 = rho1 * rhohat1 - R1
        r2 = rho2 * rhohat2 - R2
        r3 = rho3 * rhohat3 - R3

        denom = f1 * g3 - f3 * g1
        r2_refined = (g3 * r1 - g1 * r3) / denom
        v2_refined = (f3 * r1 - f1 * r3) / (f3 * g1 - f1 * g3)

        f1, g1 = f_and_g(r2_refined, v2_refined, tau1, order_four)
        f3, g3 = f_and_g(r2_refined, v2_refined, tau3, order_four)
        a1 = g3 / (f1 * g3 - f3 * g1)
        a3 = -g1 / (f1 * g3 - f3 * g1)

        mag_r2 = np.linalg.norm(r2_refined)
        diff = np.linalg.norm(r2_refined - r2_prev) / mag_r2

        if diff < tolerance:
            break

        r2_prev = r2_refined

    r2_ecliptic = equatorial_to_ecliptic(r2_refined)
    v2_ecliptic = equatorial_to_ecliptic(v2_refined)

    return r2_ecliptic, v2_ecliptic

