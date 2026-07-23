# Shyla Agrawal
# 7/7 - 7/21
# Orbit Determination Library

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# constants
AU = 149597870.7     # kilometers in one AU
DAY = 86400          # seconds in one day
GAUSSIAN_K = 0.0172020989484
C_AU = 173.14463267  # Speed of light in AU/day

# helper functions

def AUday_to_AUGday(r_vec, v_vec):
    return np.array(r_vec), np.array(v_vec) / GAUSSIAN_K

def hms_to_deg(hours, minutes, seconds):
    return (hours + minutes/60.0 + seconds/3600.0) * 15.0

def dms_to_deg(is_pos, degrees, minutes, seconds):
    res = degrees + minutes/60.0 + seconds/3600.0
    return res if is_pos else -res

def angular_momentum(r_vec, v_vec):
    r_vec, v_vec = AUday_to_AUGday(r_vec, v_vec)
    return np.cross(r_vec, v_vec)

def ecliptic_to_equatorial(vector):
    eps = np.radians(23.4374)
    
    rotation = np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(eps), -np.sin(eps)],
        [0.0, np.sin(eps), np.cos(eps)]
    ])
    
    return rotation @ vector

def equatorial_to_ecliptic(vector):
    eps = np.radians(23.4374)
    rotation = np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(eps), np.sin(eps)],
        [0.0, -np.sin(eps), np.cos(eps)]
    ])
    return rotation @ vector

def radec_to_decimal_to_rhohat(ra_str, dec_str):
    h, m, s = map(float, str(ra_str).strip().split())
    ra = np.radians(hms_to_deg(h, m, s))

    dec_str = str(dec_str).strip()
    sign = -1.0 if dec_str.startswith('-') else 1.0
    parts = dec_str.lstrip('+-').split()
    d, m, s = map(float, parts)
    dec = np.radians(sign * (d + m / 60.0 + s / 3600.0))

    return np.array([
        np.cos(ra) * np.cos(dec),
        np.sin(ra) * np.cos(dec),
        np.sin(dec)
    ])

def radec_to_decimal(ra_str, dec_str):
    h, m, s = map(float, str(ra_str).strip().split())
    ra = np.radians(hms_to_deg(h, m, s))

    dec_str = str(dec_str).strip()
    sign = -1.0 if dec_str.startswith('-') else 1.0
    parts = dec_str.lstrip('+-').split()
    d, m, s = map(float, parts)
    dec = np.radians(sign * (d + m / 60.0 + s / 3600.0))

    return ra, dec

def decimal_to_radec(ra_deg, dec_deg):
    # Ensure RA is positive in [0, 360)
    ra_deg = ra_deg % 360.0
    
    # Convert RA degrees to hours, minutes, seconds
    ra_hours_tot = ra_deg / 15.0
    h = int(ra_hours_tot)
    m = int((ra_hours_tot - h) * 60.0)
    s = ((ra_hours_tot - h) * 60.0 - m) * 60.0
    ra_str = f"{h:02d} {m:02d} {s:06.3f}"
    
    # Determine sign and convert Dec degrees to degrees, minutes, seconds
    sign = "-" if dec_deg < 0 else "+"
    abs_dec = abs(dec_deg)
    d = int(abs_dec)
    dm = int((abs_dec - d) * 60.0)
    ds = ((abs_dec - d) * 60.0 - dm) * 60.0
    dec_str = f"{sign}{d:02d} {dm:02d} {ds:05.2f}"
    
    return ra_str, dec_str

def determinent(a, b, c):
    return np.dot(np.cross(a, b), c)

# main methods

def orbital_elements(r_vec, v_vec):
    h_vec = angular_momentum(r_vec, v_vec)
    r_vec, v_vec = AUday_to_AUGday(r_vec, v_vec)
    mu = 1.0

    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    h = np.linalg.norm(h_vec)

    a = 1.0 / ((2.0/r) - v**2)
    e = np.sqrt(max(0.0, 1.0 - h**2/(mu*a)))
    i = np.arccos(np.clip(h_vec[2]/h, -1.0, 1.0))
    
    Omega = np.arctan2(h_vec[0], -h_vec[1])
    if Omega < 0:
        Omega += 2.0 * np.pi

    cos_U = (r_vec[0]*np.cos(Omega) + r_vec[1]*np.sin(Omega))/r
    sin_U = r_vec[2]/(r*np.sin(i)) if np.sin(i) != 0 else 0.0
    U = np.arctan2(sin_U, cos_U)
    if U < 0:
        U += 2.0 * np.pi

    cos_nu = np.clip((a*(1.0 - e**2)/r - 1.0)/e if e != 0 else 1.0, -1.0, 1.0)
    r_dot_v = np.dot(r_vec, v_vec)
    sin_nu = (a*(1.0 - e**2)/(e*h)) * (r_dot_v/r) if (e*h) != 0 else 0.0
    nu = np.arctan2(sin_nu, cos_nu)
    if nu < 0:
        nu += 2.0 * np.pi

    omega = U - nu
    if omega < 0:
        omega += 2.0 * np.pi

    cos_E = np.clip((1.0/e)*(1.0 - r/a) if e != 0 else 1.0, -1.0, 1.0)
    E = np.arccos(cos_E)
    if np.pi < nu < 2.0 * np.pi:
        E = 2.0 * np.pi - E

    M = E - e*np.sin(E)
    if M < 0:
        M += 2.0 * np.pi

    return np.array([a, e, np.degrees(i), np.degrees(Omega), np.degrees(omega), np.degrees(M)])

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

def f_and_g(r_vec, v_vec, tau, order_four):
    r = np.linalg.norm(r_vec)
    u = 1.0 / r**3
    z = np.dot(r_vec, v_vec) / r**2
    q = np.dot(v_vec, v_vec) / r**2 - u

    if order_four:
        f = (1.0 - 0.5*u*tau**2 + 0.5*u*z*tau**3 + (1.0/24.0)*(3.0*u*q - 15.0*u*z**2 + u**2)*tau**4)
        g = (tau - (1.0/6.0)*u*tau**3 + 0.25*u*z*tau**4)
    else:
        f = (1.0 - 0.5*u*tau**2 + 0.5*u*z*tau**3)
        g = (tau - (1.0/6.0)*u*tau**3)

    return f, g

def solve_rhos(a1, a3, rhohat1, rhohat2, rhohat3, R1, R2, R3):
    num1 = a1 * determinent(R1, rhohat2, rhohat3) - determinent(R2, rhohat2, rhohat3) + a3 * determinent(R3, rhohat2, rhohat3)
    den1 = a1 * determinent(rhohat1, rhohat2, rhohat3)
    rho1 = num1 / den1

    num2 = a1 * determinent(rhohat1, R1, rhohat3) - determinent(rhohat1, R2, rhohat3) + a3 * determinent(rhohat1, R3, rhohat3)
    den2 = -determinent(rhohat1, rhohat2, rhohat3)
    rho2 = num2 / den2

    num3 = a1 * determinent(rhohat2, R1, rhohat1) - determinent(rhohat2, R2, rhohat1) + a3 * determinent(rhohat2, R3, rhohat1)
    den3 = a3 * determinent(rhohat2, rhohat3, rhohat1)
    rho3 = num3 / den3

    return rho1, rho2, rho3

def gauss_method(t1, t2, t3, ra1, dec1, ra2, dec2, ra3, dec3, R1, R2, R3, order_four=True, tolerance=1e-10, max_iter=100):
    rhohat1 = radec_to_decimal_to_rhohat(ra1, dec1)
    rhohat2 = radec_to_decimal_to_rhohat(ra2, dec2)
    rhohat3 = radec_to_decimal_to_rhohat(ra3, dec3)

    tau1 = GAUSSIAN_K * (t1 - t2)
    tau0 = GAUSSIAN_K * (t3 - t1)
    tau3 = GAUSSIAN_K * (t3 - t2)

    a1 = tau3 / tau0
    a3 = -tau1 / tau0

    rho1, rho2, rho3 = solve_rhos(a1, a3, rhohat1, rhohat2, rhohat3, R1, R2, R3)

    r1_init = rho1 * rhohat1 - R1
    r2_init = rho2 * rhohat2 - R2
    r3_init = rho3 * rhohat3 - R3

    v12 = (r2_init - r1_init) / (t2 - t1)
    v23 = (r3_init - r2_init) / (t3 - t2)

    v2_init = ((t3 - t2) * v12 + (t2 - t1) * v23) / (t3 - t1) / GAUSSIAN_K

    r2_refined = r2_init
    v2_refined = v2_init
    r2_prev = r2_init

    for _ in range(max_iter):
        t1_corr = t1 - (rho1 / C_AU)
        t2_corr = t2 - (rho2 / C_AU)
        t3_corr = t3 - (rho3 / C_AU)

        tau1 = GAUSSIAN_K * (t1_corr - t2_corr)
        tau3 = GAUSSIAN_K * (t3_corr - t2_corr)

        f1, g1 = f_and_g(r2_refined, v2_refined, tau1, order_four)
        f3, g3 = f_and_g(r2_refined, v2_refined, tau3, order_four)

        denom_fg = f1 * g3 - f3 * g1
        a1 = g3 / denom_fg
        a3 = -g1 / denom_fg

        rho1, rho2, rho3 = solve_rhos(a1, a3, rhohat1, rhohat2, rhohat3, R1, R2, R3)

        r1 = rho1 * rhohat1 - R1
        r3 = rho3 * rhohat3 - R3

        r2_refined = (g3 * r1 - g1 * r3) / denom_fg
        v2_refined = (f3 * r1 - f1 * r3) / (f3 * g1 - f1 * g3)

        mag_r2 = np.linalg.norm(r2_refined)
        diff = np.linalg.norm(r2_refined - r2_prev) / mag_r2

        if diff < tolerance:
            break

        r2_prev = r2_refined

    v2_refined_au = v2_refined * GAUSSIAN_K
    return r2_refined, v2_refined_au

def compute_coefficients(filename):
    star_positions = pd.read_csv(filename)

    # column vectors
    x = star_positions["x_pix"]
    y = star_positions["y_pix"]
    ra = star_positions["ra_deg"]
    dec = star_positions["dec_deg"]

    # compute the xy matrix
    N = len(star_positions)
    xy_matrix = np.array([
        [N, np.sum(x), np.sum(y)],
        [np.sum(x), np.sum(x**2), np.sum(x*y)],
        [np.sum(y), np.sum(x*y), np.sum(y**2)]
    ])

    # compute LHS vector for RA and Dec
    lhs_ra = np.array([np.sum(ra), np.sum(x * ra), np.sum(y * ra)])
    lhs_dec = np.array([np.sum(dec), np.sum(x * dec), np.sum(y * dec)])

    # compute the plate solution coefficients 
    b1, a11, a12 = np.linalg.solve(xy_matrix, lhs_ra)
    b2, a21, a22 = np.linalg.solve(xy_matrix, lhs_dec)

    return b1, a11, a12, b2, a21, a22

def xy_to_radec(x_pix, y_pix, b1, a11, a12, b2, a21, a22):
    # Converts a location in an image (x_pix, y_pix) to an RA and DEC coordinate in degrees
    # Uses the 6 plate solve coefficients to do so
    
    return (b1 + a11 * x_pix + a12 * y_pix, b2 + a21 * x_pix + a22 * y_pix)

def compute_sigmas(coefficients, star_df):
    b1, a11, a12, b2, a21, a22 = coefficients
    
    x = star_df["x_pix"]
    y = star_df["y_pix"]
    ra_cat = star_df["ra_deg"]
    dec_cat = star_df["dec_deg"]
    
    N = len(star_df)
    
    ra_fit = b1 + a11 * x + a12 * y
    ra_diff = ra_fit - ra_cat
    
    dec_fit = b2 + a21 * x + a22 * y
    dec_diff = dec_fit - dec_cat
    
    sigma_ra = np.sqrt(np.sum(ra_diff**2) / (N - 3))
    sigma_dec = np.sqrt(np.sum(dec_diff**2) / (N - 3))
    
    return sigma_ra, sigma_dec
