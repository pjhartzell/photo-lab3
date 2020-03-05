import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from tabulate import tabulate


# Load fiducial coordinates ----------------------------------------------------
pxy_fid_27 = np.loadtxt("image_27_fid.txt")
pxy_fid_28 = np.loadtxt("image_28_fid.txt")


# Principal Point Correction ---------------------------------------------------
def pp_offset_corr(pxy, x_p, y_p, print_flag=False):
    p = pxy[:,0]
    x = pxy[:,1]
    y = pxy[:,2]
    x_p_vec = x_p * np.ones((1, p.shape[0]))
    y_p_vec = y_p * np.ones((1, p.shape[0]))
    x_bar = x - x_p_vec
    y_bar = y - y_p_vec

    if print_flag:
        result_table = np.hstack((pxy,
                                  x_p_vec.reshape(-1,1),
                                  y_p_vec.reshape(-1,1),
                                  x_bar.reshape(-1,1),
                                  y_bar.reshape(-1,1)))
        headers = ["#", "x", "y", "x_correction", "y_correction",
                   "x_bar", "y_bar"]
        print("\nPRINCIPAL POINT OFFSET CORRECTION (mm)")
        print(tabulate(result_table, headers, tablefmt="simple",
            floatfmt=(".0f", ".4f", ".4f", ".4f", ".4f", ".4f", ".4f")))

    return np.vstack((p, x_bar, y_bar)).T

# Calibrated principal point (mm)
x_p = -0.006
y_p = 0.006

# Image 27
pxy_bar_27 = pp_offset_corr(pxy_fid_27, x_p, y_p, print_flag=False)
# Image 28
pxy_bar_28 = pp_offset_corr(pxy_fid_28, x_p, y_p, print_flag=False)


# Radial Lens Distortion Correction --------------------------------------------
def balanced_radial_corr(pxy, k0, k1, k2, k3, print_flag=False):
    p = pxy[:,0]
    x_bar = pxy[:,1]
    y_bar = pxy[:,2]
    r = np.sqrt(x_bar**2 + y_bar**2)
    delta_x_radial = -x_bar * (k0 + k1*r**2 + k2*r**4 + k3*r**6)
    delta_y_radial = -y_bar * (k0 + k1*r**2 + k2*r**4 + k3*r**6)

    if print_flag:
        result_table = np.hstack((pxy,
                                  delta_x_radial.reshape(-1,1),
                                  delta_y_radial.reshape(-1,1)))
        headers = ["#", "x_bar", "y_bar", "delta_x_radial", "delta_y_radial"]
        print("\nRADIAL LENS DISTORTION CORRECTION (mm)")
        print(tabulate(result_table, headers, tablefmt="simple",
            floatfmt=(".0f", ".4f", ".4f", ".4f", ".4f")))

    return np.vstack((p, delta_x_radial, delta_y_radial)).T

# Radial lens distortion coefficients
k0 = 0.8878e-4
k1 = -0.1528e-7
k2 = 0.5256e-12
k3 = 0

# Image 27
pxy_radial_27 = balanced_radial_corr(pxy_bar_27, k0, k1, k2, k3,
                                     print_flag=False)
# Image 28
pxy_radial_28 = balanced_radial_corr(pxy_bar_28, k0, k1, k2, k3,
                                     print_flag=False)


# Decentering Lens Distortion Correction ---------------------------------------
def decentering_corr(pxy, p1, p2, print_flag=False):
    p = pxy[:,0]
    x_bar = pxy[:,1]
    y_bar = pxy[:,2]
    r = np.sqrt(x_bar**2 + y_bar**2)
    delta_x_dec = -(p1*(r**2 + 2*x_bar**2) + 2*p2*x_bar*y_bar)
    delta_y_dec = -(p2*(r**2 + 2*y_bar**2) + 2*p1*x_bar*y_bar)

    if print_flag:
        result_table = np.hstack((pxy,
                                  delta_x_dec.reshape(-1,1),
                                  delta_y_dec.reshape(-1,1)))
        headers = ["#", "x_bar", "y_bar",
                   "delta_x_decentering", "delta_y_decentering"]
        print("\nDECENTERING LENS DISTORTION CORRECTION (mm)")
        print(tabulate(result_table, headers, tablefmt="simple",
            floatfmt=(".0f", ".4f", ".4f", ".4f", ".4f")))

    return np.vstack((p, delta_x_dec, delta_y_dec)).T

# Decentering lens distortion coefficients
p1 = 0.1346e-6
p2 = 0.1224e-7

# Image 27
pxy_decentering_27 = decentering_corr(pxy_bar_27, p1, p2, print_flag=False)
# Image 28
pxy_decentering_28 = decentering_corr(pxy_bar_28, p1, p2, print_flag=False)


# Atmospheric Refraction Correction --------------------------------------------
def atmospheric_corr(pxy, h, H, c, print_flag=False):
    p = pxy[:,0]
    x_bar = pxy[:,1]
    y_bar = pxy[:,2]
    r = np.sqrt(x_bar**2 + y_bar**2)
    K = ((2410*H)/(H**2-6*H+250) - ((2410*h)/(h**2-6*h+250))*(h/H)) / 1000000
    delta_x_atm = -x_bar*K*(1+(r**2)/(c**2))
    delta_y_atm = -y_bar*K*(1+(r**2)/(c**2))

    if print_flag:
        result_table = np.hstack((pxy,
                                  delta_x_atm.reshape(-1,1),
                                  delta_y_atm.reshape(-1,1)))
        headers = ["#", "x_bar", "y_bar",
                   "delta_x_atmosphere", "delta_y_atmosphere"]
        print("\nATMOSPHERIC LENS DISTORTION CORRECTION (mm)")
        print(tabulate(result_table, headers, tablefmt="simple",
            floatfmt=(".0f", ".4f", ".4f", ".4f", ".4f")))

    return np.vstack((p, delta_x_atm, delta_y_atm)).T

# Parameters for atmospheric corrections
H = 1.860        # km
h = 1.100        # km
c = 153.358     # mm

# Image 27
pxy_atmosphere_27 = atmospheric_corr(pxy_bar_27, h, H, c, print_flag=False)
# Image 28
pxy_atmosphere_28 = atmospheric_corr(pxy_bar_28, h, H, c, print_flag=False)


# Final Refined Coordinates ----------------------------------------------------
def refined_coords(pxy, pxy_bar, pxy_radial, pxy_decenter, pxy_atm, 
                   print_flag=False):
    p = pxy[:,0]
    x = pxy[:,1]
    y = pxy[:,2]
    xy_refined = (pxy_bar[:,1:]
                 + pxy_radial[:,1:]
                 + pxy_decenter[:,1:]
                 + pxy_atm[:,1:])

    if print_flag:
        result_table = np.hstack((pxy, xy_refined))
        headers = ["#", "x_original", "y_original",
                   "x_refined", "y_refined"]
        print("\nFINAL REFINED COORDINATES (mm)")
        print(tabulate(result_table, headers, tablefmt="simple",
            floatfmt=(".0f", ".4f", ".4f", ".4f", ".4f")))

    return np.hstack((p.reshape(-1,1), xy_refined))

# Image 27
pxy_prime_27 = refined_coords(pxy_fid_27, pxy_bar_27, pxy_radial_27,
                              pxy_decentering_27, pxy_atmosphere_27,
                              print_flag=False)
# Image 28
pxy_prime_28 = refined_coords(pxy_fid_28, pxy_bar_28, pxy_radial_28,
                              pxy_decentering_28, pxy_atmosphere_28,
                              print_flag=False)


# Save Corrected Coords --------------------------------------------------------
np.savetxt("image_27_corrected.txt", pxy_prime_27, fmt='%0.5f')
np.savetxt("image_28_corrected.txt", pxy_prime_28, fmt='%0.5f')


