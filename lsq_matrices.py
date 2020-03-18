import numpy as np
from rotations import M_rot

def A_and_w(omega, phi, kappa, bx, by, bz, left, right, c):
    # Update rotation matrix with current estimates
    M = M_rot(omega, phi, kappa)
    # Form A and w matrices
    for i in range(left.shape[0]):
        # Left image coords
        xL = left[i,1]
        yL = left[i,2]
        # Right image coords; update first with current rotation
        right_prime = M.T.dot(np.array([[right[i,1]],
                                        [right[i,2]],
                                        [-c]]))
        xR_prime = right_prime[0]
        yR_prime = right_prime[1]
        zR_prime = right_prime[2]
        # Variables for determinants
        A = 
        B = 
        C = 
        D = 
        E = 
        F = 

