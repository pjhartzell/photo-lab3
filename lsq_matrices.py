import numpy as np
from rotations import M_rot

def A_and_w(x_hat, left, right, bx, c):
    # Unpack x_hat
    by = x_hat[0,0]
    bz = x_hat[1,0]
    omega = x_hat[2,0]
    phi = x_hat[3,0]
    kappa = x_hat[4,0]
    # Update rotation matrix with current estimates
    M = M_rot(omega, phi, kappa)
    # Prep some variables
    num_points = left.shape[0]
    A_mat = np.zeros((num_points, 5))
    w_vec = np.zeros((num_points, 1))
    # Build A and w matrics
    for i in range(num_points):
        # Left image coords
        xL = left[i,1]
        yL = left[i,2]
        # Right image coords; update first with current rotation
        right_prime = M.T.dot(np.array([[right[i,1]],
                                        [right[i,2]],
                                        [-c]]))
        xR_prime = right_prime[0,0]
        yR_prime = right_prime[1,0]
        zR_prime = right_prime[2,0]
        # Variables for determinants
        A = -yR_prime*np.sin(omega) + zR_prime*np.cos(omega)
        B = xR_prime*np.sin(omega)
        C = -xR_prime*np.cos(omega)
        D = -yR_prime*np.cos(omega)*np.cos(phi) - zR_prime*np.sin(omega)*np.cos(phi)
        E = xR_prime*np.cos(omega)*np.cos(phi) - zR_prime*np.sin(phi)
        F = xR_prime*np.sin(omega)*np.cos(phi) + yR_prime*np.sin(phi)
        # Determinants
        dby = np.linalg.det(np.array([[0,        1,        0       ],
                                      [xL,       yL,      -c       ],
                                      [xR_prime, yR_prime, zR_prime]]))
        dbz = np.linalg.det(np.array([[0,        0,        1       ],
                                      [xL,       yL,      -c       ],
                                      [xR_prime, yR_prime, zR_prime]]))
        domega = np.linalg.det(np.array([[bx, by,       bz      ],
                                         [xL, yL,      -c       ],
                                         [0, -zR_prime, yR_prime]]))
        dphi = np.linalg.det(np.array([[bx, by,  bz],
                                       [xL, yL, -c ],
                                       [A,  B,   C ]]))
        dkappa = np.linalg.det(np.array([[bx, by,  bz],
                                         [xL, yL, -c ],
                                         [D,  E,   F ]]))
        misclosure = np.linalg.det(np.array([[bx,       by,       bz      ],
                                             [xL,       yL,      -c       ],
                                             [xR_prime, yR_prime, zR_prime]]))
        # Populate current row in A and w
        A_mat[i,:] = [dby, dbz, domega, dphi, dkappa]
        w_vec[i,0] = misclosure

    return A_mat, w_vec
