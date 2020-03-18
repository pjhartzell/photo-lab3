import numpy as numpy

def R1(omega):
    sin_omega = np.sin(omega)
    cos_omega = np.cos(omega)
    return np.array([[1,  0,          0        ],
                     [0,  cos_omega, -sin_omega],
                     [0,  sin_omega,  cos_omega]])

def R2(phi):
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    return np.array([[ cos_phi, 0, sin_phi],
                     [ 0,       1, 0       ],
                     [-sin_phi, 0, cos_phi]])

def R3(kappa):
    sin_kappa = np.sin(kappa)
    cos_kappa = np.cos(kappa)
    return np.array([[cos_kappa, -sin_kappa, 0],
                     [sin_kappa,  cos_kappa, 0],
                     [0,          0,         1]])

def M_rot(omega, phi, kappa):
    return R3(kappa).dot(R2(phi)).dot(R1(omega))
