# left and right image coordinate pairs (5 pair minimum)
# bx value
# initial values for 5 unknowns

# set large delcap values
# while delcap values > threshold
#   compute misclosure with current estimates
#   compute design matrix with current estimates
#   least squares for delcap vector
#   update vector of unknowns with latest delcap

# report number of iterations
# report solved values of unknowns
# report correlation matrixx

import numpy as np
from lsq_matrices import A_and_w

class RelativeOrientation:
    # Some default values
    def __init__(self):
        self.by = 0
        self.bz = 0
        self.omega = 0
        self.phi = 0
        self.kappa = 0

    # Read image coords; order = point#, x, y
    def read_left(self, left_file):
        self.left = np.loadtxt(left_file)
    def read_right(self, right_file):
        self.right = np.loadtxt(right_file)
    
    # Match left and right image coords in case not in same order and/or non-
    # matching coordinates exist
    def match_coords(self):
        point_nums, left_indices, right_indices = np.intersect1d(
            self.left[:,0], self.right[:,0],
            assume_unique=True, return_indices=True
        )
        self.left = self.left[left_indices,:]
        self.right = self.right[right_indices,:]

    # Set bx baseline distance
    def set_bx(self, bx):
        self.bx = bx
    
    # Set flying height (for LSQ iteration termination)
    def set_agl(self, agl):
        self.agl = agl

    # Set camera constant (~f)
    def set_c(self, c):
        self.c = c

    # Least squares Relative Orientation
    def relative_orientation(self):
        # Check for minimum number of points
        if self.left.shape[0] < 5:
            print('Minimum of 5 matching points is required.')
            return 0

        # Least squares angle correction threshold in radians
        # Equivalent to 1 mm on the ground
        angle_threshold = 0.001 / self.agl
        # Least squares baseline component correction threshold in meters
        # 1 micrometer in image space
        baseline_threshold = 0.000001

        # Iterate until corrections are less than thresholds
        delta_hat = np.ones((5,1))
        while (delta_hat[0:2,0] > baseline_threshold).any() or 
                (delta_hat[2:,0] > angle_threshold).any():

            # Form A and w matrices with current parameter estimates
            A, w = A_and_w(self.omega, self.phi, self.kappa, 
                           self.bx, self.by, self.bz
                           self.left, self.right, self.c)




my_ro = RelativeOrientation()
my_ro.read_left('image_27_corrected.txt')
my_ro.read_right('image_28_corrected.txt')
my_ro.match_coords()
print(my_ro.left)
print(my_ro.right)
my_ro.set_bx(92)
my_ro.relative_orientation()