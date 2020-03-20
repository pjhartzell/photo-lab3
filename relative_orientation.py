# report number of iterations
# report solved values of unknowns
# report correlation matrixx

import numpy as np
from lsq_matrices import A_and_w

class RelativeOrientation:
    # Some default values
    def __init__(self):
        self.bx = 0.0
        self.by = 0.0
        self.bz = 0.0
        self.omega = 0.0
        self.phi = 0.0
        self.kappa = 0.0
        self.c = 0.0
        self.agl = 0.0

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

    # Least squares Relative Orientation
    def relative_orientation(self):
        # Check for minimum number of points
        if self.left.shape[0] < 5:
            print('Minimum of 5 matching points is required.')
            return 0

        # Least squares angle correction threshold (radians)
        # Equivalent to 1 mm on the ground
        angle_threshold = 0.001 / self.agl
        # Least squares baseline component correction threshold (mm)
        # 1 micrometer in image space, which is well below comparator precision
        baseline_threshold = 0.001

        # Prep some variables
        delta_hat = np.ones((5,1))
        x_hat = np.array([[self.by],
                          [self.bz],
                          [self.omega],
                          [self.phi],
                          [self.kappa]])
        # Iterate until delta_hat corrections are less than thresholds
        iter_count = 0
        while (((delta_hat[0:2,0] > baseline_threshold).any() or 
                (delta_hat[2:,0] > angle_threshold).any()) and
                 iter_count < 20):
            input()
            # Form A and w matrices with current parameter estimates
            A, w = A_and_w(x_hat, self.left, self.right, self.bx, self.c)
            print(w)
            # Least squares solution
            delta_hat = -np.linalg.inv(A.T.dot(A)).dot(A.T).dot(w)
            print("delta_hat = {}".format(delta_hat))
            # Update current cumulative solution estimate
            x_hat += delta_hat
            # print("x_hat = {}".format(x_hat))
            iter_count += 1
            # print(iter_count)


my_ro = RelativeOrientation()
my_ro.read_left('image_27_corrected.txt')
my_ro.read_right('image_28_corrected.txt')
my_ro.match_coords()
my_ro.bx = 92
my_ro.c = 153.358
my_ro.agl = 1860 - 1100
my_ro.relative_orientation()