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

class RelativeOrientation:
    # Some default values
    def __init__(self):
        self.bx = 0
        self.by = 0
        self.bz = 0
        self.omega = 0
        self.phi = 0
        self.kappa = 0

    # Read image coords; order = point#, x, y
    def read_left(self, left_file):
        self.left_coords = np.loadtxt(left_file)
    def read_right(self, right_file):
        self.right_coords = np.loadtxt(right_file)
    
    # Match left and right image coords in case not in same order and/or non-
    # matching coordinates exist
    def match_coords(self):
        point_nums, left_indices, right_indices = np.intersect1d(
            self.left_coords[:,0], self.right_coords[:,0],
            assume_unique=True, return_indices=True
        )
        self.left_coords = self.left_coords[left_indices,:]
        self.right_coords = self.right_coords[right_indices,:]
    
    # Least squares Relative Orientation
    def relative_orientation(self):
        # Check for minimum number of points
        if self.left_coords.shape[0] < 5:
            print('Minimum of 5 matching points is required.')
            return 0
        



my_ro = RelativeOrientation()
my_ro.read_left('image_27_corrected.txt')
my_ro.read_right('image_28_corrected.txt')
my_ro.match_coords()
print(my_ro.left_coords)
print(my_ro.right_coords)
my_ro.relative_orientation()