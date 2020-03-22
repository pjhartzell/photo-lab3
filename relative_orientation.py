import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
from rotations import M_rot


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
        self.iter_count = 0

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

    # View the coordinate pattern in left and right images
    def plot_coords(self):
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
        ax1.plot(self.left[:,1], self.left[:,2], 'bo')
        ax1.axis('square')
        ax1.set(xlim=(-228.6, 228.6), ylim=(-228.6, 228.6))
        ax1.set_title('Left Image')
        ax2.plot(self.right[:,1], self.right[:,2], 'bo')
        ax2.axis('square')
        ax2.set(xlim=(-228.6, 228.6), ylim=(-228.6, 228.6))
        ax2.set_title('Right Image')
        plt.show()

    # Least squares Relative Orientation
    def solve_ro(self):
        # Check for minimum number of points
        if self.left.shape[0] < 5:
            print('Minimum of 5 matching points is required.')
            return 0
        # Least squares correction thresholds set to 0.1 micrometer in image
        # (and model) space, which is an order of magnitude below comparator
        # measurement precision.
        angle_threshold = 0.0001 / self.c   # Radians
        baseline_threshold = 0.0001         # Millimeters
        # Prep some variables
        delta_hat = np.ones((5,1))
        x_hat = np.array([[self.by],
                          [self.bz],
                          [self.omega],
                          [self.phi],
                          [self.kappa]])
        # Iterate until delta_hat corrections are less than thresholds
        while (((np.absolute(delta_hat[0:2,0]) > baseline_threshold).any() or 
                (np.absolute(delta_hat[2:,0]) > angle_threshold).any())):
            # Form A and w matrices with current parameter estimates
            A, w = self.A_and_w(x_hat, self.left, self.right, self.bx, self.c)
            # Least squares solution
            delta_hat = -np.linalg.inv(A.T.dot(A)).dot(A.T).dot(w)
            # Update current cumulative solution estimate
            x_hat += delta_hat
            self.iter_count += 1
        # Update class variables with final values
        self.by = x_hat[0,0]
        self.bz = x_hat[1,0]
        self.omega = x_hat[2,0]
        self.phi = x_hat[3,0]
        self.kappa = x_hat[4,0]
        # Compute correlation matrix
        self.correlation = self.compute_correlation(A)

    # Design (A) and misclosure (w) matrices
    def A_and_w(self, x_hat, left, right, bx, c):
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

    # Correlation matrix
    def compute_correlation(self, A):
        ATA = A.T.dot(A)
        corr = np.zeros(ATA.shape)
        for i in range(ATA.shape[0]):
            for j in range(ATA.shape[1]):
                corr[i,j] = ATA[i,j] / np.sqrt(ATA[i,i]*ATA[j,j])
        return corr

    # Print solved RO parameters
    def report_ro(self):
        print('\nNumber of iterations = {}'.format(self.iter_count))
        print('by = {:.3f} mm'.format(self.by))
        print('bz = {:.3f} mm'.format(self.bz))
        print('omega = {:.5f} degrees'.format(np.rad2deg(self.omega)))
        print('phi = {:.5f} degrees'.format(np.rad2deg(self.phi)))
        print('kappa = {:.5f} degrees'.format(np.rad2deg(self.kappa)))

    # Model space coordinates
    def compute_model(self):
        # Prep some variables
        num_points = self.right.shape[0]
        lambda_scale = np.zeros((num_points, 1))
        mu_scale = np.zeros((num_points, 1))
        left_model = np.zeros((num_points, 4))
        right_model = np.zeros((num_points, 4))
        model = np.zeros((num_points, 4))
        # Rotate right image coords
        M = M_rot(self.omega, self.phi, self.kappa)
        right_coords = np.hstack(
            (self.right[:,1:], -self.c * np.ones((num_points, 1)))
        )
        right_prime = M.T.dot(right_coords.T).T
        # Compute scales and model coordinates
        for i in range(num_points):
            # Unpack some values for clarity
            bx = self.bx
            by = self.by
            bz = self.bz
            xL = self.left[i,1]
            yL = self.left[i,2]
            c = self.c
            xR_prime = right_prime[i,0]
            yR_prime = right_prime[i,1]
            zR_prime = right_prime[i,2]
            # Scales
            lambda_scale[i,0] = \
                (bx*zR_prime - bz*xR_prime) / (xL*zR_prime + c*xR_prime)
            mu_scale[i,0] = (-bx*c - bz*xL) / (xL*zR_prime + c*xR_prime)
            # Model coords
            left_model[i,0] = self.left[i,0]                    # point num
            left_model[i,1] = lambda_scale[i,0] * xL            # x
            left_model[i,2] = lambda_scale[i,0] * yL            # y
            left_model[i,3] = -lambda_scale[i,0] * c            # z
            right_model[i,0] = self.right[i,0]                  # point num
            right_model[i,1] = mu_scale[i,0] * xR_prime + bx    # x
            right_model[i,2] = mu_scale[i,0] * yR_prime + by    # y
            right_model[i,3] = mu_scale[i,0] * zR_prime + bz    # z
            # Check X and Z coords for equivalence to 0.1 micrometer
            if (np.absolute(left_model[i,[1,3]] - 
                            right_model[i,[1,3]]) > 0.0001).any():
                print('Error in model space coordinate {}.'.format(
                    self.left[i,0]))
            # Compute and store final model coordinates
            model[i,0] = left_model[i,0]
            model[i,1] = left_model[i,1]
            model[i,2] = (left_model[i,2] + right_model[i,2]) / 2
            model[i,3] = left_model[i,3]
        # Make scales and model space coordinates accessible to class
        self.lambda_scale = lambda_scale
        self.mu_scale = mu_scale
        self.left_model = left_model
        self.right_model = right_model
        self.model = model

    # Print model space coordinates
    def report_model(self):
        print('\nModel Space Coordinates')
        headers = ['#', 'x (mm)', 'y (mm)', 'z (mm)']
        print(tabulate(self.model, headers, tablefmt='simple',
            floatfmt=('.0f', '.3f', '.3f', '.3f')))

    # Print model space y-parallax
    def report_parallax(self):
        parallax = self.right_model[:,2] - self.left_model[:,2]
        parallax = np.expand_dims(parallax, axis=1)
        point_nums = np.expand_dims(self.left[:,0], axis=1)
        parallax = np.hstack((point_nums, parallax))
        print('\nY Parallax')
        headers = ['#', 'Parallax (mm)']
        print(tabulate(parallax, headers, tablefmt='simple',
            floatfmt=('.0f', '.3f')))
    
    # Plot scale factors
    def plot_scale_factors(self):
        point_nums = self.left[:,0]
        plt.plot(point_nums, self.lambda_scale, label='Left Image (27)')
        plt.plot(point_nums, self.mu_scale, label='Right Image (28)')
        plt.xlabel('Point Number')
        plt.ylabel('Scale')
        plt.title('Left and Right Image Scale Factors')
        plt.legend()
        plt.show()
    
    # Print correlation matrix
    def report_correlation(self):
        print('\nCorrelation Coefficient Matrix')
        headers = ['by', 'bz', 'omega', 'phi', 'kappa']
        print(tabulate(self.correlation, headers, tablefmt='simple',
            floatfmt=('.3f', '.3f', '.3f', '.3f', '.3f')))