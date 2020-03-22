from relative_orientation import RelativeOrientation


# Create relative orientation object
my_ro = RelativeOrientation()

# Read in image coordinates, make sure the point numbers match, visualize
my_ro.read_left('image_27_given.txt')
my_ro.read_right('image_28_given.txt')
my_ro.match_coords()
my_ro.plot_coords()

# Define the baseline x component and camera constant
my_ro.bx = 92
my_ro.c = 153.358

# Solve the relative orientation parameters and print
my_ro.solve_ro()
my_ro.report_ro()

# Compute the model space coordinates and print
my_ro.compute_model()
my_ro.report_model()

# Print the y-parallax
my_ro.report_parallax()

# Plot scale factors from model space coordinate computations
my_ro.plot_scale_factors()

# Print correlation matrix from the least squares design matrix
my_ro.report_correlation()