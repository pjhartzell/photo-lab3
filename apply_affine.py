import numpy as np


def affine_transform(a, b, c, d, t_x, t_y, xy):
    A = np.array([[a, b], [c, d]])
    T = np.array([[t_x],[t_y]])
    return np.matmul(A, xy.T) + T


# Image 27
a = 0.0118994263
b = 0.0000002998
c = -0.0000013405
d = 0.0119012647
t_x = -122.0170
t_y = 123.5343

pxy_obs = np.loadtxt("image_27_obs.txt")
p = pxy_obs[:,0].reshape(-1,1)
xy_obs = pxy_obs[:,1:]
xy_fid = affine_transform(a, b, c, d, t_x, t_y, xy_obs)
pxy_fid = np.hstack((p, xy_fid.T))
np.savetxt("image_27_fid.txt", pxy_fid, fmt='%0.5f')


# Image 28
a = 0.0119000883
b = -0.0000084564
c = 0.0000074035
d = 0.0119010331
t_x = -122.1921
t_y = 123.5180

pxy_obs = np.loadtxt("image_28_obs.txt")
p = pxy_obs[:,0].reshape(-1,1)
xy_obs = pxy_obs[:,1:]
xy_fid = affine_transform(a, b, c, d, t_x, t_y, xy_obs)
pxy_fid = np.hstack((p, xy_fid.T))
np.savetxt("image_28_fid.txt", pxy_fid, fmt='%0.5f')