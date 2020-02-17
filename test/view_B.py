import matplotlib.pyplot as pl
import numpy as np
from MagField import MagField
from mpl_toolkits.mplot3d import Axes3D

Npix_box		= 128
size_box		= 40.0 # physical size box

magfield	= MagField(size_box,Npix_box,12345)

x, y, z = np.meshgrid(np.arange(-20, 20, 2),
                      np.arange(-20, 20, 2),
                      np.arange(-20, 20, 2))
Bcube	= magfield.interp_fn((x,y,z))

fig = pl.figure()
ax = fig.gca(projection='3d')

ax.quiver(x, y, z, Bcube[:,:,:,0], Bcube[:,:,:,1], Bcube[:,:,:,2], length=1, normalize=True)

pl.show()
