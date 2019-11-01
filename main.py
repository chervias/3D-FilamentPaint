import numpy as np
import healpy as hp
import matplotlib.pyplot as pl
from src/Functions import *
from src/Sky import Sky
from src/MagField import MagField
from src/FilPop import FilPop
import multiprocessing as mp

output_tqumap	= 'test_ns512_Nf10k.fits'
nside			= 256
Npix_box		= 256
Nfil			= 10000
size_box		= 200.0 # physical size box
# Create the sky object
sky				= Sky(nside)
# Create the magnetic field object
magfield		= MagField(size_box,Npix_box,12345)
np.save('magfield_Bcube.npy',magfield.Bcube)
exit()
# Create the filament population object
population		= FilPop(Nfil,magfield)

data	= [(n,sky,population,magfield) for n in range(Nfil)]
p 		= mp.Pool(processes=2)
for x in p.imap_unordered(paint_filament_aux, data,chunksize=2):
	t,q,u		= x
	sky.Tmap 	+= t
	sky.Qmap 	+= q
	sky.Umap 	+= u

#sky.mask[pix_filament]	= 1.0
#hp.mollview(sky.mask,nest=False)
hp.mollview(sky.Tmap,nest=False)
hp.mollview(sky.Qmap,nest=False)
hp.mollview(sky.Umap,nest=False)
pl.show()

sky.save_sky(output_tqumap)
