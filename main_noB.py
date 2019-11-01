import sys
sys.path.insert(0,'src')
import numpy as np
import healpy as hp
import matplotlib.pyplot as pl
from Functions import *
from Sky import Sky
from MagField_noB import MagField_noB
from FilPop import FilPop
import multiprocessing as mp

output_tqumap	= 'test_ns512_Nf10k.fits'
nside			= 512
Npix_box		= 256
theta_LH_RMS	= 50.0 # in degrees
Nfil			= 10000
size_box		= 402.0 # physical size box
# Create the sky object
sky			= Sky(nside)
# Create the magnetic field object
magfield		= MagField_noB(size_box,'magfield_Bcube.npy')
# Create the filament population object
population		= FilPop(Nfil,theta_LH_RMS,magfield)

data	= [(n,sky,population,magfield) for n in range(Nfil)]
p 		= mp.Pool(processes=32)
for x in p.imap_unordered(paint_filament_aux, data,chunksize=20):
	t,q,u		= x
	sky.Tmap 	+= t
	sky.Qmap 	+= q
	sky.Umap 	+= u

sky.save_sky(output_tqumap)
