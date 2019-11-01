import sys
sys.path.insert(0,'src')
import numpy as np
import healpy as hp
import matplotlib.pyplot as pl
from Functions import *
from Sky import Sky
from MagField import MagField_1fil
from FilPop import FilPop_1fil
import multiprocessing as mp

nside			= 256
Npix_box		= 256
Nfil			= 1
size_box		= 100.0 # physical size box

# Create the filament population object
population	= FilPop_1fil(Nfil)

(h,w)	= pl.figaspect(3/5.)
fig		= pl.figure(figsize=(h,w))

minP 	= -0.05
maxP	= +0.05

c = 0
for direc in ['-z','45deg','+y']:
	sky			= Sky(nside)
	magfield	= MagField_1fil(size_box,Npix_box,12345,direc)
	t,q,u		= paint_filament(0,sky,population,magfield)
	sky.Tmap 	+= t
	sky.Qmap 	+= q
	sky.Umap 	+= u
	
	Tlm,Elm,Blm = hp.map2alm((sky.Tmap,sky.Qmap,sky.Umap),iter=1)
	Emap 		= hp.alm2map(Elm,nside,verbose=False)
	Bmap 		= hp.alm2map(Blm,nside,verbose=False)

	#sky.mask[pix_filament]	= 1.0
	#hp.mollview(sky.mask,nest=False)
	hp.gnomview(sky.Tmap,rot=(0,0),nest=False,reso=15,fig=1,sub=(3,5,c+1),title='T '+direc,notext=True,cbar=False,cmap='bwr')
	hp.gnomview(sky.Qmap,rot=(0,0),nest=False,reso=15,fig=1,sub=(3,5,c+2),title='Q '+direc,notext=True,cbar=False,min=minP,max=maxP,cmap='bwr')
	hp.gnomview(sky.Umap,rot=(0,0),nest=False,reso=15,fig=1,sub=(3,5,c+3),title='U '+direc,notext=True,cbar=False,min=minP,max=maxP,cmap='bwr')
	hp.gnomview(Emap,rot=(0,0),nest=False,reso=15,fig=1,sub=(3,5,c+4),title='E '+direc,notext=True,cbar=False,min=minP,max=maxP,cmap='bwr')
	hp.gnomview(Bmap,rot=(0,0),nest=False,reso=15,fig=1,sub=(3,5,c+5),title='B '+direc,notext=True,cbar=False,min=minP,max=maxP,cmap='bwr')
	c  = c + 5

pl.tight_layout()
pl.savefig('QU_sign.pdf',format='pdf')
