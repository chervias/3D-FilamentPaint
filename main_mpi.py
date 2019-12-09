import FilamentPaint
import numpy as np
import healpy as hp
from mpi4py import MPI
import sys
sys.path.insert(0,'code')
from Sky import Sky
from MagField import MagField
from FilPop import FilPop

nside			= 512
Npix_box		= 256
theta_LH_RMS	= None # in degrees
size_scale		= 0.7
size_ratio		= 0.25
Nfil			= 100
size_box		= 1500.0 # physical size box

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# rank=0 process will create the objects and distribute them
if rank==0:
	# Create the sky object
	sky				= Sky(nside)
	# Create the magnetic field object
	magfield		= MagField(size_box,Npix_box,12345)
	# Create the filament population object
	population		= FilPop(Nfil,theta_LH_RMS,size_ratio,size_scale,magfield,fixed_distance=True)
else:
	sky = None
	magfield = None
	population = None

# Broadcast the objects
comm.bcast(sky,root=0)
comm.bcast(magfield,root=0)
comm.bcast(population,root=0)

# rank=0 will create an array with size Nfil
Nscatter = Nfil // size
if comm.rank == 0:
	send_buff = np.arange(Nfil, dtype=np.int32)
else:
	send_buff = np.empty(Nfil, dtype=np.int32)
rcv_buff = np.empty(Nscatter, dtype=np.int32)
comm.Scatter(send_buff, rcv_buff, root=0)
print('Process',rank,'received the numbers',rcv_buff)

tqu_total = np.zeros((3,12*nside**2))

for n in rcv_buff:
	#print(type(n))
	tqu_total			+= FilamentPaint.Paint_Filament(n,sky,population,magfield)