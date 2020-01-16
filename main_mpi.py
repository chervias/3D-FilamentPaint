import FilamentPaint
import numpy as np
import healpy as hp
from mpi4py import MPI
import sys
sys.path.insert(0,'code')
from Sky import Sky
from MagField import MagField
from FilPop import FilPop

output_tqumap	= 'test_ns512_112k_maa45_sr0p1_sl5p0.fits'
nside			= 512
Npix_box		= 256
theta_LH_RMS	= 45. # in degrees
size_scale		= 0.7
size_ratio		= 0.1
slope			= 5.0
Nfil			= 112000
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
	population		= FilPop(Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,fixed_distance=False)
else:
	sky = None
	magfield = None
	population = None

#print(rank,population)
# Broadcast the objects
sky 		= comm.bcast(sky,root=0)
magfield 	= comm.bcast(magfield,root=0)
population	= comm.bcast(population,root=0)
#print(rank,population)
#exit()

# rank=0 will create an array with size Nfil
Nscatter = Nfil // size
if comm.rank == 0:
	send_buff = np.arange(Nfil, dtype=np.int32)
else:
	send_buff = np.empty(Nfil, dtype=np.int32)
rcv_buff = np.empty(Nscatter, dtype=np.int32)
comm.Scatter(send_buff, rcv_buff, root=0)
#print('Process',rank,'received the numbers',rcv_buff)
tqu_total = np.zeros((3,12*nside**2))

for n in rcv_buff:
	#print('Process=',rank,'is working on filament',n,end='')
	tqu_total			+= FilamentPaint.Paint_Filament(n,sky,population,magfield)

# put a barrier to make sure all processeses are finished
comm.Barrier()
if rank==0:
	# only processor 0 will actually get the data
	tqu_final = np.zeros_like(tqu_total)
else:
	tqu_final = None
# use MPI to get the totals 
comm.Reduce([tqu_total, MPI.DOUBLE],[tqu_final, MPI.DOUBLE],op = MPI.SUM,root = 0)

if rank==0:
	sky.Tmap = tqu_final[0,:]
	sky.Qmap = tqu_final[1,:]
	sky.Umap = tqu_final[2,:]
	sky.save_sky(output_tqumap)
