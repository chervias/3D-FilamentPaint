import FilamentPaint
import numpy as np
import healpy as hp
from mpi4py import MPI
import sys
sys.path.insert(0,'code')
from Sky import Sky,SkyAux
from MagField import MagField
from FilPop import FilPop

# usage: python main_mpi.py Nfils nside theta_LH size_ratio slope size_box size_scale

nside			= int(sys.argv[2])
Npix_box		= 256
theta_LH_RMS	= float(sys.argv[3]) # in degrees, if -1 all filaments point up
size_scale		= float(sys.argv[7]) # size scale
size_ratio		= float(sys.argv[4])
slope			= float(sys.argv[5])
Nfil			= int(sys.argv[1])
size_box		= float(sys.argv[6]) # physical size box

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank==0:
	output_tqumap	= 'test_ns2048/tqumap_ns%s_%s_maa%s_sr%s_sl%s_minsize%s_wLarson.fits'%(str(nside),str(Nfil),str(theta_LH_RMS).replace('.','p'),str(size_ratio).replace('.','p'),str(slope).replace('.','p'),str(size_scale).replace('.','p'))
	#print(output_tqumap)

# rank=0 process will create the objects and distribute them
if rank==0:
	# Create the sky object
	sky				= Sky(nside)
	r_unit_vectors  = sky.r_unit_vectors
	local_triad     = sky.local_triad
	skyaux          = SkyAux(nside)
	# Create the magnetic field object
	magfield		= MagField(size_box,Npix_box,12345)
	# Create the filament population object
	population		= FilPop(nside,Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,1234,fixed_distance=True)
else:
	skyaux         = SkyAux(nside)
	r_unit_vectors = np.empty((12*nside**2,3))
	local_triad    = np.empty((12*nside**2,3,3))
	magfield       = None
	population     = None

# Broadcast the objects
# The Sky object is too big to be broadcasted. I will broadcast its elements and 
comm.Bcast(r_unit_vectors,root=0)
comm.Bcast(local_triad,root=0)
skyaux.r_unit_vectors = r_unit_vectors
skyaux.local_triad    = local_triad

magfield 	= comm.bcast(magfield,root=0)
population	= comm.bcast(population,root=0)

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
	tqu_total			+= FilamentPaint.Paint_Filament(n,skyaux,population,magfield)

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
