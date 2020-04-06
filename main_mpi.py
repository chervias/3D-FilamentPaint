import FilamentPaint
import numpy as np
import healpy as hp
from mpi4py import MPI
import sys
sys.path.insert(0,'code')
from Sky import *
#from MagFieldCopy import MagField
from MagField import MagField
#from FilPopCopy import FilPop
from FilPop import FilPop

# usage: python main_mpi.py Nfils nside theta_LH size_ratio slope size_box size_scale suffix

nside			= int(sys.argv[2])
Npix_box		= 256

try:
	theta_LH_RMS	= float(sys.argv[3]) # in degrees, if -1 the filaments are perpendicular to the LOS
except ValueError:
	theta_LH_RMS	= None

size_scale		= float(sys.argv[7]) # size scale
size_ratio		= float(sys.argv[4])
slope			= float(sys.argv[5])
Nfil			= int(sys.argv[1])
size_box		= float(sys.argv[6]) # physical size box
suffix = str(sys.argv[8])

#dust_template = '/home/chervias/CMB-work/Filaments/COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits'

shared_comm = MPI.COMM_WORLD.Split_type(MPI.COMM_TYPE_SHARED)
#print("Shared comm contains: ", shared_comm.Get_size(), " processes")
size_pool = shared_comm.Get_size()
rank = shared_comm.Get_rank()

if rank==0:
	output_tqumap	= 'test_ns2048/tqumap_ns%s_%s_maa%s_sr%s_sl%s_minsize%s_%s.fits'%(str(nside),str(Nfil),str(theta_LH_RMS).replace('.','p'),str(size_ratio).replace('.','p'),str(slope).replace('.','p'),str(size_scale).replace('.','p'),suffix)
	#print(output_tqumap)

double_size = MPI.DOUBLE.Get_size()
size = (12*nside**2,3,3)
if rank==0:
	total_size = np.prod(size)
	nbytes = total_size * double_size
else:
	nbytes = 0

shared_comm.Barrier()
win = MPI.Win.Allocate_shared(nbytes, double_size, comm=shared_comm)
# Construct the array                                                                                                                     
buf, itemsize = win.Shared_query(0)
_storedZModes = np.ndarray(buffer=buf, dtype=np.double, shape=size)

win.Fence()
# rank 0 will only fill the array
if rank==0:
	r_unit_vectors  = get_r_unit_vectors(nside)
	local_triad     = get_local_triad(nside,r_unit_vectors)
	#print("RANK: ", shared_comm.Get_rank() , " is filling the array ")
	win.Put(local_triad,0,0)
	#print("RANK: ", shared_comm.Get_rank() , " SUCCESSFULLY filled the array ")
win.Fence()

shared_comm.Barrier()

if rank==0:
	# Create the magnetic field object
	magfield		= MagField(size_box,Npix_box,23456)
	# Create the filament population object
	population		= FilPop(Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,2345,fixed_distance=False)
	#population		= FilPop(nside,Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,1234,dust_template,fixed_distance=False)
else:
	magfield       = None
	population     = None

# Broadcast the objects
magfield 	= shared_comm.bcast(magfield,root=0)
population	= shared_comm.bcast(population,root=0)

if shared_comm.rank == 0:
	#test = np.arange(population.realNfil,dtype='int32')
	test = np.arange(0,population.realNfil,dtype='int32')
	split = np.array_split(test, size_pool)
	split_size = [len(split[i]) for i in range(len(split))]
	split_disp = np.insert(np.cumsum(split_size), 0, 0)[0:-1]
else:
	test = None
	split = None
	split_size = None
	split_disp = None

split_size = shared_comm.bcast(split_size, root = 0)
split_disp = shared_comm.bcast(split_disp, root = 0)
test_local = np.zeros(split_size[rank],dtype='int32')
shared_comm.Scatterv([test, split_size, split_disp, MPI.INT], test_local, root=0)
print('rank ', rank, ': ', test_local)

#output_chunk = np.zeros(np.shape(split[rank])) #Create array to receive subset of data on each core, where rank specifies the core
#shared_comm.Scatterv([test,split_sizes_input, displacements_input,MPI.INT],output_chunk,root=0)

#rcv_buff = np.empty(Nscatter, dtype=np.int32)
#shared_comm.Scatterv([send_buff], rcv_buff, root=0)
#print('Process',rank,'received the numbers',output_chunk)
tqu_total = np.zeros((3,12*nside**2))

for n in test_local:
	tqu_total			+= FilamentPaint.Paint_Filament(n,nside,_storedZModes,population,magfield)
# put a barrier to make sure all processeses are finished
shared_comm.Barrier()
if rank==0:
	# only processor 0 will actually get the data
	tqu_final = np.zeros_like(tqu_total)
else:
	tqu_final = None
# use MPI to get the totals 
shared_comm.Reduce([tqu_total, MPI.DOUBLE],[tqu_final, MPI.DOUBLE],op = MPI.SUM,root = 0)

if rank==0:
	tqu_map = [tqu_final[0,:],tqu_final[1,:],tqu_final[2,:]]
	hp.write_map(output_tqumap,tqu_map,nest=False,overwrite=True)
