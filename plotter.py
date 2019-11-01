import healpy as hp
import matplotlib.pyplot as pl

tqu_map	= hp.read_map('test_ns256_Nf10k.fits',nest=False,field=(0,1,2))

cname='hot'
c = pl.get_cmap(cname)
c.set_under(alpha=0.0)

hp.mollview(tqu_map[0],title='T',sub=311,cbar=False,cmap=c)
hp.mollview(tqu_map[1],title='Q',sub=312,cbar=False,cmap=c)
hp.mollview(tqu_map[2],title='U',sub=313,cbar=False,cmap=c)
#hp.graticule(dpar=360,dmer=180)

pl.savefig("TQUmap.pdf",  bbox_inches = 'tight')
