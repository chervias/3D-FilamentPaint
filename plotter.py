import healpy as hp
import matplotlib.pyplot as pl

tqu_map	= hp.read_map('test_ns512_Nf10k.fits',nest=False,field=(0,1,2))

cname='hot'
c = pl.get_cmap(cname)
c.set_under(alpha=0.0)

hp.mollview(tqu_map[0])
hp.mollview(tqu_map[1])
hp.mollview(tqu_map[2])
pl.show()
