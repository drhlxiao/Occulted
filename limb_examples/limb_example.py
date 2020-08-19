 #######################################
# limb_example.py
# Erica Lastufka 27/6/17 

#Description: Example of overplotting the limb using Sunpy coordinate transformations. Based on example from Stuard Mumford https://gist.github.com/Cadair/de415eeeaad7f340c33cde1a8f618898#file-aia-limb-on-stereo-ipynb
#######################################

import numpy as np
import sunpy.map
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.net import vso
from datetime import timedelta as td
from datetime import datetime as dt
from sunpy.coordinates import frames

euvmap= sunpy.map.Map('../../data/stereo-aia/20130427_105530_n4eua.fts')
EUVmap={euvmap.instrument: euvmap.submap(SkyCoord((-1100, 1100) * u.arcsec, (-1100, 1100) * u.arcsec,frame=euvmap.coordinate_frame))}

m=EUVmap['SECCHI']

current_map_time=m.date

#vc = vso.VSOClient()
#instr= vso.attrs.Instrument('AIA')
#sample = vso.attrs.Sample(24 * u.hour)
#wave = vso.attrs.Wave(19.1 * u.nm, 19.45 * u.nm)
#time = vso.attrs.Time(dt.strftime(current_map_time,'%Y-%m-%dT%H:%M:%S'),dt.strftime(current_map_time +td(seconds=2),'%Y-%m-%dT%H:%M:%S')) #should have data to within 1 s
#res=vc.query(wave, sample, time, instr)

#files = vc.get(res).wait()

aiamap= sunpy.map.Map('aia_lev1_193a_2013_04_27t10_55_30_84z_image_lev1.2.fits')
AIAmap={aiamap.instrument: aiamap.submap(SkyCoord((-1100, 1100) * u.arcsec, (-1100, 1100) * u.arcsec,frame=aiamap.coordinate_frame))}

print AIAmap[AIAmap.keys()[0]].rsun_obs
print AIAmap[AIAmap.keys()[0]].rsun_meters#.to(u.arcsec)

fig = plt.figure(figsize=[3, 5])
ax = fig.add_subplot(1, 1, 1, projection=m)
m.plot(axes=ax)
                
ax.set_autoscale_on(False)
r = AIAmap[AIAmap.keys()[0]].rsun_obs.to(u.deg)#-1*u.arcsec # remove the one arcsec so it's on disk.
th = np.linspace(-360*u.deg, -180*u.deg) #x<0
                    
x = r * np.sin(th)
y = r * np.cos(th)
coords = SkyCoord(x, y, frame=AIAmap[AIAmap.keys()[0]].coordinate_frame) #,rsun=AIAmap[AIAmap.keys()[0]].rsun_meters)

hgs = coords.transform_to(frames.HeliographicStonyhurst)
hgs.D0 = m.dsun
hgs.L0 = m.heliographic_longitude-1*u.arcsec # put back the 1 arcsec                    
hgs.B0 = m.heliographic_latitude
limbcoords = hgs.transform_to(m.coordinate_frame)
ax.plot_coord(limbcoords, color='w') 

#zindex=np.searchsorted(limbcoords.Ty,0)
print euvmap.rsun_meters
print limbcoords.Ty[23], limbcoords.Tx[23]
#fig.show()

#now do it in IDL to show the difference....

#from idlpy import *
#IDL.run("fits2map,'20130427_105530_n4eua.fts',map")
#IDL.run("loadct,9")
#IDL.run("plot_map,map,center=[-600,200],fov=[10],/log")
#IDL.run("map_oplot_limb_from_earth,map")     
#IDL.run("map_oplot_limb_from_earth,map,factor=1.01")

#try the bridge later

#import pidly

#idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
#coords=idl.run('limb_example.pro')
#idl.close()
#get limbcoords from idl

#overplot on pthyon plot
#lcidl=SkyCoord(coords,frame='heliographic_stonyhurst')
#ax.plot_coord(lcidl, color='b')

fig.show()
            
