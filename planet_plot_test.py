#planet_plot_test.py
#try and plot positions of Mercury, Sun and Earth at a given datettime

import astropy.units as u
import astropy.coordinates as c
from astropy.wcs import WCS
from datetime import datetime as dt
from astropy.time import Time

obstime=Time('2011-11-03')
earth=c.get_body('Earth',obstime)
mercury=c.get_body('Mercury',obstime)
sun=c.get_sun(obstime)

print earth.geocentrictrueecliptic
print mercury.geocentrictrueecliptic
print sun.geocentrictrueecliptic

def to_spherical_cartesian(body):
   bx=body.geocentrictrueecliptic.distance*np.cos(body.geocentrictrueecliptic.lat)*np.cos(body.geocentrictrueecliptic.lon)
   by=body.geocentrictrueecliptic.distance*np.cos(body.geocentrictrueecliptic.lat)*np.cos(body.geocentrictrueecliptic.lon)
   return bx,by

mx,my=to_spherical_cartesian(mercury)
sx,sy=to_spherical_cartesian(sun)

fig,ax=plt.subplots()
ax.plot([0,0], color='b',marker='o')
ax.plot(my,my, color='r',marker='o')
ax.plot(sx,sy, color='y',marker='o')

fig.show()
