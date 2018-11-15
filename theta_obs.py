#assume everything in the eccliptic plane!
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord


def senkret(vec):
    '''calculate vectors perpendicular to input vector'''
    m1=vec[1]/vec[0]
    m2=-1./m1
    l1=np.sqrt(vec[0]**2 + vec[1]**2)
    x2=l1/np.sqrt(1+m2**2)
    vEast=(x2, m2*x2)
    vWest=(-x2,-m2*x2)
    #d1=np.rad2deg(np.arctan2(vEast[1],vEast[0]))-np.rad2deg(np.arctan2(vec[1],vec[0]))
    #d2=np.rad2deg(np.arctan2(vWest[1],vWest[0]))-np.rad2deg(np.arctan2(vec[1],vec[0]))
    #print d1,d2
    #line on which vectors lie. Make it as long as the original one
    return vEast, vWest

def theta_obs(v1,v1E,v1W,v2,v2E,v2W):
    ''' is there an theta_obs such that v1 and v2 are included between to vn_E/W ?'''
    #which limb vector of v1 is closest to vector v2?
    tv2=np.rad2deg(np.arctan2(v2[1],v2[0]))
    tv1E=np.rad2deg(np.arctan2(v1E[1],v1E[0]))
    #print tv2,tv1E
    v1Ev2=tv2-tv1E
    tv1W=np.rad2deg(np.arctan2(v1W[1],v1W[0]))    
    v1Wv2=tv2-tv1W
    #print tv2,tv1W
    #print v1Ev2,v1Wv2
    if np.abs(v1Ev2) < np.abs(v1Wv2):
        vfov={'value':v1E, 'side':'E'}
        #print 'v1 east limb closer'
    else:
        vfov={'value':v1W, 'side':'W'}
        #print 'v1 west limb closer'
    #print closest_limbvec
    tfov=np.rad2deg(np.arctan2(vfov['value'][1],vfov['value'][0]))
    #test which limb of v2 has the same sign of the x-coordinate of v1
    #we don't want to look on the wrong side of the Sun!
    #print np.sign(v2E[0]), np.sign(v1[0])
    if np.sign(v2E[0]) == np.sign(v1[0]):
        vlimb={'value':v2E, 'side':'E'}
        #print 'v2 east limb farther'
    else: 
        vlimb={'value':v2W, 'side':'W'}
        #print 'v2 west limb farther'
    vElimb={'value':v2E, 'side':'E'}
    #v2Evfov=np.rad2deg(np.arctan2(v2E[1],v2E[0]))-tfov
    #v2Wvfov=np.rad2deg(np.arctan2(v2W[1],v2W[0]))-tfov
    #print v2Evfov,v2Wvfov
    #if v2Evfov > 180.:
    #    v2Evfov=360.-v2Evfov
    #if v2Wvfov > 180.:
    #    v2Wvfov=360.-v2Wvfov
    #print v2Evfov,v2Wvfov
    #if np.abs(v2Evfov) > np.abs(v2Wvfov): #v2E/W farthest from vfov
    #    vlimb={'value':v2E, 'side':'E'}
    #    print 'v2 east limb farther'
    #else:
    #    vlimb={'value':v2W, 'side':'W'}
    #    print 'v2 west limb farther'
    tlimb= np.rad2deg(np.arctan2(vlimb['value'][1],vlimb['value'][0]))
    #print tlimb,tfov
    theta_obs=np.max([np.abs(tlimb),np.abs(tfov)])-np.min([np.abs(tlimb),np.abs(tfov)])#angle between vlimb and vfov
    theta_ex=180.-theta_obs
    #print theta_obs,theta_ex
    return theta_obs,theta_ex,vlimb,vElimb,vfov

def stereo_loc_to_lon(stereo_loc):
    '''convert location coordinates to heliographic stonyhurst and estimate the longitudinal angle from here'''
    aa=stereo_loc.transform_to('heliographic_stonyhurst')
    #obs_lon=stero_loc.observer.lon.value #negative because it's pointing backwards at the observer (+- 2 degrees off?)
    slon = np.abs(aa.lon.value) #east limb is zero, west is 180
    return slon
    
def is_flare_observed(stereo_lon,vElimb,vfov,vlimb):
    '''stereo_lon is measured from the east limb. What is the angle from vfov to the east limb? is stereo-lon included in this angle?'''
    tvE=np.rad2deg(np.arctan2(vElimb['value'][1],vElimb['value'][0]))
    tvfov=np.rad2deg(np.arctan2(vfov['value'][1],vfov['value'][0]))
    #print tvE,tvfov
    vEfov=np.abs(tvfov-tvE)
    #print vEfov
    if vlimb['side'] == 'W': #if stereo_lon is between vfov and west limb it's observed
        tvW=np.rad2deg(np.arctan2(vlimb['value'][1],vlimb['value'][0]))
        vWfov=np.abs(tvfov-tvW)
        if stereo_lon <=vWfov: #what if it's visible from the other side though? if the farthest limb is the west limb?
            joint_obs = True
        else:
            joint_obs = False
        
    else:
        if stereo_lon <=vEfov: #what if it's visible from the other side though? if the farthest limb is the west limb?
            joint_obs = True
        else:
            joint_obs = False
    return joint_obs
    
#if __name__ == "__main__":
#    #vec1=(-179.3,.5545) 
#    #vec2=(-.052,80.941) #stereo B on 2011-08-20
#    #mx=-177.034
#    #my=-2.24
#    sx=104.27
#    sy=-.1215
#    spx=803.17 #A/B lonlat
#    spy=506.909
#    #v1=(np.deg2rad(my),np.deg2rad(mx)) 
#    v2=(np.deg2rad(sy),np.deg2rad(sx)) #stereo B on 2011-08-20
#    #vec1=(np.cos(v1[0])*np.cos(v1[1]),np.cos(v1[0])*np.sin(v1[1]))
#    vec1=(-.995,-.0908) #mlonlat
#    vec2=(np.cos(v2[0])*np.cos(v2[1]),np.cos(v2[0])*np.sin(v2[1]))
#    print vec1,vec2
#    smap=sunpy.map.Map('20100831_205530_14euA.fts')
#    stereo_loc=SkyCoord(spx*u.arcsec,spy*u.arcsec,frame=smap.coordinate_frame)
#    s1E,s1W=senkret(vec1)
#    s2E,s2W=senkret(vec2)
#    thobs,thex,vlimb,vElimb,vfov=theta_obs(vec1,s1E,s1W,vec2,s2E,s2W)
#    stereo_lon=stereo_loc_to_lon(stereo_loc)
#    joint_obs=is_flare_observed(stereo_lon,vElimb,vfov,vlimb)
#    print stereo_lon
#    print joint_obs
#    fig,ax=plt.subplots()
#    ax.scatter(vec1[0],vec1[1], c='r',marker='o',label='v1',s=50)
#    ax.scatter(0,0, c='b',marker='o',label='Sun',s=150)
#    ax.scatter(vec2[0],vec2[1] , c='g',marker='o',label='v2',s=50)
#    ax.scatter(s1E[0],s1E[1],  c='r',marker='*',label='v1E')
#    ax.scatter(s1W[0],s1W[1],  c='r',marker='v',label='v1W')
#    ax.scatter(s2E[0],s2E[1],  c='g',marker='*',label='v2E')
#    ax.scatter(s2W[0],s2W[1],  c='g',marker='v',label='v2W')
#    ax.plot([0,vec1[0]],[0,vec1[1]],'-r')
#    ax.plot([0,vec2[0]],[0,vec2[1]],'-g')
#    ax.plot([0,s1E[0]],[0,s1E[1]],'--r')
#    ax.plot([0,s1W[0]],[0,s1W[1]],'-.r')
#    ax.plot([0,s2E[0]],[0,s2E[1]],'--g')
#    ax.plot([0,s2W[0]],[0,s2W[1]],'-.g')
#    ax.plot([0,s2W[0]],[0,s2W[1]],'-.g')
#    ax.set_xlabel('x_HEE')
#    ax.set_ylabel('y_HEE')
#    ax.set_aspect('equal')
#    ax.legend()
#    fig.show()
