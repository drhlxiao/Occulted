#assume everything in the eccliptic plane!
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

#DO ALL I/O in SkyCoords!

#keep everything in sunpy coordinates so that there's no confusion and the frame is in the object!
def senkret(vec):
    '''calculate vectors perpendicular to input vector.'''
    #m1=vec[1]/vec[0]
    #m2=-1./m1
    #l1=np.sqrt(vec[0]**2 + vec[1]**2)
    #x2=l1/np.sqrt(1+m2**2) #this is always going to be positive! should be able to be negative too
    vEast=(vec[1], -vec[0]) #is this right in given coordinate system? 'left'
    vWest=(-vec[1],vec[0]) #'right'
    #d1=np.rad2deg(np.arctan2(vEast[1],vEast[0]))-np.rad2deg(np.arctan2(vec[1],vec[0]))
    #d2=np.rad2deg(np.arctan2(vWest[1],vWest[0]))-np.rad2deg(np.arctan2(vec[1],vec[0]))
    #print d1,d2
    #line on which vectors lie. Make it as long as the original one
    return vEast, vWest

def theta_obs(v1,v1E,v1W,v2,v2E,v2W,flare_lon=False, plot=False):
    ''' is there an theta_obs such that v1 and v2 are included between to vn_E/W ?'''
    #input vectors are in HEE Cartesian
    #to identify correct boundary limb vectors, convert to angle theta
    #x=cos(theta), y=sin(theta), theta from -180 to + 180
    if v1[1] >= 0.0: #need positive theta
        th_v1=np.rad2deg(np.arccos(v1[0]))
    else: #need negative theta
        th_v1=-np.rad2deg(np.arccos(v1[0])) #+ 180.

    if v2[1] >= 0.0:
        th_v2=np.rad2deg(np.arccos(v2[0]))
    else:
        th_v2=-np.rad2deg(np.arccos(v2[0])) #+ 180.

    # 0 < th_v1-th_v2 < 180: (E2, W1)

    #print v1,v2,th_v1, th_v2
    if -180.0 < th_v1-th_v2 <= 0.: #this might be wrong now
        vfov={'value':v2E, 'side':'E'} #HEE-ish
        vlimb={'value':v1W, 'side':'W'} #HEE-ish


    # 180 < th_v1 - th_v2 < 360: (W2,E1)
    else:
        vfov={'value':v2W, 'side':'W'} #HEE-ish
        vlimb={'value':v1E, 'side':'E'} #HEE-ish


    tv2=np.rad2deg(np.arctan2(v2[1],v2[0]))
    tv1E=np.rad2deg(np.arctan2(v1E[1],v1E[0]))
    #print tv2,tv1E
    v1Ev2=tv2-tv1E
    tv1W=np.rad2deg(np.arctan2(v1W[1],v1W[0]))
    v1Wv2=tv2-tv1W
    #print tv2,tv1W
    #print v1Ev2,v1Wv2
    #if np.abs(v1Ev2) < np.abs(v1Wv2): #v1 is the observing instrument, ie the one we care about
    ##does this always mean we pick v2E and v1W as the limiting limbs? do that for now

    #    vfov={'value':v1E, 'side':'E'}
        #print 'v1 east limb closer'
    #else:
    #    vfov={'value':v1W, 'side':'W'}
        #print 'v1 west limb closer'
    #print closest_limbvec
    tfov=np.rad2deg(np.arctan2(vfov['value'][1],vfov['value'][0]))
    tvlimb=np.rad2deg(np.arctan2(vlimb['value'][1],vlimb['value'][0]))
    #test which limb of v2 has the same sign of the x-coordinate of v1
    #we don't want to look on the wrong side of the Sun!
    #print np.sign(v2E[0]), np.sign(v1[0])
    #if np.sign(v2E[0]) == np.sign(v1[0]):
    #    vlimb={'value':v2E, 'side':'E'}
        #print 'v2 east limb farther'
    #else:
    #    vlimb={'value':v2W, 'side':'W'}
    #    #print 'v2 west limb farther'
    #vElimb={'value':v2E, 'side':'E'} #this is 0 in STEREO/2nd intstrument HGS frame
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
    joint_obs= False
    if type(flare_lon) == float: #calculate joint_obs. flare_lon is HGS_earth (deg)
        #tvE=np.rad2deg(np.arctan2(vElimb['value'][1],vElimb['value'][0]))
        #tvfov=np.rad2deg(np.arctan2(vfov['value'][1],vfov['value'][0]))
        #vEfov=np.abs(tvfov-tvE)
        #if vlimb['side'] == 'W': #if stereo_lon is between vfov and west limb it's observed
        #    tvW=np.rad2deg(np.arctan2(vlimb['value'][1],vlimb['value'][0]))
        #   vWfov=np.abs(tvfov-tvW)

        #convert all longitudes to something between 0 and 360? <-- this might be best. add 180 to everything > 0?, flip sign on everything else? try
        #if flare_lon >=0:
#            l360=180.+flare_lon
#         else:
#            l360=-1.*flare_lon
#         if tfov >=0:
#            fov360=180.+tfov
#         else:
#            fov360=-1.*tfov
#         if tvlimb >=0:
#            limb360=180.+tvlimb
#         else:
#            limb360=-1.*tvlimb
#         #fov360=tfov+180.
#         #limb360=tvlimb+180.

#         #find min and max based on which is the acute angle: angle between fov-limb or between limb-fov
#         if limb360-fov360 > fov360-limb360:
#             maxang=
#             minang=
#         else:
#             maxang=
#             minang=

        #print 'joint_obs?', flare_lon, np.max([tfov,tvlimb]),np.min([tfov,tvlimb]),l360,np.max([fov360,limb360]),np.min([fov360,limb360])
        #if l360 <= maxang and l360 >= minang: #what if it's visible from the other side though? if the farthest limb is the west limb?
        vflare=[np.cos(np.deg2rad(flare_lon)),np.sin(np.deg2rad(flare_lon))]
        #print vflare,'vflare'
        crossab=np.cross(vlimb['value'],vflare)
        crossac=np.cross(vlimb['value'],vfov['value'])
        #print 'cross'
        crossca=np.cross(vfov['value'],vlimb['value'])
        crosscb=np.cross(vfov['value'],vflare)
        #print 'cross print ', crossab*crossac,crossca*crosscb
        if crossab*crossac >=0. and crossca*crosscb >=0:
            joint_obs = True

        #else:
        #    if stereo_lon <=vEfov: #what if it's visible from the other side though? if the farthest limb is the west limb?
        #        joint_obs = True
        #    else:
        #        joint_obs = False

    #print 'check'

    if plot:
        fig,ax=plt.subplots()
        ax.scatter(v1[1],v1[0], c='r',marker='o',label='v1',s=50)
        ax.scatter(0,0, c='b',marker='o',label='Sun',s=150)
        ax.scatter(v2[1],v2[0] , c='g',marker='o',label='v2',s=50)
        ax.scatter(v1E[1],v1E[0],  c='r',marker='*',label='v1E')
        ax.scatter(v1W[1],v1W[0],  c='r',marker='v',label='v1W')
        ax.scatter(v2E[1],v2E[0],  c='g',marker='*',label='v2E')
        ax.scatter(v2W[1],v2W[0],  c='g',marker='v',label='v2W')
        ax.plot([0,v1[1]],[0,v1[0]],'-r')
        ax.plot([0,v2[1]],[0,v2[0]],'-g')
        ax.plot([0,v1E[1]],[0,v1E[0]],'--r')
        ax.plot([0,v1W[1]],[0,v1W[0]],'-.r')
        ax.plot([0,v2E[1]],[0,v2E[0]],'--g')
        ax.plot([0,v2W[1]],[0,v2W[0]],'-.g')
        ax.plot([0,v2W[1]],[0,v2W[0]],'-.g')
        ax.plot([0,0],[0,1],'k')
        if type(flare_lon) == float:
            #plot flare vector...
            ax.plot([0,np.sin(np.deg2rad(flare_lon))],[0,np.cos(np.deg2rad(flare_lon))],'b', label='flare')
        #shade between vfov and vlimb
        ax.fill([vfov['value'][1],vlimb['value'][1],0.],[vfov['value'][0],vlimb['value'][0],0.],alpha=0.5) #x-points, y-points, always include [0,0]
        ax.set_ylabel('x_HEE')
        ax.set_xlabel('y_HEE')
        ax.set_xlim([-1.1,1.1])
        ax.set_ylim([1.1,-1.1])
        ax.set_aspect('equal')
        ax.set_title(str(joint_obs))
        ax.legend()
        fig.show()


    return theta_obs,theta_ex,vlimb,'vElimb',vfov,joint_obs

def stereo_loc_to_lon(stereo_loc):
    '''convert location coordinates to heliographic stonyhurst and estimate the longitudinal angle from here'''
    aa=stereo_loc.transform_to('heliographic_stonyhurst')
    #obs_lon=stero_loc.observer.lon.value #negative because it's pointing backwards at the observer (+- 2 degrees off?)
    slon = np.abs(aa.lon.value) #east limb is zero, west is 180
    return slon

def is_flare_observed(stereo_lon,vElimb,vfov,vlimb,plot=False,earthframe=True):
    '''stereo_lon is measured from the east limb. What is the angle from vfov to the east limb? is stereo-lon included in this angle?'''
    #stereo_lon is in coordinate frame where vElimb ==0
    #OR it's in the frame where the Earth-Sun axis is 90 degrees. Need to be consistent!
    #make this the default
    #vElimb, vfov,vlimb are all in HEE
    #earthframe East limb in HEE: [x=0,y=1] #just assuming +1 because it would be pretty unsusual for a def. with a cosine in it to start negative?
    #so if input is in HGS earth then the HEE east limb is at (0,1)
    #easiest to compare in radial coords
    #problem is how exactly to convert to HEE? same as HEEQ just with a different x-axis definition? idk for sure
    #so maybe stay Cartesian
    #Or use sunpy transforms and trust that they work...
    #if type(stereo_lon) != SkyCoord:
        #convert skycorod
    #else:
        #check frame
        #convert to HEE from HGS_earth

    #now convert everything else? but first everything else needs to be KEPT in skycoords
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

    if plot:
        fig,ax=plt.subplots() #just plot londgitudes now?
        ax.scatter(vElimb,0, c='r',marker='o',label='vElimb',s=50)
        ax.scatter(vElim + vfov,0, c='b',marker='o',label='vWlimb',s=50)
        ax.scatter(stereo_lon,0 , c='g',marker='o',label='v_flare',s=50)
        #ax.scatter(s1E[0],s1E[1],  c='r',marker='*',label='v1E')
        #ax.scatter(s1W[0],s1W[1],  c='r',marker='v',label='v1W')
        #ax.scatter(s2E[0],s2E[1],  c='g',marker='*',label='v2E')
        #ax.scatter(s2W[0],s2W[1],  c='g',marker='v',label='v2W')
        #ax.plot([0,vec1[0]],[0,vec1[1]],'-r')
        #ax.plot([0,vec2[0]],[0,vec2[1]],'-g')
        #ax.plot([0,s1E[0]],[0,s1E[1]],'--r')
        #ax.plot([0,s1W[0]],[0,s1W[1]],'-.r')
        #ax.plot([0,s2E[0]],[0,s2E[1]],'--g')
        #ax.plot([0,s2W[0]],[0,s2W[1]],'-.g')
        #ax.plot([0,s2W[0]],[0,s2W[1]],'-.g')
        ax.set_xlabel('x_HGS (stereo)')
        ax.set_ylabel('y_HGS (stereo)')
        ax.set_aspect('equal')
        ax.legend()
        fig.show()


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
