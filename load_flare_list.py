# IPython log file

execfile('code/OCFlareList.py')
fl=OCFlareList([10110321,12062521,12071911,12091702,13042736])
fl.list[0].Files.Raw['stereo']=['20101103_124530_n4eub.fts','20101103_124530_n4eub.fts']
fl.list[1].Files.Raw['stereo']=['20120625_182030_n4eub.fts','20120625_182030_n4eub.fts']
fl.list[2].Files.Raw['stereo']=['20120719_052030_n4eub.fts','20120719_052030_n4eub.fts']
fl.list[3].Files.Raw['stereo']=['20120917_062030_n4eub.fts','20120917_062030_n4eub.fts']
fl.list[4].Files.Raw['stereo']=['20130427_105530_n4eua.fts','20130427_105530_n4eua.fts']
fl.list[0].Properties.source_pos_STEREO=[-218,-393]
fl.list[1].Properties.source_pos_STEREO=[190,-125]
fl.list[2].Properties.source_pos_STEREO=[-400,-143]
fl.list[3].Properties.source_pos_STEREO=[174,200]
fl.list[4].Properties.source_pos_STEREO=[-634,101]

#add rhessi source and limb positions to calculate the distance
fl.list[0].Properties.e3={'low':[-912.6693, -326.8458],'high':[-913.1742, -327.1038
],'limb':[-909.9458, -326.3011]}
fl.list[1].Properties.e3={'low':[-951.9592, -258.1428],'high':[-926.2449, -246.5102],'limb':[-909.1021, -250.1837]}
fl.list[2].Properties.e3={'low':[971.6326, -236.0204],'high':[982.0408, -237.857],'limb':[915.9184, -228.6735]}
fl.list[3].Properties.e3={'low':[-969.2653, 254.3061],'high':[],'limb':[-922.7347, 242.0612]}
fl.list[4].Properties.e3={'low':[966.1851, 206.9305],'high':[],'limb':[931.0151, 201.9630]}
import astropy.units as u
import scipy
#for f in fl.list:
#    lc=f.plot_stereo(AIAmap=False,maps=False,oplotlimb=True,subplots=False,zoom=True,save=False)
#    pos=f.get_closest_limbcoord(limbcoords=lc,redef=False)
#    f.calc_projection_effect()
#    dlow= scipy.spatial.distance.euclidean(f.Properties.e3['low'],f.Properties.e3['limb'])
#    if f.Properties.e3['high'] !=[]:
#        dhigh= scipy.spatial.distance.euclidean(f.Properties.e3['high'],f.Properties.e3['limb'])
#    else:
#        dhigh=0.0
#    #print 'dlow ',dlow*715.,'dhigh ',dhigh*715., 'avg ', np.mean([dlow,dhigh])*715.
#    f.Properties.e4={'dlow':dlow*715,'dhigh':dhigh*715,'avg':np.mean([dlow,dhigh])*715}
#    theight=f.Properties.e5['hkm']+f.Properties.e4['dlow']
#    theighta=f.Properties.e5['hkma']+f.Properties.e4['dlow']
#    print 'total height plus projection effect: ',theighta
#    f.Properties.e5={'theight':theight,'theighta':theighta}
    
   
#fl.list[1].Properties.source_pos_STEREO=[180,-135]
