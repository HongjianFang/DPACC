#/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : cluster.py
#
# Purpose :
#
# Creation Date : 22-06-2019
#
# Last Modified : Sun Jun 23 18:58:50 2019
#
# Created By : Hongjian Fang: hfang@mit.edu 
#
#_._._._._._._._._._._._._._._._._._._._._.*/

def clustersta(stlat,stlon,ncell=5000):
    import distance
    import numpy as np
    nsta = len(stlat)
    rlat = np.random.rand(ncell,)
    rlat = np.arccos(2*rlat-1)-np.pi/2
    rlon = np.random.rand(ncell,)
    rlon = 2*np.pi*rlon

    stlon[stlon<0]+=360
    stidx = np.zeros(nsta,)
    aslat = np.zeros(nsta,)
    aslon = np.zeros(nsta,)

    for ii in range(nsta):
        dis = distance.distance_mesh(np.radians(stlat[ii]),np.radians(stlon[ii]),rlat,rlon)
        idx = np.argmin(dis)
        stidx[ii] = idx
        aslat[ii] = rlat[idx]
        aslon[ii] = rlon[idx]
    return stidx,aslat,aslon
