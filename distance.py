#/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : distance.py
#
# Purpose :
#
# Creation Date : 15-01-2019
#
# Last Modified : Tue Jan 15 18:06:38 2019
#
# Created By : Hongjian Fang: hfang@mit.edu 
#
#_._._._._._._._._._._._._._._._._._._._._.*/
def distance_mesh(lat1,lon1,lat2,lon2):
      import numpy as np
      R = 6371

      dLat = lat2-lat1
      dLon = lon2-lon1
      a = np.sin(dLat/2) * np.sin(dLat/2) +  np.sin(dLon/2) * np.sin(dLon/2) * np.cos(lat1) * np.cos(lat2)
      c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
      d = R * c
      return d
