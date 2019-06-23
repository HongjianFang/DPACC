#/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : main.py
#
# Purpose :
#
# Creation Date : 22-06-2019
#
# Last Modified : Sun Jun 23 19:16:48 2019
#
# Created By : Hongjian Fang: hfang@mit.edu 
#
#_._._._._._._._._._._._._._._._._._._._._.*/
from collections import defaultdict
from scipy import interpolate
import random
import pickle
import numpy as np
import obspy
import distaz
from matplotlib import pyplot as plt
from obspy.core.util import AttribDict
import os
from scipy.fftpack import fft,ifft
import glob
from scipy.signal import hilbert
#from obspy.taup import TauPyModel
import pandas as pd
#from obspy.core.stream import Stream
from obspy.signal.filter import bandpass
import yaml
import cluster
#import acc
#import single_event
#import plot_event

def autocorr_td(tr,stalign_a,stalign_b,conlen=15):
        N = tr.stats.npts
        delta = tr.stats.delta
        yf = fft(tr.data)
        ayf = np.abs(yf)
        conlen = 15
        myf = np.convolve(ayf, np.ones((conlen,))/conlen, mode='same')
        yf = yf/myf
        xx = np.real(ifft(yf))
        
        stalign_a = int(stalign_a/delta)
        stalign_b = int(stalign_b/delta)
        winlen = stalign_b - stalign_a
        xxw = xx[stalign_a:stalign_b]
        acorr = np.correlate(xx,xxw,mode='full')
        acorr = acorr[winlen:]
        maxloc = np.argmax(abs(acorr))
        acorr = np.roll(acorr,-maxloc)
        
        tr.data = acorr
        return tr

def autocorr_fd(tr,conlen=15):
        N = tr.stats.npts
        delta = tr.stats.delta
        data = np.zeros(2*N,)
        data[:N] = tr.data
        yf = fft(data)
        ayf = np.abs(yf)
        #conlen = 10
        myf = np.convolve(ayf, np.ones((conlen,))/conlen, mode='same')
        yf = yf/myf
        myf = yf*np.conj(yf)
        xx = np.real(ifft(myf))
        tr.data = xx[:N]
        return tr

def single_event(ievent):
    nsta = len(stalist)
    stapos = np.zeros((nsta,2))
    idx = 0
    staposcl = np.zeros((nsta,2))
    for ista in range(nsta):
        stanet = stalist[ista].split('.')[0]
        staname = stalist[ista].split('.')[1]
        stasub = irissta[(irissta['net']==stanet) & (irissta['sta']==staname)]
        stlat = stasub.iloc[0]['lat']
        stlon = stasub.iloc[0]['lon']
        stapos[ista,0] = stlat
        stapos[ista,1] = stlon

        
    strmacc1 = obspy.core.stream.Stream()
    strmori = obspy.core.stream.Stream()

    if evdep1<30.0:
        print ( evdep1 )
        return strmacc1,strmori
          

    idx = np.arange(nsta)
    for ista in range(len(idx)):
                stlat = stapos[idx[ista],0]
                stlon = stapos[idx[ista],1]
                dis1 = distaz.DistAz(evlat1,evlon1,stlat,stlon)

                if dis1.delta<mindist or dis1.delta>maxdist:
                    continue
                #return strmacc1,strmori

                trace = '/'.join([eqdir,evname1,datatype,stalist[idx[ista]]])

                strm = obspy.read(trace)
                tr = strm[0]
                if tr.stats.npts < 50:
                    continue
                #return strmacc1,strmori
                tr.stats.coordinates=AttribDict({'latitude':stlat,'longitude':stlon,'elevation':0})
                

                if phase == 'P':
                    parr = ftelep(dis1.delta,evdep1)[0]  
                else:
                    parr = fteles(dis1.delta,evdep1)[0]
                tr.stats.distance = dis1.delta#inc_angle
                tr.stats.baz = dis1.baz
                tr.trim(evtime1+parr-trimb,evtime1+parr+trima,pad=True,fill_value=0)
                tr.resample(rsample)
                npts = tr.stats.npts
                tr.stats.starttime = 0
                tr.detrend()
                tr.taper(tpratio)
                trdata = np.zeros(npts+20*rsample,)
                trdata[10*rsample:-10*rsample] = tr.data
                trdata = bandpass(trdata,frqmin,frqmax,tr.stats.sampling_rate,2,True)
                tr.data = trdata[10*rsample:-10*rsample]
                tr.normalize()
                
                if envolope:
                    tr.data = np.abs(hilbert(tr.data))
                strmori.append(tr.copy())
                
                
                if timedomain:
                    tr = autocorr_td(tr,sig_bs,sig_es,conlen=conlen)
                else:
                    tr = autocorr_fd(tr,phaseshift=phishift)


                npts = tr.stats.npts
                tr.taper(tpratio)
                trdata = np.zeros(npts+20*rsample,)
                trdata[10*rsample:-10*rsample] = tr.data
                trdata = bandpass(trdata,frqmin,frqmax,tr.stats.sampling_rate,2,True)
                tr.data = trdata[10*rsample:-10*rsample]
                tr.normalize()
                if np.isnan(tr.data).any():
                    #return strmacc1,strmori
                    continue

                strmacc1.append(tr.copy())
                
    print ('finishing acc')
    return strmacc1,strmori

def plot_event(strmacc1,strmori):
    datasave_sub = defaultdict(list)
    strmstack = obspy.core.stream.Stream()
    accstack = obspy.core.stream.Stream()
    tr = strmacc1[0].copy()
    reftime0 = fppdp(refdismax,evdep1+mindep)
    reftime1 = fppdp(refdismax,evdep1+maxdep)
    reftime0s = fspdp(refdismax,evdep1+mindep)
    reftime1s = fspdp(refdismax,evdep1+maxdep)
       
    if len(strmacc1) < 10:
        print ('small no. of traces:',len(strmacc1))
        return
    staposcl = np.zeros((len(strmacc1),2))
    for idx in range(len(strmacc1)):
        staposcl[idx,0] = strmacc1[idx].stats.coordinates['latitude']
        staposcl[idx,1] = strmacc1[idx].stats.coordinates['longitude']
        
    stidx,aslat,aslon = cluster.clustersta(staposcl[:,0],staposcl[:,1])
    uidx,ucounts = np.unique(stidx,return_counts=True)
    idx = np.where(ucounts>5)[0]
    useidx = int(len(idx)*0.8)
    
    npts = strmacc1[0].stats.npts
    stacksub = np.zeros((npts,useidx))
    
    starttime = strmacc1[0].stats.starttime
    stacklinear = np.zeros((npts,))

    buall = np.zeros(npts,)
    phi = np.zeros((npts,),dtype=complex)
    idx = random.sample(idx,useidx)
    reftime = fppdp(refdismax,evdep1)
    for ii in idx:
        bu = np.zeros(npts,)
        data = np.zeros(npts,)
        cellidx = np.where(stidx==uidx[ii])[0]
        reflat = np.rad2deg(aslat[cellidx[0]])
        reflon = np.rad2deg(aslon[cellidx[0]])
        disref = distaz.DistAz(evlat1,evlon1,reflat,reflon)
        sidx = random.sample(cellidx,1)[0]
        refdis = strmacc1[sidx].stats.distance
        ctime = fppdp(refdis,evdep1)
        ctime = reftime-ctime
        data = np.roll(strmacc1[sidx].data[:npts],int(np.round(ctime*rsample)))
        phi = phi+np.exp(1j*np.angle(hilbert(data)))
        buall = buall+data
        for jj in range(len(cellidx)):
            refdis = strmacc1[cellidx[jj]].stats.distance
            ctime = fppdp(refdis,evdep1)
            ctime = reftime-ctime
            data = np.roll(strmacc1[cellidx[jj]].data[:npts],int(np.round(ctime*rsample)))
            phi = phi+np.exp(1j*np.angle(hilbert(data)))
            bu = bu+data
        bu = bu/np.max(abs(bu))
        tr.data = bu
        tr.stats.distance = disref.delta
        tr.stats.baz = disref.baz
        tr.normalize()
        accstack.append(tr.copy())
    nsta = len(idx)
    stacklinear = buall/np.max(abs(buall))   
    if pws:
        aphi = (np.abs(phi)/nsta)**mu
        buall = bu*aphi
        stacklinear = buall/np.max(abs(buall))*aphi

    if len(strmori) < 20:
        print ('small no. of traces:',len(strmori))
        #return
    staposcl = np.zeros((len(strmori),2))
    for idx in range(len(strmori)):
        staposcl[idx,0] = strmori[idx].stats.coordinates['latitude']
        staposcl[idx,1] = strmori[idx].stats.coordinates['longitude']
    
    stidx,aslat,aslon = cluster.clustersta(staposcl[:,0],staposcl[:,1],ncell=8000)
    uidx,ucounts = np.unique(stidx,return_counts=True)
    idx = np.where(ucounts>=5)[0]
    if len(idx) > 20:
        idx = random.sample(range(len(idx)),20)
    npts = strmori[0].stats.npts
    for ii in idx:
        data = np.zeros(npts,)
        cellidx = np.where(stidx==uidx[ii])[0]
        reflat = np.rad2deg(aslat[cellidx[0]])
        reflon = np.rad2deg(aslon[cellidx[0]])
        disref = distaz.DistAz(evlat1,evlon1,reflat,reflon)
        bu = np.zeros(npts,)
        for icell in cellidx:
            data = strmori[icell].data
            bu = bu+data[:npts]
        nsta = len(cellidx)
        aphi = (np.abs(phi)/nsta)**mu
        bu = bu#*aphi
        tr.data = bu
        tr.stats.distance = disref.delta
        tr.stats.baz = disref.baz
        tr.normalize()
        strmstack.append(tr.copy())
    trec = np.linspace(-trimb,trima,npts) 
    fig = plt.figure(figsize=(15,20))
    ax1 = fig.add_subplot(511)
    bg = 0
    ed = np.min([int(reftime1*rsample+100),npts])
    for ii in range(len(strmstack)):
        offset = strmstack[ii].stats.distance/10.0
        data = offset+strmstack[ii].data
        ax1.plot(trec[:endtime*rsample],data[:endtime*rsample])
    ax1.set_title('Origin waveforms (epi) '+str(ievent),fontsize=15)
    ax1.set_ylabel('Epi (x10)')

    ax2 = fig.add_subplot(512)
    naz = 4
    aziint = 360.0/naz
    for ii in range(len(strmstack)):
        offset = strmstack[ii].stats.baz/aziint
        data = offset+strmstack[ii].data
        ax2.plot(trec[:endtime*rsample],data[:endtime*rsample])
    ax2.set_title('Origin waveforms (azi)')
    ax2.set_ylabel('Azimuth (x90)')
    

    ddep = np.linspace(mindep,maxdep,ndep)
    nsta = len(strmacc1)
    npts = strmacc1[0].stats.npts
    depstack = np.zeros(npts-200,)
    

    stack = np.zeros((ndep,npts))
        
    evdep2 = evdep1#+ddep[idep]
    reftime = fppdp(refdismax,evdep2)
    data = strmacc1[0].data.copy()

    bu = np.zeros(npts,)
    phi = np.zeros((npts,),dtype=complex)
    for ii in range(0,nsta):
        refdis = strmacc1[ii].stats.distance
        ctime = fppdp(refdis,evdep2)
        ctime = reftime-ctime
        data = np.roll(strmacc1[ii].data[:npts],int(np.round(ctime*rsample)))
        phi = phi+np.exp(1j*np.angle(hilbert(data)))
        bu = bu+data

    aphi = (np.abs(phi)/nsta)**mu
    bu = bu#*aphi
    if pws:
        stack = bu*aphi
    else:
        stack = bu    
    
    depc = np.arange(mindep,maxdep,5)
    ndepc = len(depc)
    stackazi = np.zeros((npts,ndepc,naz))
    phiazi = np.zeros((npts,ndepc,naz),dtype=complex)
    nstaazi = np.zeros((ndepc,naz))
    for jj in range(ndepc):
        evdepc = evdep1+depc[jj]
        reftime = fppdp(refdismax,evdepc)
        for ii in range(len(strmacc1)):
            bu = np.zeros(npts,)
            phi = np.zeros((npts,),dtype=complex)        
            data = np.zeros(npts,)
            refdis = strmacc1[ii].stats.distance
            ctime = fppdp(refdis,evdepc)
            ctime = reftime-ctime
            data = np.roll(strmacc1[ii].data[:npts],int(np.round(ctime*rsample)))
            phi = np.exp(1j*np.angle(hilbert(data)))        
            stbaz = strmacc1[ii].stats.baz
            idx = int(stbaz/aziint)
            stackazi[:,jj,idx] += data
            phiazi[:,jj,idx] += phi
            nstaazi[jj,idx] += 1

    bg = np.max([int(reftime0*rsample-100),0])
    ed = np.min([int(reftime1s*rsample+100),npts])
    offset = 0
    ttime = np.linspace(bg/rsample,ed/rsample,ed-bg)
    
    axacc = fig.add_subplot(513)
    for ii in range(len(accstack)):
        offset = accstack[ii].stats.distance/10.0
        data = offset+accstack[ii].data
        axacc.plot(ttime,data[bg:ed])
    axacc.set_title('ACC (epi) '+str(ievent),fontsize=15)
    axacc.set_ylabel('Epi (x10)')
    
    ax3 = fig.add_subplot(514)
    stackp = stacklinear[bg:ed]
    stackp = stack[bg:ed]

    ax3.plot(ttime,stackp/np.max(abs(stackp))*2,'r-',linewidth=3.0)
    ax3.set_title('ACC stacking')
    
    for ii in range(naz):
        for jj in range(ndepc):
            offset = 3+ii*2#nboots*2+
            if pws:
                aphi = (np.abs(phiazi[:,jj,ii])/nstaazi[jj,ii])**mu
                stackazi[:,jj,ii] = stackazi[:,jj,ii]*aphi
            data = offset+stackazi[:,jj,ii]/np.max(abs(stackazi[:,jj,ii]))*2
            ax3.plot(ttime,data[bg:ed])
        ax3.text(bg/np.real(rsample),offset,str(int(nstaazi[jj,ii])),fontsize=15)
        
    ax3.axvline(reftime0,linestyle='--',color='magenta')
    ax3.axvline(reftime1,linestyle='--',color='magenta')
    ax3.axvline(reftime0s,linestyle='--',color='blue')
    ax3.axvline(reftime1s,linestyle='--',color='blue')


    datasave_sub['evid'] = evname1
    datasave_sub['data'] = stacklinear
    
    datasave_sub['pP0'] = reftime0*rsample
    datasave_sub['pP1'] = reftime1*rsample
    
    depstack = np.zeros(ndep,)
    data = np.zeros((naz,npts))
    for ii in range(naz):
        for jj in range(ndepc):        
            data[ii,:] += stackazi[:,jj,ii]/np.max(abs(stackazi[:,jj,ii]))
    stack = np.zeros(npts,)
    for ii in range(naz):
        if nstaazi[0,ii] >5:
            stack += np.abs(hilbert(data[ii,:]*nstaazi[0,ii]/nsta))
    datasave_sub['stack'] = stack
    ax4 = fig.add_subplot(515)
    for idep in range(ndep):
        evdep2 = evdep1+ddep[idep]
        pPtime = fppdp(refdismax,evdep2)[0]
        sPtime = fspdp(refdismax,evdep2)[0]
    
        depstack[idep] = 0.6*stack[int(sPtime*rsample)]**2+0.4*stack[int(pPtime*rsample)]**2
    ax4.plot(ddep,depstack,'rd-')
    ax4.set_xlabel('Depth')
    ymax = np.max([depstack.max(),10.0])
    ax4.set_ylim([-0.1,ymax])
    ax4.set_title('Mapping to depth_'+str(evdep1))
    fig.savefig(figdir+'/'+evname1+'_Mw'+str(evmag1)+'_td.png',dpi=200)
    plt.ioff()
    plt.close()
    datasave_sub['stackenv'] = stack
    datasave_sub['sdep'] = np.vstack([ddep,depstack])
    datasave_sub['evlat'] = evlat1
    datasave_sub['evlon'] = evlon1
    datasave_sub['evdep_cat'] = evdep1
    datasave_sub['evdep_inv'] = evdep1+ddep[np.argmax(depstack)]
    datasave_sub['dist'] = refdismax
    datasave_sub['saveid'] = ievent
    return datasave_sub

with open('parameters.in','r') as fin:
    par = yaml.load(fin)

trimb = par['trimb']
trima = par['trima']
frqmin = par['frqmin']
frqmax = par['frqmax'] 
rsample = par['rsample']
endtime = par['endtime']
tpratio = par['tpratio']
sig_bs = par['sig_bs']
sig_es = par['sig_es']
envolope = par['envolope']
timedomain = par['timedomain']
ndep = par['ndep']
mu = par['mu']
mindep = par['mindep']
maxdep = par['maxdep']
mindist = par['mindist']
maxdist = par['maxdist']
phase = par['phase'] 
comp = par['comp']
refdismax = par['refdismax']
pws = par['pws']
conlen = par['conlen'] 
datatype = par['datatype']
figdir = par['figdir']
eqdir = par['eqdir']

ndep = 301
ndis = 131
dep = np.linspace(30,330,ndep)
dis = np.linspace(30,95,ndis)
pPdP = np.load('./tables/pPdP30to330.npy')
fppdp = interpolate.interp2d(dis, dep, pPdP, kind='linear')

sPdP = np.load('./tables/sPdP30to330.npy')
fspdp = interpolate.interp2d(dis, dep, sPdP, kind='linear')

pPdP = np.load('./tables/tele30to330.npy')
ftelep = interpolate.interp2d(dis, dep, pPdP, kind='linear')

eqs = open(''.join([eqdir,'/EVENTS-INFO/event_list_pickle']))
eqs = pickle.load(eqs)
evid = []

for ii in eqs:
    evid.append(ii['event_id'])

if not os.path.isdir(figdir):
    os.makedirs(figdir)

print ('begin stacking for each station')
npoints = int(endtime*rsample+1)
#global evelist
evelist = glob.glob(eqdir+'/*a')
nevent = len(evelist) 
irissta = pd.read_table('./tables/IRISSTA9019.txt',names=('net','sta','lat','lon'),header=0,delim_whitespace=True,keep_default_na=False)

evdepinv = np.zeros(nevent)
evdep = np.zeros(nevent)
evlon = np.zeros(nevent)
sampleevent = nevent
datasave = []

for ievent in range(nevent):
    ievent1 = ievent
    evname1 = evelist[ievent1].split('/')[-1]
    evidx1 = evid.index(evname1)
    evlat1 = eqs[evidx1]['latitude']
    evlon1 = eqs[evidx1]['longitude']
    evdep1 = eqs[evidx1]['depth']
    evtime1 = eqs[evidx1]['datetime']
    evmag1 = eqs[evidx1]['magnitude']  
    stalist1 = glob.glob('/'.join([eqdir,evname1,datatype,'*.'+comp]))
    stalist = [stal.split('/')[-1] for stal in stalist1]
    print ('the',ievent,'th event of all ',nevent)
    print ('evinfo for eq:',evname1,evmag1,evlat1,evlon1,evdep1,len(stalist))
 
    strmacc1,strmori = single_event(ievent)
    datasavesub = plot_event(strmacc1,strmori)
    datasave.append(datasavesub)

if __name__ == '__main__':
    with open(eqdir+'_'+comp+'.pkl', 'wb') as fout:
        pickle.dump(datasave,fout)
    print('finshing all')
