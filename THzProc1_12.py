# THzProc 1.12

# process THz data in tvl format with various funtionalities

# 08MAR2013: change in processed files only
# 14APR2013: start version number
# 25AUG2013: fix remap on "middle points"
# 02SEP013: add peak flagging
# 17SEP2013: add depth/thickness profile line and 3D (from some parts of proc03.py)
# 19SEP2013: modified C:\Python27\Lib\site-packages\mpl_toolkits\mplot3d\axes3d.py to fix aspect ratio in 3D plot
#            per Ben Root implementation: see web pages about fixing aspect ratio in C:\f1\current\THzProc\doc,paper
# 19SEP2013: add Mayavi rendering of 3D surface plot - much faster! But all Matplotlib functionalities stop.  Furture
#            work to embed both in WxPython.
# 20SEP2013: changed font factor from default 1.5 to 0.8 (this reduces the font size) in C:\Python27\Lib\site-packages\mayavi\modules\axes.py
# 22SEP2013: fix more of remap and np.max -> np.amax
# 18NOV2013; remove baseline trend of internal reflection in near-field imaging
# 21DEC2013: add AscanOnly
# 03APR2014: add matched filter in InterAct
#   JUL2014: revise paramter settings (e.g. use list to extend gate numbers) and 3D visual sections for displaying interior of Army clear glass panel
# 26SEP2014: change to avg abs amp when SigType=2 and add pixel output in PickPt
# 27JAN2015: add handle for FSE on edge
# 08JUN2015: remove excessive 300ps waveform end amplification; scale B-scan for display 
# 20DEC2015: correct SigType>1 prog error
# 27DEC2015: allow c-scan contrast change ; turn Bscan on/off
# 28DEC2015: separate out MakeCscan base on SigType choices ; reorganize THzProc 1.10.2 into 1.11
# 29DEC2015: fix issues with PulseLen; start rewriting slicer
# 01JAN2016: add more baseline removal option
# 17APR2016: show grid (i,j) location at click on C-scan

# last update: 17APR2016

from __future__ import with_statement # this line has to be the first of all headers
import numpy as np
from scipy.fftpack import rfft,irfft
from skimage import filters
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.widgets import Cursor
from matplotlib.figure import figaspect
from mpl_toolkits.mplot3d import Axes3D
# from mayavi import mlab # "from enthought.mayavi import mlab"  moducle name changed since Mayavi 4.4.x
import struct, os, sys, copy
import pdb

# parameters to change with data file ==================================================================================================

# Follower gate is off with FollowGateOn=0.  In this case, max. pos. peak in the entire waveform is located first, 
# followed by the max. neg. peak located within +-HalfPulse. Vpp is then obtained from these two peaks and used in the C-scan.
# When follower gate is on with FollowGateOn>0, the reference peaks and Vpp are determined within the fixed gate (1st element in global 
# list gate).  Follower gates (2nd and up elements in list gate) are related to max. peak (for now) in the fixed gate.  FollowGateOn=k 
# denote the k-th follower gate (k-th+1 element in list gate) usually the k-th interface below front surface.

# Additional amplitude compensation also applies when AmpCompOn=1

# if FSEtol>0 then throw away any location whose FSE<FSEtol; set FSEtol<0 to disable this

# get 3D front surface or layer thickness profile if DepthMapOn!=0.
# FollowGatOn=0: find front surface profile; =1: find thickness between FS and 1st interface; =2: between 1st and 2nd interfaces, etc.
# DepthMapOn=1: "hole to stick out"; =-1: "hole to sink in" as it looks like           
# default uses Matplotlib routines and is much slower. If MayaviOn=1 also, then uses Mayavi to render it.
# PickPeak=1: pos peak, =2: neg. peak, =3: mid-point  (prefer 1; 3 gets errors?!)

# Sigtype=1: CscanAmpNew=Vpp within the gate
#         2: avg mag (mag = abs amp)
#         3: median mag
#         4: min mag
#         5: max mag
#         6: avg amp ~ integrated gate
#         7: median amp
#         8: min amp 
#         9: max amp

# MfOn=1 to do matched filer for IU/CRC data processing project

# dual band histogram equalization is used when VppGate=1 and settings in VppL1, VppH1, VppL2 and VppH2

# removal of baseline trend of internal reflection (near field only): 
# TrendOff= 5 - average several "blank" waveforms (i.e. not on the sample) at beginning and end of row as reference to be substracted from the sample waveforms
#               assume there are at least some (1mm or more preferred) "blank" at beginning and end of row and X spacing is sufficiently smaller (0.05mm preferred)
# TrendOff= 1 - use one waveform of each row as reference; aligned each waveform in that row with the reference at mid pt between max. and min. peaks
# TrendOff=-1 - 5 plus use one waveform of each row as reference; aligned each waveform in that row with the reference at mid pt between max. and min. peaks
#         = 2 - samea as 1 plus rescaling refwav to have same Vpp as each waveform
#         = 3 - use a external waveform file as reference; the path of the reference waveform file must be given as refname
#         = 4 - samea as 3 plus rescaling refwav to have same Vpp as each waveform

# AscanOnly=1 - only plot A-scan; no spectra plots

# PulseLen<0 enforces searcg of neg. peak within original gate

# below about slicer is out of date!!!
# when SliceOn=1, begt and endt are used as begining and ending time bins for the global time gates (within which slicing takes place)
# when FollowGateOn=1, slices also track the max. peak in follower gate
# Slhw=half time width per slice in bin unit ; Sldel=spacing between slice gate in bin unit ; SlHisEquOn=1 applies histogram equalization
# within each slice

# For 2nd gate and up, PulseLen allows variable pulse length.  if PulseLen<0, enforce searcg of neg. peak within original gate

# =======================================================================================================================================

basedir = 'F:\\RR 2017\\THz Data\\Black CMC'
filename = 'Sample 1 (res=0.25mm OD=60ps).tvl'
FollowGateOn=1
FSEtol=-0.17 ; FlagPeakOn=1 ; DepthMapOn=0 ; PeakPick=1; MayaviOn=0 ; TrendOff=0 ; AscanOnly=1
SigType=1 ; BscanOn=1 ; BscanDir=0
gate=[[100,1000],[2100, 2700]]
workf=3. ; AmpCompOn=0 ; MfOn=0 ; DepOut=0 ; skip=10 ; DepAmpTol=0.5 # threshold for depth/thickness: 0.2 for output, 0.5 for visual
AmpCor300On=1 ; AmpCor300Par=[0.,35.,5.0,1.,240.,300.,1.,4.5,4.]
Zminfac=1.0 ; Zmaxfac=1.0 ; VppGateOn=0 ; VppL1=0. ; VppH1=0.58 ; VppL2=0.58 ; VppH2=4.
SliceOn=0 ; SlHisEquOn=1 ; Slhw=50 ; Sldel=20 ; fps=15
XCorTol=0.1 ; PulseLen=3.

# other more "universal" parameters ======================================================================================================
 
#XCorTol=-0.7 # correctioin to Backlash shift in mm
YDiffTol=0.3 ; XDiffTol=0.4  # tolerance (as ratio wrt scan resolution) allowed for a scan point to deviate from its desired coord. position
tiny=1.e-4 ; CscanEdgeFac=0.05 # edge ratio around a C-scan image

# half pulse width for bracketing the follower signal (usually the front surface echo) changed 14SEP2013
HalfPulse=2. #ps    reduced from 3. 24SEP2013

#threshold for lead gate signal
fthres=0.5

# tolerance ratio for number of scan points in a X line
nlim=0.8

# index of refraction of TBC: assumed = 5
IndRfnTBC=5.

# threshold for peak detecton in peakdet 
PeakTol=0.025

# -------------------------------------------------------------------------------------------------------------------------------------
# color palettes to use: 'gray', 'Blues','Greens','Reds','jet','hsv','rainbow' (add '_r' to reverse)

clrmap=['gray', 'gray_r','Blues_r','Reds_r','Greens','jet','hsv','rainbow']
marker='s' ; marklinwid=0.01 ; marksize=7 ; figh=8. #inch

LightSpeed=0.299 #mm/ps speed of light in vacuum

# plotting colors
pltclr=['r-','g-','b-','r--','g--','b--','c-','c--','m-','m--','y-','y--','k-','k--','r-.','g-.','b-.','r:','g:','b:']

# _______________________________________________________________________________________________________________________________________________________________________

def main():

  # last update: 08JUN015
  
  dat=DataFile(filename,basedir=basedir)
  # waveform length = col_size-2 (offset the X and y coords.)
  # the actual data waveforms come after the first three waveforms which contain ref, filter, etc.
  X=dat.data[3:]['x'] ; Y=dat.data[3:]['y'] ; waveform=dat.data['waveform'][3:]
  datlen=len(waveform) ; wavlen=len(waveform[0])
  
  wavtlen,Xres,Yres,Xstep,Ystep,Xmin,Xmax,Ymin,Ymax,subsamp,Navg,Xspeed,ScanType,axis1= HeaderInfo(dat.header)
  
  delt=wavtlen/(wavlen-1) ; delf=1./(wavlen*delt)
  time=t=np.linspace(0.,wavtlen,wavlen) ; freq=np.linspace(0.,(wavlen/2)*delf,wavlen/2)
  CscanAmp=np.zeros(datlen) ; nHalfPulse=int(HalfPulse/delt) ; print "nHalfPulse=",nHalfPulse
  
  print
  print 'asn len=',wavlen, 'asn Tlen=',wavtlen, 'delt=',delt, ' delf=',delf,'subsamp=',subsamp,' #avg=',Navg,' scan type=',ScanType
  print 'X min=',Xmin, ' max=',Xmax,'Y min=',Ymin, ' max=',Ymax,' scan step X=',Xstep,' Y=',Ystep,' res X=',Xres,' Y=',Yres
  TrueXres=(Xmax-Xmin)/float(Xstep-1) ; TrueYres=(Ymax-Ymin)/float(Ystep-1)  # these are what ScanAcquire actually use 25A2013
  print 'true res X=',TrueXres,' Y=',TrueYres
  print
  
  XCor=XCorTol*Xres ; YDiff=Yres*YDiffTol ; XChkTol=10
  
  WaveformNew,CscanAmpNew,Xnew,Ynew,pos,Xstep,Ystep=  \
  ReMap(waveform,X,Y,Xmax,Xmin,Ymin,TrueXres,TrueYres,Xstep,Ystep,ScanType,axis1,wavlen,wavtlen,XCorTol,YDiffTol,XChkTol,tiny)
  # BE AWARE: Xstep, Ystep may be updated after call to ReMap

  del waveform ; del X ; del Y # release memory

  # remove excessive amplitude amplification for 300 ps waveforms
  if AmpCor300On==1 and abs(wavtlen-300.)<tiny:
    AmpCor300(AmpCor300Par[0],AmpCor300Par[1],AmpCor300Par[2],AmpCor300Par[3],AmpCor300Par[4],AmpCor300Par[5],AmpCor300Par[6], \
              AmpCor300Par[7],AmpCor300Par[8],wavlen,delt,Xstep,Ystep,WaveformNew)
  
  if FollowGateOn>0:
    BinRange=copy.deepcopy(gate) # this way for sure makes BinRange a separate copy of gate
  elif FollowGateOn==0:  
    BinRange=[[0,wavlen]]
  PeakBin=FindPeaks(WaveformNew,Xstep,Ystep,wavlen,nHalfPulse,fthres,BinRange, FollowGateOn)

  if DepthMapOn!=0:
    DepThick=FindDepThick(DepthMapOn,PeakPick,WaveformNew,PeakBin,Xstep,Ystep,wavlen,delt,IndRfnTBC,DepAmpTol)
    X3d=np.linspace(Xmin,Xmax,Xstep) ; Y3d=np.linspace(Ymin,Ymax,Ystep)
    X3d,Y3d=np.meshgrid(X3d,Y3d)
    if DepOut==1:
      OutDep(basedir,skip,X3d,Y3d,DepThick)
  else:
    DepThick=0. ; X3d=0. ; Y3d=0. # dummies

  if MayaviOn==1 and DepthMapOn!=0:
    Mayavi3D(X3d,Y3d,DepThick)
  else:
    if SliceOn==1:
      pass
      #Slicer(datlen,wavlen,wavtlen,delt,delf,time,freq,Xstep,Ystep,axis1,WaveformNew,Xnew,Ynew,CscanAmpNew,  \
      #       workf,FollowGateOn,SlHisEquOn,Slhw,Sldel,fps,SigType) # do not use! out of date!!!!!
    else:
      InterAct(datlen,wavlen,wavtlen,delt,delf,time,freq,Xstep,Ystep,axis1,WaveformNew,Xnew,Ynew,CscanAmpNew,BscanDir,  \
               workf,FollowGateOn,VppGateOn,DepthMapOn,SigType,DepAmpTol,VppL1,VppH1,VppL2,VppH2, \
               X3d,Y3d,DepThick,PeakBin,nHalfPulse)

  #  plt.show() has been moved into InterAct

# _______________________________________________________________________________________________________________________________________________________________________

def InterAct(datlen,wavlen,wavtlen,delt,delf,time,freq,Xstep,Ystep,axis1,WaveformNew,Xnew,Ynew,CscanAmpNew,BscanDir, \
               workf,FollowGateOn,VppGateOn,DepthMapOn,SigType,AmpTol,VppL1,VppH1,VppL2,VppH2, \
               X3d,Y3d,DepThick,PeakBin,nHalfPulse):
  
# use FindPeaks
# add matched filter 03APR2014
# separate out MakeCscan 28DEC2015
# last update: 28DEC2015
  gatewavt=np.zeros(wavlen)
  maxfloc=np.zeros((Ystep,Xstep),dtype=np.int16)
  maxfVpp=np.zeros((Ystep,Xstep),dtype='<f') ; maxfVppAll=0.
    
# matched filter application for IU/CRC data processing project  
  mffile="AL foil & paper on back of 5mm glass fiber plate Focus@FS then 1mm down (100ps 4096pt 1avg no purg 0.5mm res IU SP project).tvl"
  if MfOn==1 and filename==mffile:  
    CscanAmpMf=MatchedFilter(delt,k,Waveform,Xstep,Ystep,PeakBin)
    Xmax=np.amax(Xnew) ; Xmin=np.min(Xnew) ; Ymax=np.amax(Ynew) ; Ymin=np.min(Ynew)
    fig7=plt.figure('Matched filtered C-scan',figsize=(figh*(Xmax-Xmin)/(Ymax-Ymin),figh)) # set width proportionally
    ax7=fig7.add_subplot(111)
    extent=(Xmin,Xmax,Ymax,Ymin)
    #im7=ax7.imshow(CscanAmpMf,interpolation='bilinear',cmap=cm.jet,origin='upper',aspect=1.,extent=extent,picker=True)
    im7=ax7.imshow(CscanAmpMf,interpolation='bilinear',cmap=cm.jet,origin='upper',vmin=0.,vmax=0.7e-10,aspect=1.,extent=extent,picker=True)
    ax7.set_ylim(Ymax+np.abs(Ymax)*CscanEdgeFac,Ymin-np.abs(Ymin)*CscanEdgeFac) 
    ax7.set_xlim(Xmin-np.abs(Xmin)*CscanEdgeFac,Xmax+np.abs(Xmax)*CscanEdgeFac)  
    if axis1=='Turntable':
      ax7.set_xlabel("Angular scan position (deg)")
    else:
      ax7.set_xlabel("X scan position (mm)")
    ax7.set_ylabel("Y scan position (mm)")
    plt.colorbar(im7,orientation='horizontal')
    ax7.grid(True)
  
 # make C-scan depending on the choice of SigType

  MakeCscan(WaveformNew,PeakBin,FollowGateOn,AmpCompOn,VppGateOn,Xstep,Ystep,SigType,FSEtol,VppL1,VppH1,VppL2,VppH2,CscanAmpNew,maxfVpp,maxfVppAll) 
  
  Zmax=np.amax(CscanAmpNew) ; Zmin=np.min(CscanAmpNew) ; diff=Zmax-Zmin  #27DEC2015 change
  print 'new Zmax,Zmin,diff,Zminfac,Zmaxfac=',Zmax,Zmin,diff,Zminfac,Zmaxfac
        
# done processing -> plot C-scan from gated A-scan data ==============================================================================

  Xmax=np.amax(Xnew) ; Xmin=np.min(Xnew) ; Ymax=np.amax(Ynew) ; Ymin=np.min(Ynew)
  fig1=plt.figure('raw C-scan',figsize=(figh*(Xmax-Xmin)/(Ymax-Ymin),figh)) # set width proportionally
  ax1=fig1.add_subplot(111)
  
# plot the raw (interpolation='none') C-scan and set origin "upper" as the C-scan is top-down mirror to the actual image 
  extent=(Xmin,Xmax,Ymax,Ymin)
  if VppGateOn==1:
    im1=ax1.imshow(CscanAmpNew,interpolation='none',cmap=cm.gray,origin='upper',vmin=VppL1,vmax=VppH2,extent=extent,picker=True,
                   norm=True)
  else:
    # put the following two lines in to make plots equal (9/20/2016)
    #Zmin = 0.0
    #Zmax = 0.75
    # 27DEC2015 change to allow conttrast variation
    im1=ax1.imshow(CscanAmpNew,interpolation='none',cmap=cm.gray,origin='upper',vmin=Zmin*Zminfac,vmax=Zmax*Zmaxfac,aspect=1., \
                  extent=extent,picker=True)
  ax1.set_ylim(Ymax+np.abs(Ymax)*CscanEdgeFac,Ymin-np.abs(Ymin)*CscanEdgeFac) 
  ax1.set_xlim(Xmin-np.abs(Xmin)*CscanEdgeFac,Xmax+np.abs(Xmax)*CscanEdgeFac)  # correced 17FEB2013
  if axis1=='Turntable':
    ax1.set_xlabel("Angular scan position (deg)")
  else:
    ax1.set_xlabel("X scan position (mm)")
  ax1.set_ylabel("Y scan position (mm)")
  # ax1.set_title("Gated peak-peak amplitude C-scan in new grid")
  # plt.colorbar(im1,orientation='horizontal')
  plt.colorbar(im1)
  ax1.grid(True)
  
# plot enhanced Cscan for quality output
  fig4=plt.figure('Interpolated C-scan',figsize=(figh*(Xmax-Xmin)/(Ymax-Ymin),figh))
  ax4=fig4.add_subplot(111)
                                 # PW pit project: use TeraView's DataProcess to try out optimal Vmin, Vmax on the raw image: vmin=0.1972,vmax=0.9237  

                                 # make edge stay out by using sciski's sobel filter.  The data must betweeen 0 and 1
                                 # 18NOV2012 test on metal ruler not good: image distorted; edge not better resolved
                                 #CscanAmpNew2=CscanAmpNew.copy() ; CscanAmpNew2=CscanAmpNew2/np.amax(CscanAmpNew2)
                                 #CscanAmpNew2=filter.sobel(CscanAmpNew2)  
                                 #ax3.imshow(CscanAmpNew2,interpolation='None',cmap=cm.jet_r,origin='upper',extent=extent)
                                 #ax4.imshow(CscanAmpNew,interpolation='bilinear',cmap=cm.gray,origin='upper',extent=extent) # further interpolate
  if VppGateOn==1:
    im4=ax4.imshow(CscanAmpNew,interpolation='bilinear',cmap=cm.jet,origin='upper',vmin=VppL1,vmax=VppH2,extent=extent) 
  else:
     
    """    
# BEGIN special section for EPRI cable project 23SEP2013: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    
    # flip both in X (angle) and Y and make circular shift in X to match up cable #1's two halves.  Y flipped by reversing set_ylim.
    # matching X: 1st half 258 = 2nd half 8 (after flipping)
    
    temp=np.zeros((Ystep,Xstep),dtype='<f')
    temp[:,::1]= CscanAmpNew[:,::-1] # flip x
    first=258 ; second=8
    temp[:,0:Xstep+second-first]= CscanAmpNew[:,first-second:Xstep] # "backward" circular shift x also flip?!
    temp[:,Xstep-first+second:Xstep]= CscanAmpNew[:,0:first-second]
    temp[:,::1]= temp[:,::-1] # flip x again
    im4=ax4.imshow(temp,interpolation='bilinear',cmap=cm.jet,origin='upper',aspect=0.5,extent=extent) # without keyword extent, only plot upper one half?
    ax4.set_ylim(Ymin,Ymax) # this reverse y axis
    ax4.set_xlim(Xmin,Xmax)
  
# END special section for EPRI cable project 23SEP2013: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
    """
    
    # 27DEC2015 change to allow conttrast variation
    im4=ax4.imshow(CscanAmpNew,interpolation='bilinear',cmap=cm.jet,origin='upper',vmin=Zmin*Zminfac,vmax=Zmax*Zmaxfac,extent=extent)  
  ax4.set_ylim(Ymax,Ymin) ; ax4.set_xlim(Xmin,Xmax)
  
  if axis1=='Turntable':
    ax4.set_xlabel("Angular scan position (deg)")
  else:
    ax4.set_xlabel("X scan position (mm)")
  ax4.set_ylabel("Y scan position (mm)")
  # ax4.set_title("Further processed C-scan in new grid")
  # plt.colorbar(im4,orientation='horizontal')
  plt.colorbar(im4)
  ax4.grid(True)
  
# plot 3D depth/thickness profile 19SEP2013   much Slower !

  if DepthMapOn!=0: 
    fig6=plt.figure('3D depth/thickness',figsize=(figh*(Xmax-Xmin)/(Ymax-Ymin),figh))
    ax6=fig6.gca(projection='3d')
    ax6.plot_surface(X3d,Y3d,DepThick,rstride=1, cstride=1,cmap=cm.jet,linewidth=0)
    # add pbaspect feature 19SEP2013 per Ben Root implementation: see web pages about fixing aspect ratio in C:\f1\current\THzProc\doc,paper
    # keep long axis scale close to 1 and scale the others proportional to it; otherwise, the plot may get clipped   
    scale=1.1
    if (Xmax-Xmin)>(Ymax-Ymin):
      ax6.pbaspect = [scale,scale*(Ymax-Ymin)/(Xmax-Xmin), scale]
    else:
      ax6.pbaspect = [scale*(Xmax-Xmin)/(Ymax-Ymin),scale, scale]
    #fig6.set_ylim(Ymax,Ymin)  screw up the plot?!
    ax6.set_xlabel("X scan position (mm)")
    ax6.set_ylabel("Y scan position (mm)")
    ax6.set_zlabel("Front Surface Profile (mils)")
     
  cursor = Cursor(ax1, useblit=True, color='red', linewidth=1)

# plot A-scan, spectrum and processed data as picked
  PickPt(fig1,Xnew,Ynew,WaveformNew,CscanAmpNew,maxfVpp,maxfVppAll,time,wavlen,wavtlen,freq,delt,Xstep,Ystep,  \
         Xmin,Xmax,Ymax,Ymin,axis1,BscanDir,PeakBin,DepThick,FollowGateOn)
  
  plt.show() # this line has to be here; ohterwise, C-scan cross cursor will not show 

# _______________________________________________________________________________________________________________________________________________________________________


def MakeCscan(WaveformNew,PeakBin,FollowGateOn,AmpCompOn,VppGateOn,Xstep,Ystep,SigType,FSEtol,VppL1,VppH1,VppL2,VppH2,CscanAmpNew,maxfVpp,maxfVppAll):

# make C-scan based on choice of SigType:

# Sigtype=1: CscanAmpNew=Vpp within the gate
#         2: avg mag (mag = abs amp)
#         3: median mag
#         4: min mag
#         5: max mag
#         6: avg amp ~ integrated gate
#         7: median amp
#         8: min amp 
#         9: max amp

# to be filled: CscanAmpNew,maxfVpp,maxfVppAll

# last update: 29DEC2015

  # use Vpp
  if SigType==1:
    for i in range(Ystep):
      for j in range(Xstep):
        CscanAmpNew[i,j]=WaveformNew[i,j,PeakBin[0,FollowGateOn,i,j]]-WaveformNew[i,j,PeakBin[1,FollowGateOn,i,j]]  
  
  # use avg mag within the gate (abs amplitude within the gate, sum up then avg) add: 20DEC2015
  elif SigType==2: 
    for i in range(Ystep):
      for j in range(Xstep):
        L=PeakBin[3,FollowGateOn,i,j] ; R=PeakBin[4,FollowGateOn,i,j]
        CscanAmpNew[i,j]=np.sum(np.abs(WaveformNew[i,j,L:R]))/(R-L)

# median magnitude within the gate        
  elif SigType==3: 
    for i in range(Ystep):
      for j in range(Xstep):
        L=PeakBin[3,FollowGateOn,i,j] ; R=PeakBin[4,FollowGateOn,i,j]
        CscanAmpNew[i,j]=np.median(np.abs(WaveformNew[i,j,L:R]))
  
  # min magnitude within the gate        
  elif SigType==4: 
    for i in range(Ystep):
      for j in range(Xstep):
        L=PeakBin[3,FollowGateOn,i,j] ; R=PeakBin[4,FollowGateOn,i,j]
        CscanAmpNew[i,j]=np.min(np.abs(WaveformNew[i,j,L:R]))
        
  # max magnitude within the gate
  elif SigType==5: 
    for i in range(Ystep):
      for j in range(Xstep):
        L=PeakBin[3,FollowGateOn,i,j] ; R=PeakBin[4,FollowGateOn,i,j]
        CscanAmpNew[i,j]=np.amax(np.abs(WaveformNew[i,j,L:R]))

  # avg amp within the gate ~ integrated gate
  elif SigType==6: 
    for i in range(Ystep):
      for j in range(Xstep):
        L=PeakBin[3,FollowGateOn,i,j] ; R=PeakBin[4,FollowGateOn,i,j]
        CscanAmpNew[i,j]=np.sum(WaveformNew[i,j,L:R])/(R-L)

# median amp within the gate
  elif SigType==7: 
    for i in range(Ystep):
      for j in range(Xstep):
        L=PeakBin[3,FollowGateOn,i,j] ; R=PeakBin[4,FollowGateOn,i,j]
        CscanAmpNew[i,j]=np.median(WaveformNew[i,j,L:R])

  # min amp within the gate
  elif SigType==8: 
    for i in range(Ystep):
      for j in range(Xstep):
        L=PeakBin[3,FollowGateOn,i,j] ; R=PeakBin[4,FollowGateOn,i,j]
        CscanAmpNew[i,j]=np.min(WaveformNew[i,j,L:R])

  # max amp within the gate
  elif SigType==9: 
    for i in range(Ystep):
      for j in range(Xstep):
        L=PeakBin[3,FollowGateOn,i,j] ; R=PeakBin[4,FollowGateOn,i,j]
        CscanAmpNew[i,j]=np.amax(WaveformNew[i,j,L:R])        
  
  if FSEtol>0: # 27JAN2015 FSE handle: if FSE too small, throw away this location
    for i in range(Ystep):
      for j in range(Xstep):
        if WaveformNew[i,j,PeakBin[0,0,i,j]]-WaveformNew[i,j,PeakBin[1,0,i,j]]<FSEtol:
          CscanAmpNew[i,j]=0.           
          
  """
  # use avg abs amp (abs amplitude integration over time then avg)
  elif SigType==3: # use avg abs amp (abs amplitude integration over time then avg)  changed: 26SEP2014, 27JAN2015 FSE handle
    if FSEtol>0:
      for i in range(Ystep):
        for j in range(Xstep):
          if WaveformNew[i,j,PeakBin[0,0,i,j]]-WaveformNew[i,j,PeakBin[1,0,i,j]]<FSEtol:
            CscanAmpNew[i,j]=0.
          else: 
            L=PeakBin[0,FollowGateOn,i,j]-nHalfPulse ; L=(L if L>0 else 0)
            R=PeakBin[1,FollowGateOn,i,j]+nHalfPulse ; R=(R if R<=wavlen else wavlen)
            CscanAmpNew[i,j]=np.sum(np.abs(WaveformNew[i,j,L:R]))/(R-L) #changed: 26SEP2014, 27JAN2015 FSE handle
    else:
      for i in range(Ystep):
        for j in range(Xstep):
          L=PeakBin[0,FollowGateOn,i,j]-nHalfPulse ; L=(L if L>0 else 0)
          R=PeakBin[1,FollowGateOn,i,j]+nHalfPulse ; R=(R if R<=wavlen else wavlen)
          if R==L:
            print i,j,L
          CscanAmpNew[i,j]=np.sum(np.abs(WaveformNew[i,j,L:R]))/(R-L) #changed: 26SEP2014, 27JAN2015 FSE handle
 """   
# amplitude compensation algorithm to be updated !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  if FollowGateOn>0 and AmpCompOn==1:
    for i in range(Ystep):  # detemine the Vpp of the signal in follower gate, normally front surface echo
        for j in range(Xstep):
          maxfVpp[i,j]=WaveformNew[i,j,PeakBin[0,0,i,j]]-WaveformNew[i,j,PeakBin[1,0,i,j]]
    maxfVppAll=np.amax(maxfVpp[:,:])    
    for i in range(Ystep):
        for j in range(Xstep):
          CscanAmpNew[i,j]=CscanAmpNew[i,j]*maxfVppAll/(maxfVpp[i,j]+maxfVppAll)      
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if VppGateOn==1: 
    for i in range(Ystep):
      for j in range(Xstep): #set up a dual-band hist equal
       if CscanAmpNew[i,j]>VppH2 or CscanAmpNew[i,j]<VppL1 or (CscanAmpNew[i,j]>VppH1 and CscanAmpNew[i,j]<VppL2): 
         CscanAmpNew[i,j]=0.
         
# _______________________________________________________________________________________________________________________________________________________________________  

def PickPt(fig1,X,Y,waveform,CscanAmp,maxfVpp,maxfVppAll,time,wavlen,wavtlen,freq,delt,Xstep,Ystep,Xmin,Xmax,Ymax,
           Ymin,axis1,BscanDir,PeakBin,DepThick,FollowGateOn):
  
# last update: 27DEC2015
  
  nHalfPulse=int(HalfPulse/delt)
  
  if AscanOnly==1:
    fig2=plt.figure('A-scan only',figsize=(8.,4.))
    ax21=fig2.add_subplot(111)
  else:
    fig2=plt.figure('A-scan and spectra',figsize=(8.,10.))
    ax21=fig2.add_subplot(311)
    ax22=fig2.add_subplot(312)
    ax23=fig2.add_subplot(313)
    
  if BscanOn==1:  # 27DEC2015 change
    if BscanDir==0: # horizontal crossline
      Bscan=np.zeros((wavlen,Xstep))
      BscanRatio=(Xmax-Xmin)/wavtlen
    elif BscanDir==1: # vertical crossline
      Bscan=np.zeros((wavlen,Ystep))
      BscanRatio=(Ymax-Ymin)/wavtlen
    fig3=plt.figure('B-scan',figsize=(figh*1.75*BscanRatio,figh))
    ax3=fig3.add_subplot(111)
  
  if DepthMapOn!=0: 
    fig5=plt.figure('Depth/Thickness Line Profile')
    ax5=fig5.add_subplot(111)
    LineX=np.linspace(Xmin,Xmax,Xstep)
  
  wav=np.zeros(wavlen) ; wav2=np.zeros(wavlen)
  delx=(Xmax-Xmin)/float(Xstep) ; dely=(Ymax-Ymin)/float(Ystep)

# onpick1 uses FindPeaks

  def onpick1(event): # =========================================================================================================================
    
# articst=AxesImage has no ind attribute; use mouseevent x,ydata instead   
    #xx=event.mouseevent.xdata ; yy=event.mouseevent.ydata  <- does not work as xdata, ydata are float64???
    #Do this: xx=float(event.mouseevent.xdata) ; yy=float(event.mouseevent.ydata)
# convert pixel locations in physical dimension (mm) into array indices 
    xID=int((event.mouseevent.xdata-Xmin)/delx) ; yID=int((event.mouseevent.ydata-Ymin)/dely)
    
# !!! wav=waveform[i] is only referencing, i.e. waveform[i] will change if changes are made in wav
    wav[0:]=waveform[yID,xID,0:]
    
    if BscanOn==1:  # 27DEC2015 change
      if BscanDir==0:
        for i in range(Xstep):
          Bscan[0:,i]=waveform[yID,i,0:]
      elif BscanDir==1:
        for i in range(Ystep):
          Bscan[0:,i]=waveform[i,xID,0:]
      
# plot waveform
    ax21.cla()
    if AmpCompOn==1: 
      wav[:]=wav[:]*maxfVppAll/(maxfVpp[yID,xID]+maxfVppAll) # <<<<<<<<
    # take Vpp fixed gate if no following 
    if FollowGateOn==0:
      GateVpp=wav[PeakBin[0,0,yID,xID]]-wav[PeakBin[1,0,yID,xID]]
    else:
      GateVpp=wav[PeakBin[0,FollowGateOn,yID,xID]]-wav[PeakBin[1,FollowGateOn,yID,xID]]
    ax21.plot(time,wav,'-r') 
    if FollowGateOn>0:
      ax21.plot([time[PeakBin[3,0,yID,xID]],time[PeakBin[3,0,yID,xID]]],[min(wav)*1.1,np.amax(wav)*1.1],'--k') # add the two lead gates
      ax21.plot([time[PeakBin[4,0,yID,xID]],time[PeakBin[4,0,yID,xID]]],[min(wav)*1.1,np.amax(wav)*1.1],'--k')
      ax21.plot([time[PeakBin[3,FollowGateOn,yID,xID]],time[PeakBin[3,FollowGateOn,yID,xID]]],[min(wav)*1.1,np.amax(wav)*1.1],'-b') # add the two follower gates
      ax21.plot([time[PeakBin[4,FollowGateOn,yID,xID]-1],time[PeakBin[4,FollowGateOn,yID,xID]-1]],[min(wav)*1.1,np.amax(wav)*1.1],'-g')
    ax21.set_xlim(0,wavtlen)
    ax21.set_ylim(np.min(wav)*1.1,np.amax(wav)*1.1)
    ax21.set_xlabel("Time (ps)")
    ax21.set_ylabel("Amplitude")                              # increase decimal places to 5  20DEC2015
    # now also print I,J   17APR2016
    ax21.set_title("I=%d J=%d X=%.4f Y=%.4f mm, GatVpp=%.5f, PxlAmp=%.5f" % (xID,yID,event.mouseevent.xdata,event.mouseevent.ydata,GateVpp,CscanAmp[yID,xID]))
    ax21.grid(True)
    
    if FlagPeakOn==1:  # only show as many as gates
      for i in range(PeakBin.shape[1]):
        ax21.scatter(time[PeakBin[0,i,yID,xID]], wav[PeakBin[0,i,yID,xID]],color='b',marker='+') # lead gate
        ax21.scatter(time[PeakBin[1,i,yID,xID]], wav[PeakBin[1,i,yID,xID]],color='g',marker='x') # follower gate
        #ax21.scatter(time[PeakBin[1,:,yID,xID]], wav[PeakBin[1,:,yID,xID]],color='k',marker='+')
    """
     maxtab, mintab = peakdet(wav,PeakTol*np.amax(np.abs(wav)),time) 
     ax21.scatter(maxtab[:,0], maxtab[:,1],color='b',marker='+')
     ax21.scatter(mintab[:,0], mintab[:,1],color='k',marker='+')
    """
# plot full spectrum
    if AscanOnly!=1:
      ax22.cla()
      wavfc=ReFFT2(wav,delt) ; wavfm=np.abs(wavfc)
      ax22.plot(freq[1:],wavfm[1:],label='full waveform',color='r',linestyle='-') # skip DC
      ax22.set_xlim(0.,workf)
      #ax22.set_ylim(0.,0.5)
      #ax22.set_ylim(np.min(np.log10(wavfm[2:-1]))*1.1,np.amax(np.log10(wavfm[2:-1]))*1.1)
      ax22.set_xlabel("Frequency (THz)")
      ax22.set_ylabel("Amplitude")
      ax22.legend()
      ax22.grid(True)
# plot gated spectrum
      ax23.cla()
      wav2[:]=0.
      if FollowGateOn>0:
        wav2[PeakBin[3,FollowGateOn,yID,xID]:PeakBin[4,FollowGateOn,yID,xID]]=wav[PeakBin[3,FollowGateOn,yID,xID]:PeakBin[4,FollowGateOn,yID,xID]]
      else:
        wav2[PeakBin[3,0,yID,xID]:PeakBin[4,0,yID,xID]]=wav[PeakBin[3,0,yID,xID]:PeakBin[4,0,yID,xID]]
      wavgatefc=ReFFT2(wav2,delt) ; wavgatefm=np.abs(wavgatefc)
      if FollowGateOn>0:
        ax23.plot(freq[1:],wavgatefm[1:],label='within blue-green gate',color='r',linestyle='-')  # skip DC
      else:
        ax23.plot(freq[1:],wavgatefm[1:],label='within black-black gate',color='r',linestyle='-')
      ax23.set_xlim(0.,workf)
      ax23.set_xlabel("Frequency (THz)")
      ax23.set_ylabel("Amplitude")
      ax23.legend()
      ax23.grid(True)
      
    fig2.canvas.draw()
    
# plot B-scan
    if BscanOn==1:  # 27DEC2015 change
      ax3.cla()
      if BscanDir==0: 
        extent=(Xmin,Xmax,0.,wavtlen)
      elif BscanDir==1: 
        extent=(Ymin,Ymax,0.,wavtlen)
      AspRatio=(Ymax-Ymin)/(Xmax-Xmin) # properly scale to display 08JUN2015
      ax3.imshow(Bscan,interpolation='bilinear',cmap=cm.seismic,origin='lower',aspect=AspRatio,extent=extent)
      #ax3.set_ylim(Ymax*(1.+CscanEdgeFac),Ymin*(1.+CscanEdgeFac)) 
      #ax3.set_xlim(Xmin*(1.+CscanEdgeFac),Xmax*(1.+CscanEdgeFac))
      ax3.set_ylabel("Time (ps)")
      if BscanDir==0:
        if axis1=='Turntable':
          ax3.set_xlabel("Angular scan position (deg)")
        else:
          ax3.set_xlabel("X scan position (mm)")
        ax3.set_title("Line Y=%.4f mm" % (event.mouseevent.ydata))
      elif BscanDir==1:  
        ax3.set_xlabel("Y scan position (mm)")
        ax3.set_title("Line X=%.4f mm" % (event.mouseevent.xdata))
      ax3.grid(True)
      fig3.canvas.draw()
    
# plot depth/thickness line 16SEP2013
    if DepthMapOn!=0:
      ax5.cla()
      # get thcikness across the whole scan line on thich the point is picked   
      ax5.plot(LineX,DepThick[yID,:],'-r')
      ax5.set_xlabel("X Scan Position (mm)")
      if DepthMapOn==1:
        ax5.set_ylabel("TBC Thickness (mils)")
        ax5.set_title("TBC thickness across scan line Y="+str(Y[yID])+"mm")
      elif np.abs(DepthMapOn)==2:
        ax5.set_ylabel("Front Surface Profile (mils)")
        ax5.set_title("Front surface profile across scan line Y="+str(Y[yID])+"mm")
      ax5.set_xlim(Xmin,Xmax)
      #ax5.set_ylim(np.amin(DepThick[yID,:])*1.1,np.amax(DepThick[yID,:])*1.1)  does not work! indexing incorrect?!
      ax5.grid(True)
      fig5.canvas.draw()
      
  fig1.canvas.mpl_connect('pick_event', onpick1)

# _______________________________________________________________________________________________________________________________________________________________________ 

def Mayavi3D(X3d,Y3d,DepThick):

# use Mayavi to plot 3D
# lastupdate: 16JUL2014
    
  if FollowGateOn>0:
    Zlabel="Layer "+str(FollowGateOn)+"Thickness (0.01mm)"
  else:
    Zlabel="Front Surface Profile (0.01mm)"
  mlab.figure(bgcolor=(1.,1.,1.),fgcolor=(0.,0.,0.)) # the fgcolor dictates the axes color including title and labels
  surf = mlab.mesh(X3d, Y3d, DepThick*100.)
  mlab.outline(surf,color=(.7, .7, .7))
  t=mlab.axes(surf, xlabel='X scan position (mm)', ylabel='Y scan position (mm)',zlabel=Zlabel, \
  nb_labels=5,x_axis_visibility=True, y_axis_visibility=True, z_axis_visibility=True)
  # the line below does not work -> font factor can't be changed by code -> go to axes.py to hardwire it
  #t.CubeAxesActor2D_property.font_factor = 0.9   
  t.title_text_property.font_family = 'time'
  t.title_text_property.font_size = 4  # this does not work right unless font factor changed 
  #t.title_text_property.color=(0.7,0.7,0.7)  # this overwrite the fgcolor in figure function
  t.label_text_property.font_family = 'time'
  t.label_text_property.font_size = 1  # this does not work right unless font factor changed 
  #t.label_text_property.color=(0.3,0.3,0.3) # this overwrite the fgcolor in figure function
  t.label_text_property.bold = False
  # mlab.text(0., 0., '3D Geometry Rendering', z=np.amax(DepTick), width=0.1)
  mlab.show()
  
# _______________________________________________________________________________________________________________________________________________________________________

def OutDep(basedir,skip,X3d,Y3d,DepThick):

# output front surface or depth profile to a textfile
# lastupdate: 15JUL2014
  
  outfile=open(os.path.join(basedir,'DepOut.txt'), 'w')
  
  print skip,X3d.shape[0],Y3d.shape[1]
  for i in range(skip,X3d.shape[0]-skip):
    for j in range(skip,Y3d.shape[1]-skip):
      outfile.write('%14.8f %14.8f %14.8f' % (X3d[i,j],Y3d[i,j],DepThick[i,j]))
      outfile.write('\n')
  outfile.close()

# _______________________________________________________________________________________________________________________________________________________________________

def FindDepThick(Surf,PickPeak,waveform,peak,Xstep,Ystep,wavlen,delt,n,AmpTol):
  
# Surf=DepthMapOn, peak=PeakBin
  
# find front surface profile or layer thickness given PickSurf choice and gate info. from array peak

# PickPeak=1: pos peak, =2: neg. peak, =3: mid-point in peak (prefer 1; 3 gets errors?!)
# need to improve edge detection

# last update: 15JUL2014

  cs=np.sqrt(1.-(np.sin(16.6*np.pi/180.)/n)**2)
  #ref=cs*(1000./25.4)*0.299*0.5*wavlen*delt
  DepThick=np.zeros((Ystep,Xstep),dtype='<f') ; FSamp=np.zeros((Ystep,Xstep),dtype='<f')
  
  for i in range(Ystep):
    for j in range(Xstep):
      FSamp[i,j]=waveform[i,j,peak[0,0,i,j]]-waveform[i,j,peak[1,0,i,j]]
  thres=AmpTol*np.amax(FSamp)
  
  if FollowGateOn>=1: # =1: find thickness between FS and 1st interface; =2: between 1st and 2nd interfaces, etc.
    for i in range(Ystep):
      for j in range(Xstep):
        # on edges, pulses get distortion and reduced amplitude, and TOF may be shorter than it really is.
        # So set threshold on Vpp amp to reduce these artifacts of earlier TOFs.
        if FSamp[i,j]>thres:
          DepThick[i,j]=cs*0.299*0.5*(peak[PickPeak-1,FollowGateOn,i,j]-peak[PickPeak-1,FollowGateOn-1,i,j])*delt/n #mm
        else:
          DepThick[i,j]=0.
  else: # FollowGateOn=0, so find FS profile: default"hole to stick out"
    for i in range(Ystep):
      for j in range(Xstep):
        if FSamp[i,j]>thres:
          DepThick[i,j]=cs*0.299*0.5*peak[PickPeak-1,0,i,j]*delt #mm
        else:
          DepThick[i,j]=0.
    if Surf==-1:
      deep=np.amax(DepThick)
      for i in range(Ystep):
        for j in range(Xstep):
          if DepThick[i,j]>tiny:
            DepThick[i,j]=deep-DepThick[i,j] # hole to sink in, as it looks like.  The deepest hole has zero thickness
    
  return DepThick
  
# _______________________________________________________________________________________________________________________________________________________________________ 

def FindPeaks(waveform,Xstep,Ystep,wavlen,nHalfPulse,fthres,BinRange, FollowGateOn):

# given list BinRange, find peak(s) and gate(s) in wavform in terms of bin locations to PeakBin[i,npeak,Ystep,Xstep], 
# i=0: pos. peak, =1 neg., =2 half way, 3=left gate, 4=right gate

# BinRange[0][0-1] are for the leading peak, ususally FSE, BinRange[1][0-1] is for 1st layer, etc.
# assume max. and min. peaks vary within 2*nHalfPulse width and individual pulses are bound within their corresponding BinRange
# future improvement: use of peakdet

# last update: 29DEC2015

  
  npeak=len(BinRange)
  PeakBin=np.zeros((5,npeak,Ystep,Xstep),dtype=np.int16)
  
  # if FollowGateOn>0, BinRange[0][0-1] are used as fixed gate following pos. peak of leading pulse (usually FSE).
  # all trailing gates are related to this pos. peak
  
  if FollowGateOn>0: 
    
    if npeak<2: # must have at least 2 peaks (including the lead peak to be followed)
      print "incorrect BinRange setting!"
      sys.exit(1)
    for i in range(Ystep):
      for j in range(Xstep):
        PeakBin[0,0,i,j]=np.argmax(waveform[i,j,BinRange[0][0]:BinRange[0][1]+1])+BinRange[0][0] # important! need to add BinRange[0][0] to argmax result
        # check more carefully to see if FSE is smaller (due to defocusing, etc) than the largest peak in the gate
        if PeakBin[0,0,i,j]-nHalfPulse-BinRange[0][0] > 0: # in case PeakBin[0,0,i,j] too early
          itmp=(np.where(waveform[i,j,BinRange[0][0]:PeakBin[0,0,i,j]-nHalfPulse]>fthres*waveform[i,j,PeakBin[0,0,i,j]]))
          if len(itmp[0])!=0: # see if any peak above the threshold before maxfloc[i,j]
            itmp2=itmp[0][0]+BinRange[0][0] # np.where returns tuple in itmp
            itmp2=np.argmax(waveform[i,j,itmp2:itmp2+nHalfPulse])+itmp2
            PeakBin[0,0,i,j]=itmp2
        # find the neg. peak within pulse width=2*nHalfPulse
        L=PeakBin[0,0,i,j]-nHalfPulse ; L=(L if L>0 else 0)
        R=PeakBin[0,0,i,j]+nHalfPulse ; R=(R if R<=wavlen else wavlen)
        PeakBin[1,0,i,j]=np.argmin(waveform[i,j,L:R])+L
        PeakBin[2,0,i,j]=(PeakBin[0,0,i,j]+PeakBin[1,0,i,j])/2
        PeakBin[3,0,i,j]=BinRange[0][0] ; PeakBin[4,0,i,j]=BinRange[0][1]
        # all trailing gates are then shifted wrt the pos. leading peak, PeakBin[0,0,i,j]
        for k in range(1,npeak):
          L=PeakBin[0,0,i,j]+BinRange[k][0] ; R=PeakBin[0,0,i,j]+BinRange[k][1]+1
          if L<0:
            L=0
          elif L>wavlen:
            print "incorrect left gate setting in ",k,"th gate!"
            sys.exit(1)
          if R<L:  #corrected 15SEP2013
            print "incorrect right gate setting in ",k,"th gate!"
            sys.exit(1)
          elif R>wavlen:
            R=wavlen
          PeakBin[3,k,i,j]=L ; PeakBin[4,k,i,j]=R
          PeakBin[0,k,i,j]=np.argmax(waveform[i,j,L:R])+L
          if PulseLen<0.: # 29DEC2015: enforce searcg of neg. peak within original gate
            L2=L ; R2=R
          else:
            LL=int(PulseLen*nHalfPulse)
            if LL<1:  LL=1 # 29DEC2015: to prevent zero-length gate
            L2=PeakBin[0,k,i,j]-LL ; L2=(L2 if L2>L else L)
            R2=PeakBin[0,k,i,j]+LL ; R2=(R2 if R2<=wavlen else wavlen)
          PeakBin[1,k,i,j]=np.argmin(waveform[i,j,L2:R2])+L2
          PeakBin[2,k,i,j]=(PeakBin[0,k,i,j]+PeakBin[1,k,i,j])/2
          
  else: # =0: the whole waveform (the BinRange is given only [0,wavlen] in main)

    for i in range(Ystep):
      for j in range(Xstep):
        PeakBin[0,0,i,j]=np.argmax(waveform[i,j,BinRange[0][0]:BinRange[0][1]+1])+BinRange[0][0]
        L2=PeakBin[0,0,i,j]-nHalfPulse ; L2=(L2 if L2>BinRange[0][0] else BinRange[0][0])
        R2=PeakBin[0,0,i,j]+nHalfPulse ; R2=(R2 if R2<=BinRange[0][1] else BinRange[0][1])
        PeakBin[1,0,i,j]=np.argmin(waveform[i,j,L2:R2])+L2
        PeakBin[2,0,i,j]=(PeakBin[0,0,i,j]+PeakBin[1,0,i,j])/2
        PeakBin[3,0,i,j]=BinRange[0][0] ; PeakBin[4,0,i,j]=BinRange[0][1]
    
  return PeakBin

# _______________________________________________________________________________________________________________________________________________________________________

# DO NOT USE ! OUT OF DATE!

def Slicer(datlen,wavlen,wavtlen,delt,delf,time,freq,Xstep,Ystep,axis1,WaveformNew,Xnew,Ynew,CscanAmpNew,PeakBin,  \
            workf,FollowGateOn,AmpCompOn,VppGateOn,SlHisEquOn,Slhw,Sldel,fps,SigType):

# last update: 28DEC2015  


# make C-scan depending on the choice of SigType

  MakeCscan(WaveformNew,PeakBin,FollowGateOn,AmpCompOn,VppGateOn,Xstep,Ystep,SigType,FSEtol,VppL1,VppH1,VppL2,VppH2,CscanAmpNew,maxfVpp,maxfVppAll) 

  maxfloc=np.zeros((Ystep,Xstep),dtype=np.int16) # int32 does not work?!
  if FollowGateOn != 1 and FollowGateOn != 2:
    if begt<0: begt=-begt # in case forgot to set begt > 0
    for i in range(Ystep):
      for j in range(Xstep):
        CscanAmpNew[i,j]=np.amax(WaveformNew[i,j,begt:endt+1])-np.min(WaveformNew[i,j,begt:endt+1])
  else:
    tmpb=0 ; tmpe=wavlen
    for i in range(Ystep):
      for j in range(Xstep):
        maxfloc[i,j]=np.argmax(WaveformNew[i,j,fbeg:fend+1])+fbeg # important! need to add fbeg to argmax result
        
        itmp=(np.where(WaveformNew[i,j,fbeg:maxfloc[i,j]-nHalfPulse]>fthres*WaveformNew[i,j,maxfloc[i,j]]))
        if len(itmp[0])!=0: # see if any peak above the threshold before maxfloc[i,j]
          itmp3=itmp[0][0]+fbeg # np.where returns tuple in itmp
          itmp3=np.argmax(WaveformNew[i,j,itmp3:itmp3+nHalfPulse])+itmp3
          maxfloc[i,j]=itmp3
        
        L=maxfloc[i,j]+begt ; R=maxfloc[i,j]+endt+1
        if L<0:  # this will happen only if begt<0
          tmpb=np.min([tmpb,L]) # np.min or np.amax take only arrays/lists
          L=0
        elif L>wavlen:
          print "incorrect left gate setting!"
          sys.exit(1)
        if R<0:
          print "incorrect right gate setting!"
          sys.exit(1)
        elif R>wavlen:
          tmpe=np.amax([tmpe,R])
          R=wavlen
          print 'tmpe,R',tmpe,R
        CscanAmpNew[i,j]=np.amax(WaveformNew[i,j,L:R])-np.min(WaveformNew[i,j,L:R])
        
        
        
        
    begt=begt-tmpb ; endt=endt-(tmpe-wavlen) # shrinked to the inner-most bounds for the global gates
    
    
    
# compensate amplitude drop (due to curvature of front surface, etc.)  
    if FollowGateOn == 2:
      maxfVpp=np.zeros((Ystep,Xstep))
      for i in range(Ystep):  # detemine the Vpp of the signal in the lead gate, normally front surface echo
          for j in range(Xstep):
            maxfVpp[i,j]=np.amax(WaveformNew[i,j,maxfloc[i,j]-nHalfPulse:maxfloc[i,j]+nHalfPulse])  \
                          -np.min(WaveformNew[i,j,maxfloc[i,j]-nHalfPulse:maxfloc[i,j]+nHalfPulse])
      maxfVppAll=np.amax(maxfVpp[:,:])    
      for i in range(Ystep):
          for j in range(Xstep):
            CscanAmpNew[i,j]=CscanAmpNew[i,j]*maxfVppAll/maxfVpp[i,j]  # this is not the same as in InterAct
            
            
            
            
            
  print endt,begt,Slhw,Sldel
  nslice=int((endt-begt-2*Slhw)/Sldel) ; print '# slice=',nslice
  MaxCscanAmp=np.amax(CscanAmpNew) # MaxCscanAmp may be same as maxfVppAll, depending on the gate locations
  
  print 'new Zmax,Zmin,diff=', MaxCscanAmp,np.min(CscanAmpNew),np.amax(CscanAmpNew)-np.min(CscanAmpNew)

# now has the maxfloc[i,j] info., go do the slicing 
  files=[] ; mpgname='slice.mpg'
  Xmax=np.amax(Xnew) ; Xmin=np.min(Xnew) ; Ymax=np.amax(Ynew) ; Ymin=np.min(Ynew)
  fig1=plt.figure(figsize=(12.,14.)) # fig1=plt.figure(figsize=(12.,5.))for R-R CMC "dog bone"
                                      #fig1=plt.figure(figsize=(12.,14.)) this is good for R-R #2 coated sample
                                      #fig1=plt.figure(figsize=(figh*(Xmax-Xmin)/(Ymax-Ymin),figh)) # set width proportionally
                                      #ax1=fig1.add_subplot(212)
                                      #ax2=fig1.add_subplot(211)
  ax2=plt.subplot2grid((2,1),(0,0))  # try using grid control to change the sizes
  ax1=plt.subplot2grid((2,1),(1,0))
                                      #the image still only occupy the middle cell using "ax1=plt.subplot2grid((3,1),(1,0),colspan=2)" ???
                                      #  Change figsize works better
  extent=(Xmin,Xmax,Ymax,Ymin)       
  LL=np.zeros((Ystep,Xstep),np.int16) # need to separately track the Ls for each scan position per slice, in order to 
  for i in range(nslice):             # match with the one shown in Ascan
    for j in range(Ystep):
      for k in range(Xstep):
        LL[j,k]=maxfloc[j,k]+begt+i*Sldel ; R=LL[j,k]+2*Slhw
        CscanAmpNew[j,k]=np.amax(WaveformNew[j,k,LL[j,k]:R])-np.min(WaveformNew[j,k,LL[j,k]:R])
    ax2.plot(time,WaveformNew[Ystep/2,Xstep/2,:],'-r',label='right double coat area') # <--------------------------- to remove for others
    ax2.plot(time,WaveformNew[Ystep/2,Xstep/5,:],'-b',label='left single coat area')  # <--------------------------- to remove for others

# amplitude compensation algorithm to be updated !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if FollowGateOn == 2: # compensation for amplitude drop (often due to curvature of front surface)
      for j in range(Ystep): 
          for k in range(Xstep):
            CscanAmpNew[j,k]=CscanAmpNew[j,k]*maxfVppAll/maxfVpp[j,k]
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    #ax2.plot(time,WaveformNew[Ystep/2,Xstep/2,:],'-b',label='center of sample') # <---------------------------------- to remove for others
    ax2.plot([time[LL[Ystep/2,Xstep/2]],time[LL[Ystep/2,Xstep/2]]],[-5.,5.],'--k') # add the two time gates
    ax2.plot([time[LL[Ystep/2,Xstep/2]+2*Slhw],time[LL[Ystep/2,Xstep/2]+2*Slhw]],[-5.,5.],'--k')
    ax2.set_xlim(0,wavtlen)
    ax2.set_ylim(np.min(WaveformNew[Ystep/2,Xstep/2,:])*1.1,np.amax(WaveformNew[Ystep/2,Xstep/2,:])*1.1)
    ax2.set_xlabel("Time (ps)")
    ax2.set_ylabel("Amplitude (a.u)")
    ax2.legend()
    ax2.grid(True)
    if SlHisEquOn==1:  
      ax1.imshow(CscanAmpNew,interpolation='bicubic',cmap=cm.jet_r,origin='upper',extent=extent)
# plot the raw (interpolation='none') C-scan and set origin "upper" as the C-scan is top-down mirror to the actual image 
    else:
      ax1.imshow(CscanAmpNew,interpolation='bicubic',cmap=cm.jet_r,origin='upper',vmin=0.,vmax=MaxCscanAmp,extent=extent)
    ax1.set_ylim(Ymax*(1.+CscanEdgeFac),Ymin*(1.+CscanEdgeFac)) 
    ax1.set_xlim(Xmin*(1.+CscanEdgeFac),Xmax*(1.+CscanEdgeFac))
    if axis1=='Turntable':
      ax1.set_xlabel("Angular scan position (deg)")
    else:
      ax1.set_xlabel("X scan position (mm)")
    ax1.set_ylabel("Y scan position (mm)")
    ax1.grid(True)        
    fname = '_tmp%03d.png'%i
    fig1.savefig(fname)
    files.append(fname)
    ax1.cla() ; ax2.cla()
    print 'Done making frame',i
   
  cmd='mencoder mf://_tmp*.png -mf type=png:fps=%02d -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o '%fps
  os.system(cmd+mpgname)
  print 'Done making movie '+mpgname

# _______________________________________________________________________________________________________________________________________________________________________


def AmpCor300(TLL,TLR,ALL,ALR,TRL,TRR,ARL,ARR,p,wavlen,delt,Xstep,Ystep,WaveformNew):
  
# fit A=a*T**b on both ends to remove excessive amplitude amplification for 300 ps waveforms
# last update: 08JUN2015

  nLL=int(TLL/delt) ; nLR=int(TLR/delt) 
  Tl=np.linspace(nLL*delt,nLR*delt,nLR-nLL+1)
  KL=np.power(ALL/ALR,1./p)
  betaL=TLR+(TLR-TLL)/(KL-1.)
  alphaL=ALL/(betaL-TLL)**p
  Al=alphaL*(betaL-Tl)**p
  
  nRL=int(TRL/delt) ; nRR=int(TRR/delt) ; print 'nLL,nRR',nLL,nRR
  Tr=np.linspace(nRL*delt,nRR*delt,nRR-nRL+1)
  KR=np.power(ARL/ARR,1./p)
  betaR=TRR-(TRL-TRR)/(KR-1.)
  alphaR=ARL/(betaR-TRL)**p
  Ar=alphaR*(betaR-Tr)**p
  
  #plt.figure()
  #plt.plot(WaveformNew[10,20,:])
  
  if nRR==wavlen:
    nRR-=1  # prevent index out of range
  for i in range(Ystep):
    for j in range(Xstep):
       WaveformNew[i,j,nLL:nLR+1]=WaveformNew[i,j,nLL:nLR+1]/Al[0:]
       WaveformNew[i,j,nRL:nRR+1]=WaveformNew[i,j,nRL:nRR+1]/Ar[0:]
       
  #plt.plot(WaveformNew[10,20,:],'r')
  #plt.show()
  #sys.exit()

# _______________________________________________________________________________________________________________________________________________________________________

def ReMap(waveform,X,Y,Xmax,Xmin,Ymin,Xres,Yres,Xstep,Ystep,ScanType,axis1,wavlen,wavtlen,XCorTol,YDiffTol,XChkTol,tiny):

# given the non-uniform X (or S, theata) cooirdinates, remap to a uniform rectangular grid
# 1st method: remap within each line by interpolation
# I/O: transformed from 1-D array "waveform" on non-uniform grid at 1-D X and Y coord. arrays for each element (scan position) 
#      of "waveform" into 2-D "WaveformNew" on uniform grid at much shorter 1-D Xnew and Ynew arrays only designated to common 
#      coord.
# handle both bi- and uni-directional scans

# CAN'T RELY ON TVL'S Xstep !!!!! A major bug found! Tvl reported correct Xstep (most time?), but in some larger scans,
# the actual data points on a X line came out short by as much as 5% (although the scan did reach the desired physical scan 
# position). Has to do time-consuming check on every data point (the old way takes a shortcut to check only around the estimated
# ends)!  DEC 2012

# the angular coord. of rotational scan was found incorrect also.  The actual pyhysic position seemed to be right (had tested) and
# the number of scan points are also fine, but the angular coords. stored in data file are off significantly. Before Teraview fixes
# it, the temp.remedy is to map the max. and min. of inccorect coords. to the desired range and scale all coords. accordingly.
# FEB 2013

# now use true Xres and Yres.  Fix the 'middle points' at the ends. It seems the Y coord. and X spacing in the central 'stable' portion
# of a line are accurate to 0.01mm on average.  X coords. are generally off due to backlash.  These observations are based on the 6 PW
# TBC tablet data 126 X 26 @0.2mm, 30ps, 4096pt, 1avg.
# 26AUG2013

# fix flyback and rotational scan and np.max -> np.amax
# 22SEP2013

# add baseline trend removal using waveforms around edges
# 15NOV2013

# add more baseline removal option
# 01JAN2016

# last update: 01JAN2016

  # Xres and Yres are true values set by ScanAcquire
  XChkLen=int(Xstep/4) ; tol=10 ; #print Xstep
  XCor=XCorTol*Xres ; YDiff=Yres*YDiffTol ; XDiff=Xres*XDiffTol ; pos=[] ; Ynew=[] #Ynew=np.zeros(Ystep,dtype='<f')
  # ibeg=(1 if np.abs(Y[0]-Y[1])>YDiff else 0) # see if the 1st scan position is an outlier <--- this is not enough; replaced by below

  datlen=len(waveform)
# for rotational scan, just inspect Y coord. of initial points if X[0] is off too much; throw them away if they are way out of line
  if axis1=='Turntable' and np.abs(X[0]-Xmin)>5.*XDiff:
    for i in range(datlen):
      if np.abs(Y[i]-Ymin)<YDiff: 
        ibeg=i
        print 'rotational scan starting point:',i
        break
  else:
# for linear scan, first inspect both coords. of initial points; throw them away if they are way out of line   15DEC2012 add
# (e.g. staring 2 Ys of "R-R 30492E blade 1 mil purge scan 2 HD (50mm 40ps 2048pt 1avg).tvl" were way off)      
    for i in range(datlen):
      if np.abs(Y[i]-Ymin)<YDiff and np.abs(X[i]-Xmin)<XDiff:
        ibeg=i
        break
  
# next find coordinate info and make a slight corection to backlash for bi-directional scan
  i=1 ; j=ibeg ; nline=0 ; Ynow=Ymin ; Ynext=Ymin+Yres
  #print ibeg,Xstep,XChkLen,j+Xstep-XChkLen,j+Xstep+tol
  #print 'datalen', datlen
  while (j<datlen):
    for k in range(j+Xstep-XChkLen,j+Xstep+tol): # jump to check last XChkLen pts.
      if np.abs(Y[k]-Ynow)>YDiff:
        # print i,k,Y[k], Ynow,YDiff
        nline=k-ibeg
        j=k-1
        break
                          #print 'nline=',nline
    if nline<nlim*Xstep: print 'line ',i,'has ',nline,' data point; too short!' # give a warning for short line which is likely problem in the data
    # compute the avg Y of the lcentral portion of the line; will use it for the new "map"
    tmp=np.sum(Y[ibeg+XChkLen:j-XChkLen])/float(j-2*XChkLen-ibeg) ; Ynew.append(tmp)      
    if ScanType=='2D Image Scan with encoder': # make the correction, even the pulse-on position is not too bad given the backlash
      X[ibeg:j]=X[ibeg:j]+(XCor if np.mod(i,2)==0 else -XCor) 
    pos.append([nline,X[ibeg],Y[ibeg],ibeg,X[j],Y[j],j,tmp,Ynow,Ynow-tmp])
    for k in range(j+1,j+tol): # find the starting point of a new line
      if np.abs(Y[k]-Ynext)<YDiff:  # skip as many "middle points" as needed   
        ibeg=k
        j=ibeg
        nline=0
        Ynow=Ynext
        i=i+1
        Ynext=Ymin+float(i)*Yres
        break
    if i==Ystep:
      nline=datlen-j+1
      tmp=np.sum(Y[ibeg+XChkLen:datlen-XChkLen])/float(datlen-2*XChkLen-ibeg) ; Ynew.append(tmp)      
      if ScanType=='2D Image Scan with encoder': # make the correction, even the pulse-on position is not too bad given the backlash
        X[ibeg:]=X[ibeg:]+(XCor if np.mod(i,2)==0 else -XCor) 
      pos.append([nline,X[ibeg],Y[ibeg],ibeg,X[-1],Y[-1],datlen,tmp,Ynow,Ynow-tmp])
      break
            
  Ynew=np.array(Ynew,dtype='<f')
  
                    #for k in range(len(pos)):
                    #  print k,pos[k]
  
  if i !=Ystep:
    print '!!! Ystep is wrong! should be',i,'!!!'  # this is fatal error: the steps in Y found should agree with Ystep 99.99%
    Ystep=i
  Xstep=np.sum([pos[i][0] for i in range(Ystep)])/Ystep
  print 'avg Xstep=',Xstep, ' as new Xstep' # this is the new Xstep to be used from this point on
  
  print 'max. Y coord. deviation=',np.amax([pos[i][-1] for i in range(Ystep)]),'mm'
  
# then determine the X bound for a smaller, trimmed grid

  if ScanType=='2D Image Scan with encoder': #bi-directional
                 #if np.sign(pos[0][1])!=np.sign(pos[1][1]) or \
                 #  (np.sign(pos[0][1])==np.sign(pos[1][1]) and np.abs(pos[0][1]-pos[1][1])>Xres): # test bi-directional or unidirectional
                         #print [pos[i][1] for i in range(0,Ystep,2)]+[pos[i][4] for i in range(1,Ystep,2)]
                         #print [pos[i][4] for i in range(0,Ystep,2)]+[pos[i][1] for i in range(1,Ystep,2)]
    Xbeg=np.amax([pos[i][1] for i in range(0,Ystep,2)]+[pos[i][4] for i in range(1,Ystep,2)])
    Xend=np.min([pos[i][4] for i in range(0,Ystep,2)]+[pos[i][1] for i in range(1,Ystep,2)])
  elif ScanType=='Flyback 2D scan': # unidirectional
    Xbeg=np.amax([pos[i][1] for i in range(Ystep)])
    Xend=np.min([pos[i][4] for i in range(Ystep)]) # was Xbeg=np.max([pos[:][1]), etc. Got Wrong number!
  else:
    print "unknown scan type!"
    sys.exit(1)
  
# correct angular coords. of rotational scan: map the current wrong coords. closer to the (Xmin,Xmax) size 
# enlarge the size a bit (by adding tiny at each end) to allow later interpolation possible for all new positions   
  scale=(Xmax-Xmin+2.*tiny)/(Xend-Xbeg)
  if axis1=='Turntable' and (scale>1.1 or scale<0.9): 
    X[:]=(X[:]-Xbeg)*scale+Xmin
    print 'new Xbeg=Xmin, Xend=Xmax'
    Xnew=np.linspace(Xmin,Xmax,Xstep)  
  else:
    print 'new Xbeg,Xend=',Xbeg+tiny,Xend-tiny,'mm which is',(Xend-Xbeg-2.*tiny),'mm long comparing with desired ',Xmax-Xmin,'mm'
    Xnew=np.linspace(Xbeg+tiny,Xend-tiny,Xstep) # add tiny to make sure all new allocated points are inside bounds
    
  WaveformNew=np.zeros((Ystep,Xstep,wavlen),dtype='<f')
  CscanAmpNew=np.zeros((Ystep,Xstep),dtype='<f')

# remap each line by interpolation -----------------------------------------------------------------------------------------------------------------
# WARNING for rotational scan! Xs in pos list are now incorrect (still the values beofore corrections)   
  if ScanType=='2D Image Scan with encoder':
    for i in range(0,Ystep,2):
      b=pos[i][3] ; e=pos[i][6]
      for j in range(Xstep):
                            #print 'forward b,e',b,e
        for k in range(b,e):
          if Xnew[j]>=X[k] and Xnew[j]<X[k+1]:
                           #if i==78: print 'i,j,k,Xnew[j],X[k],X[k+1]',i,j,k,Xnew[j],X[k],X[k+1]
            w1=np.abs(Xnew[j]-X[k]) ; w2=np.abs(X[k+1]-Xnew[j]) ; w=w1+w2
            WaveformNew[i,j,:]=(w2/w)*waveform[k][:]+(w1/w)*waveform[k+1][:]
            b=k
            break
    for i in range(1,Ystep,2): # reverse direction
      b=pos[i][3] ; e=pos[i][6]
      for j in range(Xstep-1,-1,-1):
                          #print ' reverse b,e',b,e
        for k in range(b,e):
                          #if i==79 and j==0: print 'b,e,i,j,k,Xnew[j],X[k],X[k+1]',b,e,i,j,k,Xnew[j],X[k],X[k+1]
          if Xnew[j]<X[k] and Xnew[j]>=X[k+1]:
                          #if i==79: print 'i,j,k,Xnew[j],X[k],X[k+1]',i,j,k,Xnew[j],X[k],X[k+1]
            w1=np.abs(X[k]-Xnew[j]) ; w2=np.abs(Xnew[j]-X[k+1]) ; w=w1+w2
            WaveformNew[i,j,:]=(w2/w)*waveform[k][:]+(w1/w)*waveform[k+1][:]
            b=k
            break
    
  elif ScanType=='Flyback 2D scan':
    for i in range(Ystep):
      b=pos[i][3] ; e=pos[i][6]
      for j in range(Xstep):
        for k in range(b,e):
          if Xnew[j]>=X[k] and Xnew[j]<X[k+1]:
            w1=Xnew[j]-X[k] ; w2=X[k+1]-Xnew[j] ; w=w1+w2
            WaveformNew[i,j,:]=(w2/w)*waveform[k][:]+(w1/w)*waveform[k+1][:]
            b=k
            break
  else:
    print "unknown scan type!"
    sys.exit(1)
    
# remove baseline trend in apertured near-field -------------------------------------------------------------------------------------------------

  # average several "blank" waveforms (i.e. not on the sample) at beginning and end of row as reference to be substracted from the sample waveforms
  # assume there are at least some (1mm or more preferred) "blank" at beginning and end of row and X spacing is sufficiently smaller (0.05mm preferred)
  if TrendOff==5: 
    refwav=np.zeros(wavlen,dtype='<f')
    if Xres<0.11:
      print '!!! there are at least 1mm blank at beginning and end of row, or baseline removal may not work properly !!!'
    if Xres>0.11: #mm
      for i in range(Ystep):
        refwav[:]=(WaveformNew[i,1,:]+WaveformNew[i,-1,:])/2. # take 2nd points, as first may deviate from the correct locations
        for j in range(Xstep):    
          WaveformNew[i,j,:]=WaveformNew[i,j,:]-refwav[:]
    elif Xres<0.11 and Xres>0.051:
      for i in range(Ystep):
        refwav[:]=(WaveformNew[i,3,:]+WaveformNew[i,-3,:]+WaveformNew[i,7,:]+WaveformNew[i,-7,:])/4.
        for j in range(Xstep):    
          WaveformNew[i,j,:]=WaveformNew[i,j,:]-refwav[:]
    elif Xres<0.051:
      for i in range(Ystep):
        refwav[:]=(WaveformNew[i,5,:]+WaveformNew[i,-5,:]+WaveformNew[i,10,:]+WaveformNew[i,-10,:]+WaveformNew[i,15,:]+WaveformNew[i,-15,:])/6.
        for j in range(Xstep):    
          WaveformNew[i,j,:]=WaveformNew[i,j,:]-refwav[:]
        
  # use the avg (option 5) waveform of each row as reference; aligned each waveform in that row with the rescaled (option 2) reference
  # probably the best option of 2nd gen
  elif TrendOff==-1: 
    refwav=np.zeros(wavlen,dtype='<f')
    early=np.zeros(100,dtype=int) ; late=np.zeros(100,dtype=int) ; total=float(Xstep*Ystep)
    if Xres<0.11:
        print '!!! there are at least 1mm blank at beginning and end of row, or baseline removal may not work properly !!!'
    for i in range(Ystep):
      if Xres>0.11: #mm
        refwav[:]=(WaveformNew[i,1,:]+WaveformNew[i,-1,:])/2. # take 2nd points, as first may deviate from the correct locations
      elif Xres<0.11 and Xres>0.051:
        refwav[:]=(WaveformNew[i,3,:]+WaveformNew[i,-3,:]+WaveformNew[i,7,:]+WaveformNew[i,-7,:])/4.
      elif Xres<0.051:
        refwav[:]=(WaveformNew[i,5,:]+WaveformNew[i,-5,:]+WaveformNew[i,10,:]+WaveformNew[i,-10,:]+WaveformNew[i,15,:]+WaveformNew[i,-15,:])/6.
      H=np.argmax(refwav[0:wavlen/3]) ; L=np.argmin(refwav[0:wavlen/3]) ; refctr=(H+L)/2 ; Vpprefctr=refwav[H]-refwav[L]
      for j in range(Xstep):
        H=np.argmax(WaveformNew[i,j,0:wavlen/3]) ; L=np.argmin(WaveformNew[i,j,0:wavlen/3]) ; ctr=(H+L)/2
        Vppctr=WaveformNew[i,j,H]-WaveformNew[i,j,L] ; fac=Vppctr/Vpprefctr
        if refctr<ctr:
          shift=ctr-refctr
          late[shift]+=1
          WaveformNew[i,j,shift:]=WaveformNew[i,j,shift:]-fac*refwav[0:-shift] ; WaveformNew[i,j,0:shift]=WaveformNew[i,j,shift]
        else:
          shift=refctr-ctr
          early[shift]+=1
          WaveformNew[i,j,0:wavlen-shift]=WaveformNew[i,j,0:wavlen-shift]-fac*refwav[shift:]
          WaveformNew[i,j,wavlen-shift:]=WaveformNew[i,j,wavlen-shift-1]

  # use one waveform of each row as reference; aligned each waveform in that row with the reference at mid pt between max. and min. peaks in the first 1/3 of waveform
  # this is the best option of 1st gen
  elif TrendOff==1: 
    refwav=np.zeros(wavlen,dtype='<f')
    early=np.zeros(100,dtype=int) ; late=np.zeros(100,dtype=int) ; total=float(Xstep*Ystep)
    for i in range(Ystep):
      refwav[:]=WaveformNew[i,1,:]
      H=np.argmax(refwav[0:wavlen/3]) ; L=np.argmin(refwav[0:wavlen/3]) ; refctr=(H+L)/2
      for j in range(Xstep):
        H=np.argmax(WaveformNew[i,j,0:wavlen/3]) ; L=np.argmin(WaveformNew[i,j,0:wavlen/3]) ; ctr=(H+L)/2
        if refctr<ctr:
          shift=ctr-refctr
          late[shift]+=1
          WaveformNew[i,j,shift:]=WaveformNew[i,j,shift:]-refwav[0:-shift] ; WaveformNew[i,j,0:shift]=WaveformNew[i,j,shift]
        else:
          shift=refctr-ctr
          early[shift]+=1
          WaveformNew[i,j,0:wavlen-shift]=WaveformNew[i,j,0:wavlen-shift]-refwav[shift:]
          WaveformNew[i,j,wavlen-shift:]=WaveformNew[i,j,wavlen-shift-1]

        
  elif TrendOff==2: # also rescale refwav to have same Vpp as each waveform 
    refwav=np.zeros(wavlen,dtype='<f')
    early=np.zeros(100,dtype=int) ; late=np.zeros(100,dtype=int) ; total=float(Xstep*Ystep)
    for i in range(Ystep):
      refwav[:]=WaveformNew[i,1,:]
      H=np.argmax(refwav[0:wavlen/3]) ; L=np.argmin(refwav[0:wavlen/3])
      refctr=(H+L)/2 ; Vpprefctr=refwav[H]-refwav[L]
      for j in range(Xstep):
        H=np.argmax(WaveformNew[i,j,0:wavlen/3]) ; L=np.argmin(WaveformNew[i,j,0:wavlen/3])
        ctr=(H+L)/2 ; Vppctr=WaveformNew[i,j,H]-WaveformNew[i,j,L]
        fac=Vppctr/Vpprefctr # rescale refwav to have same Vpp as each waveform (assuming there is a system amp variation during scan)
        if refctr<ctr:
          shift=ctr-refctr
          late[shift]+=1
          WaveformNew[i,j,shift:]=WaveformNew[i,j,shift:]-fac*refwav[0:-shift] ; WaveformNew[i,j,0:shift]=WaveformNew[i,j,shift]
        else:
          shift=refctr-ctr
          early[shift]+=1
          WaveformNew[i,j,0:wavlen-shift]=WaveformNew[i,j,0:wavlen-shift]-fac*refwav[shift:]
          WaveformNew[i,j,wavlen-shift:]=WaveformNew[i,j,wavlen-shift-1]
    
  elif TrendOff==3 or TrendOff==4:
    
    time,tmp,n,delt,Tdat=open_asn_dat2(basedir,refname)
    if abs(Tdat-wavtlen)>0.2:
      print 'reference waveform time length not agree with data!'
      sys.exit(1)
    refwav=np.zeros(wavlen,dtype='<f')
    if n==wavlen:
      refwav[:]=tmp[:]
    elif n==2*wavlen:
      refwav[:]=tmp[::2]
    else:
      print 'invalid length of reference waveform for baseline removal!'
      sys.exit(1) 
    early=np.zeros(50,dtype=int) ; late=np.zeros(50,dtype=int) ; total=float(Xstep*Ystep)
    H=np.argmax(refwav[0:wavlen/3]) ; L=np.argmin(refwav[0:wavlen/3])
    refctr=(H+L)/2 ; Vpprefctr=refwav[H]-refwav[L]
    
    if TrendOff==3:
      for i in range(Ystep):
        for j in range(Xstep):
          H=np.argmax(WaveformNew[i,j,0:wavlen/3]) ; L=np.argmin(WaveformNew[i,j,0:wavlen/3]) ; ctr=(H+L)/2
          if refctr<ctr:
            shift=ctr-refctr
            late[shift]+=1
            WaveformNew[i,j,shift:]=WaveformNew[i,j,shift:]-refwav[0:-shift] ; WaveformNew[i,j,0:shift]=WaveformNew[i,j,shift]
          else:
            shift=refctr-ctr
            early[shift]+=1
            WaveformNew[i,j,0:wavlen-shift]=WaveformNew[i,j,0:wavlen-shift]-refwav[shift:]
            WaveformNew[i,j,wavlen-shift:]=WaveformNew[i,j,wavlen-shift-1]
   
    elif TrendOff==4:
      for i in range(Ystep):
        for j in range(Xstep):
          H=np.argmax(WaveformNew[i,j,0:wavlen/3]) ; L=np.argmin(WaveformNew[i,j,0:wavlen/3])
          ctr=(H+L)/2 ; Vppctr=WaveformNew[i,j,H]-WaveformNew[i,j,L]
          fac=Vppctr/Vpprefctr # rescale refwav to have same Vpp as each waveform
          if refctr<ctr:
            shift=ctr-refctr
            late[shift]+=1
            WaveformNew[i,j,shift:]=WaveformNew[i,j,shift:]-fac*refwav[0:-shift] ; WaveformNew[i,j,0:shift]=WaveformNew[i,j,shift]
          else:
            shift=refctr-ctr
            early[shift]+=1
            WaveformNew[i,j,0:wavlen-shift]=WaveformNew[i,j,0:wavlen-shift]-fac*refwav[shift:]
            WaveformNew[i,j,wavlen-shift:]=WaveformNew[i,j,wavlen-shift-1]
  
  if TrendOff!=0 and TrendOff!=5:
    print
    print 'stat of baseline removal'
    print 'shift   late     %     no/early    %'
    for i in range(50):
      if late[i]!=0 or early[i]!=0:
        print '%2d %10d %6.3f %10d %6.3f' %(i, late[i],100.*(late[i]/total),early[i],100.*(early[i]/total))

# Xstep,Ystep may be updated on return. CscanAmpNew is created but contains just zeros 
  return WaveformNew,CscanAmpNew,Xnew,Ynew,pos,Xstep,Ystep

# _______________________________________________________________________________________________________________________________________________________________________

def HeaderInfo(header):
  
# get key parameters from the header
# last update: 17FEB2013

  for item in header:
    b=item.strip()
    a=b.split()
    if a[0]=='RSDLAMP:':
      wavt=float(a[1])
    elif a[0]=='RESOLUTION:':
      xres=float(a[1])
    elif a[0]=='Y_RESOLUTION:':
      yres=float(a[1])
    elif a[0]=='XSTEPS:':
      xstep=int(float(a[1])) # can't directly int(a[1]) of string a[1] ?!
    elif a[0]=='YSTEPS:':
      ystep=int(float(a[1]))
    elif a[0]=='XMIN:':
      xmin=float(a[1])
    elif a[0]=='XMAX:':
      xmax=float(a[1])
    elif a[0]=='YMIN:':
      ymin=float(a[1])
    elif a[0]=='YMAX:':
      ymax=float(a[1])
    elif a[0]=='SUBSAMPLING:':
      subsamp=int(float(a[1]))     
    elif a[0]=='AVERAGES:':
      navg=int(float(a[1]))  
    elif a[0]=='XSPEED:':
      xspeed=float(a[1])
    elif a[0]=='SCAN_NAME:':
      scantype=b[10:].strip()
    elif a[0]=='AXIS1:':
      axis1=b[7:].strip()
      
  return wavt,xres,yres,xstep,ystep,xmin,xmax,ymin,ymax,subsamp,navg,xspeed,scantype,axis1

# _______________________________________________________________________________________________________________________________________________________________________ 

def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = np.arange(len(v))
    
    v = np.asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN
    
    lookformax = True
    
    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
    
    #return maxtab, mintab # return as list for easier processing 14JUL2013
    return np.array(maxtab), np.array(mintab) 

# _______________________________________________________________________________________________________________________________________________________________________   

def MatchedFilter(delt,k,Waveform,Xstep,Ystep,PeakBin):
  
# perform matched filter using filter templates from paer L area of GF plate  (IU/CRC 2013.4 proect)
# last update: 03APR2014
  
  # paper L in "AL foil & paper on back of 5mm glass fiber plate Focus@FS then 1mm down (100ps 4096pt 1avg no purg 0.5mm res IU SP project).tvl"
  refcord = [[22,14],[24,15],[26,16],[28,17],[30,18],[32,16],[34,15],[35,14],[36,13],[37,12]]
  nref=len(refcord)
  
  # only take data in gate; k indicate which signal
  n=PeakBin[4,k,0,0]-PeakBin[3,k,0,0]+1
  tart=np.zeros(n,dtype='<f') ; reft=np.zeros(n,dtype='<f') ; reff=np.zeros((nref,n/2),complex) ; outt=np.zeros((nref,n),dtype='<f')
  out=np.zeros(nref)
  Cscan=np.zeros((Ystep,Xstep),dtype='<f')

  # prepare filters in freq. domain
  for ii in range(nref):
    i=refcord[ii][0] ; j=refcord[ii][1]
    reft[0:]=Waveform[i,j,PeakBin[3,k,i,j]:PeakBin[4,k,i,j]+1]
    reff[ii,:]=ReFFT2(reft[:],delt)
  
  for i in range(Ystep):
    for j in range(Xstep):
      tart[0:]=Waveform[i,j,PeakBin[3,k,i,j]:PeakBin[4,k,i,j]+1]
      tarf=ReFFT2(tart,delt)
      for ii in range(nref):
        outt[ii,:]=IReFFT2(tarf * np.conj(reff[ii,:]),delt)
        out[ii]=np.amax(outt[ii,:])
      out.sort() # in-place sort
      
      Cscan[i,j]=np.prod(out[nref:nref-4:-1])
      
  return Cscan

# _______________________________________________________________________________________________________________________________________________________________________
  
def ReFFT2(Rein,delt):
  
  n=len(Rein) ; Rein=np.asarray(Rein,dtype=np.float)   # add this line 06NOV2011 as rfft in Python 2.7 does not support >f4 type
  ff=np.zeros(n) ; fout=np.zeros(n/2,complex) ; Reout=np.zeros(n/2) ; Imout=np.zeros(n/2)
  ff=rfft(Rein)
  
  Reout[1:]=ff[range(1,n-2,2)] ; Reout[0]=ff[0]
  Imout[1:]=ff[range(2,n-1,2)] ; Imout[0]=ff[-1]
  fout.real=Reout*delt ; fout.imag=Imout*delt
  return fout
  
# _______________________________________________________________________________________________________________________________________________________________________  
  
def IReFFT2(fin,delt):
    
  # rewritten 20JUN2010, 06NOV2011
  
  n=len(fin)*2
  Reout2=np.zeros(n)
  
  Reout2[range(1,n-2,2)]=fin.real[1:] ; Reout2[0]=fin.real[0]
  Reout2[range(2,n-1,2)]=fin.imag[1:] ; Reout2[-1]=fin.imag[0]
  Reout2=np.asarray(Reout2,dtype=np.float)   # add this line 06NOV2011 as rfft in Python 2.7 does not support >f4 type
  
  Reout2=irfft(Reout2)/delt
  
  return Reout2

# _______________________________________________________________________________________________________________________________________________________________________

def moving_average(x, n, type='simple',rescale='yes'):
    """
    compute an n period moving average.

    type is 'simple' | 'exponential'
    recale is "yes' | 'no'
    
    modified from an example in Matplotlib 16APR2011 
    update: 05JUL2011

    """
    x = np.asarray(x)
    if type=='simple':
        weights = np.ones(n)
    else:
        weights = np.exp(np.linspace(-1., 0., n))

    weights /= weights.sum()

    nn=len(x)
    a =  np.convolve(x, weights, mode='full')[:nn]
    # shift back the delay and let the end nn point be tha same as the last valid point
    a[0:nn-n/2] = a[n/2:nn] ; a[nn-n/2+1:]=a[nn-n/2]  
    
    if rescale=='yes': a=a*(np.amax(x)-np.min(x))/(np.amax(a)-np.min(a))
        
    return a

# _______________________________________________________________________________________________________________________________________________________________________

class DataFile(object):
  def __init__(self, fname, basedir=None):
    if basedir is not None:
      fname = os.path.join(basedir, fname)
    with open(fname, 'rb') as fobj:
      self.header = list(iter(fobj.readline, "ENDOFHEADER\r\n"))
      # the first 4-byte word after the header gives the length of each row
      firstword = fobj.read(4)
      # convert the 4-byte string to a float
      col_size = struct.unpack(">f", firstword)[0]
      # move the read point back so we can read the first row in it's entirety
      fobj.seek(-4,1)
            
      # define a compound datatype for the data: the two coordinate values
      # followed by the THz waveform
      bin_dtype = [("x", ">f"),("y", ">f"), ("waveform", ">f", col_size-2)]
            
      # read the data into an array
      self.data = np.fromfile(fobj, bin_dtype)
  
# _______________________________________________________________________________________________________________________________________________________________________   

def open_asn_dat2(basedir,fname):  # modified 19NOV2013
    
  fmt=fname[-3:] # assume standard 3-character extension
  #basedir = os.getcwd()
  infile=open(os.path.join(basedir,fname), 'r')
    
  x=[] ; y=[]
  count=0
  if fmt=='asn':
    for line in infile:
      if count>0:  # skip the first line which contains the column headers
        xval,yval,t1,t2,t3,t4 = line.split()
        x.append(float(xval)) ; y.append(float(yval))
      count+=1
  elif fmt=='dat':
    for line in infile:
      xval,yval=line.split()
      x.append(float(xval)) ; y.append(float(yval))   
    
  x = np.array(x) ; y = np.array(y) ; n=len(y) ;  totalT=x[-1]-x[0] ; delt=totalT/(n-1)
  infile.close()
  return x,y,n,delt,totalT
  
# _______________________________________________________________________________________________________________________________________________________________________


if __name__ == '__main__':
  main()