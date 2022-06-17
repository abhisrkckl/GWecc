#!/usr/bin/python3

import numpy as np
from enterprise_GWecc import GWecc, eccentric_cw_delay
from enterprise.pulsar import Pulsar
import matplotlib.pyplot as plt

year = 365.25*24*3600
day = 24*3600
parsec = 102927125.0
MSun = 4.92703806e-6
ns = 1e-9

def hms_to_rad(hh,mm,ss):
    sgn = np.sign(hh)
    return sgn * (sgn*hh + mm/60 + ss/3600) * np.pi/12

def dms_to_rad(dd,mm,ss):
    sgn = np.sign(dd)
    return sgn * (sgn*dd + mm/60 + ss/3600) * np.pi/180

M = 1e8
q = 1
Pb0 = 5    # years
Omega = 0
i = 0
l0 = gamma0 = 0
z = 0

psr = Pulsar("J0437-4715.IPTADR2/J0437-4715.fit.par", "J0437-4715.IPTADR2/J0437-4715.IPTADR2.tim")

RA_P  = psr.phi                 # rad
DEC_P = np.pi/2 - psr.theta     # rad
D_P = psr.pdist[0]*1000         # pc
     
RA_GW = hms_to_rad(  4,0,0)
DEC_GW = dms_to_rad(-45,0,0)
D_GW = 1e9 # pc 

toas = psr.toas / (24*3600)      # MJD
t0   = min(toas)                 # MJD

for idx,e0 in enumerate([0.1,0.5,0.8]):

    ax = plt.subplot(311+idx)
    
    # Calling EccentricResiduals function through the enterprise interface
    ecc_gw = eccentric_cw_delay( cos_gwtheta = np.sin(DEC_GW), 
                                 gwphi = RA_GW,
                                 log10_dist = np.log10(D_GW/1e6),
                                 log10_h=None,
                                 psi = Omega, 
                                 cos_inc = np.cos(i),
                                 log10_M = np.log10(M), 
                                 q = q,
                                 log10_F = np.log10(2/(Pb0*year)), 
                                 e0 = e0,  
                                 gamma0 = gamma0, 
                                 l0 = l0, 
                                 tref = max(psr.toas),
                                 z = z,
                                 p_dist=0,
                                 psrTerm=True,
                                 evolve=True)
    ecc_gw_fn = ecc_gw("ecc_gw", psr=psr)
    res2 = ecc_gw_fn()
    
    #res = np.asarray(res) - np.mean(res)
    res2 = np.asarray(res2) - np.mean(res2)
    
    #plt.plot(toas/365.25, res/ns, marker='x', ls='dotted')
    plt.plot(toas/365.25, res2/ns, marker='x', ls='', c='red')
    plt.ylabel("$R(t)$  (ns)", fontsize=14)
    plt.text(t0/365.25,0.8*ax.yaxis.get_data_interval()[1],"e="+str(e0),size=15, ha="center", va="center", bbox=dict(boxstyle="round",facecolor='cyan',alpha=0.1))
    if idx<2:
        plt.tick_params(labelbottom=False,labelsize=12)
    else:
        plt.tick_params(labelsize=12)
        plt.xlabel("$t$ (years)", fontsize=14)

plt.suptitle("PSR J0437-4715")

plt.show()
