import numpy as np
from enterprise_GWecc import GWecc
import matplotlib.pyplot as plt
import faulthandler; faulthandler.enable()

year = 365.25*24*3600
parsec = 102927125.0
MSun = 4.92703806e-6

def PN_periastron_advance(M, q, n, e):
	ipow = np.power
	cbrt = np.cbrt
	eta = np.abs(1-q)/(1+q)
	xi = cbrt( ipow(M*n, 2) )
	OTS = np.sqrt(1-e**2)
	M_PI = np.pi
	return     (3*xi/ipow(OTS,2) 
                    + ((78 + ipow(e,2)*(51 - 26*eta) - 28*eta)*ipow(xi,2))/(4.*ipow(OTS,4)) 
                    + ((18240 - 25376*eta + 896*ipow(eta,2) + ipow(e,4)*(2496 - 1760*eta + 1040*ipow(eta,2)) 
                        + (1920 + ipow(e,2)*(3840 - 1536*eta) - 768*eta)*OTS + 492*eta*ipow(M_PI,2) 
                        + ipow(e,2)*(28128 - 27840*eta + 5120*ipow(eta,2) 
                        + 123*eta*ipow(M_PI,2)))*ipow(xi,3))/(128.*ipow(OTS,6))	)

M = 1e9*MSun
q = 1
Pb = 1.5*year
n0 = 2*np.pi/Pb
DGW = 1e9*parsec
Omega = 0
i = 0
t0 = 0
l0 = gamma0 = 0
ts = np.linspace(0,10*year,5000)

hpx0 = np.zeros([5000,2])
hpx1 = np.zeros([5000,2])

e0 = 0.1
k = PN_periastron_advance(M,q,n0,e0)
for idx,t in enumerate(ts):
	l = n0*t
	gamma = k*n0*t
	
	hp,hx = GWecc.EccentricWaveform_fn(M, q, Omega, i, n0, e0, l, gamma, DGW)
	
	hpx0[idx] = hp,hx
	

hp1, hx1 = GWecc.EccentricWaveform_px(M/MSun, q, Omega, i, t0, Pb/year, e0, l0, gamma0, DGW, 0, ts/(24*3600))
ax = plt.subplot(321)
plt.plot(ts/year, hp1, label="$h_+$ (Conservative+Reactive)")
plt.plot(ts/year, hpx0[:,0], 'r--', label="$h_+$ (Conservative only)")
plt.grid()
plt.xlim([0,10])
#plt.ylim([-6e-11,6e-11])
#plt.legend(loc='upper right',fontsize=11)
#plt.xlabel("$t$ (years)")
plt.ylabel("$h_{+}$", fontsize=14)
plt.text(1,0.8*ax.yaxis.get_data_interval()[1],"e="+str(e0),size=15, ha="center", va="center", bbox=dict(boxstyle="round",facecolor='cyan',alpha=0.1))
plt.tick_params(labelbottom=False,labelsize=12)

plt.subplot(322)
plt.plot(ts/year, hx1, label="$h_\\times$ (Conservative+Reactive)")
plt.plot(ts/year, hpx0[:,1], 'r--', label="$h_\\times$ (Conservative only)")
plt.ylabel("$h_{\\times}$", fontsize=14)
plt.grid()
plt.xlim([0,10])
#plt.ylim([-6e-11,6e-11])
#plt.legend(loc='upper right',fontsize=11)
plt.tick_params(labelbottom=False,labelsize=12)
#plt.xlabel("$t$ (years)")
#plt.ylabel("$h_{+,\\times}$")

#################################################################################################

e0 = 0.5
k = PN_periastron_advance(M,q,n0,e0)
for idx,t in enumerate(ts):
	l = n0*t
	gamma = k*n0*t
	
	hp,hx = GWecc.EccentricWaveform_fn(M, q, Omega, i, n0, e0, l, gamma, DGW)
	
	hpx0[idx] = hp,hx
hp1, hx1 = GWecc.EccentricWaveform_px(M/MSun, q, Omega, i, t0, Pb/year, e0, l0, gamma0, DGW, 0, ts/(24*3600))
ax = plt.subplot(323)
plt.plot(ts/year, hp1, label="$h_+$ (Conservative+Reactive)")
plt.plot(ts/year, hpx0[:,0], 'r--', label="$h_+$ (Conservative only)")
plt.grid()
plt.xlim([0,10])
#plt.ylim([-6e-11,6e-11])
#plt.legend(loc='upper right',fontsize=11)
#plt.xlabel("$t$ (years)")
plt.ylabel("$h_{+}$", fontsize=14)
plt.text(1,0.8*ax.yaxis.get_data_interval()[1],"e="+str(e0),size=15, ha="center", va="center", bbox=dict(boxstyle="round",facecolor='cyan',alpha=0.1))
plt.tick_params(labelbottom=False,labelsize=12)

plt.subplot(324)
plt.plot(ts/year, hx1, label="$h_\\times$ (Conservative+Reactive)")
plt.plot(ts/year, hpx0[:,1], 'r--', label="$h_\\times$ (Conservative only)")
plt.ylabel("$h_{\\times}$", fontsize=14)
plt.grid()
plt.xlim([0,10])
#plt.ylim([-7e-11,7e-11])
#plt.legend(loc='upper right',fontsize=11)
plt.tick_params(labelbottom=False,labelsize=12)

e0 = 0.8
k = PN_periastron_advance(M,q,n0,e0)
for idx,t in enumerate(ts):
	l = n0*t
	gamma = k*n0*t
	
	hp,hx = GWecc.EccentricWaveform_fn(M, q, Omega, i, n0, e0, l, gamma, DGW)
	
	hpx0[idx] = hp,hx
hp1, hx1 = GWecc.EccentricWaveform_px(M/MSun, q, Omega, i, t0, Pb/year, e0, l0, gamma0, DGW, 0, ts/(24*3600))
ax = plt.subplot(325)
plt.plot(ts/year, hp1, label="$h_+$ (Conservative+Reactive)")
plt.plot(ts/year, hpx0[:,0], 'r--', label="$h_+$ (Conservative only)")
plt.grid()
plt.xlim([0,10])
#plt.ylim([-6e-11,6e-11])
#plt.legend(loc='upper right',fontsize=11)
#plt.xlabel("$t$ (years)")
plt.ylabel("$h_{+}$", fontsize=14)
plt.xlabel("$t$ (years)", fontsize=13)
plt.text(1,0.8*ax.yaxis.get_data_interval()[1],"e="+str(e0),size=15, ha="center", va="center", bbox=dict(boxstyle="round",facecolor='cyan',alpha=0.1))
plt.tick_params(labelsize=12)

plt.subplot(326)
plt.plot(ts/year, hx1, label="Conservative+Reactive")
plt.plot(ts/year, hpx0[:,1], 'r--', label="Conservative only")
plt.grid()
plt.xlim([0,10])
#plt.ylim([-6e-11,6e-11])
plt.legend(loc='upper left',fontsize=11)
plt.xlabel("$t$ (years)", fontsize=13)
plt.ylabel("$h_{\\times}$", fontsize=14)
plt.tick_params(labelsize=12)

plt.show()
