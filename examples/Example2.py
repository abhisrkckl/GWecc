#!/usr/bin/python3

import enterprise_GWecc as GWecc
import numpy as np
import matplotlib.pyplot as plt

year = 365*24*3600
ns = 1e-9

toas = year*np.linspace(0,15,1000)

RA_P = 0.65
DEC_P = -1.12
D_P = 700
	 
RA_GW = np.pi/3
DEC_GW = -np.pi/3.5
D_GW = 1e9

Omega = 0
i = np.pi/3

M = 1e10
q = 1

Pb = 5*year
n0 = 2*np.pi/Pb
e0 = 0.5
gamma0 = 0
l0 = 0
t0 = 0

residuals = GWecc.GWecc.EccentricResiduals(M, q,
                                           Omega, i,
                                           t0, n0, e0, l0, gamma0,
                                           D_GW, RA_GW, DEC_GW, 
                    					   D_P, RA_P, DEC_P, 
					                       GWecc.GWecc.ResidualsMethod_Num,
                    					   GWecc.GWecc.ResidualsTerms_Earth,
					                       toas);

Rp,Rx = GWecc.GWecc.EccentricResiduals_px(M, q,
                                          Omega, i,
                                          t0, n0, e0, l0, gamma0,
                                          D_GW,
                                          GWecc.GWecc.ResidualsMethod_Num,
                                          toas)

cosmu,Fp,Fx = GWecc.GWecc.AntennaPattern(RA_GW, DEC_GW, RA_P, DEC_P)
residuals1 = np.array(Rp)*Fp + np.array(Rx)*Fx;

"""residuals = GWecc.eccentric_cw_delay(toas, 
	     	       		     RA_psr, DEC_psr, D_psr,
	     	       		     RA_GW, DEC_GW, log10_D_GW,
	     	       		     psi, i,
	     	       		     log10_M, q,
	     	       		     log10_f0_GW, e0, gamma0, l0, t0,
	     	       		     psrTerm=False)"""

#plt.plot(toas/Pb, np.asarray(residuals)/ns)
#plt.plot(toas/Pb, np.asarray(residuals1)/ns)
plt.plot(toas/Pb, np.asarray(residuals)/np.asarray(residuals1))
plt.xlabel("t/Pb")
plt.ylabel("$\Delta_{GW}$ (ns)")
plt.grid()
plt.show()

