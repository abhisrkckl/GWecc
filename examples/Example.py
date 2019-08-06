#!/usr/bin/python3

import enterprise_GWecc as GWecc
import numpy as np
import matplotlib.pyplot as plt

year = 365.25*24*3600
ns = 1e-9

toas = year*np.linspace(0,15,1000)

RA_psr = 0.65
DEC_psr = -1.12
D_psr = 700
	 
RA_GW = np.pi/3
DEC_GW = np.pi/3.5
log10_D_GW = 9

psi = 0
i = np.pi/3

log10_M = 10
q = 1

log10_f0_GW = -8
e0 = 0.01
gamma0 = 0
l0 = 0
t0 = 0

residuals = GWecc.eccentric_cw_delay(toas, 
	             	       		     RA_psr, DEC_psr, D_psr,
	             	       		     RA_GW, DEC_GW, log10_D_GW,
	             	       		     psi, i,
	             	       		     log10_M, q,
	             	       		     log10_f0_GW, e0, gamma0, l0, t0,
	             	       		     psrTerm=False)

plt.plot(toas/year, np.asarray(residuals)/ns)
plt.xlabel("t (yr)")
plt.ylabel("$\Delta_{GW}$ (ns)")
plt.show()

