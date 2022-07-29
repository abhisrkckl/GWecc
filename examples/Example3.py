#!/usr/bin/python3

import enterprise_GWecc as GWecc
import numpy as np
import matplotlib.pyplot as plt

year = 365.25 * 24 * 3600
ns = 1e-9
# parsec = 102927125.0

toas = 365.25 * np.linspace(0, 4, 1000)

RA_P = 1.65
DEC_P = -1.12
D_P = 1000

RA_GW = 3.1
DEC_GW = -np.pi / 3.5
D_GW = 1e9

Omega = 0
i = np.pi / 3

M = 5e9
q = 1

Pb0 = 2
e0 = 0.01
gamma0 = 0
l0 = 0
t0 = 0

z = 0.0

cosmu, Fp, Fx = GWecc.GWecc.AntennaPattern(RA_GW, DEC_GW, RA_P, DEC_P)
delay = -D_P * (1 - cosmu) / (1 + z)
res_px = GWecc.GWecc.EccentricResiduals_px(
    M,
    q,
    Omega,
    i,
    t0,
    Pb0,
    e0,
    l0,
    gamma0,
    D_GW,
    delay,
    z,
    GWecc.GWecc.ResidualsMethod_Num,
    GWecc.GWecc.ResidualsTerms_Earth,
    toas,
)


res_from_px = np.array(res_px[0]) * Fp + np.array(res_px[1]) * Fx
res = GWecc.GWecc.EccentricResiduals(
    M,
    q,
    Omega,
    i,
    t0,
    Pb0,
    e0,
    l0,
    gamma0,
    D_GW,
    RA_GW,
    DEC_GW,
    D_P,
    RA_P,
    DEC_P,
    z,
    GWecc.GWecc.ResidualsMethod_Num,
    GWecc.GWecc.ResidualsTerms_Earth,
    toas,
)

plt.subplot(211)

plt.plot(toas / 365.25, np.asarray(res_px[0]) / ns - np.mean(res_px[0]) / ns, label="+")
plt.plot(toas / 365.25, np.asarray(res_px[1]) / ns - np.mean(res_px[1]) / ns, label="x")
plt.xlabel("t (year)")
plt.ylabel("$\Delta_{+,\\times}$ (ns)")
plt.legend()
plt.grid()

plt.subplot(212)
plt.plot(toas / 365.25, np.asarray(res_from_px) / ns - np.mean(res_from_px) / ns)
plt.plot(toas / 365.25, np.asarray(res) / ns - np.mean(res) / ns)
plt.xlabel("t (year)")
plt.ylabel("$\Delta_{GW}$ (ns)")
plt.grid()

plt.show()
