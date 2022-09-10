"""Fe statistic function plots
Singh et al (In prep)
"""

import enterprise_GWecc as GWecc
import numpy as np
import matplotlib.pyplot as plt

year = 365.25 * 24 * 3600
ns = 1e-9

ntoas = 1000
toas = 365.25 * np.linspace(0, 15, ntoas)

RA_P = 1.65
DEC_P = -1.12
D_P = 1000

RA_GW = 3.1
DEC_GW = -np.pi / 3.5
D_GW = 1e9

M = 5e9
q = 1

Omega = 0
i = np.pi / 3

Pb0 = 2
e0 = 0.2
gamma0 = 0
l0 = 0
tref = 0

z = 0.0

res_N = GWecc.GWecc.EccentricResiduals(
    M,
    q,
    Omega,
    i,
    tref,
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

As1 = GWecc.GWecc.FeStatFuncs(
    M, q, tref, Pb0, e0, l0, D_GW, RA_GW, DEC_GW, RA_P, DEC_P, z, toas, GWecc.GWecc.ResidualsMethod_Num
)

As2 = GWecc.GWecc.FeStatFuncs(
    M, q, tref, Pb0, e0, l0, D_GW, RA_GW, DEC_GW, RA_P, DEC_P, z, toas, GWecc.GWecc.ResidualsMethod_Adb
)

plt.subplot(311)
plt.suptitle("Functions needed for Fe statistics")
res_N = np.asarray(res_N) - np.average(res_N)
plt.plot(toas / 365.25, res_N / ns, label="PTA signal (Earth term) (ns)")
plt.legend()

plt.subplot(312)
for idx, Ai in enumerate(As1):
    Ai -= np.mean(Ai)
    plt.plot(toas / 365.25, Ai, label="$A_{%d}$ (Earth term) (Num)" % (idx + 1))

plt.subplot(313)
for idx, Ai in enumerate(As2):
    plt.plot(toas / 365.25, Ai, label="$A_{%d}$ (Earth term) (Adb)" % (idx + 1))

plt.legend()
plt.xlabel("t (yr)")

plt.show()
