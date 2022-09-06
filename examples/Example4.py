"""Execution time vs ntoas
"""

import enterprise_GWecc as GWecc
import numpy as np
import matplotlib.pyplot as plt

year = 365.25 * 24 * 3600
ns = 1e-9

ntoas = 10
toas = 365.25 * np.linspace(0, 15, ntoas)

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
e0 = 0.5
gamma0 = 0
l0 = 0
t0 = 0

z = 0.0

plt.subplot(221)
plt.suptitle("Comparing numerical and analytic residuals")

import time


res_N = GWecc.GWecc.EccentricResiduals(
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

runtime_num = []
runtime_adb = []
runtime_anl = []
ntoass = [10, 50, 100, 500, 1000, 5000]
for ntoas in ntoass:
    toas = 365.25 * np.linspace(0, 15, ntoas)
    
    start = time.time()
    res_N = [
        GWecc.GWecc.EccentricResiduals(
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
        for i in range(10000 // ntoas)
    ]
    end = time.time()
    runtime_num.append((end - start) / 10000)
    
    start = time.time()
    res_N = [
        GWecc.GWecc.EccentricResiduals(
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
            GWecc.GWecc.ResidualsMethod_Adb,
            GWecc.GWecc.ResidualsTerms_Earth,
            toas,
        )
        for i in range(10000 // ntoas)
    ]
    end = time.time()
    runtime_adb.append((end - start) / 10000)

    start = time.time()
    res_N = [
        GWecc.GWecc.EccentricResiduals(
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
            GWecc.GWecc.ResidualsMethod_Anl,
            GWecc.GWecc.ResidualsTerms_Earth,
            toas,
        )
        for i in range(10000 // ntoas)
    ]
    end = time.time()
    runtime_anl.append((end - start) / 10000)
    
plt.loglog(ntoass, runtime_adb, label="adb")
plt.loglog(ntoass, runtime_anl, label="anl")
plt.loglog(ntoass, runtime_num, label="num")
plt.legend()
plt.show()