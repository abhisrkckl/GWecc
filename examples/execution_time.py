"""Execution time vs ntoas 
"""

import enterprise_GWecc as GWecc
import numpy as np
import matplotlib.pyplot as plt
import time

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
#e0 = 0.5
gamma0 = 0
l0 = 0
t0 = 0

z = 0.0


res_N = GWecc.GWecc.EccentricResiduals(
    M,
    q,
    Omega,
    i,
    t0,
    Pb0,
    0.5,
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
runtime_pm = []

ntoas = 10000
nreps = 2

start = time.time()
res_N = [
    GWecc.GWecc.EccentricResiduals(
        M,
        q,
        Omega,
        i,
        t0,
        Pb0,
        0,
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
    for i in range(nreps)
]
end = time.time()
runtime_circ = (end-start)  / (nreps*ntoas)

es = np.linspace(0.001, 0.8, 10)
for e0 in es:
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
        for i in range(nreps)
    ]
    end = time.time()
    runtime_num.append((end - start) / (nreps*ntoas))

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
        for i in range(nreps)
    ]
    end = time.time()
    runtime_adb.append((end - start) / (nreps*ntoas))

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
        for i in range(nreps)
    ]
    end = time.time()
    runtime_anl.append((end - start) / (nreps*ntoas))

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
            GWecc.GWecc.ResidualsMethod_PM,
            GWecc.GWecc.ResidualsTerms_Earth,
            toas,
        )
        for i in range(nreps)
    ]
    end = time.time()
    runtime_pm.append((end - start) / (nreps*ntoas))

plt.plot(es, np.array(runtime_adb)/runtime_circ, marker='+', ls='--', label="Analytic")
plt.plot(es, np.array(runtime_pm)/runtime_circ, marker='x', ls='--', label="Fourier")
plt.plot(es[es<=0.3], np.array(runtime_anl)[es<=0.3]/runtime_circ, marker='o', ls='--', label="Post-circular")
plt.plot(es, np.array(runtime_num)/runtime_circ, marker='^', ls='--', label="Numerical")
plt.yscale('log')
plt.xlabel("$e_{t0}$", fontsize=12)
plt.ylabel("Execution time per TOA\n(Normalized by execution time for circular)", fontsize=12)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=11)
plt.show()
