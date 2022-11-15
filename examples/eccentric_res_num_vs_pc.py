"""Comparing numerical and post-circular residuals (High eccentricity)
    -- needs fixing.
"""

from enterprise_GWecc.GWecc import (
    eccentric_residuals,
    ResidualsMethod_PC,
    ResidualsMethod_Num,
    ResidualsTerms_Earth,
    ResidualsTerms_Pulsar,
)
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

res_N = eccentric_residuals(
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
    ResidualsMethod_Num,
    ResidualsTerms_Earth,
    toas,
)
res_A = eccentric_residuals(
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
    ResidualsMethod_PC,
    ResidualsTerms_Earth,
    toas,
)
res_N = np.array(res_N) + res_A[0]
plt.plot(toas / 365.25, np.asarray(res_A) / ns, label="PC")
plt.plot(toas / 365.25, np.asarray(res_N) / ns, label="Num")
plt.xlabel("t (year)")
plt.ylabel("$\Delta_{GW}$ (ns)")
plt.legend()
plt.grid()

plt.subplot(222)
res_N = eccentric_residuals(
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
    ResidualsMethod_Num,
    ResidualsTerms_Pulsar,
    toas,
)
res_A = eccentric_residuals(
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
    ResidualsMethod_PC,
    ResidualsTerms_Pulsar,
    toas,
)
res_N = np.array(res_N) + res_A[0]
plt.plot(toas / 365.25, np.asarray(res_A) / ns)
plt.plot(toas / 365.25, np.asarray(res_N) / ns)
plt.xlabel("t (year)")
plt.ylabel("$\Delta_{GW}$ (ns)")
plt.grid()

plt.subplot(223)
e0 = 0.5
plt.suptitle("Comparing numerical and post-circular residuals (e=0.5)")
res_N = eccentric_residuals(
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
    ResidualsMethod_Num,
    ResidualsTerms_Earth,
    toas,
)
res_A = eccentric_residuals(
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
    ResidualsMethod_PC,
    ResidualsTerms_Earth,
    toas,
)
res_N = np.array(res_N) + res_A[0]
plt.plot(toas / 365.25, np.asarray(res_A) / ns)
plt.plot(toas / 365.25, np.asarray(res_N) / ns)
plt.xlabel("t (year)")
plt.ylabel("$\Delta_{GW}$ (ns)")
plt.grid()

plt.subplot(224)
res_N = eccentric_residuals(
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
    ResidualsMethod_Num,
    ResidualsTerms_Pulsar,
    toas,
)
res_A = eccentric_residuals(
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
    ResidualsMethod_PC,
    ResidualsTerms_Pulsar,
    toas,
)
res_N = np.array(res_N) + res_A[0]
plt.plot(toas / 365.25, np.asarray(res_A) / ns)
plt.plot(toas / 365.25, np.asarray(res_N) / ns)
plt.xlabel("t (year)")
plt.ylabel("$\Delta_{GW}$ (ns)")
plt.grid()

plt.show()
