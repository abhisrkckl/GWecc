"""Fe statistic function plots
Singh et al (In prep)
"""

from turtle import color
import enterprise_GWecc as GWecc
import numpy as np
import matplotlib.pyplot as plt

year = 365.25 * 24 * 3600
ns = 1e-9

ntoas = 1000
toas = 365.25 * np.linspace(0, 10, ntoas)

RA_P = np.pi / 2
DEC_P = -np.pi / 3
D_P = 1000

RA_GW = np.pi
DEC_GW = -np.pi / 4
D_GW = 1e9

M = 5e9
q = 1

Omega = 0
i = np.pi / 3

Pb0 = 5
e0 = 0.2
gamma0 = 0
l0 = 0
tref = 0

z = 0.0

e0s = [0.1, 0.3, 0.6]
for jdx, e0 in enumerate(e0s):
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

    As = GWecc.GWecc.FeStatFuncs(
        M,
        q,
        tref,
        Pb0,
        e0,
        l0,
        D_GW,
        RA_GW,
        DEC_GW,
        RA_P,
        DEC_P,
        z,
        toas,
        GWecc.GWecc.ResidualsMethod_Adb,
    )

    plt.subplot(331 + jdx)
    res_N = np.asarray(res_N) - np.average(res_N)
    plt.plot(toas / 365.25, res_N / ns)
    plt.title(f"e0 = {e0}", fontsize=14)
    plt.axhline(0, color="grey", ls="--")
    if jdx == 0:
        plt.ylabel("PTA signal (Earth term) (ns)", fontsize=14)
    plt.tick_params(axis="both", labelsize=11)

    plt.subplot(334 + jdx)
    for idx, Ai in enumerate(As[:3]):
        plt.plot(toas / 365.25, Ai, label="$A_{%d}$" % (idx + 1))
    if jdx == 0:
        plt.ylabel("$A_i$ (Earth term)", fontsize=14)
        plt.legend(framealpha=0.3)
    plt.axhline(0, color="grey", ls="--")
    plt.tick_params(axis="both", labelsize=11)

    plt.subplot(337 + jdx)
    for idx, Ai in enumerate(As[3:]):
        plt.plot(toas / 365.25, Ai, label="$A_{%d}$" % (idx + 4))
    if jdx == 0:
        plt.ylabel("$A_i$ (Earth term)", fontsize=14)
        plt.legend(framealpha=0.3)
    plt.xlabel("t (yr)", fontsize=14)
    plt.axhline(0, color="grey", ls="--")
    plt.tick_params(axis="both", labelsize=11)

plt.show()
