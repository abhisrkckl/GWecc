"""Figure 5 of Susobhanan+ 2020
+/x polarizations of Earth & Pulsar terms (adiabatic)
"""

import numpy as np
from enterprise_GWecc import GWecc
import matplotlib.pyplot as plt

year = 365.25 * 24 * 3600
day = 24 * 3600
parsec = 102927125.0
MSun = 4.92703806e-6
ns = 1e-9


def hms_to_rad(hh, mm, ss):
    sgn = np.sign(hh)
    return sgn * (sgn * hh + mm / 60 + ss / 3600) * np.pi / 12


def dms_to_rad(dd, mm, ss):
    sgn = np.sign(dd)
    return sgn * (sgn * dd + mm / 60 + ss / 3600) * np.pi / 180


M = 1e9
q = 1
Pb0 = 1.5  # years
# n0 = 2*np.pi/Pb
Omega = 0
i = 0
t0 = 0
l0 = gamma0 = 0
z = 0

RA_P = hms_to_rad(4, 37, 15.81476)
DEC_P = dms_to_rad(-47, 15, 8.6242)
D_P = 156.3  # pc

RA_GW = hms_to_rad(4, 0, 0)
DEC_GW = dms_to_rad(-45, 0, 0)
D_GW = 1e9  # pc

ntoas = 5000
toas = 365.25 * np.linspace(0, 10, ntoas)  # days

cosmu, Fp, Fx = GWecc.AntennaPattern(RA_GW, DEC_GW, RA_P, DEC_P)
delay = -1000  # -D_P*(1-cosmu) / (1+z)

for idx, e0 in enumerate([0.1, 0.5, 0.8]):

    for jdx, term in enumerate(
        [GWecc.ResidualsTerms_Earth, GWecc.ResidualsTerms_Pulsar]
    ):

        ax = plt.subplot(321 + idx * 2 + jdx)

        res_px = GWecc.EccentricResiduals_px(
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
            GWecc.ResidualsMethod_Adb,
            term,
            toas,
        )
        res_p = np.asarray(res_px[0]) - np.mean(res_px[0])
        res_x = np.asarray(res_px[1]) - np.mean(res_px[1])
        plt.plot(toas / 365.25, res_p, label="$s_+$", color="b")
        plt.plot(toas / 365.25, res_x, "--", label="$s_\\times$", color="g")
        plt.grid()
        plt.xlim([0, 10])

        if term == GWecc.ResidualsTerms_Earth:
            plt.ylabel("$s_{+,\\times}(t)$  (ns)", fontsize=14)
            plt.text(
                1,
                0.8 * ax.yaxis.get_data_interval()[1],
                "e=" + str(e0),
                size=15,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round", facecolor="cyan", alpha=0.5),
            )

        if idx == 0:
            title = (
                "Earth term" if term == GWecc.ResidualsTerms_Earth else "Pulsar term"
            )
            plt.title(title)

        if idx == 2 and jdx == 1:
            plt.legend(loc="lower right", fontsize=12)

        if idx < 2:
            plt.tick_params(labelbottom=False, labelsize=12)
        else:
            plt.tick_params(labelsize=12)
            plt.xlabel("$t=t_E$ (years)", fontsize=14)

plt.show()
