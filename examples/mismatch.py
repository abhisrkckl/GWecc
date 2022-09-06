from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

from enterprise_GWecc.GWecc import (
    ResidualsTerms_Earth, 
    EccentricResiduals,
    ResidualsMethod_Adb,
    ResidualsMethod_Num
)

def hms_to_rad(hh, mm, ss):
    sgn = np.sign(hh)
    return sgn * (sgn * hh + mm / 60 + ss / 3600) * np.pi / 12


def dms_to_rad(dd, mm, ss):
    sgn = np.sign(dd)
    return sgn * (sgn * dd + mm / 60 + ss / 3600) * np.pi / 180

def pulsar_fit_matrix(ts):
    """ts in days
    """
    yr = 365.25

    nts = len(ts)

    I = np.eye(nts)

    M = np.transpose([np.ones_like(ts), ts, ts**2, np.cos(2*np.pi*ts/yr), np.sin(2*np.pi*ts/yr)])

    return I - np.dot(M, np.linalg.solve(M.T @ M, M.T))

def mismatch_function(ts):
    P = pulsar_fit_matrix(ts)

    def mismatch(s1, s2):
        s1_s2 = s1 @ P @ s2
        s1_s1 = s1 @ P @ s1
        s2_s2 = s2 @ P @ s2

        return 1 - s1_s2 / np.sqrt(s1_s1*s2_s2)
    
    return mismatch


q = 1
Pb0 = 1.5  # years
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

T = 30*365.25
ntoas = int(T/10) # 10-day cadence
ts = np.linspace(0, T, ntoas)  # days

delay = 1000

term = ResidualsTerms_Earth

mismatch = mismatch_function(ts)

es = np.linspace(0.01, 0.85, 20)
Ms = 10**np.linspace(7, 9, 3)
Pb0s = [0.5, 1.5, 5]

for idx,Pb0 in enumerate(Pb0s):
    ax = plt.subplot(311+idx)
    for M in Ms:
        print(f"M = {M:e}", )
        mismatches = []
        for e0 in es:
            res_adb = EccentricResiduals(
                M,
                q,
                Omega,
                i,
                T,
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
                ResidualsMethod_Adb,
                ResidualsTerms_Earth,
                ts,
            )

            res_num = EccentricResiduals(
                M,
                q,
                Omega,
                i,
                T,
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
                ts,
            )

            res_adb = np.asarray(res_adb) - np.mean(res_adb)

            res_num = np.asarray(res_num) - np.mean(res_num)

            mismatches.append( mismatch(res_adb, res_num) )

        plt.plot(es, mismatches, label=f"M={M:.0e} MSun")
        plt.yscale('log')
        plt.legend(fontsize=14)
        plt.xlabel("e", fontsize=14)
        plt.ylabel("mismatch", fontsize=14)
        plt.tick_params(axis='both', labelsize=12)
    
    plt.text(
        0.3,
        0.6 * ax.yaxis.get_data_interval()[1],
        f"Pb0={Pb0} yr",
        size=14,
        ha="center",
        va="center",
        bbox=dict(boxstyle="round", facecolor="grey", alpha=0.4),
    )
    
plt.show()