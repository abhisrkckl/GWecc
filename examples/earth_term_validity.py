import numpy as np
import matplotlib.pyplot as plt

import astropy.units as au
import astropy.constants as ac


def Nfunc1(e0, eta):
    return (
        (
            20368
            - 14784 * eta
            + e0**2 * (219880 - 159600 * eta)
            + e0**4 * (197022 - 141708 * eta)
            + e0**6 * (11717 - 8288 * eta)
        )
        * eta
    ) / (280.0 * (1 - e0**2) ** 4.5)


def fmax_3(M, e0, eta, T):
    N1 = Nfunc1(e0, eta)
    M = M * ac.GM_sun / ac.c**3
    T = au.Quantity(T, "year")
    return (
        (((np.pi / 25 / T**2 / N1 / M ** (7 / 3))) ** (3 / 13) / np.pi).to("Hz").value
    )


T = 30
eta = 0.25

"""
plt.subplot(121)
e0s = np.linspace(0, 0.9, 5)
Ms = 10 ** np.linspace(7, 9, 3)
for e0 in e0s:
    fmaxs = fmax_3(Ms, e0, eta, T) * 1e9
    plt.loglog(Ms, fmaxs, label=f"e0={e0:0.1f}")
plt.legend(framealpha=0.4)
plt.xlabel("M (Msun)")
plt.ylabel("Fmax (nHz)")
plt.title(f"eta={eta}")


plt.subplot(122)
e0s = np.linspace(0, 0.9, 20)
for M in Ms:
    fmaxs = fmax_3(M, e0s, eta, T) * 1e9
    plt.plot(e0s, fmaxs, label=f"M={M:0.0e} Msun")
plt.yscale("log")
plt.legend(framealpha=0.4)
plt.xlabel("e0")

plt.show()"""

e0s = np.linspace(0, 0.9, 100)
Ms = 10 ** np.linspace(7, 9, 40)
e0s, Ms = np.meshgrid(e0s, Ms)
Fmaxs = fmax_3(Ms, e0s, eta, T)
cf = plt.contourf(e0s, np.log10(Ms), np.log10(Fmaxs)+9)
plt.colorbar(cf)
plt.xlabel("$e_0$", fontsize=13)
plt.ylabel("$log_{10}[M]$ $(M_{sun})$", fontsize=13)
plt.title("$log_{10}[F_{\\max}]$ (nHz)")
plt.show()