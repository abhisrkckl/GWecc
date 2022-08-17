import numpy as np
from enterprise_GWecc.enterprise_GWecc_cosmoz import (
    eccentric_cw_delay_Planck18,
    eccentric_cw_waveform_Planck18,
)
from enterprise.signals import signal_base


def spline_coeffs(t12, s12, h12):
    t1, t2 = t12
    s1, s2 = s12
    h1, h2 = h12

    T = np.array(
        [
            [t1**3, t1**2, t1, 1],
            [3 * t1**2, 2 * t1, 1, 0],
            [t2**3, t2**2, t2, 1],
            [3 * t2**2, 2 * t2, 1, 0],
        ]
    )
    S = np.array([s1, h1, s2, h2])
    return np.linalg.solve(T, S)


class PTASignalSpline:
    def __init__(self, ptasignal_func, waveform_func, toas, Nspline):
        self.toas = toas
        self.tmax = np.max(toas)
        self.tmin = np.min(toas)
        self.tspan = self.tmax - self.tmin

        self.ptasignal_func = ptasignal_func
        self.waveform_func = waveform_func

        self.Nspline = Nspline
        self.ts = np.linspace(self.tmin, self.tmax, Nspline + 1)
        self.ss = self.ptasignal_func(self.ts)
        self.hs = self.waveform_func(self.ts)

        self.coeffss = [
            spline_coeffs(
                (self.ts[i], self.ts[i + 1]),
                (self.ss[i], self.ss[i + 1]),
                (self.hs[i], self.hs[i + 1]),
            )
            for i in range(Nspline)
        ]

        self.eval = np.vectorize(self.evaluate)

    def evaluate(self, t):
        i = min(int((t - self.tmin) / self.tspan * self.Nspline), self.Nspline - 1)

        A = self.coeffss[i]
        T = np.array([t**3, t**2, t, 1])
        return np.dot(A, T)


@signal_base.function
def eccentric_cw_delay_Planck18_spline(
    toas,
    theta,
    phi,
    pdist,
    cos_gwtheta,
    gwphi,
    psi,
    cos_inc,
    log10_M,
    q,
    log10_F,
    e0,
    gamma0,
    l0,
    tref,
    log10_zc,
    zp,
    p_dist=1.0,
    psrTerm=False,
    evolve=True,
    xspline=1
):
    gwecc_s_func = lambda ts: eccentric_cw_delay_Planck18(
        ts,
        theta,
        phi,
        pdist,
        cos_gwtheta,
        gwphi,
        psi,
        cos_inc,
        log10_M,
        q,
        log10_F,
        e0,
        gamma0,
        l0,
        tref,
        log10_zc,
        zp,
        p_dist=p_dist,
        psrTerm=psrTerm,
        evolve=evolve,
    )
    gwecc_h_func = lambda ts: eccentric_cw_waveform_Planck18(
        ts,
        theta,
        phi,
        pdist,
        cos_gwtheta,
        gwphi,
        psi,
        cos_inc,
        log10_M,
        q,
        log10_F,
        e0,
        gamma0,
        l0,
        tref,
        log10_zc,
        zp,
        p_dist=p_dist,
        psrTerm=psrTerm,
    )

    nharms = 18.64801851 * (1 - e0**2) ** -1.5 - 14.04695398
    Nspline = int(nharms * xspline * (np.max(toas) - np.min(toas)) * (10**log10_F) + 1)

    if Nspline < len(toas)/2:
        gwecc_spline = PTASignalSpline(gwecc_s_func, gwecc_h_func, toas, Nspline=Nspline)
        return gwecc_spline.eval(toas)
    else:
        return gwecc_s_func(toas)
