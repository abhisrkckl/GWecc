import numpy as np
from enterprise_GWecc.GWecc import (
    ResidualsMethod_Num,
    ResidualsMethod_Adb,
    ResidualsMethod_PC,
    ResidualsMethod_PM,
    ResidualsTerms_Both,
    ResidualsTerms_Earth,
    ResidualsTerms_Pulsar,
)

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
Pb0 = 5  # years
e0 = 0.3
psi = 0
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

delay = -1000  # -D_P*(1-cosmu) / (1+z)

terms = [ResidualsTerms_Both, ResidualsTerms_Earth, ResidualsTerms_Pulsar]
methods = [
    ResidualsMethod_Num,
    ResidualsMethod_Adb,
    ResidualsMethod_PC,
    ResidualsMethod_PM,
]
