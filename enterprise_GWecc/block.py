from tkinter import N
import numpy as np

from enterprise.signals.parameter import Uniform, Constant
from enterprise_GWecc import eccentric_cw_delay_Planck18
from enterprise_extensions.deterministic import CWSignal

def gwecc_block(tref, skyloc=None, log10_F=None, ecc=None, redshift=None, psrTerm=False, name='gwecc'):
    """enterprise_extensios-style convinience function for creating a Signal for GWecc"""
    if skyloc is None:
        cos_gwtheta = Uniform(-1, 1)(f"{name}_cos_gwtheta")
        gwphi = Uniform(0, 2*np.pi)(f"{name}_gwphi")
    else:
        cos_gwtheta = Constant(skyloc[0])(f"{name}_cos_gwtheta")
        gwphi = Constant(skyloc[1])(f"{name}_gwphi")

    if log10_F is None:
        log10_F = Uniform(-9.0, -7.0)(f'{name}_log10_F')
    else:
        log10_F = Constant(log10_F)(f'{name}_log10_F')

    if ecc is None:
        e0 = Uniform(0.0, 0.8)(f'{name}_e0')
    else:
        e0 = Constant(ecc)(f'{name}_e0')

    psi = Uniform(-np.pi/2, np.pi)(f'{name}_psi')
    cos_inc = Uniform(-1, 1)(f'{name}_cos_inc')

    log10_M = Uniform(6.0, 9.1)(f'{name}_log10_M')
    q = Uniform(0.01, 1)(f'{name}_q')

    gamma0 = Uniform(-np.pi/2, np.pi)(f'{name}_psi')
    l0 = Uniform(-np.pi, np.pi)(f'{name}_psi')
    
    tref = Constant(tref)(f'{name}_tref')
    
    if redshift is None:
        z = Uniform(0, 1)(f'{name}_q')
    else:
        z = Constant(redshift)(f'{name}_q')
    
    wf = eccentric_cw_delay_Planck18(cos_gwtheta=cos_gwtheta, gwphi=gwphi, 
                                     psi=psi, cos_inc=cos_inc,
                                     log10_M=log10_M, q=q,
                                     log10_F=log10_F, e0=e0, gamma0=gamma0, l0=l0, tref=tref,
                                     z=z,
                                     psrTerm=psrTerm)
    
    # I am setting ecc=False to avoid eccess parameters.
    ecw = CWSignal(wf, ecc=False, psrTerm=psrTerm, name=name)

    return ecw

