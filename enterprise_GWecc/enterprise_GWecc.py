from .GWecc import EccentricResiduals, ResidualsMethod_Num, ResidualsTerms_Both, ResidualsTerms_Earth
#from enterprise.signals import signal_base
import numpy as np

year_to_s   = 365.25*24*3600

#@signal_base.function
def eccentric_cw_delay(toas,
                       RA_psr, DEC_psr, D_psr,
                       RA_GW, DEC_GW, log10_D_GW,
                       psi, i,
                       log10_M, q,
                       log10_f0_GW, e0, gamma0, l0, t0,
                       z,
                       psrTerm=False,
                       evolve=True):
    """
    toas        are Pulsar TOAs                 in MJD
    
    RA_PSR      is  RA of pulsar                in rad
    DEC_PSR     is  DEC of pulsar               in rad
    D_PSR       is  distance to pulsar          in parsec
    
    RA_GW       is  RA of GW source             in rad
    DEC_GW      is  DEC of GW source            in rad
    log10_D_GW  is  log10 distance to GW source in parsec
    
    psi         is  Polarization angle          in rad
    i           is  Inclination of GW source    in rad
    
    log10_M     is  log10 Total Mass of GW src  in SolarMass
    q           is  Mass ratio of GW source
    
    log10_f0_GW is  log10 GW frequency at t=t0  in Hz
    e0          is  Eccentricity at t=t0
    gamma0      is  Periastron angle at t=t0    in rad
    l0          is  Mean anomaly at t=t0        in rad
    t0          is  Reference time              in MJD
    
    z           is  redshift
    
    psrTerm     is  [boolean] Whether to add pulsar term 
       
    evolve      is  [boolean] Whether to evolve phase exactly   ---   Not implemented (Always True)
    
    
    Returns:    
            TOA delays due to GWs from eccentric binary sources (in s)
    """
        
    D_GW = 10.**log10_D_GW   
    M = 10.**log10_M
    
    n0 = np.pi*(10.**log10_f0_GW)    # GW frequency is twice the orbital frequency.
    Pb0 = 2*np.pi/n0 / year_to_s

    residuals_method = ResidualsMethod_Num
    residuals_terms = ResidualsTerms_Both if psrTerm else ResidualsTerms_Earth
    
    return np.asarray(
               EccentricResiduals(M, q,
                                  psi, i,
                                  t0, Pb0, e0, l0, gamma0,
                                  D_GW, RA_GW, DEC_GW,
                                  D_psr,  RA_psr,  DEC_psr,
                                  z,
                                  residuals_method,
                                  residuals_terms,
                                  toas)    
            )
