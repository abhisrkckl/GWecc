from .GWecc import EccentricResiduals, ResidualsMethod_Num, ResidualsTerms_Both, ResidualsTerms_Earth
from enterprise.signals import signal_base
import numpy as np

year_to_s   = 365.25*24*3600

@signal_base.function
def eccentric_cw_delay(toas, 
                       theta, phi, pdist, 
                       cos_gwtheta, gwphi, log10_dist,
                       psi, cos_inc,
                       log10_M, q,
                       log10_f0_GW, e0, gamma0, l0, t0,
                       z,
                       p_dist=1.0,
                       psrTerm=False,
                       evolve=True):
    """
    ======================================================================================================================
    PARAM       DESCRIPTION                             UNIT        EXPRESSION              COMMENTS
    ======================================================================================================================
    toas        are Pulsar TOAs                         in s                                (Taken from the Pulsar object)
    
    theta       is  Polar angle of pulsar               in rad      = pi/2 - DEC_psr        (Taken from the Pulsar object)      
    phi         is  Azimuthal angle of pulsar           in rad      = RA_psr                (Taken from the Pulsar object)
    pdist       is  distance to pulsar (value,error)    in kpc                              (Taken from the Pulsar object)
    
    p_dist      is  Pulsar distance prior               in kpc                              (Parametrizes the departure from the measured pulsar distance given in pdist.)
    
    cos_gwtheta is  Cos of Polar angle of GW source     in rad      = cos(pi/2 - DEC_GW)     
    gwphi       is  Azimithal angle of GW source        in rad      = RA_GW
    log10_dist  is  log10 distance to GW source         in Mpc      = log10(D_GW)
    
    psi         is  Polarization angle                  in rad      
    cos_inc     is  Cos of Inclination of GW source     in rad      = cos(i)
    
    log10_M     is  log10 Total Mass of GW src          in MSun     = log10(M)              (Not to be confused with the chirp mass.)
    q           is  Mass ratio of GW source                         = m1/m2
    
    log10_f0_GW is  log10 GW frequency at t=t0          in Hz       = log10(2/Pb0)
    e0          is  Eccentricity at t=t0
    gamma0      is  Periastron angle at t=t0            in rad
    l0          is  Mean anomaly at t=t0                in rad                              (Should be set to zero if t0 is a free parameter.)
    t0          is  Reference time                      in s                                (Should be set to a fixed epoch if l0 is a free parameter.)
    
    z           is  redshift
    
    psrTerm     is  [boolean] Whether to add pulsar term 
       
    evolve      is  [boolean] Whether to evolve phase exactly                               Not implemented (Always True)
    
    ======================================================================================================================
    
    Returns:    
            TOA delays due to GWs from eccentric binary sources (in s)
    """
    
    toas = toas / (24*3600)
    t0   = t0 / (24*3600)
    
    # Check the precision of these parameters later. 
    RA_psr  = phi
    DEC_psr = np.pi/2 - theta
    
    #converts from kpc to pc
    if type(pdist) is not tuple:
        D_psr = pdist*1000 
    else:
        D_psr = (pdist[0]+pdist[1]*p_dist)*1000 
    
    gwtheta = np.arccos(cos_gwtheta)
    DEC_GW = np.pi/2 - gwtheta
    RA_GW  = gwphi
    
    D_GW = 1e6 * (10.**log10_dist)   
    M = 10.**log10_M
    
    i = np.arccos(cos_inc)
    
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
