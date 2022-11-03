import pytest
import numpy as np
from astropy.cosmology import Planck18
from enterprise.pulsar import Pulsar
from enterprise_GWecc.GWecc import EccentricResiduals
from enterprise_GWecc import eccentric_cw_delay, eccentric_cw_delay_Planck18
from default_test_params import (
    M,
    q,
    Omega,
    i,
    Pb0,
    e0,
    l0,
    gamma0,
    RA_GW,
    DEC_GW,
    methods,
    year,
    day,
    ResidualsTerms_Earth,
)


@pytest.fixture
def testpsr():
    psr = Pulsar(
        "J0437-4715.IPTADR2/J0437-4715.fit.par",
        "J0437-4715.IPTADR2/J0437-4715.IPTADR2.tim",
    )
    return psr


z = 0.01
D_GW = 1e6 * Planck18.luminosity_distance(z).value

@pytest.mark.parametrize("method", methods)
def test_eccentric_cw_delay(testpsr, method):
    toas = testpsr.toas / (24 * 3600)  # MJD
    tref = max(toas)  # MJD
    
    RA_P = testpsr.phi
    DEC_P = np.pi/2 - testpsr.theta
    D_P = testpsr.pdist[0] * 1000  # pc

    res = EccentricResiduals(
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
        method,
        ResidualsTerms_Earth,
        toas,
    )

    ecc_gw = eccentric_cw_delay(
        cos_gwtheta=np.sin(DEC_GW),
        gwphi=RA_GW,
        log10_dist=np.log10(D_GW / 1e6),
        log10_h=None,
        psi=Omega,
        cos_inc=np.cos(i),
        log10_M=np.log10(M),
        q=q,
        log10_F=np.log10(2 / (Pb0 * year)),
        e0=e0,
        gamma0=gamma0,
        l0=l0,
        tref=tref * day,
        z=z,
        p_dist=0,
        psrTerm=False,
        evolve=True,
    )
    ecc_gw_fn = ecc_gw("ecc_gw", psr=testpsr)
    res2 = ecc_gw_fn()

    res = np.asarray(res) - np.mean(res)
    res2 = np.asarray(res2) - np.mean(res2)

    assert np.allclose(res, res2)


@pytest.mark.parametrize("method", methods)
def test_eccentric_cw_delay_Planck18(testpsr, method):
    toas = testpsr.toas / (24 * 3600)  # MJD
    tref = max(toas)  # MJD

    ecc_gw_1 = eccentric_cw_delay(
        cos_gwtheta=np.sin(DEC_GW),
        gwphi=RA_GW,
        log10_dist=np.log10(D_GW / 1e6),
        log10_h=None,
        psi=Omega,
        cos_inc=np.cos(i),
        log10_M=np.log10(M),
        q=q,
        log10_F=np.log10(2 / (Pb0 * year)),
        e0=e0,
        gamma0=gamma0,
        l0=l0,
        tref=tref * day,
        z=z,
        p_dist=0,
        psrTerm=False,
        evolve=True,
    )
    ecc_gw_1_fn = ecc_gw_1("ecc_gw_1", psr=testpsr)
    res_1 = ecc_gw_1_fn()

    ecc_gw_2 = eccentric_cw_delay_Planck18(
        cos_gwtheta=np.sin(DEC_GW),
        gwphi=RA_GW,
        psi=Omega,
        cos_inc=np.cos(i),
        log10_M=np.log10(M),
        q=q,
        log10_F=np.log10(2 / (Pb0 * year)),
        e0=e0,
        gamma0=gamma0,
        l0=l0,
        tref=tref * day,
        log10_zc=np.log10(z),
        zp=0,
        p_dist=0,
        psrTerm=False,
        evolve=True,
    )
    ecc_gw_2_fn = ecc_gw_2("ecc_gw_2", psr=testpsr)
    res_2 = ecc_gw_2_fn()

    res_1 -= np.mean(res_1)
    res_2 -= np.mean(res_2)

    assert np.allclose(res_1, res_2)
