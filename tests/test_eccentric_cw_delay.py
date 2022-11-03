import pytest
import numpy as np
from enterprise.pulsar import Pulsar
from enterprise_GWecc.GWecc import EccentricResiduals
from enterprise_GWecc import eccentric_cw_delay
from default_test_params import (
    M,
    q,
    Omega,
    i,
    Pb0,
    e0,
    l0,
    gamma0,
    D_GW,
    RA_GW,
    DEC_GW,
    RA_P,
    DEC_P,
    z,
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


@pytest.mark.parametrize("method", methods)
def test_eccentric_cw_delay(testpsr, method):
    toas = testpsr.toas / (24 * 3600)  # MJD
    tref = max(toas)  # MJD
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
