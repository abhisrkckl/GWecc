import pytest
import numpy as np
from enterprise_GWecc.GWecc import eccentric_waveform, eccentric_waveform_px
from default_test_params import (
    M,
    q,
    psi,
    i,
    t0,
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
    terms,
    toas,
    ntoas,
)


@pytest.mark.parametrize("term", terms)
def test_eccentric_waveform(term):
    h = eccentric_waveform(
        M,
        q,
        psi,
        i,
        t0,
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
        term,
        toas,
    )

    assert np.all(np.isfinite(h)) and not np.all(h == 0) and len(h) == ntoas


def test_eccentric_waveform_px():
    hp, hx = eccentric_waveform_px(
        M,
        q,
        psi,
        i,
        t0,
        Pb0,
        e0,
        l0,
        gamma0,
        D_GW,
        0,
        toas,
    )

    assert np.all(np.isfinite(hp)) and np.all(np.isfinite(hx))
    assert not np.all(hp == 0) and not np.all(hx == 0)
    assert len(hp) == ntoas and len(hx) == ntoas
