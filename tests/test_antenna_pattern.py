import pytest
import numpy as np
from itertools import product as outer_product
from enterprise_GWecc.GWecc import (
    antenna_pattern,
    eccentric_residuals,
    eccentric_residuals_px,
)
from default_test_params import (
    RA_GW,
    DEC_GW,
    RA_P,
    DEC_P,
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
    D_P,
    delay,
    z,
    methods,
    terms,
    toas,
)


@pytest.mark.parametrize("method, term", outer_product(methods, terms))
def test_antenna_pattern(method, term):
    cosmu, Fp, Fx = antenna_pattern(RA_GW, DEC_GW, RA_P, DEC_P)

    assert np.all(np.isfinite([cosmu, Fp, Fx]))

    res_px = eccentric_residuals_px(
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
        delay,
        z,
        method,
        term,
        toas,
    )
    res_from_px = np.array(res_px[0]) * Fp + np.array(res_px[1]) * Fx

    res = eccentric_residuals(
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
        method,
        term,
        toas,
    )

    assert np.allclose(res, res_from_px)
