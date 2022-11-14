import pytest
import numpy as np
from itertools import product as outer_product
from enterprise_GWecc.GWecc import eccentric_residuals, eccentric_residuals_px
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
    methods,
    terms,
    toas,
    ntoas,
    delay,
)


@pytest.mark.parametrize("method, term", outer_product(methods, terms))
def test_eccentric_residuals(method, term):
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

    assert np.all(np.isfinite(res)) and not np.all(res == 0) and len(res) == ntoas


@pytest.mark.parametrize("method, term", outer_product(methods, terms))
def test_eccentric_residuals_px(method, term):

    res_p, res_x = eccentric_residuals_px(
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

    assert np.all(np.isfinite(res_p)) and np.all(np.isfinite(res_x))
    assert not np.all(res_p == 0) and not np.all(res_x == 0)
    assert len(res_p) == ntoas and len(res_x) == ntoas
