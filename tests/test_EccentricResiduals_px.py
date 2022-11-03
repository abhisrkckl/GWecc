import pytest
import numpy as np
from itertools import product as outer_product
from enterprise_GWecc.GWecc import (
    EccentricResiduals_px,
    ResidualsMethod_Num,
    ResidualsMethod_Adb,
    ResidualsMethod_Anl,
    ResidualsMethod_PM,
    ResidualsTerms_Both,
    ResidualsTerms_Earth,
    ResidualsTerms_Pulsar,
)
from default_test_params import *

methods = [
    ResidualsMethod_Num,
    ResidualsMethod_Adb,
    ResidualsMethod_Anl,
    ResidualsMethod_PM,
]
terms = [ResidualsTerms_Both, ResidualsTerms_Earth, ResidualsTerms_Pulsar]


@pytest.mark.parametrize("method, term", outer_product(methods, terms))
def test_EccentricResiduals_px(method, term):

    delay = -1000  # -D_P*(1-cosmu) / (1+z)

    res_p, res_x = EccentricResiduals_px(
        M,
        q,
        Omega,
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
