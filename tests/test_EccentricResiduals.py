import pytest
import numpy as np
from itertools import product as outer_product
from enterprise_GWecc.GWecc import (
    EccentricResiduals,
    ResidualsMethod_Num,
    ResidualsMethod_Adb,
    ResidualsMethod_Anl,
    ResidualsMethod_PM,
    ResidualsTerms_Both,
    ResidualsTerms_Earth,
    ResidualsTerms_Pulsar
)
from default_test_params import *

methods = [ResidualsMethod_Num, ResidualsMethod_Adb, ResidualsMethod_Anl, ResidualsMethod_PM]
terms = [ResidualsTerms_Both, ResidualsTerms_Earth, ResidualsTerms_Pulsar]

@pytest.mark.parametrize("method, term", outer_product(methods, terms))
def test_EccentricResiduals(method, term):
    res = EccentricResiduals(
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

    assert np.all(np.isfinite(res)) and len(res) == ntoas