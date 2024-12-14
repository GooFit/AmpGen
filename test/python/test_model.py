from ampgen import make_ampgen_model 
import pytest
import numpy as np 
import json

def test_amplitudes(model): 
    (m, refFileName, k, evt) = model
    ref = json.load( open(refFileName, 'r') )
    refv = ref[ f"{k}_{evt}"]
    amplitudes = m.amplitudes( m.phsp_event(which=evt))
    delta = np.subtract(amplitudes[k], refv)
    assert np.isclose(amplitudes[k], refv, equal_nan=True, rtol=5e-6, atol=5e-6).all(), f"{k:<80}/{evt} doesn't match reference,  Δ: {delta} rel: {delta/refv}" 

