import numpy as np
from helperScripts.calculateB import calculateB_linear

def test_B_value_for_zero_length():
    result = calculateB_linear(np.array([100]), np.array([0]))
    assert np.allclose(result, 1.0), "B should be 1.0 when length_of_element is zero"