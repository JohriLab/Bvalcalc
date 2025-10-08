import numpy as np
import pytest
from pathlib import Path
from Bvalcalc import get_params, calculateB_linear, calculateB_unlinked, calculateB_hri


class TestCalculateBLinear:
    """Test suite for calculateB_linear function."""
    
    def setup_method(self):
        """Set up test parameters."""
        self.params = get_params(str(Path(__file__).parent / "testparams" / "nogcBasicParams.py"))
    
    def test_B_value_for_zero_length(self):
        """Test that B = 1.0 when length_of_element is zero."""
        result = calculateB_linear(distance_to_element=100, length_of_element=0, params=self.params)
        assert np.allclose(result, 1.0), "B should be 1.0 when length_of_element is zero"
    
    def test_B_value_for_zero_distance(self):
        """Test B calculation when distance_to_element is zero."""
        result = calculateB_linear(distance_to_element=0, length_of_element=1000, params=self.params)
        assert result < 1.0, "B should be less than 1.0 when element is adjacent"
        assert result > 0.0, "B should be greater than 0.0"
    
    def test_B_value_scalar_inputs(self):
        """Test B calculation with scalar inputs."""
        result = calculateB_linear(distance_to_element=500, length_of_element=10000, params=self.params)
        assert isinstance(result, (float, np.ndarray)), "Result should be numeric"
        assert result < 1.0, "B should be less than 1.0 for linked selection"
        assert result > 0.0, "B should be greater than 0.0"
    
    def test_B_value_array_inputs(self):
        """Test B calculation with numpy array inputs."""
        distances = np.array([100, 500, 1000])
        lengths = np.array([1000, 5000, 10000])
        result = calculateB_linear(distance_to_element=distances, length_of_element=lengths, params=self.params)
        assert isinstance(result, np.ndarray), "Result should be numpy array"
        assert result.shape == distances.shape, "Result shape should match input shape"
        assert np.all(result < 1.0), "All B values should be less than 1.0"
        assert np.all(result > 0.0), "All B values should be greater than 0.0"
    
    def test_B_value_with_gene_conversion(self):
        """Test B calculation with gene conversion enabled."""
        gc_params = get_params(str(Path(__file__).parent / "testparams" / "gcBasicParams.py"))
        result = calculateB_linear(distance_to_element=500, length_of_element=10000, params=gc_params)
        assert isinstance(result, (float, np.ndarray)), "Result should be numeric"
        assert result < 1.0, "B should be less than 1.0 for linked selection"
        assert result > 0.0, "B should be greater than 0.0"
    
    def test_B_value_monotonicity(self):
        """Test that B decreases as distance decreases (stronger selection)."""
        distances = np.array([100, 500, 1000])
        length = 10000
        results = calculateB_linear(distance_to_element=distances, length_of_element=length, params=self.params)
        # B should decrease as distance decreases (stronger selection)
        assert results[0] <= results[1] <= results[2], "B should decrease with decreasing distance"


class TestCalculateBUnlinked:
    """Test suite for calculateB_unlinked function."""
    
    def setup_method(self):
        """Set up test parameters."""
        self.params = get_params(str(Path(__file__).parent / "testparams" / "nogcBasicParams.py"))
    
    def test_B_value_for_zero_unlinked_L(self):
        """Test that B = 1.0 when unlinked_L is zero."""
        result = calculateB_unlinked(unlinked_L=0, params=self.params)
        assert np.allclose(result, 1.0), "B should be 1.0 when unlinked_L is zero"
    
    def test_B_value_scalar_input(self):
        """Test B calculation with scalar input."""
        result = calculateB_unlinked(unlinked_L=200000, params=self.params)
        assert isinstance(result, (float, np.ndarray)), "Result should be numeric"
        assert result < 1.0, "B should be less than 1.0 for unlinked selection"
        assert result > 0.0, "B should be greater than 0.0"
    
    def test_B_value_array_input(self):
        """Test B calculation with numpy array input."""
        unlinked_Ls = np.array([100000, 200000, 500000])
        result = calculateB_unlinked(unlinked_L=unlinked_Ls, params=self.params)
        assert isinstance(result, np.ndarray), "Result should be numpy array"
        assert result.shape == unlinked_Ls.shape, "Result shape should match input shape"
        assert np.all(result < 1.0), "All B values should be less than 1.0"
        assert np.all(result > 0.0), "All B values should be greater than 0.0"
    
    def test_B_value_monotonicity(self):
        """Test that B decreases as unlinked_L increases (more selection)."""
        unlinked_Ls = np.array([100000, 200000, 500000])
        results = calculateB_unlinked(unlinked_L=unlinked_Ls, params=self.params)
        # B should decrease as unlinked_L increases (more selection)
        assert results[0] >= results[1] >= results[2], "B should decrease with increasing unlinked_L"


class TestCalculateBHri:
    """Test suite for calculateB_hri function."""
    
    def setup_method(self):
        """Set up test parameters."""
        self.params = get_params(str(Path(__file__).parent / "testparams" / "nogcBasicParams.py"))
    
    def test_B_value_for_zero_interfering_L(self):
        """Test that B' = distant_B when interfering_L is zero."""
        distant_B = 0.9
        result = calculateB_hri(distant_B=distant_B, interfering_L=0, params=self.params)
        assert np.allclose(result, distant_B), "B' should equal distant_B when interfering_L is zero"
    
    def test_B_value_scalar_inputs(self):
        """Test B' calculation with scalar inputs."""
        distant_B = 0.95
        interfering_L = 50000
        result = calculateB_hri(distant_B=distant_B, interfering_L=interfering_L, params=self.params)
        assert isinstance(result, (float, np.ndarray)), "Result should be numeric"
        assert result < distant_B, "B' should be less than distant_B due to interference"
        assert result > 0.0, "B' should be greater than 0.0"
    
    def test_B_value_array_inputs(self):
        """Test B' calculation with numpy array inputs."""
        distant_Bs = np.array([0.9, 0.95, 0.99])
        interfering_Ls = np.array([10000, 50000, 100000])
        result = calculateB_hri(distant_B=distant_Bs, interfering_L=interfering_Ls, params=self.params)
        assert isinstance(result, np.ndarray), "Result should be numpy array"
        assert result.shape == distant_Bs.shape, "Result shape should match input shape"
        assert np.all(result < distant_Bs), "All B' values should be less than distant_B values"
        assert np.all(result > 0.0), "All B' values should be greater than 0.0"
    
    def test_B_value_monotonicity_with_interfering_L(self):
        """Test that B' decreases as interfering_L increases."""
        distant_B = 0.95
        # Test with individual calls since mixing scalar and array returns scalar
        result1 = calculateB_hri(distant_B=distant_B, interfering_L=10000, params=self.params)
        result2 = calculateB_hri(distant_B=distant_B, interfering_L=50000, params=self.params)
        result3 = calculateB_hri(distant_B=distant_B, interfering_L=100000, params=self.params)
        # B' should decrease as interfering_L increases (more interference)
        assert result1 >= result2 >= result3, "B' should decrease with increasing interfering_L"
    
    def test_B_value_monotonicity_with_distant_B(self):
        """Test that B' behavior with different distant_B values."""
        # Test with individual calls since the relationship is complex
        result1 = calculateB_hri(distant_B=0.8, interfering_L=50000, params=self.params)
        result2 = calculateB_hri(distant_B=0.9, interfering_L=50000, params=self.params)
        result3 = calculateB_hri(distant_B=0.95, interfering_L=50000, params=self.params)
        # All results should be valid B values
        assert all(0 < r < 1 for r in [result1, result2, result3]), "All B' values should be between 0 and 1"


class TestApiIntegration:
    """Test suite for API integration and parameter handling."""
    
    def test_get_params_with_file(self):
        """Test get_params with a specific parameter file."""
        params = get_params(str(Path(__file__).parent / "testparams" / "nogcBasicParams.py"))
        assert isinstance(params, dict), "get_params should return a dictionary"
        required_keys = ["r", "u", "g", "k", "t1", "t1half", "t2", "t3", "t4", "f1", "f2", "f3", "f0", "t_constant"]
        for key in required_keys:
            assert key in params, f"Parameter {key} should be present in params"
    
    def test_functions_without_params(self):
        """Test that functions work when params is None (uses environment variable)."""
        # This test relies on the conftest.py setting BCALC_params
        result1 = calculateB_unlinked(unlinked_L=100000)
        assert isinstance(result1, (float, np.ndarray)), "Should work without explicit params"
        
        result2 = calculateB_linear(distance_to_element=500, length_of_element=10000)
        assert isinstance(result2, (float, np.ndarray)), "Should work without explicit params"
        
        result3 = calculateB_hri(distant_B=0.95, interfering_L=50000)
        assert isinstance(result3, (float, np.ndarray)), "Should work without explicit params"
    
    def test_consistent_results_with_same_params(self):
        """Test that functions give consistent results with the same parameters."""
        params = get_params(str(Path(__file__).parent / "testparams" / "nogcBasicParams.py"))
        
        # Test calculateB_unlinked
        result1a = calculateB_unlinked(unlinked_L=200000, params=params)
        result1b = calculateB_unlinked(unlinked_L=200000, params=params)
        assert np.allclose(result1a, result1b), "calculateB_unlinked should give consistent results"
        
        # Test calculateB_linear
        result2a = calculateB_linear(distance_to_element=500, length_of_element=10000, params=params)
        result2b = calculateB_linear(distance_to_element=500, length_of_element=10000, params=params)
        assert np.allclose(result2a, result2b), "calculateB_linear should give consistent results"
        
        # Test calculateB_hri
        result3a = calculateB_hri(distant_B=0.95, interfering_L=50000, params=params)
        result3b = calculateB_hri(distant_B=0.95, interfering_L=50000, params=params)
        assert np.allclose(result3a, result3b), "calculateB_hri should give consistent results"


class TestExtendedOutput:
    """Test suite for extended output functionality from manuscript Table 1 and SA.1/SA.2."""
    
    def test_unlinked_constant_dfe(self):
        """Test calculateB_unlinked with constant DFE (Table 1 scenario)."""
        # Load UnlinkedParams with constant DFE
        params = get_params(str(Path(__file__).parent / "testparams" / "UnlinkedParams.py"), constant_dfe=True)
        
        # Test with the scenario from Table 1
        result = calculateB_unlinked(unlinked_L=200000, params=params)
        
        # Verify the exact expected value from manuscript
        expected = 0.9166085033148563
        assert np.allclose(result, expected, rtol=1e-10), f"Expected {expected}, got {result}"
    
    def test_unlinked_variable_dfe(self):
        """Test calculateB_unlinked with variable DFE (Table 1 scenario)."""
        # Load UnlinkedParams with variable DFE
        params = get_params(str(Path(__file__).parent / "testparams" / "UnlinkedParams.py"), constant_dfe=False)
        
        # Test with the scenario from Table 1
        result = calculateB_unlinked(unlinked_L=200000, params=params)
        
        # Verify the exact expected value from manuscript
        expected = 0.7570197890708311
        assert np.allclose(result, expected, rtol=1e-10), f"Expected {expected}, got {result}"
    
    def test_hri_constant_dfe(self):
        """Test calculateB_hri with constant DFE (Table SA.1 scenario)."""
        # Load HriParams with constant DFE
        params = get_params(str(Path(__file__).parent / "testparams" / "HriParams.py"), constant_dfe=True)
        
        # Test with the scenario from Table SA.1
        result_linear = calculateB_linear(distance_to_element=0, length_of_element=10000, params=params)
        result_hri = calculateB_hri(distant_B=1.0, interfering_L=10000, params=params)
        
        # Verify the exact expected values from manuscript
        expected_linear = 0.4723665527410146
        expected_hri = 0.4793601357930167
        
        assert np.allclose(result_linear, expected_linear, rtol=1e-10), f"Expected linear {expected_linear}, got {result_linear}"
        assert np.allclose(result_hri, expected_hri, rtol=1e-10), f"Expected hri {expected_hri}, got {result_hri}"
    
    def test_hri_variable_dfe(self):
        """Test calculateB_hri with variable DFE (Table SA.2 scenario)."""
        # Load HriParams with variable DFE
        params = get_params(str(Path(__file__).parent / "testparams" / "HriParams.py"), constant_dfe=False)
        
        # Test with the scenario from Table SA.2
        result_linear = calculateB_linear(distance_to_element=0, length_of_element=10000, params=params)
        result_hri = calculateB_hri(distant_B=1.0, interfering_L=10000, params=params)
        
        # Verify the exact expected values from manuscript
        expected_linear = 0.06190177776976752
        expected_hri = 0.310727648509544
        
        assert np.allclose(result_linear, expected_linear, rtol=1e-10), f"Expected linear {expected_linear}, got {result_linear}"
        assert np.allclose(result_hri, expected_hri, rtol=1e-10), f"Expected hri {expected_hri}, got {result_hri}" 