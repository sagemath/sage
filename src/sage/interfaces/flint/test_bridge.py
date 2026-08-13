"""
Test module for SageMath-FLINT bridge (one-way conversion)
"""
import unittest
import flint
import mpmath
from sage.all import RR, CC, matrix, pi, I, RealField, ComplexField
from sage.interfaces.flint.bridge import (
    sage_to_flint_arb, 
    sage_to_flint_acb, 
    sage_matrix_to_flint
)

class TestFlintSageBridge(unittest.TestCase):
    """
    Test suite for the SageMath-to-FLINT bridge functions.
    """
    
    def setUp(self):
        """Set up test cases."""
        # Real number examples at different precisions
        self.real_std = RR(1.234)
        self.real_high = RealField(200)(pi)
        self.real_exact = RR(1)
        
        # Complex number examples
        self.complex_std = CC(1.234 + 5.678*I)
        self.complex_high = ComplexField(200)(pi + I*pi)
        self.complex_exact = CC(1 + I)
        
        # Matrix examples
        self.matrix_real = matrix(RR, 2, 2, [1.1, 2.2, 3.3, 4.4])
        self.matrix_complex = matrix(CC, 2, 2, [1+I, 2-I, 3+2*I, 4])

    def test_real_conversion_basic(self):
        """Test basic conversion of SageMath RealNumber to FLINT arb."""
        arb_num = sage_to_flint_arb(self.real_std)
        self.assertIsInstance(arb_num, flint.arb)
        self.assertAlmostEqual(float(arb_num), float(self.real_std))

    def test_real_conversion_exact(self):
        """Test exact number conversion of SageMath RealNumber to FLINT arb."""
        arb_num = sage_to_flint_arb(self.real_exact)
        self.assertIsInstance(arb_num, flint.arb)
        self.assertEqual(float(arb_num), float(self.real_exact))

    def test_real_conversion_high_precision(self):
        """Test high precision conversion of SageMath RealNumber to FLINT arb."""
        arb_num = sage_to_flint_arb(self.real_high)
        self.assertEqual(str(arb_num).replace('[', '').replace(']', '')[:15], str(self.real_high)[:15])

    def test_complex_conversion_basic(self):
        """Test basic conversion of SageMath ComplexNumber to FLINT acb."""
        acb_num = sage_to_flint_acb(self.complex_std)
        self.assertIsInstance(acb_num, flint.acb)
        self.assertAlmostEqual(complex(acb_num), complex(self.complex_std))

    def test_complex_conversion_high_precision(self):
        """Test high precision conversion of SageMath ComplexNumber to FLINT acb."""
        acb_num = sage_to_flint_acb(self.complex_high)
        
        flint_real_str = str(acb_num.real).replace('[', '').replace(']', '')
        flint_imag_str = str(acb_num.imag).replace('[', '').replace(']', '')
        
        self.assertEqual(flint_real_str[:15], str(self.complex_high.real())[:15])
        self.assertEqual(flint_imag_str[:15], str(self.complex_high.imag())[:15])

    def test_complex_conversion_exact(self):
        """Test exact number conversion of SageMath ComplexNumber to FLINT acb."""
        acb_num = sage_to_flint_acb(self.complex_exact)
        self.assertIsInstance(acb_num, flint.acb)
        self.assertEqual(complex(acb_num), complex(self.complex_exact))

    def test_matrix_conversion_real(self):
        """Test conversion of SageMath real matrix to FLINT arb_mat."""
        arb_mat = sage_matrix_to_flint(self.matrix_real)
        self.assertIsInstance(arb_mat, flint.arb_mat)
        self.assertEqual(arb_mat.nrows(), self.matrix_real.nrows())
        self.assertEqual(arb_mat.ncols(), self.matrix_real.ncols())
        # Check elements
        for i in range(self.matrix_real.nrows()):
            for j in range(self.matrix_real.ncols()):
                self.assertAlmostEqual(float(arb_mat[i,j]), float(self.matrix_real[i,j]))

    def test_matrix_conversion_complex(self):
        """Test conversion of SageMath complex matrix to FLINT acb_mat."""
        acb_mat = sage_matrix_to_flint(self.matrix_complex)
        self.assertIsInstance(acb_mat, flint.acb_mat)
        self.assertEqual(acb_mat.nrows(), self.matrix_complex.nrows())
        self.assertEqual(acb_mat.ncols(), self.matrix_complex.ncols())
        # Check elements
        for i in range(self.matrix_complex.nrows()):
            for j in range(self.matrix_complex.ncols()):
                self.assertAlmostEqual(
                    complex(acb_mat[i,j]), 
                    complex(self.matrix_complex[i,j])
                )

    def test_error_handling(self):
        """Test error handling for invalid input types."""
        # Test with invalid input for sage_to_flint_arb
        with self.assertRaises(Exception):
            sage_to_flint_arb("not a number")
        
        # Test with invalid input for sage_to_flint_acb
        with self.assertRaises(Exception):
            sage_to_flint_acb("not a number")

if __name__ == "__main__":
    unittest.main()