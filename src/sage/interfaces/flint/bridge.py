"""
flint_sage_bridge.py - Bridge between SageMath and Python-FLINT
Solves GitHub issue #39836
This module provides one-way conversion from SageMath numerical types
to Python-FLINT arbitrary precision types.
"""
import flint
import mpmath
from sage.all import RR, CC

def sage_to_flint_arb(sage_real):
    """
    Convert SageMath RealNumber to FLINT arb type.
    
    This preserves more precision than a simple float conversion by using
    mpmath as an intermediary.
    
    Parameters:
    -----------
    sage_real : sage.rings.real_mpfr.RealNumber
        A SageMath real number
    
    Returns:
    --------
    flint.arb
        The corresponding FLINT arb number
    
    Examples:
    ---------
    >>> from sage.all import RR, pi
    >>> sage_to_flint_arb(RR(pi))
    [3.1415926535897932384626433832795028842 +/- ...]
    """
    try:
        # Get the precision from the SageMath object if available
        prec = sage_real.prec() if hasattr(sage_real, 'prec') else 53
        
        # Get exact string representation from SageMath
        if hasattr(sage_real, 'str'):
            # Use the exact string representation if available
            str_val = sage_real.str(no_sci=True)
        else:
            # Fall back to regular string representation
            str_val = str(sage_real).replace('[', '').replace(']', '')
            
        # Convert using mpmath for precision preservation
        with mpmath.workdps(prec):
            mp_val = mpmath.mpf(str_val)
            result = flint.arb(mp_val, prec=prec)
            
        return result
    except (TypeError, ValueError):
        # Fallback method
        return flint.arb(float(sage_real))

def sage_to_flint_acb(sage_complex):
    """
    Convert SageMath ComplexNumber to FLINT acb type.
    
    This preserves more precision than a simple complex conversion by using
    mpmath as an intermediary.
    
    Parameters:
    -----------
    sage_complex : sage.rings.complex_mpfr.ComplexNumber
        A SageMath complex number
    
    Returns:
    --------
    flint.acb
        The corresponding FLINT acb complex number
    
    Examples:
    ---------
    >>> from sage.all import CC, pi, I
    >>> sage_to_flint_acb(CC(pi + I*pi))
    ([3.1415926535897932384626433832795028842 +/- ...] + [3.1415926535897932384626433832795028842 +/- ...]*I)
    """
    try:
        # Get the precision from the SageMath object if available
        prec = sage_complex.prec() if hasattr(sage_complex, 'prec') else 53
        
        # Get string representations of real and imaginary parts
        real_part = sage_complex.real()
        imag_part = sage_complex.imag()
        
        if hasattr(real_part, 'str'):
            real_str = real_part.str(no_sci=True)
            imag_str = imag_part.str(no_sci=True)
        else:
            real_str = str(real_part).replace('[', '').replace(']', '')
            imag_str = str(imag_part).replace('[', '').replace(']', '')
        
        # Convert using mpmath with appropriate precision
        with mpmath.workdps(prec):
            mp_val = mpmath.mpc(real_str, imag_str)
            result = flint.acb(mp_val, prec=prec)
            
        return result
    except (TypeError, ValueError):
        # Fallback method
        return flint.acb(complex(sage_complex))

def sage_matrix_to_flint(sage_matrix):
    """
    Convert a SageMath matrix with RR/CC entries to FLINT matrix.
    
    Parameters:
    -----------
    sage_matrix : sage.matrix
        A SageMath matrix with real or complex entries
    
    Returns:
    --------
    flint.arb_mat or flint.acb_mat
        A FLINT matrix with corresponding entries
    
    Examples:
    ---------
    >>> from sage.all import matrix, RR
    >>> M = matrix(RR, 2, 2, [1.1, 2.2, 3.3, 4.4])
    >>> sage_matrix_to_flint(M)
    [1.10000000000000 2.20000000000000]
    [3.30000000000000 4.40000000000000]
    """
    rows, cols = sage_matrix.nrows(), sage_matrix.ncols()
    
    # Detect if we're working with real or complex entries
    is_complex = False
    for i in range(rows):
        for j in range(cols):
            if sage_matrix[i,j].parent() == CC or (hasattr(sage_matrix[i,j], 'imag') and not sage_matrix[i,j].imag().is_zero()):
                is_complex = True
                break
        if is_complex:
            break
    
    # Create the appropriate matrix type
    if is_complex:
        result = flint.acb_mat(rows, cols)
        for i in range(rows):
            for j in range(cols):
                result[i,j] = sage_to_flint_acb(sage_matrix[i,j])
    else:
        result = flint.arb_mat(rows, cols)
        for i in range(rows):
            for j in range(cols):
                result[i,j] = sage_to_flint_arb(sage_matrix[i,j])
                
    return result