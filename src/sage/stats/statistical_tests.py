'''Statstical tests designed to be used in doctests.

Many SAGE functions that produce random elements do so with an implicit understanding that the elements are coming from
a certain distribution. For discrete objects, uniform sampling is often expected.

This module is intended to hold such statistical tests.
'''
import warnings

def sigmas_from_uniform(list_of_counts):
    '''Approximation for the likelihood of the observations given a uniform distribution.

    Given $n$ possible outcomes and $k$ samples, the expected number of occurrences of each outcome is $k/n$.
    Let $x_i$ be the number of occurences for outcome $i$.

    For sufficiently large $k,n$:
    chi_square = $sum_i (x_i-k/n)**2/(k/n)$ is approximately distributed as a chi squared with $k-1$ degrees of freedom.
    sigmas = $(chi_square_stat-(k-1))/(2k-2)**.5$ is approximately distributed as Gaussian with mean 0 and variance 1.

    Negative sigmas indicates that the counts are *more* uniform than you would expect at random.
    Positive sigmas indicates that the counts are *less* uniform than you would expect at random.
    
    EXAMPLES::
    
        sage: mylist = [58, 46, 50, 63, 48, 59, 53, 58, 50, 65, 48, 55, 61, 59, 40, 49, 60, 48, 56, 49, 60, 52, 54, 57, 41, 55, 55, 50, 54, 46, 57, 54, 53, 63, 61, 45, 54, 48, 53, 48, 54, 53, 57, 43, 59, 75, 53, 43, 55, 44, 56, 47, 65, 63, 61, 56, 59, 47, 33, 54, 53, 54, 53, 54, 54, 65, 62, 65, 59, 50, 63, 60, 46, 36, 45, 48, 50, 53, 52, 54, 53, 59, 61, 51, 57, 65, 50, 60, 43, 65, 38, 49, 66, 56, 42, 61, 66, 49, 49, 63, 49, 44, 57, 42, 50, 61, 46, 37, 42, 54, 49, 47, 45, 60, 49, 50, 48, 56, 51, 62, 57, 54, 47, 47, 63, 61, 53, 43, 50, 50, 50, 51, 51, 56, 51, 47, 71, 51, 49, 57, 61, 66, 45, 60, 48, 56, 55, 46, 52, 66, 52, 56, 42, 47, 59, 56, 47, 57, 52, 56, 49, 61, 64, 61, 68, 57, 51, 61]
        sage: abs(sigmas_from_uniform(mylist)) < 2 # within two standard deviations of the mean
        True
        sage: mylist = [100] * 50
        sage: sigmas_from_uniform(mylist) < -2 # this is *more* uniform than expected
        True
        sage: mylist = [100] * 10
        sage: sigmas_from_uniform(mylist) < -2 # this is *more* uniform than expected 
        doctest:warning
            ...
        UserWarning: Asymptotic approximations are suspect with a small number of observations.
        True
        
    TESTS::
    
        sage: abs(sigmas_from_uniform([]))
        Traceback (most recent call last):
            ...
        ZeroDivisionError: list_of_counts must contain at least two nonzero elements. 
        sage: abs(sigmas_from_uniform([100]))
        Traceback (most recent call last):
            ...
        ZeroDivisionError: list_of_counts must contain at least two nonzero elements. 
        sage: abs(sigmas_from_uniform([0,0,0]))
        Traceback (most recent call last):
            ...
        ZeroDivisionError: list_of_counts must contain at least two nonzero elements. 
    
    AUTHORS:

    '''
    num_samples = sum(list_of_counts)
    num_counts = len(list_of_counts)
    if num_counts < 50 or num_samples < 20:
        warnings.warn('Asymptotic approximations are suspect with a small number of observations.')
    try:
        expected = num_samples / num_counts
        if expected < 5:
            warnings.warn('Chi squared test works best with at least 5 x [number of samples].')
        chi_square_stat = sum((expected - count)**2 / expected for count in list_of_counts)  # chisquare(k-1)
        sigmas = (chi_square_stat - (num_counts-1)) * (2*num_counts-2)**-.5  # N(0,1)
    except (ZeroDivisionError):
        raise ZeroDivisionError("list_of_counts must contain at least two nonzero elements.")
    return sigmas
