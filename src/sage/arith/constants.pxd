# sage_setup: distribution = sagemath-categories
#
# Various mathematical constants, represented as double precision hex
# float. We use these instead of decimal constants like 3.1415...
# because the hex floats are exactly representable as "double", so there
# shouldn't be any rounding issues in the compiler.
# See https://github.com/sagemath/sage/issues/23919#comment:15
#
# Hex floats are standardized in C99, but GCC accepts them
# unconditionally, also for C++ code.
#
# Note: if we ever need a lot of constants like this, we should
# consider auto-generating this list.
#

cdef extern from *:
    double M_PI "0x3.243f6a8885a3p+0"         # π
    double M_EULER "0x9.3c467e37db0c8p-4"     # γ
    double M_LN2 "0xb.17217f7d1cf78p-4"       # log(2)
    double M_LN10 "0x2.4d763776aaa2cp+0"      # log(10)
    double M_LNPI "0x1.250d048e7a1bdp+0"      # log(π)
    double M_1_LN2 "0x1.71547652b82fep+0"     # 1/log(2)
    double M_1_LN10 "0x6.f2dec549b9438p-4"    # 1/log(10)
    double M_1_LNPI "0xd.fa22fdd8cfd98p-4"    # 1/log(π)
    double M_LN2_LN10 "0x4.d104d427de7fcp-4"  # log(2)/log(10)

    double LOG_TEN_TWO_PLUS_EPSILON "0x3.5269e12f346e4p+0"  # log(10,2) rounded up
