"""
Tests for the fast univariate series expansion in Pynac
-------------------------------------------------------

This is a transcription of part of a testsuite that was
originally written for Symengine. Instead of matching
the first few terms of the expansion we pick specific
coefficients and compare to known values.

AUTHORS
-------
(c) 2016 Ralf Stephan <ralf@ark.in-berlin.de>

See https://github.com/symengine/symengine/blob/master/symengine/tests/basic/

TESTS::

    sage: def test(ex, coeff, value):
    ....:     assert(ex.series(x, coeff+3).coefficient(x,coeff) == value)

    sage: test(sin(x), 9, 1/362880)
    sage: test(sin(x)+cos(x), 8, 1/40320)
    sage: test(sin(x)*cos(x), 11, -4/155925)
    sage: test(sin(atan(x)), 27, -1300075/8388608)
    sage: test(cos(x/(1-x)), 11, -125929/362880)
    sage: test(1/(1-x), 99, 1)
    sage: test(x/(1-x-x^2), 35, 9227465)
    sage: test(x^3/(1-2*x^2), 49, 8388608)
    sage: test(1/(1-sin(x)), 10, 1382/14175)
    sage: test(1/x, -1, 1)
    sage: test(1/x/(1-x), 0, 1)
    sage: test(sqrt(4-x), 6, -21/2097152)
    sage: test((1-x)^(-2/3), 10, 1621477/4782969)
    sage: test(sqrt(1-x), 10, -2431/262144)
    sage: test(sqrt(cos(x)), 8, -559/645120)
    sage: test(cos(x)^(-2/3), 10, 701/127575)
    sage: test(log(1+x), 98, -1/98)
    sage: test(log(cos(x)), 12, -691/935550)
    sage: test(log(1/(1-x)), 48, 1/48)
    sage: test(exp(x), 9, 1/362880)
    sage: test(exp(log(1+x)), 0, 1)
    sage: test(exp(log(1+x)), 1, 1)
    sage: test(log(exp(x)), 1, 1)
    sage: test(exp(sin(x)), 10, -2951/3628800)
    sage: test(cos(x)^sin(x), 16, 1381/2661120)
    sage: test(tan(x), 13, 21844/6081075)
    sage: test(tan(x/(1-x)), 12, 1303712/14175)
    sage: test(cot(x), 13, -4/18243225)
    sage: test(cot(log(1+x)), 8, 4399/41472)
    sage: test(sec(x), 14, 199360981/87178291200)
    sage: test(sec(x*sqrt(1-x)), 11, -21463/103680)
    sage: test(csc(x), 13, 8191/37362124800)
    sage: test(csc(x/(1-x)), 14, 355857510913/100590336000)
    sage: test(asin(x), 15, 143/10240)
    sage: test(asin(x/(1-x)), 16, 1259743/2048)
    sage: test(atan(x), 19, -1/19)
    sage: test(atan(x/(1-x)), 33, 65536/33)
    sage: test(sinh(x), 9, 1/362880)
    sage: test(sinh(x/(1-x)), 10, 325249/40320)
    sage: test(cosh(x), 10, 1/3628800)
    sage: test(cosh(x/(1-x)), 11, 3756889/362880)
    sage: test(tanh(x), 13, 21844/6081075)
    sage: test(tanh(x/(1-x)), 14, 225979/66825)
    sage: test(coth(x), 13, 4/18243225)
    sage: test(coth(x/(1-x)), 14, -3651803/16372125)
    sage: test(sech(x), 16, 3878302429/4184557977600)
    sage: test(sech(x/(1-x)), 16, -34746888589811/4184557977600)
    sage: test(csch(x), 13, -8191/37362124800)
    sage: test(csch(x/(1-x)), 14, 2431542527/11176704000)
    sage: test(atanh(x), 99, 1/99)
    sage: test(atanh(x/(1-x)), 16, 2048)
    sage: test(asinh(x), 15, -143/10240)
    sage: test(asinh(x/(1-x)), 16, -3179/2048)
"""
