.. -*- coding: utf-8 -*-

.. linkall

.. _prep-quickstart-numerical-analysis:

Sage Quickstart for Numerical Analysis
======================================

This `Sage <http://www.sagemath.org/>`_ quickstart tutorial was
developed for the MAA PREP Workshop "Sage: Using Open\-Source
Mathematics Software with Undergraduates" (funding provided by NSF DUE
0817071).  It is licensed under the Creative Commons
Attribution\-ShareAlike 3.0 license (`CC BY\-SA
<http://creativecommons.org/licenses/by-sa/3.0/>`_).

Sage includes many tools for numerical analysis investigations.

The place to begin is the default decimal number types in Sage.

Basic Analysis
--------------

- The ``RealField`` class using arbitrary precision (implemented with
  `MPFR <http://www.mpfr.org/>`_).

- The default real numbers (``RR``) is ``RealField(53)`` (i.e., 53 bits
  of precision).

- But you can make them live at whatever precision you wish.

::

    sage: ring=RealField(3)

To print the actual number (without rounding off the last few imprecise
digits to only display correct digits), call the ``.str()`` method::

    sage: print(ring('1').nextabove())
    1.2

::

    sage: print(ring('1').nextabove().str())
    1.2
    sage: print(ring('1').nextbelow().str())
    0.88

Let's change our precision.

::

    sage: ring=RealField(20)
    sage: print(ring('1').nextabove().str())
    1.0000019
    sage: print(ring('1').nextbelow().str())
    0.99999905

You can also specify the rounding mode.

::

    sage: ringup=RealField(3,rnd='RNDU')
    sage: ringdown=RealField(3,rnd='RNDD')

::

    sage: ring(1/9).str()
    '0.11111116'

::

    sage: ringup(1/9).str()
    '0.13'

::

    sage: ringdown(1/9).str()
    '0.10'

::

    sage: ring(1/9).str(base=2)
    '0.00011100011100011100100'

Let's see exactly what the fraction is.

::

    sage: ring(1/9).exact_rational()
    233017/2097152

::

    sage: n(233017/2097152)
    0.111111164093018

Converting to floating point binary (IEEE format)
-------------------------------------------------

::

    sage: x=13.28125

Here is ``x`` as a binary floating point number.

::

    sage: x.str(base=2)
    '1101.0100100000000000000000000000000000000000000000000'

From this binary floating point representation, we can construct ``x`` again.

::

    sage: xbin=2^3 +2^2 + 2^0+2^-2+2^-5; xbin
    425/32

::

    sage: xbin.n()
    13.2812500000000

To check, let's ask Sage what ``x`` is as an exact fraction (no rounding
involved).

::

    sage: x.exact_rational()
    425/32

Now let's convert ``x`` to the IEEE 754 binary format that is commonly
used in computers.  For `IEEE 754 <http://grouper.ieee.org/groups/754/>`_,
the first step in getting the binary format is to normalize the number,
or express the number as a number between 1 and 2 multiplied by a power of 2.
For our ``x`` above, we multiply by `2^{-3}` to get a number between 1 and 2.

::

    sage: (x*2^(-3)).str(base=2)
    '1.1010100100000000000000000000000000000000000000000000'

The fraction part of the normalized number is called the *mantissa* (or
*significand* ).

::

    sage: f=(x*2^(-3)).frac() # .frac() gets the fraction part, the part after the decimal
    sage: mantissa=f.str(base=2)
    sage: mantissa
    '0.10101001000000000000000000000000000000000000000000000'

Since we multiplied by :math:`2^{-3}` to normalize the number, we need
to store this information.  We add 1023 in order to not have to worry
about storing negative numbers.  In the below, ``c`` will be our
exponent.

::

    sage: c=(3+1023).str(base=2)
    sage: c
    '10000000010'

Note that ``c`` has 11 bits, which is exactly what we want.

::

    sage: len(c)
    11

Evaluating ``mantissa[2:54]`` will give
the first 52 binary digits after the decimal point of the
mantissa.  Note that we don't need to store the leading 1 before the
decimal point because it will always be there from the way we normalized
things.  This lets us get 53\-bit precision using only 52 bits of
storage.

::

    sage: len(mantissa[2:54])
    52

Since the original number was positive, our sign bit is zero.

::

    sage: sign='0'

So here is our 64\-bit double\-precision floating point number.

::

    sage: sign+' '+c+' '+mantissa[2:54] # the [2:] just chops off the '0.', since we just need to store the digits after the decimal point
    '0 10000000010 1010100100000000000000000000000000000000000000000000'

::

    sage: len(sign+c+mantissa[2:54]) # it's 64 bits!
    64

Here we convert back to our original number from the floating point
representation that we constructed.

::

    sage: ((-1)^(int(sign)) * 2^(int(c,base=2)-1023)*(1+RR(mantissa[:54], base=2)))
    13.2812500000000

::

    sage: x
    13.2812500000000

So they agree!

Sage uses a cutting\-edge numerical library, MPFR, to carry out precise
floating point arithmetic using any precision a user specifies.  MPFR
has a slightly different convention for normalization.  In MPFR, we
normalize by multiplying by an appropriate power of 2 to make the
mantissa an integer, instead of a binary fraction.  This allows us to
use big integer libraries and sophisticated techniques to carry out
calculations at an arbitrary precision.

::

    sage: x.sign_mantissa_exponent()
    (1, 7476679068876800, -49)

::

    sage: 7476679068876800*2^(-49)
    425/32

Note that the mantissa here has the same zero/nonzero bits as the
mantissa above (before we chopped off the leading 1 above).

::

    sage: 7476679068876800.str(base=2)
    '11010100100000000000000000000000000000000000000000000'

Interval Arithmetic
-------------------

Sage also lets you compute using intervals to keep track of error
bounds.  These basically use the round up and round down features shown
above.

::

    sage: ring=RealIntervalField(10)
    sage: a=ring(1/9)
    sage: a
    0.112?

The question mark notation means that the number is contained in the
interval found by incrementing and decrementing the last digit of the
number.  See the `documentation for real interval fields
<http://doc.sagemath.org/html/en/reference/sage/rings/real_mpfi.html>`_ for
details.  In the above case, Sage is saying that 1/9 is somewhere
between 0.111 and 0.113.  Below, we see that ``1/a`` is somewhere
between 8.9 and 9.1.

::

    sage: 1/a
    9.0?

We can get a more precise estimate of the interval if we explicitly
print out the interval.

::

    sage: print((1/a).str(style='brackets'))
    [8.9843 .. 9.0157]

Included Software
-----------------

Scipy (included in Sage) has a lot of numerical algorithms.  See `the
Scipy docs <http://docs.scipy.org/doc/scipy/reference/>`_.

Mpmath is also included in Sage, and contains a huge amount of numerical
stuff.  See `the mpmath codebase <https://github.com/fredrik-johansson/mpmath/>`_.

The `Decimal python module
<http://docs.python.org/library/decimal.html>`_ has also been useful for
textbook exercises which involved rounding in base 10.

Plotting with precision
-----------------------

Sometimes plotting involves some rather bad rounding errors because
plotting calculations are done with machine\-precision floating point
numbers.

::

    sage: f(x)=x^2*(sqrt(x^4+16)-x^2)
    sage: plot(f,(x,0,2e4))
    Graphics object consisting of 1 graphics primitive

We can instead make a function that specifically evaluates all
intermediate steps to 100 bits of precision using the ``fast_callable``
system.

::

    sage: R=RealField(100) # 100 bits
    sage: g=fast_callable(f, vars=[x], domain=R)
    sage: plot(g,(x,0,2e4))
    Graphics object consisting of 1 graphics primitive

