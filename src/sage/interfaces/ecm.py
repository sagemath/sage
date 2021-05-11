r"""
The Elliptic Curve Factorization Method

The elliptic curve factorization method (ECM) is the fastest way to
factor a **known composite** integer if one of the factors is
relatively small (up to approximately 80 bits / 25 decimal digits). To
factor an arbitrary integer it must be combined with a primality
test. The :meth:`ECM.factor` method is an example for how to combine
ECM with a primality test to compute the prime factorization of integers.

Sage includes GMP-ECM, which is a highly optimized implementation of
Lenstra's elliptic curve factorization method.  See
http://ecm.gforge.inria.fr for more about GMP-ECM.

AUTHORS:

These people wrote GMP-ECM:
Pierrick Gaudry, Jim Fougeron,
Laurent Fousse, Alexander Kruppa,
Dave Newman, Paul Zimmermann

BUGS:

Output from ecm is non-deterministic. Doctests should set the random
seed, but currently there is no facility to do so.

TESTS:

Check that the issues from :trac:`27199` are fixed::

    sage: n = 16262093986406371
    sage: ecm = ECM()
    sage: ecm.factor(n, B1=10)
    [1009, 1009, 1733, 3023, 3049]

    sage: n = 1308301 * (10^499 + 153)
    sage: ECM(B1=600).one_curve(n, c=1, sigma=10)
    [1308301, 100...00153]
"""

###############################################################################
#       Copyright (C) 2006, William Stein <wstein@gmail.com>
#       Copyright (C) 2006, Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2013, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 3 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################

import re
import subprocess

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ


class ECM(SageObject):

    def __init__(self, B1=10, B2=None, **kwds):
        r"""
        Create an interface to the GMP-ECM elliptic curve method
        factorization program.

        See http://ecm.gforge.inria.fr

        INPUT:

        - ``B1`` -- integer. Stage 1 bound

        - ``B2`` -- integer. Stage 2 bound (or interval B2min-B2max)

        In addition the following keyword arguments can be used:

        - ``x0`` -- integer `x`. use `x` as initial point

        - ``sigma`` -- integer `s`. Use s as curve generator [ecm]

        - ``A`` -- integer `a`. Use a as curve parameter [ecm]

        - ``k`` -- integer `n`. Perform `>= n` steps in stage 2

        - ``power`` -- integer `n`. Use `x^n` for Brent-Suyama's
          extension

        - ``dickson`` -- integer `n`. Use `n`-th Dickson's polynomial
          for Brent-Suyama's extension

        - ``c`` -- integer `n`. Perform `n` runs for each input

        - ``pm1`` --  boolean. perform P-1 instead of ECM

        - ``pp1`` --  boolean. perform P+1 instead of ECM

        - ``q`` -- boolean. quiet mode

        - ``v`` -- boolean. verbose mode

        - ``timestamp`` --  boolean. print a time stamp with each number

        - ``mpzmod`` -- boolean. use GMP's mpz_mod for mod reduction

        - ``modmuln`` -- boolean. use Montgomery's MODMULN for mod reduction

        - ``redc`` -- boolean. use Montgomery's REDC for mod reduction

        - ``nobase2`` -- boolean. Disable special base-2 code

        - ``base2`` -- integer `n`. Force base 2 mode with 2^n+1 (n>0)
          or 2^n-1 (n<0)

        - ``save`` -- string filename. Save residues at end of stage 1
          to file

        - ``savea`` -- string filename. Like -save, appends to
          existing files

        - ``resume`` -- string filename. Resume residues from file,
          reads from stdin if file is "-"

        - ``primetest`` -- boolean. Perform a primality test on input

        - ``treefile`` -- string. Store product tree of F in files f.0
          f.1 ...

        - ``i`` -- integer. increment B1 by this constant on each run

        - ``I`` -- integer `f`. auto-calculated increment for B1
          multiplied by `f` scale factor.

        - ``inp`` -- string. Use file as input (instead of redirecting
          stdin)

        - ``b`` -- boolean. Use breadth-first mode of file processing

        - ``d`` -- boolean. Use depth-first mode of file processing
          (default)

        - ``one`` -- boolean. Stop processing a candidate if a factor
          is found (looping mode )

        - ``n`` -- boolean. Run ecm in 'nice' mode (below normal
          priority)

        - ``nn`` -- boolean. Run ecm in 'very nice' mode (idle
          priority)

        - ``t`` -- integer `n`. Trial divide candidates before P-1,
          P+1 or ECM up to `n`.

        - ``ve`` -- integer `n`. Verbosely show short (`< n`
          character) expressions on each loop

        - ``B2scale`` -- integer. Multiplies the default B2 value

        - ``go`` -- integer. Preload with group order val, which can
          be a simple expression, or can use N as a placeholder for
          the number being factored.

        - ``prp`` -- string. use shell command cmd to do large
          primality tests

        - ``prplen`` -- integer.  only candidates longer than this
          number of digits are 'large'

        - ``prpval`` -- integer. value>=0 which indicates the prp
          command foundnumber to be PRP.

        - ``prptmp`` -- file. outputs n value to temp file prior to
          running (NB. gets deleted)

        - ``prplog`` -- file. otherwise get PRP results from this file
          (NB. gets deleted)

        - ``prpyes`` -- string. Literal string found in prplog file
          when number is PRP

        - ``prpno`` -- string. Literal string found in prplog file
          when number is composite
        """
        self._cmd = self._make_cmd(B1, B2, kwds)

    def _make_cmd(self, B1, B2, kwds):
        ecm = ['ecm']
        options = []
        for x, v in kwds.items():
            if v is False:
                continue
            options.append('-{0}'.format(x))
            if (v is not True) and (v != ''):
                options.append(str(v))
        if B2 is None:
            args = [str(B1)]
        else:
            args = [str(B1), str(B2)]
        return ecm + options + args

    def _run_ecm(self, cmd, n):
        """
        Run ECM and return output as string.

        INPUT:

        - ``cmd`` -- list of strings. The command.

        - ``n`` -- integer suitable for ECM. No argument checking is
          performed.

        OUTPUT:

        String.

        EXAMPLES::

            sage: ecm._run_ecm(['cat'], 1234)
            '1234'
        """
        from subprocess import Popen, PIPE

        # Under normal usage this program only returns ASCII; anything
        # else mixed is garbage and an error
        # So just accept latin-1 without encoding errors, and let the
        # output parser deal with the rest
        p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE, encoding='latin-1')
        out, err = p.communicate(input=str(n))
        if err != '':
            raise ValueError(err)
        return out

    def __call__(self, n):
        """
        Call syntax.

        INPUT:

        - ``n`` -- integer.

        OUTPUT:

        String. The ECM output.

        EXAMPLES::

            sage: print(ecm(3))    # random output
            GMP-ECM 6.4.4 [configured with MPIR 2.6.0, --enable-asm-redc] [ECM]
            Input number is 3 (1 digits)
            ********** Factor found in step 1: 3
            Found input number N
        """
        n = self._validate(n)
        return self._run_ecm(self._cmd, n)

    def interact(self):
        """
        Interactively interact with the ECM program.

        EXAMPLES::

            sage: ecm.interact()    # not tested
        """
        print("Enter numbers to run ECM on them.")
        print("Press control-D to exit.")
        subprocess.call(self._cmd)

    # Recommended settings from
    # http://www.mersennewiki.org/index.php/Elliptic_Curve_Method
    _recommended_B1_list = {15: 2000,
                            20: 11000,
                            25: 50000,
                            30: 250000,
                            35: 1000000,
                            40: 3000000,
                            45: 11000000,
                            50: 44000000,
                            55: 110000000,
                            60: 260000000,
                            65: 850000000,
                            70: 2900000000}

    def _B1_table_value(self, factor_digits, min=15, max=70):
        """
        Return key in ``_recommended_B1_list``.

        INPUT:

        - ``factor_digits`` -- integer. Number of digits.

        - ``min``, ``max`` -- integer. Min and max values.

        OUTPUT:

        Integer. A key in _recommended_B1_list.

        EXAMPLES::

            sage: ecm._B1_table_value(33)
            35
        """
        if factor_digits < min:
            factor_digits = min
        if factor_digits > max:
            raise ValueError('too many digits')
        step = 5
        return ((factor_digits + step - 1) // step) * step

    def recommended_B1(self, factor_digits):
        r"""
        Return recommended ``B1`` setting.

        INPUT:

        - ``factor_digits`` -- integer. Number of digits.

        OUTPUT:

        Integer. Recommended settings from
        http://www.mersennewiki.org/index.php/Elliptic_Curve_Method

        EXAMPLES::

            sage: ecm.recommended_B1(33)
            1000000
        """
        return self._recommended_B1_list[self._B1_table_value(factor_digits)]

    _parse_status_re = re.compile(
        r'Using B1=(\d+), B2=(\d+), polynomial ([^,]+), sigma=(\d+)')

    _found_input_re = re.compile('Found input number N')

    _found_factor_re = re.compile(
        r'Found (?P<primality>.*) factor of [\s]*(?P<digits>\d+) digits: (?P<factor>\d+)')

    _found_cofactor_re = re.compile(
        r'(?P<primality>.*) cofactor (?P<cofactor>\d+) has [\s]*(?P<digits>\d+) digits')

    def _parse_output(self, n, out):
        """
        Parse the ECM output

        INPUT:

        - ``n`` -- integer. The ECM input number.

        - ``out`` -- string. The stdout from the ECM invocation.

        OUTPUT:

        List of pairs ``(integer, bool)`` consisting of factors of the
        ECM input and whether they are deemed to be probable
        prime. Note that ECM is not a good primality test, and there
        is a sizeable probability that the "probable prime" is
        actually composite.

        EXAMPLES::

            sage: out = '\n'.join([
            ....:   'GMP-ECM 6.4.4 [configured with MPIR 2.6.0, --enable-asm-redc] [ECM]',
            ....:   'Input number is 1632143 (7 digits)',
            ....:   'Using B1=40, B2=480, polynomial x^1, sigma=3145777366',
            ....:   'Step 1 took 0ms',
            ....:   'Step 2 took 0ms',
            ....:   'Run 2 out of 1000:',
            ....:   'Using B1=40, B2=480, polynomial x^1, sigma=2101568373',
            ....:   'Step 1 took 0ms',
            ....:   'Step 2 took 0ms',
            ....:   '********** Factor found in step 2: 1632143',
            ....:   'Found input number N'])
            sage: ecm._parse_output(1632143, out)
            [(1632143, True)]
            sage: ecm.get_last_params()['sigma']
            '2101568373'

            sage: from sage.interfaces.ecm import TEST_ECM_OUTPUT_1
            sage: n1 = 508021860739623467191080372196682785441177798407961
            sage: ecm._parse_output(n1, TEST_ECM_OUTPUT_1)
            [(79792266297612017, True), (6366805760909027985741435139224233, True)]

            sage: from sage.interfaces.ecm import TEST_ECM_OUTPUT_2
            sage: ecm._parse_output(32193213281156929, TEST_ECM_OUTPUT_2)
            [(179424673, True), (179424673, True)]

            sage: from sage.interfaces.ecm import TEST_ECM_OUTPUT_3, TEST_ECM_OUTPUT_4
            sage: n3 = 66955751844124594814248420514215108438425124740949701470891
            sage: ecm._parse_output(n3, TEST_ECM_OUTPUT_3)
            [(197002597249, True),
             (339872432034468861533158743041639097889948066859, False)]
            sage: ecm._parse_output(n3, TEST_ECM_OUTPUT_4)
            [(265748496095531068869578877937, False),
             (251951573867253012259144010843, True)]
        """
        out_lines = out.lstrip().splitlines()
        if not out_lines[0].startswith('GMP-ECM'):
            raise ValueError('invalid output')
        result = []
        for line in out_lines:
            # print('parsing line >>{0}<<'.format(line))
            m = self._parse_status_re.match(line)
            if m is not None:
                group = m.groups()
                self._last_params = {'B1': group[0], 'B2': group[1],
                                     'poly': group[2], 'sigma': group[3]}
                continue
            m = self._found_input_re.match(line)
            if m is not None:
                return [(n, True)]
            m = self._found_factor_re.match(line)
            if m is not None:
                factor = m.group('factor')
                primality = m.group('primality')
                assert primality in ['prime', 'composite', 'probable prime']
                result += [(ZZ(factor), primality != 'composite')]
                continue  # cofactor on the next line
            m = self._found_cofactor_re.match(line)
            if m is not None:
                cofactor = m.group('cofactor')
                primality = m.group('primality')
                assert primality in ['Prime', 'Composite', 'Probable prime']
                result += [(ZZ(cofactor), primality != 'Composite')]
                # assert len(result) == 2
                return result
        raise ValueError('failed to parse ECM output')

    def one_curve(self, n, factor_digits=None, B1=2000, algorithm="ECM", **kwds):
        """
        Run one single ECM (or P-1/P+1) curve on input n.

        Note that trying a single curve is not particularly useful by
        itself. One typically needs to run over thousands of trial
        curves to factor `n`.

        INPUT:

        - ``n`` -- a positive integer

        - ``factor_digits`` -- integer. Decimal digits estimate of the
          wanted factor.

        - ``B1`` -- integer. Stage 1 bound (default 2000)

        - ``algorithm`` -- either "ECM" (default), "P-1" or "P+1"

        OUTPUT:

        a list ``[p, q]`` where p and q are integers and n = p * q.
        If no factor was found, then p = 1 and q = n.

        .. WARNING::

            Neither p nor q in the output is guaranteed to be prime.

        EXAMPLES::

            sage: f = ECM()
            sage: n = 508021860739623467191080372196682785441177798407961
            sage: f.one_curve(n, B1=10000, sigma=11)
            [1, 508021860739623467191080372196682785441177798407961]
            sage: f.one_curve(n, B1=10000, sigma=1022170541)
            [79792266297612017, 6366805760909027985741435139224233]
            sage: n = 432132887883903108009802143314445113500016816977037257
            sage: f.one_curve(n, B1=500000, algorithm="P-1")
            [67872792749091946529, 6366805760909027985741435139224233]
            sage: n = 2088352670731726262548647919416588631875815083
            sage: f.one_curve(n, B1=2000, algorithm="P+1", x0=5)
            [328006342451, 6366805760909027985741435139224233]
        """
        n = self._validate(n)
        if factor_digits is not None:
            B1 = self.recommended_B1(factor_digits)
        if algorithm == "P-1":
            kwds['pm1'] = ''
        elif algorithm == "P+1":
            kwds['pp1'] = ''
        elif algorithm == "ECM":
            pass
        else:
            raise ValueError('unknown algorithm')
        cmd = self._make_cmd(B1, None, kwds)
        out = self._run_ecm(cmd, n)
        try:
            factors = self._parse_output(n, out)
            return [factors[0][0], factors[1][0]]
        except (ValueError, IndexError):
            # output does not end in factorization (ValueError)
            # or factors has only one element above (IndexError)
            return [ZZ(1), n]

    def _find_factor(self, n, factor_digits, B1, **kwds):
        """
        Helper for :meth:`find_factor`.

        INPUT:

        See :meth:`find_factor`.

        OUTPUT:

        List of pairs ``(integer, bool)`` consisting of factors of the
        ECM input and whether they are probable prime. Note that ECM
        is not a good primality test and there is a sizeable chance
        that a "probable prime" is actually composite.

        EXAMPLES::

            sage: f = ECM()
            sage: n = 508021860739623467191080372196682785441177798407961
            sage: sorted(f._find_factor(n, None, 2000))
            [(79792266297612017, True),
             (6366805760909027985741435139224233, True)]
        """
        n = self._validate(n)
        kwds.setdefault('c', 1000000000)
        kwds.setdefault('I', 1)
        if factor_digits is not None:
            B1 = self.recommended_B1(factor_digits)
        kwds['one'] = True
        cmd = self._make_cmd(B1, None, kwds)
        out = self._run_ecm(cmd, n)
        return self._parse_output(n, out)

    def find_factor(self, n, factor_digits=None, B1=2000, **kwds):
        """
        Return a factor of n.

        See also :meth:`factor` if you want a prime factorization of
        `n`.

        INPUT:

        - ``n`` -- a positive integer,

        - ``factor_digits`` -- integer or ``None`` (default). Decimal
          digits estimate of the wanted factor.

        - ``B1`` -- integer. Stage 1 bound (default 2000). This is
          used as bound if ``factor_digits`` is not specified.

        - ``kwds`` -- optional keyword parameters.

        OUTPUT:

        List of integers whose product is n. For certain lengths of
        the factor, this is the best algorithm to find a
        factor.

        .. NOTE::

            ECM is not a good primality test. Not finding a
            factorization is only weak evidence for `n` being
            prime. You should run a **good** primality test before
            calling this function.

        EXAMPLES::

            sage: f = ECM()
            sage: n = 508021860739623467191080372196682785441177798407961
            sage: f.find_factor(n)
            [79792266297612017, 6366805760909027985741435139224233]

        Note that the input number cannot have more than 4095 digits::

            sage: f = 2^2^14+1
            sage: ecm.find_factor(f)
            Traceback (most recent call last):
            ...
            ValueError: n must have at most 4095 digits
        """
        factors = self._find_factor(n, factor_digits, B1, **kwds)
        return [factor[0] for factor in factors]

    def factor(self, n, factor_digits=None, B1=2000, proof=False, **kwds):
        """
        Return a probable prime factorization of `n`.

        Combines GMP-ECM with a primality test, see
        :meth:`~sage.rings.integer.Integer.is_prime`. The primality
        test is provable or probabilistic depending on the `proof`
        flag.

        Moreover, for small `n` PARI is used directly.

        .. WARNING::

            There is no mathematical guarantee that the factors
            returned are actually prime if ``proof=False``
            (default). It is extremely likely, though. Currently,
            there are no known examples where this fails.

        INPUT:

        - ``n`` -- a positive integer

        - ``factor_digits`` -- integer or ``None`` (default). Optional
          guess at how many digits are in the smallest factor.

        - ``B1`` -- initial lower bound, defaults to 2000 (15 digit
          factors). Used if ``factor_digits`` is not specified.

        - ``proof`` -- boolean (default: ``False``). Whether to prove
          that the factors are prime.

        - ``kwds`` -- keyword arguments to pass to ecm-gmp. See help
          for :class:`ECM` for more details.

        OUTPUT:

        A list of integers whose product is n.

        .. NOTE::

            Trial division should typically be performed, but this is
            not implemented (yet) in this method.

            If you suspect that n is the product of two
            similarly-sized primes, other methods (such as a quadratic
            sieve -- use the qsieve command) will usually be faster.

            The best known algorithm for factoring in the case where
            all factors are large is the general number field
            sieve. This is not implemented in Sage; You probably want
            to use a cluster for problems of this size.

        EXAMPLES::

            sage: ecm.factor(602400691612422154516282778947806249229526581)
            [45949729863572179, 13109994191499930367061460439]
            sage: ecm.factor((2^197 + 1)/3)  # long time
            [197002597249, 1348959352853811313, 251951573867253012259144010843]
            sage: ecm.factor(179427217^13) == [179427217] * 13
            True
        """
        n = self._validate(n)
        factors = [n]                 # factors that need to be factorized futher
        probable_prime_factors = []   # output prime factors
        while factors:
            n = factors.pop()

            # Step 0: Primality test
            if n.is_prime(proof=proof):
                probable_prime_factors.append(n)
                continue

            # Step 1: Use PARI directly for small primes
            if n.ndigits() < 15:
                for p, e in n.factor(algorithm='pari'):
                    probable_prime_factors.extend([p] * e)
                continue

            # Step 2: Deal with small factors efficiently
            # Step 2+1/3: Determine if N is a perfect power
            if n.is_perfect_power():
                base, exp = n.perfect_power()
                factors.extend([base] * exp)
                continue

            # Step 2+2/3: Do trial division to remove small prime
            # factors, and maybe some other factorization algorithms
            # that perform well on small ranges. This all depends on
            # the kind of number you are trying to factor (todo)

            # Step 3: Call find_factor until a factorization is found
            n_factorization = [n]
            while len(n_factorization) == 1:
                n_factorization = self.find_factor(n,B1=B1)
            factors.extend(n_factorization)

        return sorted(probable_prime_factors)

    def get_last_params(self):
        """
        Return the parameters (including the curve) of the last ecm run.

        In the case that the number was factored successfully, this
        will return the parameters that yielded the factorization.

        OUTPUT:

        A dictionary containing the parameters for the most recent
        factorization.

        EXAMPLES::

            sage: ecm.factor((2^197 + 1)/3)             # long time
            [197002597249, 1348959352853811313, 251951573867253012259144010843]
            sage: ecm.get_last_params()                 # random output
            {'poly': 'x^1', 'sigma': '1785694449', 'B1': '8885', 'B2': '1002846'}
        """
        return self._last_params

    def time(self, n, factor_digits, verbose=False):
        """
        Print a runtime estimate.

        BUGS:

        This method should really return something and not just print
        stuff on the screen.

        INPUT:

        - ``n`` -- a positive integer

        - ``factor_digits`` -- the (estimated) number of digits of the
          smallest factor

        OUTPUT:

        An approximation for the amount of time it will take to find a
        factor of size factor_digits in a single process on the
        current computer.  This estimate is provided by GMP-ECM's
        verbose option on a single run of a curve.

        EXAMPLES::

            sage: n = next_prime(11^23)*next_prime(11^37)
            sage: ecm.time(n, 35)                  # random output
            Expected curves: 910, Expected time: 23.95m

            sage: ecm.time(n, 30, verbose=True)     # random output
            GMP-ECM 6.4.4 [configured with MPIR 2.6.0, --enable-asm-redc] [ECM]
            Running on localhost.localdomain
            Input number is 304481639541418099574459496544854621998616257489887231115912293 (63 digits)
            Using MODMULN [mulredc:0, sqrredc:0]
            Using B1=250000, B2=128992510, polynomial Dickson(3), sigma=3244548117
            dF=2048, k=3, d=19110, d2=11, i0=3
            Expected number of curves to find a factor of n digits:
            35  40  45  50  55  60  65  70  75  80
            4911  70940  1226976  2.5e+07  5.8e+08  1.6e+10  2.7e+13  4e+18  5.4e+23  Inf
            Step 1 took 230ms
            Using 10 small primes for NTT
            Estimated memory usage: 4040K
            Initializing tables of differences for F took 0ms
            Computing roots of F took 9ms
            Building F from its roots took 16ms
            Computing 1/F took 9ms
            Initializing table of differences for G took 0ms
            Computing roots of G took 8ms
            Building G from its roots took 16ms
            Computing roots of G took 7ms
            Building G from its roots took 16ms
            Computing G * H took 6ms
            Reducing  G * H mod F took 5ms
            Computing roots of G took 7ms
            Building G from its roots took 17ms
            Computing G * H took 5ms
            Reducing  G * H mod F took 5ms
            Computing polyeval(F,G) took 34ms
            Computing product of all F(g_i) took 0ms
            Step 2 took 164ms
            Expected time to find a factor of n digits:
            35  40  45  50  55  60  65  70  75  80
            32.25m  7.76h  5.60d  114.21d  7.27y  196.42y  337811y  5e+10y  7e+15y  Inf
            <BLANKLINE>
            Expected curves: 4911, Expected time: 32.25m
        """
        title_curves = 'Expected number of curves to find a factor of n digits:'
        title_time = 'Expected time to find a factor of n digits:'
        n = self._validate(n)
        B1 = self.recommended_B1(factor_digits)
        cmd = self._make_cmd(B1, None, {'v': True})
        out = self._run_ecm(cmd, n)
        if verbose:
            print(out)
        if title_time not in out:
            print('Unable to compute timing, factorized immediately')
            return

        out_lines = iter(out.splitlines())
        while next(out_lines) != title_curves:
            pass
        header_curves = next(out_lines)
        curve_count_table = next(out_lines)

        while next(out_lines) != title_time:
            pass
        header_time = next(out_lines)
        time_table = next(out_lines)

        assert header_curves == header_time
        assert header_curves.split() == [
            '35', '40', '45', '50', '55', '60', '65', '70', '75', '80']
        h_min = 35
        h_max = 80
        offset = (self._B1_table_value(factor_digits, h_min, h_max) - h_min) // 5
        print('offset', offset)
        curve_count = curve_count_table.split()[offset]
        time = time_table.split()[offset]
        print('Expected curves: {0}, Expected time: {1}'.format(curve_count, time))

    def _validate(self, n):
        """
        Verify that n is positive and has at most 4095 digits.

        INPUT:

        - ``n`` -- integer.

        OUTPUT:

        The integer as a Sage integer.  This function raises a
        ValueError if the two conditions listed above are not both
        satisfied.  It is here because GMP-ECM silently ignores all
        digits of input after the 4095th!

        EXAMPLES::

            sage: ecm = ECM()
            sage: ecm._validate(3)
            3
            sage: ecm._validate(0)
            Traceback (most recent call last):
            ...
            ValueError: n must be positive
            sage: ecm._validate(10^5000)
            Traceback (most recent call last):
            ...
            ValueError: n must have at most 4095 digits
        """
        n = ZZ(n)
        if n <= 0:
            raise ValueError("n must be positive")
        if n.ndigits() > 4095:
            raise ValueError("n must have at most 4095 digits")
        return n


# unique instance
ecm = ECM()


# Tests
TEST_ECM_OUTPUT_1 = """
GMP-ECM 6.4.4 [configured with MPIR 2.6.0, --enable-asm-redc] [ECM]
Input number is 508021860739623467191080372196682785441177798407961 (51 digits)
Using B1=2000, B2=147396, polynomial x^1, sigma=2005325688
Step 1 took 1ms
Step 2 took 2ms
Run 2 out of 1000000000:
Using B1=2399, B2=2399-186156, polynomial x^1, sigma=3689070339
Step 1 took 3ms
Step 2 took 2ms
[...]
Run 29 out of 1000000000:
Using B1=16578, B2=16578-3162402, polynomial x^1, sigma=2617498039
Step 1 took 12ms
Step 2 took 17ms
********** Factor found in step 2: 79792266297612017
Found prime factor of 17 digits: 79792266297612017
Prime cofactor 6366805760909027985741435139224233 has 34 digits
"""

TEST_ECM_OUTPUT_2 = """
GMP-ECM 6.4.4 [configured with MPIR 2.6.0, --enable-asm-redc] [ECM]
Input number is 32193213281156929 (17 digits)
Using B1=2000, B2=147396, polynomial x^1, sigma=434130265
Step 1 took 2ms
Step 2 took 3ms
********** Factor found in step 2: 179424673
Found prime factor of  9 digits: 179424673
Prime cofactor 179424673 has 9 digits
"""

TEST_ECM_OUTPUT_3 = """
GMP-ECM 6.4.4 [configured with MPIR 2.6.0, --enable-asm-redc] [ECM]
Input number is 66955751844124594814248420514215108438425124740949701470891 (59 digits)
Using B1=2000, B2=147396, polynomial x^1, sigma=553262339
Step 1 took 3ms
Step 2 took 4ms
Run 2 out of 1000000000:
Using B1=2399, B2=2399-186156, polynomial x^1, sigma=557154369
Step 1 took 5ms
Step 2 took 4ms
Run 3 out of 1000000000:
Using B1=2806, B2=2806-224406, polynomial x^1, sigma=478195111
Step 1 took 5ms
Step 2 took 4ms
********** Factor found in step 2: 197002597249
Found prime factor of 12 digits: 197002597249
Composite cofactor 339872432034468861533158743041639097889948066859 has 48 digits
"""

TEST_ECM_OUTPUT_4 = """
GMP-ECM 6.4.4 [configured with MPIR 2.6.0, --enable-asm-redc] [ECM]
Input number is 66955751844124594814248420514215108438425124740949701470891 (59 digits)
Using B1=2000, B2=147396, polynomial x^1, sigma=1881424010\n
Step 1 took 4ms
Step 2 took 2ms
********** Factor found in step 2: 265748496095531068869578877937
Found composite factor of 30 digits: 265748496095531068869578877937
Prime cofactor 251951573867253012259144010843 has 30 digits
"""
