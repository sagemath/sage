r"""
Watkins Symmetric Power `L`-function Calculator

SYMPOW is a package to compute special values of symmetric power
elliptic curve L-functions. It can compute up to about 64 digits of
precision. This interface provides complete access to sympow, which
is a standard part of Sage (and includes the extra data files).

.. note::

    Each call to ``sympow`` runs a complete
    ``sympow`` process. This incurs about 0.2 seconds
    overhead.

AUTHORS:

- Mark Watkins (2005-2006): wrote and released sympow

- William Stein (2006-03-05): wrote Sage interface

ACKNOWLEDGEMENT (from sympow readme):

-  The quad-double package was modified from David Bailey's
   package: http://crd.lbl.gov/~dhbailey/mpdist/

-  The ``squfof`` implementation was modified from
   Allan Steel's version of Arjen Lenstra's original LIP-based code.

-  The ``ec_ap`` code was originally written for the
   kernel of MAGMA, but was modified to use small integers when
   possible.

-  SYMPOW was originally developed using PARI, but due to licensing
   difficulties, this was eliminated. SYMPOW also does not use the
   standard math libraries unless Configure is run with the -lm
   option. SYMPOW still uses GP to compute the meshes of inverse
   Mellin transforms (this is done when a new symmetric power is added
   to datafiles).
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################

import os

from sage.structure.sage_object import SageObject
from sage.misc.all import pager
from sage.misc.verbose import verbose
import sage.rings.all


class Sympow(SageObject):
    r"""
    Watkins Symmetric Power `L`-function Calculator

    Type ``sympow.[tab]`` for a list of useful commands
    that are implemented using the command line interface, but return
    objects that make sense in Sage.

    You can also use the complete command-line interface of sympow via
    this class. Type ``sympow.help()`` for a list of
    commands and how to call them.
    """
    def _repr_(self):
        """
        Return a string describing this calculator module
        """
        return "Watkins Symmetric Power L-function Calculator"

    def __call__(self, args):
        """
        Used to call sympow with given args
        """
        cmd = 'sympow %s' % args
        with os.popen(cmd) as f:
            v = f.read().strip()
        verbose(v, level=2)
        return v

    def _fix_err(self, err):
        w = err
        j = w.rfind('./sympow')
        if j != -1:
            w = w[:j - 1] + "sympow('" + w[j + 9:] + ')'
        return w

    def _curve_str(self, E):
        return '-curve "%s"' % (str(list(E.minimal_model().a_invariants())).replace(' ', ''))

    def L(self, E, n, prec):
        r"""
        Return `L(\mathrm{Sym}^{(n)}(E, \text{edge}))` to prec digits of
        precision, where edge is the *right* edge. Here `n` must be
        even.

        INPUT:


        -  ``E`` - elliptic curve

        -  ``n`` - even integer

        -  ``prec`` - integer


        OUTPUT:


        -  ``string`` - real number to prec digits of precision
           as a string.


        .. note::

           Before using this function for the first time for a given
           `n`, you may have to type ``sympow('-new_data n')``,
           where ``n`` is replaced by your value of `n`.

        If you would like to see the extensive output sympow prints when
        running this function, just type ``set_verbose(2)``.

        EXAMPLES:

        These examples only work if you run ``sympow -new_data 2`` in a
        Sage shell first. Alternatively, within Sage, execute::

            sage: sympow('-new_data 2')  # not tested

        This command precomputes some data needed for the following
        examples. ::

            sage: a = sympow.L(EllipticCurve('11a'), 2, 16)  # not tested
            sage: a                                          # not tested
            '1.057599244590958E+00'
            sage: RR(a)                                      # not tested
            1.05759924459096
        """
        if n % 2 == 1:
            raise ValueError("n (=%s) must be even" % n)
        if prec > 64:
            raise ValueError("prec (=%s) must be at most 64" % prec)
        if prec < 1:
            raise ValueError("prec (=%s) must be at least 1" % prec)
        v = self('-sp %sp%s %s' % (n, prec, self._curve_str(E)))
        i = v.rfind(': ')
        if i == -1:
            print(self._fix_err(v))
            raise RuntimeError("failed to compute symmetric power")
        x = v[i + 2:]
        return x

    def Lderivs(self, E, n, prec, d):
        r"""
        Return `0^{th}` to `d^{th}` derivatives of
        `L(\mathrm{Sym}^{(n)}(E,s)` to prec digits of precision, where
        `s` is the right edge if `n` is even and the center
        if `n` is odd.

        INPUT:


        -  ``E`` - elliptic curve

        -  ``n`` - integer (even or odd)

        -  ``prec`` - integer

        -  ``d`` - integer


        OUTPUT: a string, exactly as output by sympow

        .. note::

           To use this function you may have to run a few commands
           like ``sympow('-new_data 1d2')``, each which takes a
           few minutes. If this function fails it will indicate what commands
           have to be run.

        EXAMPLES::

            sage: print(sympow.Lderivs(EllipticCurve('11a'), 1, 16, 2))  # not tested
            ...
             1n0: 2.538418608559107E-01
             1w0: 2.538418608559108E-01
             1n1: 1.032321840884568E-01
             1w1: 1.059251499158892E-01
             1n2: 3.238743180659171E-02
             1w2: 3.414818600982502E-02
        """
        if prec > 64:
            raise ValueError("prec (=%s) must be at most 64" % prec)
        if prec < 1:
            raise ValueError("prec (=%s) must be at least 1" % prec)
        v = self('-sp %sp%sd%s %s' % (n, prec, d, self._curve_str(E)))
        return self._fix_err(v)

    def modular_degree(self, E):
        """
        Return the modular degree of the elliptic curve E, assuming the
        Stevens conjecture.

        INPUT:


        -  ``E`` - elliptic curve over Q


        OUTPUT:


        -  ``integer`` - modular degree


        EXAMPLES: We compute the modular degrees of the lowest known
        conductor curves of the first few ranks::

            sage: sympow.modular_degree(EllipticCurve('11a'))
            1
            sage: sympow.modular_degree(EllipticCurve('37a'))
            2
            sage: sympow.modular_degree(EllipticCurve('389a'))
            40
            sage: sympow.modular_degree(EllipticCurve('5077a'))
            1984
            sage: sympow.modular_degree(EllipticCurve([1, -1, 0, -79, 289]))
            334976
        """
        v = self('%s -moddeg' % self._curve_str(E))
        s = 'Modular Degree is '
        i = v.find(s)
        if i == -1:
            print(self._fix_err(v))
            raise RuntimeError("failed to compute modular degree")
        return sage.rings.all.Integer(v[i + len(s):])

    def analytic_rank(self, E):
        r"""
        Return the analytic rank and leading `L`-value of the
        elliptic curve `E`.

        INPUT:


        -  ``E`` - elliptic curve over Q


        OUTPUT:


        -  ``integer`` - analytic rank

        -  ``string`` - leading coefficient (as string)


        .. note::

           The analytic rank is *not* computed provably correctly in
           general.

        .. note::

           In computing the analytic rank we consider
           `L^{(r)}(E,1)` to be `0` if
           `L^{(r)}(E,1)/\Omega_E > 0.0001`.

        EXAMPLES: We compute the analytic ranks of the lowest known
        conductor curves of the first few ranks::

            sage: sympow.analytic_rank(EllipticCurve('11a'))
            (0, '2.53842e-01')
            sage: sympow.analytic_rank(EllipticCurve('37a'))
            (1, '3.06000e-01')
            sage: sympow.analytic_rank(EllipticCurve('389a'))
            (2, '7.59317e-01')
            sage: sympow.analytic_rank(EllipticCurve('5077a'))
            (3, '1.73185e+00')
            sage: sympow.analytic_rank(EllipticCurve([1, -1, 0, -79, 289]))
            (4, '8.94385e+00')
            sage: sympow.analytic_rank(EllipticCurve([0, 0, 1, -79, 342]))  # long time
            (5, '3.02857e+01')
            sage: sympow.analytic_rank(EllipticCurve([1, 1, 0, -2582, 48720]))  # long time
            (6, '3.20781e+02')
            sage: sympow.analytic_rank(EllipticCurve([0, 0, 0, -10012, 346900]))  # long time
            (7, '1.32517e+03')
        """
        v = self('%s -analrank' % self._curve_str(E))
        s = 'Analytic Rank is '
        i = v.rfind(s)
        if i == -1:
            print(self._fix_err(v))
            raise RuntimeError("failed to compute analytic rank")
        j = v.rfind(':')
        r = sage.rings.all.Integer(v[i + len(s):j])
        i = v.rfind(' ')
        L = v[i + 1:]
        return r, L

    def new_data(self, n):
        """
        Pre-compute data files needed for computation of n-th symmetric
        powers.
        """
        print(self('-new_data %s' % n))

    def help(self):
        h = """
     sympow('-sp 2p16 -curve "[1,2,3,4,5]"')
will compute L(Sym^2 E,edge) for E=[1,2,3,4,5] to 16 digits of precision
The result
 2n0: 8.370510845377639E-01
 2w0: 8.370510845377586E-01
consists of two calculations using different parameters to test the
functional equation (to see if these are sufficiently close).

     sympow('-sp 3p12d2,4p8 -curve "[1,2,3,4,5]"')
will compute the 0th-2nd derivatives of L(Sym^3 E,center) to 12 digits
and L(Sym^4 E,edge) to 8 digits

Special cases: When a curve has CM, Hecke power can be used instead

     sympow('-sp 7p12d1 -hecke -curve "[0,0,1,-16758,835805]"')

will compute the 0th-1st derivatives of L(Sym^7 psi,special) to 12 digits.

Bloch-Kato numbers can be obtained for powers not congruent to 0 mod 4:

     sympow('-sp 2bp16 -curve "[1,2,3,4,5]"')

should return
 2n0: 4.640000000000006E+02
 2w0: 4.639999999999976E+02
which can be seen to be very close to the integer 464.

Modular degrees can be computed with the -moddeg option.

     sympow('-curve "[1,2,3,4,5]" -moddeg')

should return
Modular Degree is 464

Analytic ranks can be computed with the -analrank option.

     sympow('-curve "[1,2,3,4,5]" -analrank')

should return
Analytic Rank is 1 : L'-value 3.51873e+00

and (if the mesh file for the fifth derivative is present)

     sympow('-curve "[0,0,1,-79,342]" -analrank')

should return
Analytic Rank is 5 : leading L-term 3.02857e+01

========================================================================

Adding new symmetric powers:

SYMPOW keeps data for symmetric powers in the directory datafiles
If a pre-computed mesh of inverse Mellin transform values is not
available for a given symmetric power, SYMPOW will fail. The command

     sympow('-new_data 2')

will add the data the 2nd symmetric power, while

     sympow('-new_data 3d2')

will add the data for the 2nd derivative and 3rd symmetric power,

     sympow('-new_data 6d0h')

will add the data for the 0th derivative of the 6th Hecke power, and

     sympow('-new_data 4c')

will add data for the 4th symmetric power for curves with CM
(these need to be done separately for powers divisible by 4).

The mesh files are stored in binary form, and thus endian-ness is a
worry when moving from one platform to another.

To enable modular degree computations, the 2nd symmetric power must
be extant, and analytic rank requires derivatives of the 1st power.

===================================================================

Output of "!sympow -help":

 -bound #      an upper BOUND for how many ap to compute
 -help         print the help message and exit
 -info [] []   only report local information for primes/sympows
               1st argument is prime range, 2nd is sympow range
 -local        only report local information (bad primes)
 -curve []     input a curve in [a1,a2,a3,a4,a6] form
 -label []     label the given curve
 -quiet        turn off some messages
 -verbose      turn on some messages
 -rootno #     compute the root number of the #th symmetric power
 -moddeg       compute the modular degree
 -analrank     compute the analytic rank
 -sloppy []    for use with -analrank; have X sloppy digits
 -nocm         abort if curve has complex multiplication
 -noqt         ignore even powers of non-minimal quad twists
 -hecke        compute Hecke symmetric powers for a CM curve
 -maxtable     set the max size of factor tables: 2^27 default
 -sp  []       argument to specify which powers
               this is a comma separated list
               in each entry, the 1st datum is the sympow
               then could come b which turns Bloch-Kato on
               then could come w# which specifies how many tests
               then could come s# which says # sloppy digits
               then must come p# which specifies the precision
                    or P# which says ignore BOUND for this power
               then must come d# which says the derivative bound
                    or D# which says do only this derivative
                    (neither need be indicated for even powers)
               default is 2w3s1p32,3bp16d1,4p8
 -new_data []  will compute inverse Mellin transform mesh for
               the given data: the format is [sp]d[dv]{h,c}
               sp is the symmetric power, dv is the derivative,
               h indicates Hecke powers, and c indicates CM case
               d[dv] is given only for odd or Hecke powers
               Examples: 1d3 2 2d1h 3d2 4 4c 5d0 6 7d0h 11d1 12c
               NOTE: new_data runs a shell script that uses GP
Other options are used internally/recursively by -new_data
"""
        pager()(h)


# An instance
sympow = Sympow()
