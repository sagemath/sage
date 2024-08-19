r"""
Interface to Frobby for fast computations on monomial ideals.

The software package Frobby provides a number of computations on
monomial ideals. The current main feature is the socle of a monomial
ideal, which is largely equivalent to computing the maximal standard
monomials, the Alexander dual or the irreducible decomposition.

Operations on monomial ideals are much faster than algorithms designed
for ideals in general, which is what makes a specialized library for
these operations on monomial ideals useful.

AUTHORS:

- Bjarke Hammersholt Roune (2008-04-25): Wrote the Frobby C++
  program and the initial version of the Python interface.

.. NOTE::

    The official source for Frobby is <https://www.broune.com/frobby>,
    which also has documentation and papers describing the algorithms used.
"""

from subprocess import Popen, PIPE
from sage.misc.misc_c import prod

from sage.cpython.string import bytes_to_str, str_to_bytes


class Frobby:
    def __call__(self, action, input=None, options=[], verbose=False):
        r"""
        This function calls Frobby as a command line program using streams
        for input and output. Strings passed as part of the command get
        broken up at whitespace. This is not done to the data passed
        via the streams.

        INPUT:

        - ``action`` -- string telling Frobby what to do
        - ``input`` -- ``None`` or string that is passed to Frobby as standard in
        - ``options`` -- list of options without the dash in front
        - ``verbose`` -- boolean (default: ``False``); print detailed information

        OUTPUT: string; what Frobby wrote to the standard output stream

        EXAMPLES:

        We compute the lcm of an ideal provided in Monos format. ::

            sage: frobby("analyze", input="vars x,y,z;[x^2,x*y];", # optional - frobby
            ....:     options=["lcm", "iformat monos", "oformat 4ti2"])
            ' 2 1 0\n\n2 generators\n3 variables\n'


        We get an exception if frobby reports an error. ::

            sage: frobby("do_dishes") # optional - frobby
            Traceback (most recent call last):
            ...
            RuntimeError: Frobby reported an error:
            ERROR: No action has the prefix "do_dishes".

        AUTHOR:

        - Bjarke Hammersholt Roune (2008-04-27)
        """
        command = ['frobby'] + action.split()
        for option in options:
            command += ('-' + option.strip()).split()

        if verbose:
            print("Frobby action: ", action)
            print("Frobby options: ", repr(options))
            print("Frobby command: ", repr(command))
            print("Frobby input:\n", input)

        process = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        if input:
            frinput = str_to_bytes(input)
        else:
            frinput = None
        output, err = process.communicate(input=frinput)
        output = bytes_to_str(output)
        err = bytes_to_str(err)

        if verbose:
            print("Frobby output:\n", output)
            print("Frobby error:\n", err)
        if process.poll() != 0:
            raise RuntimeError("Frobby reported an error:\n" + err)

        return output

    def alexander_dual(self, monomial_ideal):
        r"""
        This function computes the Alexander dual of the passed-in
        monomial ideal. This ideal is the one corresponding to the
        simplicial complex whose faces are the complements of the
        nonfaces of the simplicial complex corresponding to the input
        ideal.

        INPUT:

        - ``monomial_ideal`` -- the monomial ideal to decompose

        OUTPUT: the monomial corresponding to the Alexander dual

        EXAMPLES:

        This is a simple example of computing irreducible decomposition. ::

            sage: # optional - frobby
            sage: (a, b, c, d) = QQ['a,b,c,d'].gens()
            sage: id = ideal(a * b, b * c, c * d, d * a)
            sage: alexander_dual = frobby.alexander_dual(id)
            sage: true_alexander_dual = ideal(b * d, a * c)
            sage: alexander_dual == true_alexander_dual # use sets to ignore order
            True

        We see how it is much faster to compute this with frobby than the built-in
        procedure for simplicial complexes::

            sage: # optional - frobby
            sage: t=simplicial_complexes.PoincareHomologyThreeSphere()
            sage: R=PolynomialRing(QQ,16,'x')
            sage: I=R.ideal([prod([R.gen(i-1) for i in a]) for a in t.facets()])
            sage: len(frobby.alexander_dual(I).gens())
            643
        """
        frobby_input = self._ideal_to_string(monomial_ideal)
        frobby_output = self('alexdual', input=frobby_input)
        return self._parse_ideals(frobby_output, monomial_ideal.ring())[0]

    def hilbert(self, monomial_ideal):
        r"""
        Compute the multigraded Hilbert-Poincaré series of the input
        ideal. Use the -univariate option to get the univariate series.

        The Hilbert-Poincaré series of a monomial ideal is the sum of all
        monomials not in the ideal. This sum can be written as a (finite)
        rational function with `(x_1-1)(x_2-1)...(x_n-1)` in the denominator,
        assuming the variables of the ring are `x_1,x2,...,x_n`. This action
        computes the polynomial in the numerator of this fraction.

        INPUT:

        - ``monomial_ideal`` -- a monomial ideal

        OUTPUT:

        A polynomial in the same ring as the ideal.

        EXAMPLES::

            sage: R.<d,b,c>=QQ[] # optional - frobby
            sage: I=[d*b*c,b^2*c,b^10,d^10]*R # optional - frobby
            sage: frobby.hilbert(I) # optional - frobby
            d^10*b^10*c + d^10*b^10 + d^10*b*c + b^10*c + d^10 + b^10 + d*b^2*c + d*b*c + b^2*c + 1
        """
        frobby_input = self._ideal_to_string(monomial_ideal)
        frobby_output = self('hilbert', input=frobby_input)
        ring = monomial_ideal.ring()
        lines = frobby_output.split('\n')
        if lines[-1] == '':
            lines.pop(-1)
        if lines[-1] == '(coefficient)':
            lines.pop(-1)
        lines.pop(0)
        resul = 0
        for l in lines:
            lis = [int(_) for _ in l.split()]
            resul += lis[0]+prod([ring.gen(i)**lis[i+1] for i in range(len(lis)-1)])
        return resul

    def associated_primes(self, monomial_ideal):
        r"""
        This function computes the associated primes of the passed-in
        monomial ideal.

        INPUT:

        - ``monomial_ideal`` -- the monomial ideal to decompose

        OUTPUT:

        A list of the associated primes of the monomial ideal. These ideals
        are constructed in the same ring as monomial_ideal is.

        EXAMPLES::

            sage: R.<d,b,c>=QQ[] # optional - frobby
            sage: I=[d*b*c,b^2*c,b^10,d^10]*R # optional - frobby
            sage: frobby.associated_primes(I)   # optional - frobby
            [Ideal (d, b) of Multivariate Polynomial Ring in d, b, c over Rational Field,
            Ideal (d, b, c) of Multivariate Polynomial Ring in d, b, c over Rational Field]
        """
        frobby_input = self._ideal_to_string(monomial_ideal)
        frobby_output = self('assoprimes', input=frobby_input)
        lines = frobby_output.split('\n')
        lines.pop(0)
        if lines[-1] == '':
            lines.pop(-1)
        lists = [[int(_) for _ in a.split()] for a in lines]

        def to_monomial(exps):
            return [v ** e for v, e in zip(monomial_ideal.ring().gens(), exps) if e != 0]
        return [monomial_ideal.ring().ideal(to_monomial(a)) for a in lists]

    def dimension(self, monomial_ideal):
        r"""
        This function computes the dimension of the passed-in
        monomial ideal.

        INPUT:

        - ``monomial_ideal`` -- the monomial ideal to decompose

        OUTPUT: the dimension of the zero set of the ideal

        EXAMPLES::

            sage: R.<d,b,c>=QQ[] # optional - frobby
            sage: I=[d*b*c,b^2*c,b^10,d^10]*R # optional - frobby
            sage: frobby.dimension(I)   # optional - frobby
            1
        """
        frobby_input = self._ideal_to_string(monomial_ideal)
        frobby_output = self('dimension', input=frobby_input)
        return int(frobby_output)

    def irreducible_decomposition(self, monomial_ideal):
        r"""
        This function computes the irreducible decomposition of the passed-in
        monomial ideal. I.e. it computes the unique minimal list of
        irreducible monomial ideals whose intersection equals monomial_ideal.

        INPUT:

        - ``monomial_ideal`` -- the monomial ideal to decompose

        OUTPUT:

        A list of the unique irredundant irreducible components of
        monomial_ideal. These ideals are constructed in the same ring
        as monomial_ideal is.

        EXAMPLES:

        This is a simple example of computing irreducible decomposition. ::

            sage: # optional - frobby
            sage: (x, y, z) = QQ['x,y,z'].gens()
            sage: id = ideal(x ** 2, y ** 2, x * z, y * z)
            sage: decom = frobby.irreducible_decomposition(id)
            sage: true_decom = [ideal(x, y), ideal(x ** 2, y ** 2, z)]
            sage: set(decom) == set(true_decom) # use sets to ignore order
            True

        We now try the special case of the zero ideal in different rings.

        We should also try PolynomialRing(QQ, names=[]), but it has a bug
        which makes that impossible (see :issue:`3028`). ::

            sage: # optional - frobby
            sage: rings = [ZZ['x'], CC['x,y']]
            sage: allOK = True
            sage: for ring in rings:
            ....:     id0 = ring.ideal(0)
            ....:     decom0 = frobby.irreducible_decomposition(id0)
            ....:     allOK = allOK and decom0 == [id0]
            sage: allOK
            True

        Finally, we try the ideal that is all of the ring in different
        rings. ::

            sage: # optional - frobby
            sage: rings = [ZZ['x'], CC['x,y']]
            sage: allOK = True
            sage: for ring in rings:
            ....:     id1 = ring.ideal(1)
            ....:     decom1 = frobby.irreducible_decomposition(id1)
            ....:     allOK = allOK and decom1 == [id1]
            sage: allOK
            True
        """
        frobby_input = self._ideal_to_string(monomial_ideal)
        frobby_output = self('irrdecom', input=frobby_input)
        return self._parse_ideals(frobby_output, monomial_ideal.ring())

    def _parse_ideals(self, string, ring):
        r"""
        This function parses a list of irreducible monomial ideals in 4ti2
        format and constructs them within the passed-in ring.

        INPUT:

        - ``string`` -- the string to be parsed
        - ``ring`` -- the ring within which to construct the irreducible
          monomial ideals within

        OUTPUT:

        A list of the monomial ideals in the order they are listed
        in the string.

        EXAMPLES::

            sage: # optional - frobby
            sage: ring = QQ['x,y,z']
            sage: (x, y, z) = ring.gens()
            sage: string = '2 3\n1 2 3\n0 5 0\n2 3\n1 2 3\n4 5 6'
            sage: frobby._parse_ideals(string, ring)
            [Ideal (x*y^2*z^3, y^5) of Multivariate Polynomial Ring in x, y, z over Rational Field,
            Ideal (x*y^2*z^3, x^4*y^5*z^6) of Multivariate Polynomial Ring in x, y, z over Rational Field]
        """
        lines = string.split('\n')
        if lines[-1] == '':
            lines.pop(-1)
        matrices = []
        while lines:
            if lines[0].split()[1] == 'ring':
                lines.pop(0)
                lines.pop(0)
                matrices.append('1 '+str(ring.ngens())+'\n'+'0 '*ring.ngens()+'\n')
            else:
                nrows = int(lines[0].split()[0])
                nmatrix = lines.pop(0)+'\n'
                for i in range(nrows):
                    nmatrix += lines.pop(0)+'\n'
                matrices.append(nmatrix)

        def to_ideal(exps):
            if len(exps) == 0:
                return ring.zero_ideal()
            gens = [prod([v ** e for v, e in zip(ring.gens(), expo) if e != 0]) for expo in exps]
            return ring.ideal(gens or ring(1))
        return [to_ideal(self._parse_4ti2_matrix(a)) for a in matrices] or [ring.ideal()]

    def _parse_4ti2_matrix(self, string):
        r"""
        Parses a string of a matrix in 4ti2 format into a nested list
        representation.

        INPUT:

        - ``string`` -- the string to be parsed

        OUTPUT:

        A list of rows of the matrix, where each row is represented as
        a list of integers.

        EXAMPLES:

        The format is straight-forward, as this example shows. ::

            sage: # optional - frobby
            sage: string = '2 3\n1 2  3\n  0 5 0'
            sage: parsed_matrix = frobby._parse_4ti2_matrix(string)
            sage: reference_matrix = [[1, 2, 3], [0, 5, 0]]
            sage: parsed_matrix == reference_matrix
            True

        A number of syntax errors lead to exceptions. ::

            sage: string = '1 1\n one' # optional - frobby
            sage: frobby._parse_4ti2_matrix(string) # optional - frobby
            Traceback (most recent call last):
            ...
            RuntimeError: Format error: encountered non-number.
        """
        try:
            ints = [int(_) for _ in string.split()]
        except ValueError:
            raise RuntimeError("Format error: encountered non-number.")
        if len(ints) < 2:
            raise RuntimeError("Format error: " +
                               "matrix dimensions not specified.")

        term_count = ints[0]
        var_count = ints[1]
        ints = ints[2:]

        if term_count * var_count != len(ints):
            raise RuntimeError("Format error: incorrect matrix dimensions.")

        exponents = []
        for i in range(term_count):
            exponents.append(ints[:var_count])
            ints = ints[var_count:]
        return exponents

    def _ideal_to_string(self, monomial_ideal):
        r"""
        This function formats the passed-in monomial ideal in 4ti2 format.

        INPUT:

        - ``monomial_ideal`` -- the monomial ideal to be formatted as a string

        OUTPUT: string in ``4ti2`` format representing the ideal

        EXAMPLES::

            sage: # optional - frobby
            sage: ring = QQ['x,y,z']
            sage: (x, y, z) = ring.gens()
            sage: id = ring.ideal(x ** 2, x * y * z)
            sage: frobby._ideal_to_string(id) == "2 3\n2 0 0\n1 1 1\n"
            True
        """
        # There is no exponent vector that represents zero as a generator, so
        # we take care that the zero ideal gets represented correctly in the
        # 4ti2 format; as an ideal with no generators.
        if monomial_ideal.is_zero():
            gens = []
        else:
            gens = monomial_ideal.gens()
        var_count = monomial_ideal.ring().ngens()
        first_row = str(len(gens)) + ' ' + str(var_count) + '\n'
        rows = [self._monomial_to_string(_) for _ in gens]
        return first_row + "".join(rows)

    def _monomial_to_string(self, monomial):
        r"""
        This function formats the exponent vector of a monomial as a string
        containing a space-delimited list of integers.

        INPUT:

        - ``monomial`` -- the monomial whose exponent vector is to be formatted

        OUTPUT: string representing the exponent vector of monomial

        EXAMPLES::

            sage: # optional - frobby
            sage: ring = QQ['x,y,z']
            sage: (x, y, z) = ring.gens()
            sage: monomial = x * x * z
            sage: frobby._monomial_to_string(monomial) == '2 0 1\n'
            True
        """
        exponents = monomial.exponents()
        if len(exponents) != 1:
            raise RuntimeError(str(monomial) + " is not a monomial.")
        exponents = exponents[0]

        # for a multivariate ring exponents will be an ETuple, while
        # for a univariate ring exponents will be just an int. To get
        # this to work we make the two cases look alike.
        if isinstance(exponents, int):
            exponents = [exponents]
        strings = [str(exponents[var]) for var in range(len(exponents))]
        return ' '.join(strings) + '\n'


# This singleton instance is what should be used outside this file.
frobby = Frobby()
