r"""
Interface to CHomP

CHomP stands for "Computation Homology Program", and is good at
computing homology of simplicial complexes, cubical complexes, and
chain complexes.  It can also compute homomorphisms induced on
homology by maps.  See the CHomP web page http://chomp.rutgers.edu/
for more information.

AUTHOR:

- John H. Palmieri
"""

import re

_have_chomp = {}
def have_chomp(program='homsimpl'):
    """
    Return True if this computer has ``program`` installed.

    The first time it is run, this function caches its result in the
    variable ``_have_chomp`` -- a dictionary indexed by program name
    -- and any subsequent time, it just checks the value of the
    variable.

    This program is used in the routine CHomP.__call__.

    If this computer doesn't have CHomP installed, you may obtain it
    from http://chomp.rutgers.edu/.

    EXAMPLES::

        sage: from sage.interfaces.chomp import have_chomp
        sage: have_chomp() # random -- depends on whether CHomP is installed
        True
        sage: 'homsimpl' in sage.interfaces.chomp._have_chomp
        True
        sage: sage.interfaces.chomp._have_chomp['homsimpl'] == have_chomp()
        True
    """
    global _have_chomp
    if program not in _have_chomp:
        from sage.misc.sage_ostools import have_program
        _have_chomp[program] = have_program(program)
    return _have_chomp[program]


class CHomP:
    r"""
    Interface to the CHomP package.

    :param program: which CHomP program to use
    :type program: string
    :param complex: a simplicial or cubical complex
    :param subcomplex: a subcomplex of ``complex`` or None (the default)
    :param base_ring: ring over which to perform computations -- must be `\ZZ` or `\GF{p}`.
    :type base_ring: ring; optional, default `\ZZ`
    :param generators: if True, also return list of generators
    :type generators: boolean; optional, default False
    :param verbose: if True, print helpful messages as the computation
       progresses
    :type verbose: boolean; optional, default False
    :param extra_opts: options passed directly to ``program``
    :type extra_opts: string
    :return: homology groups as a dictionary indexed by dimension

    The programs ``homsimpl``, ``homcubes``, and ``homchain`` are
    available through this interface.  ``homsimpl`` computes the
    relative or absolute homology groups of simplicial complexes.
    ``homcubes`` computes the relative or absolute homology groups of
    cubical complexes.  ``homchain`` computes the homology groups of
    chain complexes.  For consistency with Sage's other homology
    computations, the answers produced by ``homsimpl`` and
    ``homcubes`` in the absolute case are converted to reduced
    homology.

    Note also that CHomP can only compute over the integers or
    `\GF{p}`.  CHomP is fast enough, though, that if you want
    rational information, you should consider using CHomP with integer
    coefficients, or with mod `p` coefficients for a sufficiently
    large `p`, rather than using Sage's built-in homology algorithms.

    See also the documentation for the functions :func:`homchain`,
    :func:`homcubes`, and :func:`homsimpl` for more examples,
    including illustrations of some of the optional parameters.

    EXAMPLES::

        sage: from sage.interfaces.chomp import CHomP
        sage: T = cubical_complexes.Torus()
        sage: CHomP()('homcubes', T) # optional - CHomP
        {0: 0, 1: Z x Z, 2: Z}

    Relative homology of a segment relative to its endpoints::

        sage: edge = simplicial_complexes.Simplex(1)
        sage: ends = edge.n_skeleton(0)
        sage: CHomP()('homsimpl', edge)  # optional - CHomP
        {0: 0}
        sage: CHomP()('homsimpl', edge, ends)  # optional - CHomP
        {0: 0, 1: Z}

    Homology of a chain complex::

        sage: C = ChainComplex({3: 2 * identity_matrix(ZZ, 2)}, degree=-1)
        sage: CHomP()('homchain', C)  # optional - CHomP
        {2: C2 x C2}
    """
    def __repr__(self):
        """
        String representation.

        EXAMPLES::

            sage: from sage.interfaces.chomp import CHomP
            sage: CHomP()           # indirect doctest
            CHomP interface
        """
        return "CHomP interface"

    def __call__(self, program, complex, subcomplex=None, **kwds):
        """
        Call a CHomP program to compute the homology of a chain
        complex, simplicial complex, or cubical complex.

        See :class:`CHomP` for full documentation.

        EXAMPLES::

            sage: from sage.interfaces.chomp import CHomP
            sage: T = cubical_complexes.Torus()
            sage: CHomP()('homcubes', T) # indirect doctest, optional - CHomP
            {0: 0, 1: Z x Z, 2: Z}
        """
        from sage.misc.temporary_file import tmp_filename
        from sage.topology.cubical_complex import CubicalComplex, cubical_complexes
        from sage.topology.simplicial_complex import SimplicialComplex, Simplex
        from sage.homology.chain_complex import HomologyGroup
        from subprocess import Popen, PIPE
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.modules.all import VectorSpace, vector
        from sage.combinat.free_module import CombinatorialFreeModule

        if not have_chomp(program):
            raise OSError("Program %s not found" % program)

        verbose = kwds.get('verbose', False)
        generators = kwds.get('generators', False)
        extra_opts = kwds.get('extra_opts', '')
        base_ring = kwds.get('base_ring', ZZ)

        if extra_opts:
            extra_opts = extra_opts.split()
        else:
            extra_opts = []

        # type of complex:
        cubical = False
        simplicial = False
        chain = False
        # CHomP seems to have problems with cubical complexes if the
        # first interval in the first cube defining the complex is
        # degenerate.  So replace the complex X with [0,1] x X.
        if isinstance(complex, CubicalComplex):
            cubical = True
            edge = cubical_complexes.Cube(1)
            original_complex = complex
            complex = edge.product(complex)
            if verbose:
                print("Cubical complex")
        elif isinstance(complex, SimplicialComplex):
            simplicial = True
            if verbose:
                print("Simplicial complex")
        else:
            chain = True
            base_ring = kwds.get('base_ring', complex.base_ring())
            if verbose:
                print("Chain complex over %s" % base_ring)

        if base_ring == QQ:
            raise ValueError("CHomP doesn't compute over the rationals, only over Z or F_p.")
        if base_ring.is_prime_field():
            p = base_ring.characteristic()
            extra_opts.append('-p%s' % p)
            mod_p = True
        else:
            mod_p = False

        #
        #    complex
        #
        try:
            data = complex._chomp_repr_()
        except AttributeError:
            raise AttributeError("Complex cannot be converted to use with CHomP.")

        datafile = tmp_filename()
        with open(datafile, 'w') as f:
            f.write(data)

        #
        #    subcomplex
        #
        if subcomplex is None:
            if cubical:
                subcomplex = CubicalComplex([complex.n_cells(0)[0]])
            elif simplicial:
                m = re.search(r'\(([^,]*),', data)
                v = int(m.group(1))
                subcomplex = SimplicialComplex([[v]])
        else:
            # replace subcomplex with [0,1] x subcomplex.
            if cubical:
                subcomplex = edge.product(subcomplex)
        #
        #    generators
        #
        if generators:
            genfile = tmp_filename()
            extra_opts.append('-g%s' % genfile)

        #
        #    call program
        #
        if subcomplex is not None:
            try:
                sub = subcomplex._chomp_repr_()
            except AttributeError:
                raise AttributeError("Subcomplex cannot be converted to use with CHomP.")
            subfile = tmp_filename()
            with open(subfile, 'w') as f:
                f.write(sub)
        else:
            subfile = ''
        if verbose:
            print("Popen called with arguments", end="")
            print([program, datafile, subfile] + extra_opts)
            print("")
            print("CHomP output:")
            print("")
        # output = Popen([program, datafile, subfile, extra_opts],
        cmd = [program, datafile]
        if subfile:
            cmd.append(subfile)
        if extra_opts:
            cmd.extend(extra_opts)
        output = Popen(cmd, stdout=PIPE).communicate()[0]
        if verbose:
            print(output)
            print("End of CHomP output")
            print("")
        if generators:
            with open(genfile, 'r') as f:
                gens = f.read()
            if verbose:
                print("Generators:")
                print(gens)
        #
        #    process output
        #
        if output.find('ERROR') != -1:
            raise RuntimeError('error inside CHomP')
        # output contains substrings of one of the forms
        # "H_1 = Z", "H_1 = Z_2 + Z", "H_1 = Z_2 + Z^2",
        # "H_1 = Z + Z_2 + Z"
        if output.find('trivial') != -1:
            if mod_p:
                return {0: VectorSpace(base_ring, 0)}
            else:
                return {0: HomologyGroup(0, ZZ)}
        d = {}
        h = re.compile("^H_([0-9]*) = (.*)$", re.M)
        tors = re.compile("Z_([0-9]*)")
        #
        #    homology groups
        #
        for m in h.finditer(output):
            if verbose:
                print(m.groups())
            # dim is the dimension of the homology group
            dim = int(m.group(1))
            # hom_str is the right side of the equation "H_n = Z^r + Z_k + ..."
            hom_str = m.group(2)
            # need to read off number of summands and their invariants
            if hom_str.find("0") == 0:
                if mod_p:
                    hom = VectorSpace(base_ring, 0)
                else:
                    hom = HomologyGroup(0, ZZ)
            else:
                rk = 0
                if hom_str.find("^") != -1:
                    rk_srch = re.search(r'\^([0-9]*)\s?', hom_str)
                    rk = int(rk_srch.group(1))
                rk += len(re.findall(r"(Z$)|(Z\s)", hom_str))
                if mod_p:
                    rk = rk if rk != 0 else 1
                    if verbose:
                        print("dimension = %s, rank of homology = %s" % (dim, rk))
                    hom = VectorSpace(base_ring, rk)
                else:
                    n = rk
                    invts = []
                    for t in tors.finditer(hom_str):
                        n += 1
                        invts.append(int(t.group(1)))
                    for i in range(rk):
                        invts.append(0)
                    if verbose:
                        print("dimension = %s, number of factors = %s, invariants = %s" % (dim, n, invts))
                    hom = HomologyGroup(n, ZZ, invts)

            #
            #    generators
            #
            if generators:
                if cubical:
                    g = process_generators_cubical(gens, dim)
                    if verbose:
                        print("raw generators: %s" % g)
                    if g:
                        module = CombinatorialFreeModule(base_ring,
                                                         original_complex.n_cells(dim),
                                                         prefix="",
                                                         bracket=True)
                        basis = module.basis()
                        output = []
                        for x in g:
                            v = module(0)
                            for term in x:
                                v += term[0] * basis[term[1]]
                            output.append(v)
                        g = output
                elif simplicial:
                    g = process_generators_simplicial(gens, dim, complex)
                    if verbose:
                        print("raw generators: %s" % gens)
                    if g:
                        module = CombinatorialFreeModule(base_ring,
                                                         complex.n_cells(dim),
                                                         prefix="",
                                                         bracket=False)
                        basis = module.basis()
                        output = []
                        for x in g:
                            v = module(0)
                            for term in x:
                                if complex._is_numeric():
                                    v += term[0] * basis[term[1]]
                                else:
                                    translate = complex._translation_from_numeric()
                                    simplex = Simplex([translate[a] for a in term[1]])
                                    v += term[0] * basis[simplex]
                            output.append(v)
                        g = output
                elif chain:
                    g = process_generators_chain(gens, dim, base_ring)
                    if verbose:
                        print("raw generators: %s" % gens)
                if g:
                    if not mod_p:
                        # sort generators to match up with corresponding invariant
                        g = [_[1] for _ in sorted(zip(invts, g), key=lambda x: x[0])]
                    d[dim] = (hom, g)
                else:
                    d[dim] = hom
            else:
                d[dim] = hom

        if chain:
            new_d = {}
            diff = complex.differential()
            if len(diff) == 0:
                return {}
            bottom = min(diff)
            top = max(diff)
            for dim in d:
                if complex._degree_of_differential == -1:  # chain complex
                    new_dim = bottom + dim
                else: # cochain complex
                    new_dim = top - dim
                if isinstance(d[dim], tuple):
                    # generators included.
                    group = d[dim][0]
                    gens = d[dim][1]
                    new_gens = []
                    dimension = complex.differential(new_dim).ncols()
                    # make sure that each vector is embedded in the
                    # correct ambient space: pad with a zero if
                    # necessary.
                    for v in gens:
                        v_dict = v.dict()
                        if dimension - 1 not in v.dict():
                            v_dict[dimension - 1] = 0
                            new_gens.append(vector(base_ring, v_dict))
                        else:
                            new_gens.append(v)
                    new_d[new_dim] = (group, new_gens)
                else:
                    new_d[new_dim] = d[dim]
            d = new_d
        return d

    def help(self, program):
        """
        Print a help message for ``program``, a program from the CHomP suite.

        :param program: which CHomP program to use
        :type program: string
        :return: nothing -- just print a message

        EXAMPLES::

            sage: from sage.interfaces.chomp import CHomP
            sage: CHomP().help('homcubes')   # optional - CHomP
            HOMCUBES, ver. ... Copyright (C) ... by Pawel Pilarczyk...
        """
        from subprocess import Popen, PIPE
        print(Popen([program, '-h'], stdout=PIPE).communicate()[0])


def homsimpl(complex=None, subcomplex=None, **kwds):
    r"""
    Compute the homology of a simplicial complex using the CHomP
    program ``homsimpl``.  If the argument ``subcomplex`` is present,
    compute homology of ``complex`` relative to ``subcomplex``.

    :param complex: a simplicial complex
    :param subcomplex: a subcomplex of ``complex`` or None (the default)
    :param base_ring: ring over which to perform computations -- must be `\ZZ` or `\GF{p}`.
    :type base_ring: ring; optional, default `\ZZ`
    :param generators: if True, also return list of generators
    :type generators: boolean; optional, default False
    :param verbose: if True, print helpful messages as the computation
       progresses
    :type verbose: boolean; optional, default False
    :param help: if True, just print a help message and exit
    :type help: boolean; optional, default False
    :param extra_opts: options passed directly to ``program``
    :type extra_opts: string
    :return: homology groups as a dictionary indexed by dimension

    EXAMPLES::

        sage: from sage.interfaces.chomp import homsimpl
        sage: T = simplicial_complexes.Torus()
        sage: M8 = simplicial_complexes.MooreSpace(8)
        sage: M4 = simplicial_complexes.MooreSpace(4)
        sage: X = T.disjoint_union(T).disjoint_union(T).disjoint_union(M8).disjoint_union(M4)
        sage: homsimpl(X)[1]  # optional - CHomP
        Z^6 x C4 x C8

    Relative homology::

        sage: S = simplicial_complexes.Simplex(3)
        sage: bdry = S.n_skeleton(2)
        sage: homsimpl(S, bdry)[3]   # optional - CHomP
        Z

    Generators: these are given as a list after the homology group.
    Each generator is specified as a linear combination of simplices::

        sage: homsimpl(S, bdry, generators=True)[3]   # optional - CHomP
        (Z, [(0, 1, 2, 3)])

        sage: homsimpl(simplicial_complexes.Sphere(1), generators=True)   # optional - CHomP
        {0: 0, 1: (Z, [(0, 1) - (0, 2) + (1, 2)])}

    TESTS:

    Generators for a simplicial complex whose vertices are not integers::

        sage: S1 = simplicial_complexes.Sphere(1)
        sage: homsimpl(S1.join(S1), generators=True, base_ring=GF(2))[3][1]  # optional - CHomP
        [('L0', 'L1', 'R0', 'R1') + ('L0', 'L1', 'R0', 'R2') + ('L0', 'L1', 'R1', 'R2') + ('L0', 'L2', 'R0', 'R1') + ('L0', 'L2', 'R0', 'R2') + ('L0', 'L2', 'R1', 'R2') + ('L1', 'L2', 'R0', 'R1') + ('L1', 'L2', 'R0', 'R2') + ('L1', 'L2', 'R1', 'R2')]
    """
    from sage.topology.simplicial_complex import SimplicialComplex
    help = kwds.get('help', False)
    if help:
        return CHomP().help('homsimpl')
    # Check types here, because if you pass homsimpl a cubical
    # complex, it tries to compute its homology as if it were a
    # simplicial complex and gets terribly wrong answers.
    if (isinstance(complex, SimplicialComplex)
        and (subcomplex is None or isinstance(subcomplex, SimplicialComplex))):
        return CHomP()('homsimpl', complex, subcomplex=subcomplex, **kwds)
    else:
        raise TypeError("Complex and/or subcomplex are not simplicial complexes.")

def homcubes(complex=None, subcomplex=None, **kwds):
    r"""
    Compute the homology of a cubical complex using the CHomP program
    ``homcubes``.  If the argument ``subcomplex`` is present, compute
    homology of ``complex`` relative to ``subcomplex``.

    :param complex: a cubical complex
    :param subcomplex: a subcomplex of ``complex`` or None (the default)
    :param base_ring: ring over which to perform computations -- must be `\ZZ` or `\GF{p}`.
    :type base_ring: ring; optional, default `\ZZ`
    :param generators: if True, also return list of generators
    :type generators: boolean; optional, default False
    :param verbose: if True, print helpful messages as the computation progresses
    :type verbose: boolean; optional, default False
    :param help: if True, just print a help message and exit
    :type help: boolean; optional, default False
    :param extra_opts: options passed directly to ``homcubes``
    :type extra_opts: string
    :return: homology groups as a dictionary indexed by dimension

    EXAMPLES::

        sage: from sage.interfaces.chomp import homcubes
        sage: S = cubical_complexes.Sphere(3)
        sage: homcubes(S)[3]   # optional - CHomP
        Z

    Relative homology::

        sage: C3 = cubical_complexes.Cube(3)
        sage: bdry = C3.n_skeleton(2)
        sage: homcubes(C3, bdry)   # optional - CHomP
        {0: 0, 1: 0, 2: 0, 3: Z}

    Generators: these are given as a list after the homology group.
    Each generator is specified as a linear combination of cubes::

        sage: homcubes(cubical_complexes.Sphere(1), generators=True, base_ring=GF(2))[1][1]   # optional - CHomP
        [[[1,1] x [0,1]] + [[0,1] x [1,1]] + [[0,1] x [0,0]] + [[0,0] x [0,1]]]
    """
    from sage.topology.cubical_complex import CubicalComplex
    help = kwds.get('help', False)
    if help:
        return CHomP().help('homcubes')
    # Type-checking is here for the same reason as in homsimpl above.
    if (isinstance(complex, CubicalComplex)
        and (subcomplex is None or isinstance(subcomplex, CubicalComplex))):
        return CHomP()('homcubes', complex, subcomplex=subcomplex, **kwds)
    else:
        raise TypeError("Complex and/or subcomplex are not cubical complexes.")


def homchain(complex=None, **kwds):
    r"""
    Compute the homology of a chain complex using the CHomP program
    ``homchain``.

    :param complex: a chain complex
    :param generators: if True, also return list of generators
    :type generators: boolean; optional, default False
    :param verbose: if True, print helpful messages as the computation progresses
    :type verbose: boolean; optional, default False
    :param help: if True, just print a help message and exit
    :type help: boolean; optional, default False
    :param extra_opts: options passed directly to ``homchain``
    :type extra_opts: string
    :return: homology groups as a dictionary indexed by dimension

    EXAMPLES::

        sage: from sage.interfaces.chomp import homchain
        sage: C = cubical_complexes.Sphere(3).chain_complex()
        sage: homchain(C)[3]   # optional - CHomP
        Z

    Generators: these are given as a list after the homology group.
    Each generator is specified as a cycle, an element in the
    appropriate free module over the base ring::

        sage: C2 = delta_complexes.Sphere(2).chain_complex()
        sage: homchain(C2, generators=True)[2]  # optional - CHomP
        (Z, [(1, -1)])
        sage: homchain(C2, generators=True, base_ring=GF(2))[2]  # optional - CHomP
        (Vector space of dimension 1 over Finite Field of size 2, [(1, 1)])

    TESTS:

    Chain complexes concentrated in negative dimensions, cochain complexes, etc.::

        sage: C = ChainComplex({-5: 4 * identity_matrix(ZZ, 2)}, degree=-1)
        sage: homchain(C)   # optional - CHomP
        {-6: C4 x C4}
        sage: C = ChainComplex({-5: 4 * identity_matrix(ZZ, 2)}, degree=1)
        sage: homchain(C, generators=True)   # optional - CHomP
        {-4: (C4 x C4, [(1, 0), (0, 1)])}
    """
    from sage.homology.chain_complex import ChainComplex_class
    help = kwds.get('help', False)
    if help:
        return CHomP().help('homchain')
    # Type-checking just in case.
    if isinstance(complex, ChainComplex_class):
        return CHomP()('homchain', complex, **kwds)
    else:
        raise TypeError("Complex is not a chain complex.")


def process_generators_cubical(gen_string, dim):
    r"""
    Process CHomP generator information for cubical complexes.

    :param gen_string: generator output from CHomP
    :type gen_string: string
    :param dim: dimension in which to find generators
    :type dim: integer
    :return: list of generators in each dimension, as described below

    ``gen_string`` has the form ::

        The 2 generators of H_1 follow:
        generator 1
        -1 * [(0,0,0,0,0)(0,0,0,0,1)]
        1 * [(0,0,0,0,0)(0,0,1,0,0)]
        ...
        generator 2
        -1 * [(0,1,0,1,1)(1,1,0,1,1)]
        -1 * [(0,1,0,0,1)(0,1,0,1,1)]
        ...

    Each line consists of a coefficient multiplied by a cube; the cube
    is specified by its "bottom left" and "upper right" corners.

    For technical reasons, we remove the first coordinate of each
    tuple, and using regular expressions, the remaining parts get
    converted from a string to a pair ``(coefficient, Cube)``, with
    the cube represented as a product of tuples.  For example, the
    first line in "generator 1" gets turned into ::

        (-1, [0,0] x [0,0] x [0,0] x [0,1])

    representing an element in the free abelian group with basis given
    by cubes.  Each generator is a list of such pairs, representing
    the sum of such elements.  These are reassembled in
    :meth:`CHomP.__call__` to actual elements in the free module
    generated by the cubes of the cubical complex in the appropriate
    dimension.

    Therefore the return value is a list of lists of pairs, one list
    of pairs for each generator.

    EXAMPLES::

        sage: from sage.interfaces.chomp import process_generators_cubical
        sage: s = "The 2 generators of H_1 follow:\ngenerator 1:\n-1 * [(0,0,0,0,0)(0,0,0,0,1)]\n1 * [(0,0,0,0,0)(0,0,1,0,0)]"
        sage: process_generators_cubical(s, 1)
        [[(-1, [0,0] x [0,0] x [0,0] x [0,1]), (1, [0,0] x [0,1] x [0,0] x [0,0])]]
        sage: len(process_generators_cubical(s, 1))  # only one generator
        1
    """
    from sage.topology.cubical_complex import Cube
    # each dim in gen_string starts with "The generator for
    # H_3 follows:".  So search for "T" to find the
    # end of the current list of generators.
    g_srch = re.compile(r'H_%s\sfollow[s]?:\n([^T]*)(?:T|$)' % dim)
    g = g_srch.search(gen_string)
    output = []
    cubes_list = []
    if g:
        g = g.group(1)
    if g:
        # project g to one end of the cylinder [0,1] x complex:
        #
        # drop the first coordinate and eliminate duplicates, at least
        # in positive dimensions, drop any line containing a
        # degenerate cube
        g = re.sub(r'\([01],', '(', g)
        if dim > 0:
            lines = g.splitlines()
            newlines = []
            for l in lines:
                x = re.findall(r'(\([0-9,]*\))', l)
                if x:
                    left, right = x
                    left = [int(a) for a in left.strip('()').split(',')]
                    right = [int(a) for a in right.strip('()').split(',')]
                    if sum([xx - yy for (xx, yy) in zip(right, left)]) == dim:
                        newlines.append(l)
                else:  # line like "generator 2"
                    newlines.append(l)
            g = newlines
        cubes_list = []
        for l in g:
            cubes = re.search(r'([+-]?)\s?([0-9]+)?\s?[*]?\s?\[\(([-0-9,]+)\)\(([-0-9,]+)\)\]', l)
            if cubes:
                if cubes.group(1) and re.search("-", cubes.group(1)):
                    sign = -1
                else:
                    sign = 1
                if cubes.group(2) and len(cubes.group(2)) > 0:
                    coeff = sign * int(cubes.group(2))
                else:
                    coeff = sign * 1
                cube1 = [int(a) for a in cubes.group(3).split(',')]
                cube2 = [int(a) for a in cubes.group(4).split(',')]
                cube = Cube(zip(cube1, cube2))
                cubes_list.append((coeff, cube))
            else:
                if cubes_list:
                    output.append(cubes_list)
                    cubes_list = []
        if cubes_list:
            output.append(cubes_list)
        return output
    else:
        return None

def process_generators_simplicial(gen_string, dim, complex):
    r"""
    Process CHomP generator information for simplicial complexes.

    :param gen_string: generator output from CHomP
    :type gen_string: string
    :param dim: dimension in which to find generators
    :type dim: integer
    :param complex: simplicial complex under consideration
    :return: list of generators in each dimension, as described below

    ``gen_string`` has the form ::

        The 2 generators of H_1 follow:
        generator 1
        -1 * (1,6)
        1 * (1,4)
        ...
        generator 2
        -1 * (1,6)
        1 * (1,4)
        ...

    where each line contains a coefficient and a simplex.  Each line
    is converted, using regular expressions, from a string to a pair
    ``(coefficient, Simplex)``, like ::

        (-1, (1,6))

    representing an element in the free abelian group with basis given
    by simplices.  Each generator is a list of such pairs,
    representing the sum of such elements.  These are reassembled in
    :meth:`CHomP.__call__` to actual elements in the free module
    generated by the simplices of the simplicial complex in the
    appropriate dimension.

    Therefore the return value is a list of lists of pairs, one list
    of pairs for each generator.

    EXAMPLES::

        sage: from sage.interfaces.chomp import process_generators_simplicial
        sage: s = "The 2 generators of H_1 follow:\ngenerator 1:\n-1 * (1,6)\n1 * (1,4)"
        sage: process_generators_simplicial(s, 1, simplicial_complexes.Torus())
        [[(-1, (1, 6)), (1, (1, 4))]]
    """
    from sage.topology.simplicial_complex import Simplex
    # each dim in gen_string starts with "The generator for H_3
    # follows:".  So search for "T" to find the end of the current
    # list of generators.
    g_srch = re.compile(r'H_%s\sfollow[s]?:\n([^T]*)(?:T|$)' % dim)
    g = g_srch.search(gen_string)
    output = []
    simplex_list = []
    if g:
        g = g.group(1)
    if g:
        lines = g.splitlines()
        for l in lines:
            simplex = re.search(r'([+-]?)\s?([0-9]+)?\s?[*]?\s?(\([0-9,]*\))', l)
            if simplex:
                if simplex.group(1) and re.search("-", simplex.group(1)):
                    sign = -1
                else:
                    sign = 1
                if simplex.group(2) and len(simplex.group(2)) > 0:
                    coeff = sign * int(simplex.group(2))
                else:
                    coeff = sign * 1
                simplex = Simplex([int(a) for a in simplex.group(3).strip('()').split(',')])
                simplex_list.append((coeff, simplex))
            else:
                if simplex_list:
                    output.append(simplex_list)
                    simplex_list = []
        if simplex_list:
            output.append(simplex_list)
        return output
    else:
        return None

def process_generators_chain(gen_string, dim, base_ring=None):
    r"""
    Process CHomP generator information for simplicial complexes.

    :param gen_string: generator output from CHomP
    :type gen_string: string
    :param dim: dimension in which to find generators
    :type dim: integer
    :param base_ring: base ring over which to do the computations
    :type base_ring: optional, default ZZ
    :return: list of generators in each dimension, as described below

    ``gen_string`` has the form ::

        [H_0]
        a1

        [H_1]
        a2
        a3

        [H_2]
        a1 - a2

    For each homology group, each line lists a homology generator as a
    linear combination of generators ``ai`` of the group of chains in
    the appropriate dimension.  The elements ``ai`` are indexed
    starting with `i=1`.  Each generator is converted, using regular
    expressions, from a string to a vector (an element in the free
    module over ``base_ring``), with ``ai`` representing the unit
    vector in coordinate `i-1`.  For example, the string ``a1 - a2``
    gets converted to the vector ``(1, -1)``.

    Therefore the return value is a list of vectors.

    EXAMPLES::

        sage: from sage.interfaces.chomp import process_generators_chain
        sage: s = "[H_0]\na1\n\n[H_1]\na2\na3\n"
        sage: process_generators_chain(s, 1)
        [(0, 1), (0, 0, 1)]
        sage: s = "[H_0]\na1\n\n[H_1]\n5 * a2 - a1\na3\n"
        sage: process_generators_chain(s, 1, base_ring=ZZ)
        [(-1, 5), (0, 0, 1)]
        sage: process_generators_chain(s, 1, base_ring=GF(2))
        [(1, 1), (0, 0, 1)]
    """
    from sage.modules.free_module_element import vector
    from sage.rings.integer_ring import ZZ
    if base_ring is None:
        base_ring = ZZ
    # each dim in gens starts with a string like
    # "[H_3]".  So search for "[" to find the end of
    # the current list of generators.
    g_srch = re.compile(r'\[H_%s\]\n([^]]*)(?:\[|$)' % dim)
    g = g_srch.search(gen_string)
    if g:
        g = g.group(1)
    if g:
        # each line in the string g is a linear
        # combination of things like "a2", "a31", etc.
        # indexing on the a's starts at 1.
        lines = g.splitlines()
        new_gens = []
        for l in lines:
            gen = re.compile(r"([+-]?)\s?([0-9]+)?\s?[*]?\s?a([0-9]*)(?:\s|$)")
            v = {}
            for term in gen.finditer(l):
                if term.group(1) and re.search("-", term.group(1)):
                    sign = -1
                else:
                    sign = 1
                if term.group(2) and len(term.group(2)) > 0:
                    coeff = sign * int(term.group(2))
                else:
                    coeff = sign * 1
                idx = int(term.group(3))
                v[idx-1] = coeff
            if v:
                new_gens.append(vector(base_ring, v))
        g = new_gens
    return g
