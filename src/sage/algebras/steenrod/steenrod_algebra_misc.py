# sage_setup: distribution = sagemath-combinat
"""
Miscellaneous functions for the Steenrod algebra and its elements

AUTHORS:

- John H. Palmieri (2008-07-30): initial version (as the file
  steenrod_algebra_element.py)

- John H. Palmieri (2010-06-30): initial version of steenrod_misc.py.
  Implemented profile functions.  Moved most of the methods for
  elements to the ``Element`` subclass of
  :class:`sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra_generic`.

The main functions here are

- :func:`get_basis_name`.  This function takes a string as input and
  attempts to interpret it as the name of a basis for the Steenrod
  algebra; it returns the canonical name attached to that basis.  This
  allows for the use of synonyms when defining bases, while the
  resulting algebras will be identical.

- :func:`normalize_profile`.  This function returns the canonical (and
  hashable) description of any profile function.  See
  :mod:`sage.algebras.steenrod.steenrod_algebra` and
  :func:`SteenrodAlgebra <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra>`
  for information on profile functions.

- functions named ``*_mono_to_string`` where ``*`` is a basis name
  (:func:`milnor_mono_to_string`, etc.).  These convert tuples
  representing basis elements to strings, for _repr_ and _latex_
  methods.
"""

# ****************************************************************************
#  Copyright (C) 2008-2010 John H. Palmieri <palmieri@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
# ****************************************************************************

######################################################
# basis names

_steenrod_milnor_basis_names = ['milnor']
_steenrod_serre_cartan_basis_names = ['serre_cartan', 'serre-cartan', 'sc',
                                      'adem', 'admissible']


def get_basis_name(basis, p, generic=None):
    """
    Return canonical basis named by string basis at the prime p.

    INPUT:

    - ``basis`` -- string

    - ``p`` -- positive prime number

    - ``generic`` -- boolean (default: 'None')

    OUTPUT: ``basis_name`` -- string

    Specify the names of the implemented bases.  The input is
    converted to lower-case, then processed to return the canonical
    name for the basis.

    For the Milnor and Serre-Cartan bases, use the list of synonyms
    defined by the variables :data:`_steenrod_milnor_basis_names` and
    :data:`_steenrod_serre_cartan_basis_names`.  Their canonical names
    are 'milnor' and 'serre-cartan', respectively.

    For the other bases, use pattern-matching rather than a list of
    synonyms:

    - Search for 'wood' and 'y' or 'wood' and 'z' to get the Wood
      bases.  Canonical names 'woody', 'woodz'.

    - Search for 'arnon' and 'c' for the Arnon C basis.  Canonical
      name: 'arnonc'.

    - Search for 'arnon' (and no 'c') for the Arnon A basis.  Also see
      if 'long' is present, for the long form of the basis.  Canonical
      names: 'arnona', 'arnona_long'.

    - Search for 'wall' for the Wall basis. Also see if 'long' is
      present.  Canonical names: 'wall', 'wall_long'.

    - Search for 'pst' for P^s_t bases, then search for the order
      type: 'rlex', 'llex', 'deg', 'revz'.  Canonical names:
      'pst_rlex', 'pst_llex', 'pst_deg', 'pst_revz'.

    - For commutator types, search for 'comm', an order type, and also
      check to see if 'long' is present.  Canonical names:
      'comm_rlex', 'comm_llex', 'comm_deg', 'comm_revz',
      'comm_rlex_long', 'comm_llex_long', 'comm_deg_long',
      'comm_revz_long'.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import get_basis_name
        sage: get_basis_name('adem', 2)
        'serre-cartan'
        sage: get_basis_name('milnor', 2)
        'milnor'
        sage: get_basis_name('MiLNoR', 5)
        'milnor'
        sage: get_basis_name('pst-llex', 2)
        'pst_llex'
        sage: get_basis_name('wood_abcdedfg_y', 2)
        'woody'
        sage: get_basis_name('wood', 2)
        Traceback (most recent call last):
        ...
        ValueError: wood is not a recognized basis at the prime 2
        sage: get_basis_name('arnon--hello--long', 2)
        'arnona_long'
        sage: get_basis_name('arnona_long', p=5)
        Traceback (most recent call last):
        ...
        ValueError: arnona_long is not a recognized basis at the prime 5
        sage: get_basis_name('NOT_A_BASIS', 2)
        Traceback (most recent call last):
        ...
        ValueError: not_a_basis is not a recognized basis at the prime 2
        sage: get_basis_name('woody', 2, generic=True)
        Traceback (most recent call last):
        ...
        ValueError: woody is not a recognized basis for the generic Steenrod algebra at the prime 2
    """
    if generic is None:
        generic = p != 2
    basis = basis.lower()
    if basis in _steenrod_milnor_basis_names:
        result = 'milnor'
    elif basis in _steenrod_serre_cartan_basis_names:
        result = 'serre-cartan'
    elif basis.find('pst') >= 0:
        if basis.find('rlex') >= 0:
            result = 'pst_rlex'
        elif basis.find('llex') >= 0:
            result = 'pst_llex'
        elif basis.find('deg') >= 0:
            result = 'pst_deg'
        elif basis.find('revz') >= 0:
            result = 'pst_revz'
        else:
            result = 'pst_revz'
    elif basis.find('comm') >= 0:
        if basis.find('rlex') >= 0:
            result = 'comm_rlex'
        elif basis.find('llex') >= 0:
            result = 'comm_llex'
        elif basis.find('deg') >= 0:
            result = 'comm_deg'
        elif basis.find('revz') >= 0:
            result = 'comm_revz'
        else:
            result = 'comm_revz'
        if basis.find('long') >= 0:
            result = result + '_long'
    elif not generic and basis.find('wood') >= 0:
        if basis.find('y') >= 0:
            result = 'woody'
        elif basis.find('z') >= 0:
            result = 'woodz'
        else:
            raise ValueError("%s is not a recognized basis at the prime %s" % (basis, p))
    elif not generic and basis.find('arnon') >= 0:
        if basis.find('c') >= 0:
            result = 'arnonc'
        else:
            result = 'arnona'
            if basis.find('long') >= 0:
                result = result + '_long'
    elif not generic and basis.find('wall') >= 0:
        result = 'wall'
        if basis.find('long') >= 0:
            result = result + '_long'
    else:
        gencase = " for the generic Steenrod algebra" if p == 2 and generic else ""
        raise ValueError("%s is not a recognized basis%s at the prime %s" % (basis, gencase, p))
    return result

######################################################
# profile functions


def is_valid_profile(profile, truncation_type, p=2, generic=None):
    r"""
    Return ``True`` if ``profile``, together with ``truncation_type``, is a valid
    profile at the prime `p`.

    INPUT:

    - ``profile`` -- when `p=2`, a tuple or list of numbers; when `p`
      is odd, a pair of such lists

    - ``truncation_type`` -- either 0 or `\infty`

    - ``p`` -- prime number (default: 2)

    - ``generic`` -- boolean (default: ``None``)

    OUTPUT: ``True`` if the profile function is valid, ``False`` otherwise

    See the documentation for :mod:`sage.algebras.steenrod.steenrod_algebra`
    for descriptions of profile functions and how they correspond to
    sub-Hopf algebras of the Steenrod algebra.  Briefly: at the prime
    2, a profile function `e` is valid if it satisfies the condition

    - `e(r) \geq \min( e(r-i) - i, e(i))` for all `0 < i < r`.

    At odd primes, a pair of profile functions `e` and `k` are valid
    if they satisfy

    - `e(r) \geq \min( e(r-i) - i, e(i))` for all `0 < i < r`.

    - if `k(i+j) = 1`, then either `e(i) \leq j` or `k(j) = 1` for all
      `i \geq 1`, `j \geq 0`.

    In this function, profile functions are lists or tuples, and
    ``truncation_type`` is appended as the last element of the list
    `e` before testing.

    EXAMPLES:

    `p=2`::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import is_valid_profile
        sage: is_valid_profile([3,2,1], 0)
        True
        sage: is_valid_profile([3,2,1], Infinity)
        True
        sage: is_valid_profile([1,2,3], 0)
        False
        sage: is_valid_profile([6,2,0], Infinity)
        False
        sage: is_valid_profile([0,3], 0)
        False
        sage: is_valid_profile([0,0,4], 0)
        False
        sage: is_valid_profile([0,0,0,4,0], 0)
        True

    Odd primes::

        sage: is_valid_profile(([0,0,0], [2,1,1,1,2,2]), 0, p=3)
        True
        sage: is_valid_profile(([1], [2,2]), 0, p=3)
        True
        sage: is_valid_profile(([1], [2]), 0, p=7)
        False
        sage: is_valid_profile(([1,2,1], []), 0, p=7)
        True
        sage: is_valid_profile(([0,0,0], [2,1,1,1,2,2]), 0, p=2, generic=True)
        True
    """
    from sage.rings.infinity import Infinity
    if generic is None:
        generic = p != 2
    if not generic:
        pro = list(profile) + [truncation_type]*len(profile)
        r = 0
        for pro_r in pro:
            r += 1 # index of pro_r
            if pro_r < Infinity:
                for i in range(1,r):
                    if pro_r < min(pro[r-i-1] - i, pro[i-1]):
                        return False
    else:
        # p odd:
        e = list(profile[0]) + [truncation_type]*len(profile[0])
        k = list(profile[1])
        if not set(k).issubset({1,2}):
            return False
        if truncation_type > 0:
            k = k + [2]
        else:
            k = k + [1]*len(profile[0])
        if len(k) > len(e):
            e = e + [truncation_type] * (len(k) - len(e))
        r = 0
        for e_r in e:
            r += 1 # index of e_r
            if e_r < Infinity:
                for i in range(1,r):
                    if e_r < min(e[r-i-1] - i, e[i-1]):
                        return False
        r = -1
        for k_r in k:
            r += 1 # index of k_r
            if k_r == 1:
                for j in range(r):
                    i = r-j
                    if e[i-1] > j and k[j] == 2:
                        return False
    return True


def normalize_profile(profile, precision=None, truncation_type='auto', p=2, generic=None):
    r"""
    Given a profile function and related data, return it in a standard form,
    suitable for hashing and caching as data defining a sub-Hopf
    algebra of the Steenrod algebra.

    INPUT:

    - ``profile`` -- a profile function in form specified below
    - ``precision`` -- integer or ``None`` (default: ``None``)
    - ``truncation_type`` -- 0 or `\infty` or 'auto' (default: 'auto')
    - ``p`` -- prime (default: 2)
    - ``generic`` -- boolean (default: ``None``)

    OUTPUT:

    a triple ``profile, precision, truncation_type``, in
    standard form as described below.

    The "standard form" is as follows: ``profile`` should be a tuple
    of integers (or `\infty`) with no trailing zeroes when `p=2`, or a
    pair of such when `p` is odd or `generic` is ``True``.  ``precision``
    should be a positive integer.  ``truncation_type`` should be 0 or `\infty`.
    Furthermore, this must be a valid profile, as determined by the
    function :func:`is_valid_profile`.  See also the documentation for
    the module :mod:`sage.algebras.steenrod.steenrod_algebra` for information
    about profile functions.

    For the inputs: when `p=2`, ``profile`` should be a valid profile
    function, and it may be entered in any of the following forms:

    - a list or tuple, e.g., ``[3,2,1,1]``
    - a function from positive integers to nonnegative integers (and
      `\infty`), e.g., ``lambda n: n+2``.  This corresponds to the
      list ``[3, 4, 5, ...]``.
    - ``None`` or ``Infinity`` -- use this for the profile function for
      the whole Steenrod algebra.  This corresponds to the list
      ``[Infinity, Infinity, Infinity, ...]``

    To make this hashable, it gets turned into a tuple.  In the first
    case it is clear how to do this; also in this case, ``precision``
    is set to be one more than the length of this tuple.  In the
    second case, construct a tuple of length one less than
    ``precision`` (default: 100).  In the last case, the empty
    tuple is returned and ``precision`` is set to 1.

    Once a sub-Hopf algebra of the Steenrod algebra has been defined
    using such a profile function, if the code requires any remaining
    terms (say, terms after the 100th), then they are given by
    ``truncation_type`` if that is 0 or `\infty`.  If
    ``truncation_type`` is 'auto', then in the case of a tuple, it
    gets set to 0, while for the other cases it gets set to `\infty`.

    See the examples below.

    When `p` is odd, ``profile`` is a pair of "functions", so it may
    have the following forms:

    - a pair of lists or tuples, the second of which takes values in
      the set `\{1,2\}`, e.g., ``([3,2,1,1], [1,1,2,2,1])``.

    - a pair of functions, one (called `e`) from positive integers to
      nonnegative integers (and `\infty`), one (called `k`) from
      nonnegative integers to the set `\{1,2\}`, e.g.,
      ``(lambda n: n+2, lambda n: 1)``.  This corresponds to the
      pair ``([3, 4, 5, ...], [1, 1, 1, ...])``.

    - ``None`` or ``Infinity`` -- use this for the profile function for
      the whole Steenrod algebra.  This corresponds to the pair
      ``([Infinity, Infinity, Infinity, ...], [2, 2, 2, ...])``.

    You can also mix and match the first two, passing a pair with
    first entry a list and second entry a function, for instance.  The
    values of ``precision`` and ``truncation_type`` are determined by
    the first entry.

    EXAMPLES:

    `p=2`::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import normalize_profile
        sage: normalize_profile([1,2,1,0,0])
        ((1, 2, 1), 0)

    The full mod 2 Steenrod algebra::

        sage: normalize_profile(Infinity)
        ((), +Infinity)
        sage: normalize_profile(None)
        ((), +Infinity)
        sage: normalize_profile(lambda n: Infinity)
        ((), +Infinity)

    The ``precision`` argument has no effect when the first argument
    is a list or tuple::

        sage: normalize_profile([1,2,1,0,0], precision=12)
        ((1, 2, 1), 0)

    If the first argument is a function, then construct a list of
    length one less than ``precision``, by plugging in the numbers 1,
    2, ..., ``precision`` - 1::

        sage: normalize_profile(lambda n: 4-n, precision=4)
        ((3, 2, 1), +Infinity)
        sage: normalize_profile(lambda n: 4-n, precision=4, truncation_type=0)
        ((3, 2, 1), 0)

    Negative numbers in profile functions are turned into zeroes::

        sage: normalize_profile(lambda n: 4-n, precision=6)
        ((3, 2, 1, 0, 0), +Infinity)

    If it doesn't give a valid profile, an error is raised::

        sage: normalize_profile(lambda n: 3, precision=4, truncation_type=0)
        Traceback (most recent call last):
        ...
        ValueError: Invalid profile
        sage: normalize_profile(lambda n: 3, precision=4, truncation_type = Infinity)
        ((3, 3, 3), +Infinity)

    When `p` is odd, the behavior is similar::

        sage: normalize_profile(([2,1], [2,2,2]), p=13)
        (((2, 1), (2, 2, 2)), 0)

    The full mod `p` Steenrod algebra::

        sage: normalize_profile(None, p=7)
        (((), ()), +Infinity)
        sage: normalize_profile(Infinity, p=11)
        (((), ()), +Infinity)
        sage: normalize_profile((lambda n: Infinity, lambda n: 2), p=17)
        (((), ()), +Infinity)

    Note that as at the prime 2, the ``precision`` argument has no
    effect on a list or tuple in either entry of ``profile``.  If
    ``truncation_type`` is 'auto', then it gets converted to either
    ``0`` or ``+Infinity`` depending on the *first* entry of
    ``profile``::

        sage: normalize_profile(([2,1], [2,2,2]), precision=84, p=13)
        (((2, 1), (2, 2, 2)), 0)
        sage: normalize_profile((lambda n: 0, lambda n: 2), precision=4, p=11)
        (((0, 0, 0), ()), +Infinity)
        sage: normalize_profile((lambda n: 0, (1,1,1,1,1,1,1)), precision=4, p=11)
        (((0, 0, 0), (1, 1, 1, 1, 1, 1, 1)), +Infinity)
        sage: normalize_profile(((4,3,2,1), lambda n: 2), precision=6, p=11)
        (((4, 3, 2, 1), (2, 2, 2, 2, 2)), 0)
        sage: normalize_profile(((4,3,2,1), lambda n: 1), precision=3, p=11, truncation_type=Infinity)
        (((4, 3, 2, 1), (1, 1)), +Infinity)

    As at the prime 2, negative numbers in the first component are
    converted to zeroes.  Numbers in the second component must be
    either 1 and 2, or else an error is raised::

        sage: normalize_profile((lambda n: -n, lambda n: 1), precision=4, p=11)
        (((0, 0, 0), (1, 1, 1)), +Infinity)
        sage: normalize_profile([[0,0,0], [1,2,3,2,1]], p=11)
        Traceback (most recent call last):
        ...
        ValueError: Invalid profile
    """
    from sage.rings.infinity import Infinity
    if truncation_type == 'zero':
        truncation_type = 0
    if truncation_type == 'infinity':
        truncation_type = Infinity
    if generic is None:
        generic = p != 2
    if not generic:
        if profile is None or profile == Infinity:
            # no specified profile or infinite profile: return profile
            # for the entire Steenrod algebra
            new_profile = ()
            truncation_type = Infinity
        elif isinstance(profile, (list, tuple)):
            # profile is a list or tuple: use it as is.  if
            # truncation_type not specified, set it to 'zero'. remove
            # trailing zeroes if truncation_type is 'auto' or 'zero'.
            if truncation_type == 'auto':
                truncation_type = 0
            # remove trailing zeroes or Infinitys
            while profile and profile[-1] == truncation_type:
                profile = profile[:-1]
            new_profile = tuple(profile)
        elif callable(profile):
            # profile is a function: turn it into a tuple.  if
            # truncation_type not specified, set it to 'infinity' if
            # the function is ever infinite; otherwise set it to
            # 0.  remove trailing zeroes if truncation_type is
            # 0, trailing Infinitys if truncation_type is oo.
            if precision is None:
                precision = 100
            if truncation_type == 'auto':
                truncation_type = Infinity
            new_profile = [max(0, profile(i)) for i in range(1, precision)]
            # remove trailing zeroes or Infinitys:
            while new_profile and new_profile[-1] == truncation_type:
                del new_profile[-1]
            new_profile = tuple(new_profile)
        if is_valid_profile(new_profile, truncation_type, p):
            return new_profile, truncation_type
        else:
            raise ValueError("Invalid profile")
    else: # p odd
        if profile is None or profile == Infinity:
            # no specified profile or infinite profile: return profile
            # for the entire Steenrod algebra
            new_profile = ((), ())
            truncation_type = Infinity
        else: # profile should be a list or tuple of length 2
            assert isinstance(profile, (list, tuple)) and len(profile) == 2, \
                "Invalid form for profile"
            e = profile[0]
            k = profile[1]
            if isinstance(e, (list, tuple)):
                # e is a list or tuple: use it as is.  if
                # truncation_type not specified, set it to 0. remove
                # appropriate trailing terms.
                if truncation_type == 'auto':
                    truncation_type = 0
                # remove trailing terms
                while e and e[-1] == truncation_type:
                    e = e[:-1]
                e = tuple(e)
            elif callable(e):
                # e is a function: turn it into a tuple.  if
                # truncation_type not specified, set it to 'infinity'
                # if the function is ever infinite; otherwise set it
                # to 0.  remove appropriate trailing terms.
                if precision is None:
                    e_precision = 100
                else:
                    e_precision = precision
                if truncation_type == 'auto':
                    truncation_type = Infinity
                e = [max(0, e(i)) for i in range(1, e_precision)]
                # remove trailing terms
                while e and e[-1] == truncation_type:
                    del e[-1]
                e = tuple(e)
            if isinstance(k, (list, tuple)):
                # k is a list or tuple: use it as is.
                k = tuple(k)
            elif callable(k):
                # k is a function: turn it into a tuple.
                if precision is None:
                    k_precision = 100
                else:
                    k_precision = precision
                k = tuple([k(i) for i in range(k_precision-1)])
            # Remove trailing ones from k if truncation_type is 'zero',
            # remove trailing twos if truncation_type is 'Infinity'.
            if truncation_type == 0:
                while k and k[-1] == 1:
                    k = k[:-1]
            else:
                while k and k[-1] == 2:
                    k = k[:-1]
            new_profile = (e, k)
        if is_valid_profile(new_profile, truncation_type, p, generic=True):
            return new_profile, truncation_type
        else:
            raise ValueError("Invalid profile")

######################################################
# string representations for elements


def milnor_mono_to_string(mono, latex=False, generic=False):
    """
    String representation of element of the Milnor basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:

    - ``mono`` -- if `generic=False`, tuple of nonnegative integers (a,b,c,...);
      if `generic=True`, pair of tuples of nonnegative integers ((e0, e1, e2,
      ...), (r1, r2, ...))

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    - ``generic`` -- whether to format generically, or for the prime 2 (default)

    OUTPUT: ``rep`` -- string

    This returns a string like ``Sq(a,b,c,...)`` when `generic=False`, or a string
    like ``Q_e0 Q_e1 Q_e2 ... P(r1, r2, ...)`` when `generic=True`.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import milnor_mono_to_string
        sage: milnor_mono_to_string((1,2,3,4))
        'Sq(1,2,3,4)'
        sage: milnor_mono_to_string((1,2,3,4),latex=True)
        '\\text{Sq}(1,2,3,4)'
        sage: milnor_mono_to_string(((1,0), (2,3,1)), generic=True)
        'Q_{1} Q_{0} P(2,3,1)'
        sage: milnor_mono_to_string(((1,0), (2,3,1)), latex=True, generic=True)
        'Q_{1} Q_{0} \\mathcal{P}(2,3,1)'

    The empty tuple represents the unit element::

        sage: milnor_mono_to_string(())
        '1'
        sage: milnor_mono_to_string((), generic=True)
        '1'
    """
    if latex:
        if not generic:
            sq = "\\text{Sq}"
            P = "\\text{Sq}"
        else:
            P = "\\mathcal{P}"
    else:
        if not generic:
            sq = "Sq"
            P = "Sq"
        else:
            P = "P"
    if mono == () or mono == (0,) or (generic and len(mono[0]) + len(mono[1]) == 0):
        return "1"
    else:
        if not generic:
            string = sq + "(" + str(mono[0])
            for n in mono[1:]:
                string = string + "," + str(n)
            string = string + ")"
        else:
            string = ""
            if len(mono[0]) > 0:
                for e in mono[0]:
                    string = string + "Q_{" + str(e) + "} "
            if len(mono[1]) > 0:
                string = string + P + "(" + str(mono[1][0])
                for n in mono[1][1:]:
                    string = string + "," + str(n)
                string = string + ")"
        return string.strip(" ")


def serre_cartan_mono_to_string(mono, latex=False, generic=False):
    r"""
    String representation of element of the Serre-Cartan basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:

    - ``mono`` -- tuple of positive integers (a,b,c,...)  when `generic=False`,
      or tuple (e0, n1, e1, n2, ...) when `generic=True`, where each ei is 0 or
      1, and each ni is positive

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    - ``generic`` -- whether to format generically, or for the prime 2 (default)

    OUTPUT: ``rep`` -- string

    This returns a string like ``Sq^{a} Sq^{b} Sq^{c} ...`` when
    `generic=False`, or a string like
    ``\beta^{e0} P^{n1} \beta^{e1} P^{n2} ...`` when `generic=True`.
    is odd.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import serre_cartan_mono_to_string
        sage: serre_cartan_mono_to_string((1,2,3,4))
        'Sq^{1} Sq^{2} Sq^{3} Sq^{4}'
        sage: serre_cartan_mono_to_string((1,2,3,4),latex=True)
        '\\text{Sq}^{1} \\text{Sq}^{2} \\text{Sq}^{3} \\text{Sq}^{4}'
        sage: serre_cartan_mono_to_string((0,5,1,1,0), generic=True)
        'P^{5} beta P^{1}'
        sage: serre_cartan_mono_to_string((0,5,1,1,0), generic=True, latex=True)
        '\\mathcal{P}^{5} \\beta \\mathcal{P}^{1}'

    The empty tuple represents the unit element 1::

        sage: serre_cartan_mono_to_string(())
        '1'
        sage: serre_cartan_mono_to_string((), generic=True)
        '1'
    """
    if latex:
        if not generic:
            sq = "\\text{Sq}"
            P = "\\text{Sq}"
        else:
            P = "\\mathcal{P}"
    else:
        if not generic:
            sq = "Sq"
            P = "Sq"
        else:
            P = "P"
    if len(mono) == 0 or mono == (0,):
        return "1"
    else:
        if not generic:
            string = ""
            for n in mono:
                string = string + sq + "^{" + str(n) + "} "
        else:
            string = ""
            index = 0
            for n in mono:
                from sage.misc.functional import is_even
                if is_even(index):
                    if n == 1:
                        if latex:
                            string = string + "\\beta "
                        else:
                            string = string + "beta "
                else:
                    string = string + P + "^{" + str(n) + "} "
                index += 1
        return string.strip(" ")


def wood_mono_to_string(mono, latex=False):
    """
    String representation of element of Wood's Y and Z bases.

    This is used by the _repr_ and _latex_ methods.

    INPUT:

    - ``mono`` -- tuple of pairs of nonnegative integers (s,t)

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    OUTPUT:

    ``string`` -- concatenation of strings of the form
    ``Sq^{2^s (2^{t+1}-1)}`` for each pair (s,t)

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import wood_mono_to_string
        sage: wood_mono_to_string(((1,2),(3,0)))
        'Sq^{14} Sq^{8}'
        sage: wood_mono_to_string(((1,2),(3,0)),latex=True)
        '\\text{Sq}^{14} \\text{Sq}^{8}'

    The empty tuple represents the unit element::

        sage: wood_mono_to_string(())
        '1'
    """
    if latex:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return "1"
    else:
        string = ""
        for (s,t) in mono:
            string = string + sq + "^{" + \
                str(2**s * (2**(t+1)-1)) + "} "
        return string.strip(" ")


def wall_mono_to_string(mono, latex=False):
    """
    String representation of element of Wall's basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:

    - ``mono`` -- tuple of pairs of nonnegative integers (m,k) with `m
      >= k`

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    OUTPUT:

    ``string`` -- concatenation of strings ``Q^{m}_{k}`` for each pair (m,k)

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import wall_mono_to_string
        sage: wall_mono_to_string(((1,2),(3,0)))
        'Q^{1}_{2} Q^{3}_{0}'
        sage: wall_mono_to_string(((1,2),(3,0)),latex=True)
        'Q^{1}_{2} Q^{3}_{0}'

    The empty tuple represents the unit element::

        sage: wall_mono_to_string(())
        '1'
    """
    if len(mono) == 0:
        return "1"
    else:
        string = ""
        for (m,k) in mono:
            string = string + "Q^{" + str(m) + "}_{" \
                + str(k) + "} "
        return string.strip(" ")


def wall_long_mono_to_string(mono, latex=False):
    """
    Alternate string representation of element of Wall's basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:

    - ``mono`` -- tuple of pairs of nonnegative integers (m,k) with `m
      >= k`

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    OUTPUT:

    ``string`` -- concatenation of strings of the form ``Sq^(2^m)``

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import wall_long_mono_to_string
        sage: wall_long_mono_to_string(((1,2),(3,0)))
        'Sq^{1} Sq^{2} Sq^{4} Sq^{8}'
        sage: wall_long_mono_to_string(((1,2),(3,0)),latex=True)
        '\\text{Sq}^{1} \\text{Sq}^{2} \\text{Sq}^{4} \\text{Sq}^{8}'

    The empty tuple represents the unit element::

        sage: wall_long_mono_to_string(())
        '1'
    """
    if latex:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return "1"
    else:
        string = ""
        for (m,k) in mono:
            for i in range(k,m+1):
                string = string + sq + "^{" + str(2**i) + "} "
        return string.strip(" ")


def arnonA_mono_to_string(mono, latex=False, p=2):
    """
    String representation of element of Arnon's A basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:

    - ``mono`` -- tuple of pairs of nonnegative integers
       (m,k) with `m >= k`

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    OUTPUT: concatenation of strings of the form ``X^{m}_{k}`` for each pair
    (m,k)

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import arnonA_mono_to_string
        sage: arnonA_mono_to_string(((1,2),(3,0)))
        'X^{1}_{2} X^{3}_{0}'
        sage: arnonA_mono_to_string(((1,2),(3,0)),latex=True)
        'X^{1}_{2} X^{3}_{0}'

    The empty tuple represents the unit element::

        sage: arnonA_mono_to_string(())
        '1'
    """
    if len(mono) == 0:
        return "1"
    else:
        string = ""
        for (m,k) in mono:
            string = string + "X^{" + str(m) + "}_{" \
                + str(k) + "} "
        return string.strip(" ")


def arnonA_long_mono_to_string(mono, latex=False, p=2):
    """
    Alternate string representation of element of Arnon's A basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:

    - ``mono`` -- tuple of pairs of nonnegative integers (m,k) with `m
      >= k`

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    OUTPUT:

    ``string`` -- concatenation of strings of the form ``Sq(2^m)``

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import arnonA_long_mono_to_string
        sage: arnonA_long_mono_to_string(((1,2),(3,0)))
        'Sq^{8} Sq^{4} Sq^{2} Sq^{1}'
        sage: arnonA_long_mono_to_string(((1,2),(3,0)),latex=True)
        '\\text{Sq}^{8} \\text{Sq}^{4} \\text{Sq}^{2} \\text{Sq}^{1}'

    The empty tuple represents the unit element::

        sage: arnonA_long_mono_to_string(())
        '1'
    """
    if latex:
        sq = "\\text{Sq}"
    else:
        sq = "Sq"
    if len(mono) == 0:
        return "1"
    else:
        string = ""
        for (m,k) in mono:
            for i in range(m,k-1,-1):
                string = string + sq + "^{" + str(2**i) + "} "
        return string.strip(" ")


def pst_mono_to_string(mono, latex=False, generic=False):
    r"""
    String representation of element of a `P^s_t`-basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:

    - ``mono`` -- tuple of pairs of integers (s,t) with `s >= 0`, `t >
      0`

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    - ``generic`` -- whether to format generically, or for the prime 2 (default)

    OUTPUT:

    ``string`` -- concatenation of strings of the form ``P^{s}_{t}``
    for each pair (s,t)

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import pst_mono_to_string
        sage: pst_mono_to_string(((1,2),(0,3)), generic=False)
        'P^{1}_{2} P^{0}_{3}'
        sage: pst_mono_to_string(((1,2),(0,3)),latex=True, generic=False)
        'P^{1}_{2} P^{0}_{3}'
        sage: pst_mono_to_string(((1,4), (((1,2), 1),((0,3), 2))), generic=True)
        'Q_{1} Q_{4} P^{1}_{2} (P^{0}_{3})^2'
        sage: pst_mono_to_string(((1,4), (((1,2), 1),((0,3), 2))), latex=True, generic=True)
        'Q_{1} Q_{4} P^{1}_{2} (P^{0}_{3})^{2}'

    The empty tuple represents the unit element::

        sage: pst_mono_to_string(())
        '1'
    """
    if len(mono) == 0:
        return "1"
    else:
        string = ""
        if not generic:
            for (s,t) in mono:
                string = string + "P^{" + str(s) + "}_{" \
                    + str(t) + "} "
        else:
            for e in mono[0]:
                string = string + "Q_{" + str(e) + "} "
            for ((s,t), n) in mono[1]:
                if n == 1:
                    string = string + "P^{" + str(s) + "}_{" \
                        + str(t) + "} "
                else:
                    if latex:
                        pow = "{%s}" % n
                    else:
                        pow = str(n)
                    string = string + "(P^{" + str(s) + "}_{" \
                        + str(t) + "})^" + pow + " "
        return string.strip(" ")


def comm_mono_to_string(mono, latex=False, generic=False):
    r"""
    String representation of element of a commutator basis.

    This is used by the _repr_ and _latex_ methods.

    INPUT:

    - ``mono`` -- tuple of pairs of integers (s,t) with `s >= 0`, `t >
      0`

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    - ``generic`` -- whether to format generically, or for the prime 2 (default)

    OUTPUT:

    ``string`` -- concatenation of strings of the form ``c_{s,t}``
    for each pair (s,t)

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import comm_mono_to_string
        sage: comm_mono_to_string(((1,2),(0,3)), generic=False)
        'c_{1,2} c_{0,3}'
        sage: comm_mono_to_string(((1,2),(0,3)), latex=True)
        'c_{1,2} c_{0,3}'
        sage: comm_mono_to_string(((1, 4), (((1,2), 1),((0,3), 2))), generic=True)
        'Q_{1} Q_{4} c_{1,2} c_{0,3}^2'
        sage: comm_mono_to_string(((1, 4), (((1,2), 1),((0,3), 2))), latex=True, generic=True)
        'Q_{1} Q_{4} c_{1,2} c_{0,3}^{2}'

    The empty tuple represents the unit element::

        sage: comm_mono_to_string(())
        '1'
    """
    if len(mono) == 0:
        return "1"
    else:
        string = ""
        if not generic:
            for (s,t) in mono:
                string = string + "c_{" + str(s) + "," \
                    + str(t) + "} "
        else:
            for e in mono[0]:
                string = string + "Q_{" + str(e) + "} "
            for ((s,t), n) in mono[1]:
                string = string + "c_{" + str(s) + "," \
                    + str(t) + "}"
                if n > 1:
                    if latex:
                        pow = "^{%s}" % n
                    else:
                        pow = "^%s" % n
                    string = string + pow
                string = string + " "
        return string.strip(" ")


def comm_long_mono_to_string(mono, p, latex=False, generic=False):
    r"""
    Alternate string representation of element of a commutator basis.

    Okay in low dimensions, but gets unwieldy as the dimension
    increases.

    INPUT:

    - ``mono`` -- tuple of pairs of integers (s,t) with `s >= 0`, `t >
      0`

    - ``latex`` -- boolean (default: ``False``); if ``True``, output
      LaTeX string

    - ``generic`` -- whether to format generically, or for the prime 2 (default)

    OUTPUT:

    ``string`` -- concatenation of strings of the form ``s_{2^s... 2^(s+t-1)}``
    for each pair (s,t)

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_misc import comm_long_mono_to_string
        sage: comm_long_mono_to_string(((1,2),(0,3)), 2)
        's_{24} s_{124}'
        sage: comm_long_mono_to_string(((1,2),(0,3)), 2, latex=True)
        's_{24} s_{124}'
        sage: comm_long_mono_to_string(((1, 4), (((1,2), 1),((0,3), 2))), 5, generic=True)
        'Q_{1} Q_{4} s_{5,25} s_{1,5,25}^2'
        sage: comm_long_mono_to_string(((1, 4), (((1,2), 1),((0,3), 2))), 3, latex=True, generic=True)
        'Q_{1} Q_{4} s_{3,9} s_{1,3,9}^{2}'

    The empty tuple represents the unit element::

        sage: comm_long_mono_to_string((), p=2)
        '1'
    """
    if len(mono) == 0:
        return "1"
    else:
        string = ""
        if not generic:
            for (s,t) in mono:
                if s + t > 4:
                    comma = ","
                else:
                    comma = ""
                string = string + "s_{"
                for i in range(t):
                    string = string + str(2**(s+i)) + comma
                string = string.strip(",") + "} "
        else:
            for e in mono[0]:
                string = string + "Q_{" + str(e) + "} "
            for ((s,t), n) in mono[1]:
                string = string + "s_{"
                for i in range(t):
                    string = string + str(p**(s+i)) + ","
                string = string.strip(",") + "}"
                if n > 1:
                    if latex:
                        pow = "^{%s}" % n
                    else:
                        pow = "^%s" % n
                    string = string + pow
                string = string + " "
        return string.strip(" ")

# miscellany:


def convert_perm(m):
    """
    Convert tuple m of nonnegative integers to a permutation in
    one-line form.

    INPUT:

    - ``m`` -- tuple of nonnegative integers with no repetitions

    OUTPUT:

    ``list`` -- conversion of ``m`` to a permutation of the set
    1,2,...,len(m)

    If ``m=(3,7,4)``, then one can view ``m`` as representing the
    permutation of the set `(3,4,7)` sending 3 to 3, 4 to 7, and 7 to
    4. This function converts ``m`` to the list ``[1,3,2]``, which
    represents essentially the same permutation, but of the set
    `(1,2,3)`. This list can then be passed to :func:`Permutation
    <sage.combinat.permutation.Permutation>`, and its signature can be
    computed.

    EXAMPLES::

        sage: sage.algebras.steenrod.steenrod_algebra_misc.convert_perm((3,7,4))
        [1, 3, 2]
        sage: sage.algebras.steenrod.steenrod_algebra_misc.convert_perm((5,0,6,3))
        [3, 1, 4, 2]
    """
    m2 = sorted(m)
    return [list(m2).index(x)+1 for x in m]
