r"""
Factory

This file contains the constructor classes and functions for `p`-adic rings and fields.

AUTHORS:

- David Roe

TESTS::

    sage: R = ZpLC(2)
    doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
    See http://trac.sagemath.org/23505 for details.
    sage: R = ZpLF(2)
    sage: R = QpLC(2)
    sage: R = QpLF(2)
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import experimental

from sage.structure.factory import UniqueFactory
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.structure.factorization import Factorization
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.structure.element import is_Element
from .padic_base_leaves import (pAdicRingCappedRelative,
                                pAdicRingCappedAbsolute,
                                pAdicRingFixedMod,
                                pAdicRingFloatingPoint,
                                pAdicRingLattice,
                                pAdicRingRelaxed,
                                pAdicFieldCappedRelative,
                                pAdicFieldFloatingPoint,
                                pAdicFieldLattice,
                                pAdicFieldRelaxed)
from . import padic_printing

######################################################
# ext_table --
# This dictionary controls what class is created by the extension
# factory when it finds a given class in the ground ring of the tower.
######################################################

from .padic_extension_leaves import (EisensteinExtensionFieldCappedRelative,
                                     EisensteinExtensionRingFixedMod,
                                     EisensteinExtensionRingCappedAbsolute,
                                     EisensteinExtensionRingCappedRelative,
                                     UnramifiedExtensionFieldCappedRelative,
                                     UnramifiedExtensionRingCappedRelative,
                                     UnramifiedExtensionRingCappedAbsolute,
                                     UnramifiedExtensionRingFixedMod,
                                     UnramifiedExtensionFieldFloatingPoint,
                                     UnramifiedExtensionRingFloatingPoint)
from .relative_extension_leaves import \
        (RelativeRamifiedExtensionRingFixedMod,
         RelativeRamifiedExtensionRingCappedAbsolute,
         RelativeRamifiedExtensionRingCappedRelative,
         RelativeRamifiedExtensionFieldCappedRelative,
         RelativeRamifiedExtensionRingFloatingPoint,
         RelativeRamifiedExtensionFieldFloatingPoint)
from functools import reduce
#This imports all of the classes used in the ext_table below.

ext_table = {}
ext_table['e', pAdicFieldCappedRelative] = EisensteinExtensionFieldCappedRelative
ext_table['e', pAdicRingCappedAbsolute] = EisensteinExtensionRingCappedAbsolute
ext_table['e', pAdicRingCappedRelative] = EisensteinExtensionRingCappedRelative
ext_table['e', pAdicRingFixedMod] = EisensteinExtensionRingFixedMod
#ext_table['e', pAdicRingFloatingPoint] = EisensteinExtensionRingFloatingPoint
#ext_table['e', pAdicFieldFloatingPoint] = EisensteinExtensionFieldFloatingPoint
#ext_table['p', pAdicFieldCappedRelative] = pAdicGeneralExtensionFieldCappedRelative
#ext_table['p', pAdicRingCappedAbsolute] = pAdicGeneralExtensionRingCappedAbsolute
#ext_table['p', pAdicRingCappedRelative] = pAdicGeneralExtensionRingCappedRelative
#ext_table['p', pAdicRingFixedMod] = pAdicGeneralExtensionRingFixedMod
ext_table['u', pAdicFieldCappedRelative] = UnramifiedExtensionFieldCappedRelative
ext_table['u', pAdicRingCappedAbsolute] = UnramifiedExtensionRingCappedAbsolute
ext_table['u', pAdicRingCappedRelative] = UnramifiedExtensionRingCappedRelative
ext_table['u', pAdicRingFixedMod] = UnramifiedExtensionRingFixedMod
ext_table['u', pAdicRingFloatingPoint] = UnramifiedExtensionRingFloatingPoint
ext_table['u', pAdicFieldFloatingPoint] = UnramifiedExtensionFieldFloatingPoint
ext_table['re', pAdicRingFixedMod] = RelativeRamifiedExtensionRingFixedMod
ext_table['re', pAdicRingCappedAbsolute] = RelativeRamifiedExtensionRingCappedAbsolute
ext_table['re', pAdicRingCappedRelative] = RelativeRamifiedExtensionRingCappedRelative
ext_table['re', pAdicFieldCappedRelative] = RelativeRamifiedExtensionFieldCappedRelative
ext_table['re', pAdicRingFloatingPoint] = RelativeRamifiedExtensionRingFloatingPoint
ext_table['re', pAdicFieldFloatingPoint] = RelativeRamifiedExtensionFieldFloatingPoint

def _canonicalize_show_prec(type, print_mode, show_prec=None):
    r"""
    Return a canonical string value for show_prec depending of the type,
    the print_mode and the given value.

    INPUT:

    - ``type`` -- a string: ``'capped-rel'``, ``'capped-abs'``, ``'fixed-mod'`` or ``'floating-point'``,
      ``'lattice-cap'`` or ``'lattice-float'``

    - ``print_mode`` -- a string: ``'series'``, ``'terse'``, ``'val-unit'``, ``'digits'``, ``'bars'``

    - ``show_prec`` -- a boolean, string or ``None``

    OUTPUT:

    A string, either ``'bigoh'``, ``'dots'`` or ``'none'``

    EXAMPLES::

        sage: from sage.rings.padics.factory import _canonicalize_show_prec
        sage: _canonicalize_show_prec('floating-point', 'series')
        'none'

        sage: _canonicalize_show_prec('capped-rel', 'series')
        'bigoh'
        sage: _canonicalize_show_prec('capped-rel', 'series', False)
        'none'

        sage: _canonicalize_show_prec('capped-abs', 'digits')
        'dots'
        sage: _canonicalize_show_prec('capped-abs', 'digits', 'bigoh')
        'bigoh'

    TESTS::

        sage: _canonicalize_show_prec('capped-abs', 'digits', 'my_precision')
        Traceback (most recent call last):
        ...
        ValueError: show_prec must be either a boolean, 'none', 'bigoh' or 'dots' when printing mode is digits
    """
    # Note that None means "choose the default for this ring", while 'none' means "don't print precision".
    if show_prec is None:
        show_prec = type not in ('floating-point', 'fixed-mod')
    if show_prec is False:
        return "none"
    if show_prec is True:
        if print_mode in ('series', 'terse', 'val-unit'):
            return "bigoh"
        else:
            return "dots"
    if print_mode in ('series', 'terse', 'val-unit'):
        if show_prec not in ('none', 'bigoh'):
            raise ValueError("show_prec must be either a boolean, 'none' or 'bigoh' when printing mode is %s" % print_mode)
    else:
        if show_prec not in ('none', 'bigoh', 'dots'):
            raise ValueError("show_prec must be either a boolean, 'none', 'bigoh' or 'dots' when printing mode is %s" % print_mode)
    return show_prec


def get_key_base(p, prec, type, print_mode, names, ram_name, print_pos, print_sep, print_alphabet, print_max_terms, show_prec, check, valid_types, label=None):
    r"""
    This implements create_key for Zp and Qp: moving it here prevents code duplication.

    It fills in unspecified values and checks for contradictions in the input.  It also standardizes irrelevant options so that duplicate parents are not created.

    EXAMPLES::

        sage: from sage.rings.padics.factory import get_key_base
        sage: get_key_base(11, 5, 'capped-rel', None, None, None, None, ':', None, None, False, True, ['capped-rel'])
        (11, 5, 'capped-rel', 'series', '11', True, '|', (), -1, 'none', None)
        sage: get_key_base(12, 5, 'capped-rel', 'digits', None, None, None, None, None, None, True, False, ['capped-rel'])
        (12,
         5,
         'capped-rel',
         'digits',
         '12',
         True,
         '|',
         ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B'),
         -1,
         'dots',
         None)
    """
    if check:
        if not isinstance(p, Integer):
            p = Integer(p)
        if not p.is_prime():
            raise ValueError("p must be prime")
        if type == 'lattice-cap':
            relative_cap = absolute_cap = None
            if prec is not None:
                try:
                    relative_cap, absolute_cap = prec
                except (ValueError, TypeError):
                    relative_cap = prec
            if relative_cap is not None:
                if relative_cap is not Infinity:
                    try:
                        relative_cap = Integer(relative_cap)
                    except TypeError:
                        raise TypeError("relative cap must be either a positive integer or infinity")
                    if relative_cap <= 0:
                        raise ValueError("relative cap must be positive")
            if absolute_cap is not None:
                try:
                    absolute_cap = Integer(absolute_cap)
                except TypeError:
                    raise TypeError("absolute cap must be an integer")
            if relative_cap is None and absolute_cap is None:
                relative_cap = DEFAULT_PREC
                absolute_cap = 2 * DEFAULT_PREC
            elif relative_cap is None:
                relative_cap = Infinity
            elif absolute_cap is None:
                absolute_cap = 2 * relative_cap
            prec = (relative_cap, absolute_cap)
        elif type == 'relaxed':
            default_prec = halting_prec = None
            secure = False
            if isinstance(prec, (list, tuple)):
                if len(prec) == 1:
                    default_prec = prec
                elif len(prec) == 2:
                    default_prec, halting_prec = prec
                else:
                    default_prec = prec[0]
                    halting_prec = prec[1]
                    secure = prec[2]
            else:
                default_prec = prec
            if default_prec is None:
                default_prec = DEFAULT_PREC
            if halting_prec is None:
                halting_prec = 2 * default_prec
            halting_prec = max(default_prec, halting_prec)
            prec = (default_prec, halting_prec, secure)
        else:
            if prec is not None:
                prec = Integer(prec)
    if prec is None:
        if type == 'lattice-cap':
            prec = (DEFAULT_PREC, 2*DEFAULT_PREC)
        else:
            prec = DEFAULT_PREC
    print_ram_name = ram_name
    if isinstance(print_mode, dict):
        if 'pos' in print_mode:
            print_pos = print_mode['pos']
        if 'ram_name' in print_mode:
            print_ram_name = print_mode['ram_name']
        if 'unram_name' in print_mode:
            # print_unram_name = print_mode['unram_name']
            pass
        if 'sep' in print_mode:
            print_sep = print_mode['sep']
        if 'alphabet' in print_mode:
            print_alphabet = print_mode['alphabet']
        if 'max_ram_terms' in print_mode:
            print_max_terms = print_mode['max_ram_terms']
        if 'max_terms' in print_mode:
            print_max_terms = print_mode['max_terms']
        if 'show_prec' in print_mode:
            show_prec = print_mode['show_prec']
        if 'mode' in print_mode:
            print_mode = print_mode['mode']
        else:
            print_mode = None
    if print_mode is None:
        print_mode = padic_printing._printer_defaults.mode()
    if print_pos is None:
        print_pos = not padic_printing._printer_defaults.allow_negatives()
    if print_sep is None:
        print_sep = padic_printing._printer_defaults.sep()
    if print_alphabet is None:
        print_alphabet = padic_printing._printer_defaults.alphabet()
    if print_max_terms is None:
        print_max_terms = padic_printing._printer_defaults.max_series_terms()

    # We eliminate irrelevant print options (e.g. print_pos if p = 2)
    if p == 2 or print_mode == 'digits':
        print_pos = True # we want this hard-coded so that we don't get duplicate parents if the keys differ.
    if print_mode == 'digits':
        print_ram_name = None
        print_alphabet = print_alphabet[:p]
    else:
        print_alphabet = []
    if print_mode != 'bars':
        print_sep = '|'
    if print_mode in ['terse', 'val-unit']:
        print_max_terms = -1

    if isinstance(names, tuple):
        names = names[0]
    if names is None and print_ram_name is None:
        name = str(p)
    elif names is not None and print_ram_name is not None:
        if not isinstance(names, str):
            names = str(names)
        if not isinstance(print_ram_name, str):
            print_ram_name = str(print_ram_name)
        if names != print_ram_name:
            raise ValueError("If both names (%s) and print_ram_name (%s) are specified, they must agree"%(names, print_ram_name))
        name = names
    else:
        if names is None:
            names = print_ram_name
        if isinstance(names, str):
            name = names
        else:
            name = str(names)
    if type not in valid_types:
        raise ValueError("type must be %s"%(", ".join(valid_types)))
    show_prec = _canonicalize_show_prec(type, print_mode, show_prec)
    key = (p, prec, type, print_mode, name, print_pos, print_sep, tuple(print_alphabet), print_max_terms, show_prec, label)
    return key

#######################################################################################################
#
#  p-Adic Fields
#  Qp -- base field
#  Qq -- unramified extension field of Qp
#  QpCR, QpLC, QpLF, QqCR -- shortcuts for capped relative and lattice versions of Qp and Qq
#
#######################################################################################################


padic_field_cache = {}
DEFAULT_PREC = Integer(20)

class Qp_class(UniqueFactory):
    r"""
    A creation function for `p`-adic fields.

    INPUT:

    - ``p`` -- integer: the `p` in `\QQ_p`

    - ``prec`` -- integer (default: ``20``) the precision cap of the field.
      In the lattice capped case, ``prec`` can either be a
      pair (``relative_cap``, ``absolute_cap``) or an integer
      (understood at relative cap).
      In the relaxed case, ``prec`` can be either a
      pair (``default_prec``, ``halting_prec``) or an integer
      (understood at default precision).
      Except in the floating point case, individual elements keep track of
      their own precision.  See TYPES and PRECISION below.

    - ``type`` -- string (default: ``'capped-rel'``) Valid types are
      ``'capped-rel'``, ``'floating-point'``, ``'lattice-cap'``, ``'lattice-float'``.
      See TYPES and PRECISION below

    - ``print_mode`` -- string (default: ``None``).  Valid modes are 'series',
      'val-unit', 'terse', 'digits', and 'bars'. See PRINTING below

    - ``names`` -- string or tuple (defaults to a string representation of
      `p`).  What to use whenever `p` is printed.

    - ``ram_name`` -- string.  Another way to specify the name; for
      consistency with the ``Qq`` and ``Zq`` and extension functions.

    - ``print_pos`` -- bool (default ``None``) Whether to only use positive
      integers in the representations of elements. See PRINTING below.

    - ``print_sep`` -- string (default ``None``) The separator character used
      in the ``'bars'`` mode. See PRINTING below.

    - ``print_alphabet`` -- tuple (default ``None``) The encoding into digits
      for use in the 'digits' mode. See PRINTING below.

    - ``print_max_terms`` -- integer (default ``None``) The maximum number of
      terms shown.  See PRINTING below.

    - ``show_prec`` -- a boolean or a string (default ``None``) Specify how
      the precision is printed. See PRINTING below.

    - ``check`` -- bool (default ``True``) whether to check if `p` is prime.
      Non-prime input may cause seg-faults (but can also be useful for
      base n expansions for example)

    - ``label`` -- string (default ``None``) used for lattice precision to
      create parents with different lattices.

    OUTPUT:

    - The corresponding `p`-adic field.

    TYPES AND PRECISION:

    There are two main types of precision for a `p`-adic element.
    The first is relative precision, which gives the number of known
    `p`-adic digits::

        sage: R = Qp(5, 20, 'capped-rel', 'series'); a = R(675); a
        2*5^2 + 5^4 + O(5^22)
        sage: a.precision_relative()
        20

    The second type of precision is absolute precision, which gives
    the power of `p` that this element is defined modulo::

        sage: a.precision_absolute()
        22

    There are several types of `p`-adic fields, depending on the methods
    used for tracking precision. Namely, we have:

    - capped relative fields (``type='capped-rel'``)

    - capped absolute fields (``type='capped-abs'``)

    - fixed modulus fields (``type='fixed-mod'``)

    - floating point fields (``type='floating-point'``)

    - lattice precision fields (``type='lattice-cap'`` or ``type='lattice-float'``)

    - exact fields with relaxed arithmetics (``type='relaxed'``)

    In the capped relative case, the relative precision of an element
    is restricted to be at most a certain value, specified at the
    creation of the field.  Individual elements also store their own
    precision, so the effect of various arithmetic operations on
    precision is tracked.  When you cast an exact element into a
    capped relative field, it truncates it to the precision cap of the
    field. ::

        sage: R = Qp(5, 5, 'capped-rel', 'series'); a = R(4006); a
        1 + 5 + 2*5^3 + 5^4 + O(5^5)
        sage: b = R(4025); b
        5^2 + 2*5^3 + 5^4 + 5^5 + O(5^7)
        sage: a + b
        1 + 5 + 5^2 + 4*5^3 + 2*5^4 + O(5^5)

    In the floating point case, elements do not track their
    precision, but the relative precision of elements is truncated
    during arithmetic to the precision cap of the field.

    In the lattice case, precision on elements is tracked by a global
    lattice that is updated after every operation, yielding better
    precision behavior at the cost of higher memory and runtime usage.
    We refer to the documentation of the function :func:`ZpLC` for a
    small demonstration of the capabilities of this precision model.

    Finally, the model for relaxed `p`-adics is quite different from any of
    the other types. In addition to storing a finite approximation, one
    also stores a method for increasing the precision.
    A quite interesting feature with relaxed `p`-adics is the possibility to
    create (in some cases) self-referent numbers, that are numbers whose
    `n`-th digit is defined by the previous ones.
    We refer to the documentation of the function :func:`ZpL` for a
    small demonstration of the capabilities of this precision model.

    PRINTING:

    There are many different ways to print `p`-adic elements.  The way
    elements of a given field print is controlled by options passed in
    at the creation of the field.  There are five basic printing modes
    (series, val-unit, terse, digits and bars), as well as various
    options that either hide some information in the print
    representation or sometimes make print representations more
    compact.  Note that the printing options affect whether different
    `p`-adic fields are considered equal.

    1. **series**: elements are displayed as series in `p`. ::

        sage: R = Qp(5, print_mode='series'); a = R(70700); a
        3*5^2 + 3*5^4 + 2*5^5 + 4*5^6 + O(5^22)
        sage: b = R(-70700); b
        2*5^2 + 4*5^3 + 5^4 + 2*5^5 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + 4*5^20 + 4*5^21 + O(5^22)

      *print_pos* controls whether negatives can be used in the
      coefficients of powers of `p`. ::

        sage: S = Qp(5, print_mode='series', print_pos=False); a = S(70700); a
        -2*5^2 + 5^3 - 2*5^4 - 2*5^5 + 5^7 + O(5^22)
        sage: b = S(-70700); b
        2*5^2 - 5^3 + 2*5^4 + 2*5^5 - 5^7 + O(5^22)

      *print_max_terms* limits the number of terms that appear. ::

        sage: T = Qp(5, print_mode='series', print_max_terms=4); b = R(-70700); repr(b)
        '2*5^2 + 4*5^3 + 5^4 + 2*5^5 + ... + O(5^22)'

      *names* affects how the prime is printed. ::

        sage: U.<p> = Qp(5); p
        p + O(p^21)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``) or 'bigoh'.
      The default is ``False`` for the ``'floating-point'`` type
      and ``True`` for all other types. ::

        sage: Qp(5)(6)
        1 + 5 + O(5^20)
        sage: Qp(5, show_prec='none')(6)
        1 + 5

        sage: QpFP(5)(6)
        1 + 5

      *print_sep* and *print_alphabet* have no effect in series mode.

      Note that print options affect equality::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    2. **val-unit**: elements are displayed as p^k*u::

        sage: R = Qp(5, print_mode='val-unit'); a = R(70700); a
        5^2 * 2828 + O(5^22)
        sage: b = R(-707/5); b
        5^-1 * 95367431639918 + O(5^19)

      *print_pos* controls whether to use a balanced representation or
      not. ::

        sage: S = Qp(5, print_mode='val-unit', print_pos=False); b = S(-70700); b
        5^2 * (-2828) + O(5^22)

      *names* affects how the prime is printed. ::

        sage: T = Qp(5, print_mode='val-unit', names='pi'); a = T(70700); a
        pi^2 * 2828 + O(pi^22)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``) or 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` type
      and ``True`` for all other types. ::

        sage: Qp(5, print_mode='val-unit', show_prec=False)(30)
        5 * 6

      *print_max_terms*, *print_sep* and *print_alphabet* have no effect.

      Equality again depends on the printing options::

        sage: R == S, R == T, S == T
        (False, False, False)

    3. **terse**: elements are displayed as an integer in base 10 or the
       quotient of an integer by a power of `p` (still in base 10)::

        sage: R = Qp(5, print_mode='terse'); a = R(70700); a
        70700 + O(5^22)
        sage: b = R(-70700); b
        2384185790944925 + O(5^22)
        sage: c = R(-707/5); c
        95367431639918/5 + O(5^19)

      The denominator, as of version 3.3, is always printed
      explicitly as a power of `p`, for predictability. ::

        sage: d = R(707/5^2); d
        707/5^2 + O(5^18)

      *print_pos* controls whether to use a balanced representation or not. ::

        sage: S = Qp(5, print_mode='terse', print_pos=False); b = S(-70700); b
        -70700 + O(5^22)
        sage: c = S(-707/5); c
        -707/5 + O(5^19)

      *name* affects how the name is printed. ::

        sage: T.<unif> = Qp(5, print_mode='terse'); c = T(-707/5); c
        95367431639918/unif + O(unif^19)
        sage: d = T(-707/5^10); d
        95367431639918/unif^10 + O(unif^10)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``) or 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` type
      and ``True`` for all other types. ::

        sage: Qp(5, print_mode='terse', show_prec=False)(6)
        6

      *print_max_terms*, *print_sep* and *print_alphabet* have no effect.

      Equality depends on printing options::

        sage: R == S, R == T, S == T
        (False, False, False)

    4. **digits**: elements are displayed as a string of base `p` digits

      Restriction: you can only use the digits printing mode for
      small primes.  Namely, `p` must be less than the length of the
      alphabet tuple (default alphabet has length 62). ::

        sage: R = Qp(5, print_mode='digits'); a = R(70700); repr(a)
        '...0000000000000004230300'
        sage: b = R(-70700); repr(b)
        '...4444444444444440214200'
        sage: c = R(-707/5); repr(c)
        '...4444444444444443413.3'
        sage: d = R(-707/5^2); repr(d)
        '...444444444444444341.33'

      Observe that the significant 0's are printed even if they are
      located in front of the number. On the contrary, unknown digits
      located after the comma appears as question marks.
      The precision can therefore be read in this mode as well.
      Here are more examples::

        sage: p = 7
        sage: K = Qp(p, prec=10, print_mode='digits')
        sage: repr(K(1))
        '...0000000001'
        sage: repr(K(p^2))
        '...000000000100'
        sage: repr(K(p^-5))
        '...00000.00001'
        sage: repr(K(p^-20))
        '...?.??????????0000000001'

      *print_max_terms* limits the number of digits that are printed.
      Note that if the valuation of the element is very negative, more
      digits will be printed. ::

        sage: S = Qp(5, print_max_terms=4); S(-70700)
        2*5^2 + 4*5^3 + 5^4 + 2*5^5 + ... + O(5^22)
        sage: S(-707/5^2)
        3*5^-2 + 3*5^-1 + 1 + 4*5 + ... + O(5^18)
        sage: S(-707/5^6)
        3*5^-6 + 3*5^-5 + 5^-4 + 4*5^-3 + ... + O(5^14)
        sage: S(-707/5^6,absprec=-2)
        3*5^-6 + 3*5^-5 + 5^-4 + 4*5^-3 + O(5^-2)
        sage: S(-707/5^4)
        3*5^-4 + 3*5^-3 + 5^-2 + 4*5^-1 + ... + O(5^16)

      *print_alphabet* controls the symbols used to substitute for digits
      greater than 9.

      Defaults to ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z')::

        sage: T = Qp(5, print_mode='digits', print_alphabet=('1','2','3','4','5')); repr(T(-70700))
        '...5555555555555551325311'

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'dots'
      (or equivalently ``True``) or 'bigoh'.
      The default is ``False`` for the ``'floating-point'`` type
      and ``True`` for all other types. ::

        sage: repr(Zp(5, print_mode='digits', show_prec=True)(6))
        '...00000000000000000011'

        sage: repr(Zp(5, print_mode='digits', show_prec='bigoh')(6))
        '11 + O(5^20)'

      *print_pos*, *name* and *print_sep* have no effect.

      Equality depends on printing options::

        sage: R == S, R == T, S == T
        (False, False, False)

    5. **bars**: elements are displayed as a string of base `p` digits
       with separators::

        sage: R = Qp(5, print_mode='bars'); a = R(70700); repr(a)
        '...4|2|3|0|3|0|0'
        sage: b = R(-70700); repr(b)
        '...4|4|4|4|4|4|4|4|4|4|4|4|4|4|4|0|2|1|4|2|0|0'
        sage: d = R(-707/5^2); repr(d)
        '...4|4|4|4|4|4|4|4|4|4|4|4|4|4|4|3|4|1|.|3|3'

      Again, note that it's not possible to read off the precision from the representation in this mode.

      *print_pos* controls whether the digits can be negative. ::

        sage: S = Qp(5, print_mode='bars',print_pos=False); b = S(-70700); repr(b)
        '...-1|0|2|2|-1|2|0|0'

      *print_max_terms* limits the number of digits that are printed.
      Note that if the valuation of the element is very negative, more
      digits will be printed. ::

        sage: T = Qp(5, print_max_terms=4); T(-70700)
        2*5^2 + 4*5^3 + 5^4 + 2*5^5 + ... + O(5^22)
        sage: T(-707/5^2)
        3*5^-2 + 3*5^-1 + 1 + 4*5 + ... + O(5^18)
        sage: T(-707/5^6)
        3*5^-6 + 3*5^-5 + 5^-4 + 4*5^-3 + ... + O(5^14)
        sage: T(-707/5^6,absprec=-2)
        3*5^-6 + 3*5^-5 + 5^-4 + 4*5^-3 + O(5^-2)
        sage: T(-707/5^4)
        3*5^-4 + 3*5^-3 + 5^-2 + 4*5^-1 + ... + O(5^16)

      *print_sep* controls the separation character. ::

        sage: U = Qp(5, print_mode='bars', print_sep=']['); a = U(70700); repr(a)
        '...4][2][3][0][3][0][0'

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'dots'
      (or equivalently ``True``) or 'bigoh'
      The default is ``False`` for the ``'floating-point'`` type
      and ``True`` for all other types. ::

        sage: repr(Qp(5, print_mode='bars', show_prec='bigoh')(6))
        '...1|1 + O(5^20)'

      *name* and *print_alphabet* have no effect.

      Equality depends on printing options::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    EXAMPLES::

        sage: K = Qp(15, check=False); a = K(999); a
        9 + 6*15 + 4*15^2 + O(15^20)
    """
    def create_key(self, p, prec = None, type = 'capped-rel', print_mode = None,
                   names = None, ram_name = None, print_pos = None,
                   print_sep = None, print_alphabet = None, print_max_terms = None, show_prec = None, check = True,
                   label = None):   # specific to Lattice precision
        r"""
        Creates a key from input parameters for ``Qp``.

        See the documentation for ``Qp`` for more information.

        TESTS::

            sage: Qp.create_key(5,40)
            (5, 40, 'capped-rel', 'series', '5', True, '|', (), -1, 'bigoh', None)
        """
        if isinstance(names, (int, Integer)):
            # old pickle; names is what used to be halt.
            names = ram_name
            ram_name = print_pos
            print_pos = print_sep
            print_alphabet = print_max_terms
            print_max_terms = check
            check = True
        if label is not None and type not in ['lattice-cap','lattice-float']:
            raise ValueError("label keyword only supported for lattice precision")
        return get_key_base(p, prec, type, print_mode, names, ram_name, print_pos, print_sep, print_alphabet, print_max_terms, show_prec, check, ['capped-rel', 'floating-point', 'lattice-cap', 'lattice-float', 'relaxed'], label)

    def create_object(self, version, key):
        r"""
        Creates an object using a given key.

        See the documentation for ``Qp`` for more information.

        TESTS::

            sage: Qp.create_object((3,4,2),(5, 41, 'capped-rel', 'series', '5', True, '|', (), -1))
            5-adic Field with capped relative precision 41
        """
        if version[0] < 3 or (version[0] == 3 and version[1] < 2) or (version[0] == 3 and version[1] == 2 and version[2] < 3):
            p, prec, type, print_mode, name = key
            print_pos, print_sep, print_alphabet, print_max_terms = None, None, None, None
        elif version[0] < 8:
            p, prec, type, print_mode, name, print_pos, print_sep, print_alphabet, print_max_terms = key
            show_prec = None
            label = None
        else:
            p, prec, type, print_mode, name, print_pos, print_sep, print_alphabet, print_max_terms, show_prec, label = key
        if (version[0] < 4 or (len(version) > 1 and version[0] == 4 and version[1] < 5) or
            (len(version) > 2 and version[0] == 4 and version[1] == 5 and version[2] < 3)):
            # keys changed in order to reduce irrelevant duplications: e.g. two Qps with print_mode 'series'
            # that are identical except for different 'print_alphabet' now return the same object.
            key = get_key_base(p, prec, type, print_mode, name, None, print_pos, print_sep, print_alphabet,
                               print_max_terms, None, False, ['capped-rel', 'fixed-mod', 'capped-abs'])
            try:
                obj = self._cache[version, key]()
                if obj is not None:
                    return obj
            except KeyError:
                pass
            p, prec, type, print_mode, name, print_pos, print_sep, print_alphabet, print_max_terms, show_prec, label = key

        if type == 'capped-rel':
            if print_mode == 'terse':
                return pAdicFieldCappedRelative(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                          'ram_name': name, 'max_terse_terms': print_max_terms, 'show_prec': show_prec}, name)
            else:
                return pAdicFieldCappedRelative(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                          'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name)
        elif type == 'floating-point':
            if print_mode == 'terse':
                return pAdicFieldFloatingPoint(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                          'ram_name': name, 'max_terse_terms': print_max_terms, 'show_prec': show_prec}, name)
            else:
                return pAdicFieldFloatingPoint(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                         'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name)
        elif type == 'relaxed':
            if print_mode == 'terse':
                return pAdicFieldRelaxed(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                'ram_name': name, 'max_terse_terms': print_max_terms, 'show_prec': show_prec}, name)
            else:
                return pAdicFieldRelaxed(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name)
        elif type[:8] == 'lattice-':
            subtype = type[8:]
            if print_mode == 'terse':
                return pAdicFieldLattice(p, prec, subtype, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                            'ram_name': name, 'max_terse_terms': print_max_terms, 'show_prec': show_prec}, name, label)
            else:
                return pAdicFieldLattice(p, prec, subtype, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                            'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name, label)
        else:
            raise ValueError("unexpected type")

Qp = Qp_class("Qp")


######################################################
# Qq -- unramified extensions
######################################################

def Qq(q, prec = None, type = 'capped-rel', modulus = None, names=None,
          print_mode=None, ram_name = None, res_name = None, print_pos = None,
       print_sep = None, print_max_ram_terms = None,
       print_max_unram_terms = None, print_max_terse_terms = None, show_prec=None, check = True, implementation = 'FLINT'):
    r"""
    Given a prime power `q = p^n`, return the unique unramified
    extension of `\QQ_p` of degree `n`.

    INPUT:

    - ``q`` -- integer, list, tuple or ``Factorization`` object. If ``q`` is an
      integer, it is the prime power `q` in `\QQ_q`. If ``q`` is a
      ``Factorization`` object, it is the factorization of the prime power `q`.
      As a tuple it is the pair ``(p, n)``, and as a list it is a single
      element list ``[(p, n)]``.

    - ``prec`` -- integer (default: ``20``) the precision cap of the field.
      Individual elements keep track of their own precision.  See
      TYPES and PRECISION below.

    - ``type`` -- string (default: ``'capped-rel'``) Valid types are
      ``'capped-rel'``, ``'floating-point'``, ``'lattice-cap'``
      and ``'lattice-float'``.  See TYPES and PRECISION below

    - ``modulus`` -- polynomial (default ``None``) A polynomial defining an
      unramified extension of `\QQ_p`.  See MODULUS below.

    - ``names`` -- string or tuple (``None`` is only allowed when `q=p`).  The
      name of the generator, reducing to a generator of the residue
      field.

    - ``print_mode`` -- string (default: ``None``).  Valid modes are ``'series'``,
      ``'val-unit'``, ``'terse'``, and ``'bars'``. See PRINTING below.

    - ``ram_name`` -- string (defaults to string representation of `p` if
      None).  ``ram_name`` controls how the prime is printed. See PRINTING
      below.

    - ``res_name`` -- string (defaults to ``None``, which corresponds to
      adding a ``'0'`` to the end of the name).  Controls how elements of
      the reside field print.

    - ``print_pos`` -- bool (default ``None``) Whether to only use positive
      integers in the representations of elements. See PRINTING below.

    - ``print_sep`` -- string (default ``None``) The separator character used
      in the ``'bars'`` mode. See PRINTING below.

    - ``print_max_ram_terms`` -- integer (default ``None``) The maximum number
      of powers of `p` shown.  See PRINTING below.

    - ``print_max_unram_terms`` -- integer (default ``None``) The maximum
      number of entries shown in a coefficient of `p`.  See PRINTING
      below.

    - ``print_max_terse_terms`` -- integer (default ``None``) The maximum
      number of terms in the polynomial representation of an element
      (using ``'terse'``).  See PRINTING below.

    - ``show_prec`` -- bool (default ``None``) whether to show the precision
      for elements.  See PRINTING below.

    - ``check`` -- bool (default ``True``) whether to check inputs.

    OUTPUT:

    - The corresponding unramified `p`-adic field.

    TYPES AND PRECISION:

    There are two types of precision for a `p`-adic element.  The first
    is relative precision, which gives the number of known `p`-adic
    digits::

        sage: R.<a> = Qq(25, 20, 'capped-rel', print_mode='series'); b = 25*a; b
        a*5^2 + O(5^22)
        sage: b.precision_relative()
        20

    The second type of precision is absolute precision, which gives
    the power of `p` that this element is defined modulo::

        sage: b.precision_absolute()
        22

    There are two types of unramified `p`-adic fields: capped relative
    fields, floating point fields.

    In the capped relative case, the relative precision of an element
    is restricted to be at most a certain value, specified at the
    creation of the field.  Individual elements also store their own
    precision, so the effect of various arithmetic operations on
    precision is tracked.  When you cast an exact element into a
    capped relative field, it truncates it to the precision cap of the
    field. ::

        sage: R.<a> = Qq(9, 5, 'capped-rel', print_mode='series'); b = (1+2*a)^4; b
        2 + (2*a + 2)*3 + (2*a + 1)*3^2 + O(3^5)
        sage: c = R(3249); c
        3^2 + 3^4 + 3^5 + 3^6 + O(3^7)
        sage: b + c
        2 + (2*a + 2)*3 + (2*a + 2)*3^2 + 3^4 + O(3^5)

    In the floating point case, elements do not track their
    precision, but the relative precision of elements is truncated
    during arithmetic to the precision cap of the field.

    MODULUS:

    The modulus needs to define an unramified extension of `\QQ_p`: when it
    is reduced to a polynomial over `\GF{p}` it should be irreducible.

    The modulus can be given in a number of forms.

    1. A **polynomial**.

      The base ring can be `\ZZ`, `\QQ`, `\ZZ_p`, `\QQ_p`, `\GF{p}`. ::

        sage: P.<x> = ZZ[]
        sage: R.<a> = Qq(27, modulus = x^3 + 2*x + 1); R.modulus()
        (1 + O(3^20))*x^3 + O(3^20)*x^2 + (2 + O(3^20))*x + 1 + O(3^20)
        sage: P.<x> = QQ[]
        sage: S.<a> = Qq(27, modulus = x^3 + 2*x + 1)
        sage: P.<x> = Zp(3)[]
        sage: T.<a> = Qq(27, modulus = x^3 + 2*x + 1)
        sage: P.<x> = Qp(3)[]
        sage: U.<a> = Qq(27, modulus = x^3 + 2*x + 1)
        sage: P.<x> = GF(3)[]
        sage: V.<a> = Qq(27, modulus = x^3 + 2*x + 1)

      Which form the modulus is given in has no effect on the unramified
      extension produced::

        sage: R == S, S == T, T == U, U == V
        (True, True, True, False)

      unless the precision of the modulus differs.  In the case of V,
      the modulus is only given to precision 1, so the resulting field
      has a precision cap of 1. ::

        sage: V.precision_cap()
        1
        sage: U.precision_cap()
        20
        sage: P.<x> = Qp(3)[]
        sage: modulus = x^3 + (2 + O(3^7))*x + (1 + O(3^10))
        sage: modulus
        (1 + O(3^20))*x^3 + (2 + O(3^7))*x + 1 + O(3^10)
        sage: W.<a> = Qq(27, modulus = modulus); W.precision_cap()
        7

    2. The modulus can also be given as a **symbolic expression**. ::

        sage: x = var('x')
        sage: X.<a> = Qq(27, modulus = x^3 + 2*x + 1); X.modulus()
        (1 + O(3^20))*x^3 + O(3^20)*x^2 + (2 + O(3^20))*x + 1 + O(3^20)
        sage: X == R
        True

      By default, the polynomial chosen is the standard lift of the
      generator chosen for `\GF{q}`. ::

        sage: GF(125, 'a').modulus()
        x^3 + 3*x + 3
        sage: Y.<a> = Qq(125); Y.modulus()
        (1 + O(5^20))*x^3 + O(5^20)*x^2 + (3 + O(5^20))*x + 3 + O(5^20)

      However, you can choose another polynomial if desired (as long as
      the reduction to `\GF{p}[x]` is irreducible). ::

        sage: P.<x> = ZZ[]
        sage: Z.<a> = Qq(125, modulus = x^3 + 3*x^2 + x + 1); Z.modulus()
        (1 + O(5^20))*x^3 + (3 + O(5^20))*x^2 + (1 + O(5^20))*x + 1 + O(5^20)
        sage: Y == Z
        False

    PRINTING:

    There are many different ways to print `p`-adic elements.  The way
    elements of a given field print is controlled by options passed in
    at the creation of the field.  There are four basic printing modes
    (``'series'``, ``'val-unit'``, ``'terse'`` and ``'bars'``; ``'digits'`` is not available), as
    well as various options that either hide some information in the
    print representation or sometimes make print representations more
    compact.  Note that the printing options affect whether different
    `p`-adic fields are considered equal.

    1. **series**: elements are displayed as series in `p`. ::

        sage: R.<a> = Qq(9, 20, 'capped-rel', print_mode='series'); (1+2*a)^4
        2 + (2*a + 2)*3 + (2*a + 1)*3^2 + O(3^20)
        sage: -3*(1+2*a)^4
        3 + a*3^2 + 3^3 + (2*a + 2)*3^4 + (2*a + 2)*3^5 + (2*a + 2)*3^6 + (2*a + 2)*3^7 + (2*a + 2)*3^8 + (2*a + 2)*3^9 + (2*a + 2)*3^10 + (2*a + 2)*3^11 + (2*a + 2)*3^12 + (2*a + 2)*3^13 + (2*a + 2)*3^14 + (2*a + 2)*3^15 + (2*a + 2)*3^16 + (2*a + 2)*3^17 + (2*a + 2)*3^18 + (2*a + 2)*3^19 + (2*a + 2)*3^20 + O(3^21)
        sage: ~(3*a+18)
        (a + 2)*3^-1 + 1 + 2*3 + (a + 1)*3^2 + 3^3 + 2*3^4 + (a + 1)*3^5 + 3^6 + 2*3^7 + (a + 1)*3^8 + 3^9 + 2*3^10 + (a + 1)*3^11 + 3^12 + 2*3^13 + (a + 1)*3^14 + 3^15 + 2*3^16 + (a + 1)*3^17 + 3^18 + O(3^19)

      *print_pos* controls whether negatives can be used in the
      coefficients of powers of `p`. ::

        sage: S.<b> = Qq(9, print_mode='series', print_pos=False); (1+2*b)^4
        -1 - b*3 - 3^2 + (b + 1)*3^3 + O(3^20)
        sage: -3*(1+2*b)^4
        3 + b*3^2 + 3^3 + (-b - 1)*3^4 + O(3^21)

      *ram_name* controls how the prime is printed. ::

        sage: T.<d> = Qq(9, print_mode='series', ram_name='p'); 3*(1+2*d)^4
        2*p + (2*d + 2)*p^2 + (2*d + 1)*p^3 + O(p^21)

      *print_max_ram_terms* limits the number of powers of `p` that appear. ::

        sage: U.<e> = Qq(9, print_mode='series', print_max_ram_terms=4); repr(-3*(1+2*e)^4)
        '3 + e*3^2 + 3^3 + (2*e + 2)*3^4 + ... + O(3^21)'

      *print_max_unram_terms* limits the number of terms that appear in a
      coefficient of a power of `p`. ::

        sage: V.<f> = Qq(128, prec = 8, print_mode='series'); repr((1+f)^9)
        '(f^3 + 1) + (f^5 + f^4 + f^3 + f^2)*2 + (f^6 + f^5 + f^4 + f + 1)*2^2 + (f^5 + f^4 + f^2 + f + 1)*2^3 + (f^6 + f^5 + f^4 + f^3 + f^2 + f + 1)*2^4 + (f^5 + f^4)*2^5 + (f^6 + f^5 + f^4 + f^3 + f + 1)*2^6 + (f + 1)*2^7 + O(2^8)'
        sage: V.<f> = Qq(128, prec = 8, print_mode='series', print_max_unram_terms = 3); repr((1+f)^9)
        '(f^3 + 1) + (f^5 + f^4 + ... + f^2)*2 + (f^6 + f^5 + ... + 1)*2^2 + (f^5 + f^4 + ... + 1)*2^3 + (f^6 + f^5 + ... + 1)*2^4 + (f^5 + f^4)*2^5 + (f^6 + f^5 + ... + 1)*2^6 + (f + 1)*2^7 + O(2^8)'
        sage: V.<f> = Qq(128, prec = 8, print_mode='series', print_max_unram_terms = 2); repr((1+f)^9)
        '(f^3 + 1) + (f^5 + ... + f^2)*2 + (f^6 + ... + 1)*2^2 + (f^5 + ... + 1)*2^3 + (f^6 + ... + 1)*2^4 + (f^5 + f^4)*2^5 + (f^6 + ... + 1)*2^6 + (f + 1)*2^7 + O(2^8)'
        sage: V.<f> = Qq(128, prec = 8, print_mode='series', print_max_unram_terms = 1); repr((1+f)^9)
        '(f^3 + ...) + (f^5 + ...)*2 + (f^6 + ...)*2^2 + (f^5 + ...)*2^3 + (f^6 + ...)*2^4 + (f^5 + ...)*2^5 + (f^6 + ...)*2^6 + (f + ...)*2^7 + O(2^8)'
        sage: V.<f> = Qq(128, prec = 8, print_mode='series', print_max_unram_terms = 0); repr((1+f)^9 - 1 - f^3)
        '(...)*2 + (...)*2^2 + (...)*2^3 + (...)*2^4 + (...)*2^5 + (...)*2^6 + (...)*2^7 + O(2^8)'

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` type
      and ``True`` for all other types. ::

        sage: U.<e> = Qq(9, 2, show_prec=False); repr(-3*(1+2*e)^4)
        '3 + e*3^2'

      *print_sep* and *print_max_terse_terms* have no effect.

      Note that print options affect equality::

        sage: R == S, R == T, R == U, R == V, S == T, S == U, S == V, T == U, T == V, U == V
        (False, False, False, False, False, False, False, False, False, False)

    2. **val-unit**: elements are displayed as `p^k u`::

        sage: R.<a> = Qq(9, 7, print_mode='val-unit'); b = (1+3*a)^9 - 1; b
        3^3 * (15 + 64*a) + O(3^7)
        sage: ~b
        3^-3 * (41 + a) + O(3)

      *print_pos* controls whether to use a balanced representation or
      not. ::

        sage: S.<a> = Qq(9, 7, print_mode='val-unit', print_pos=False); b = (1+3*a)^9 - 1; b
        3^3 * (15 - 17*a) + O(3^7)
        sage: ~b
        3^-3 * (-40 + a) + O(3)

      *ram_name* affects how the prime is printed. ::

        sage: A.<x> = Qp(next_prime(10^6), print_mode='val-unit')[]
        sage: T.<a> = Qq(next_prime(10^6)^3, 4, print_mode='val-unit', ram_name='p', modulus=x^3+385831*x^2+106556*x+321036); b = ~(next_prime(10^6)^2*(a^2 + a - 4)); b
        p^-2 * (503009563508519137754940 + 704413692798200940253892*a + 968097057817740999537581*a^2) + O(p^2)
        sage: b * (a^2 + a - 4)
        p^-2 * 1 + O(p^2)

      *print_max_terse_terms* controls how many terms of the polynomial
      appear in the unit part. ::

        sage: U.<a> = Qq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3); b = ~(17*(a^3-a+14)); b
        17^-1 * (22110411 + 11317400*a + 20656972*a^2 + ...) + O(17^5)
        sage: b*17*(a^3-a+14)
        1 + O(17^6)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` type
      and ``True`` for all other types. ::

        sage: U.<e> = Qq(9, 2, print_mode='val-unit', show_prec=False); repr(-3*(1+2*e)^4)
        '3 * (1 + 3*e)'

      *print_sep*, *print_max_ram_terms* and *print_max_unram_terms* have no
      effect.

      Equality again depends on the printing options::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    3. **terse**: elements are displayed as a polynomial of degree less
       than the degree of the extension. ::

        sage: R.<a> = Qq(125, print_mode='terse')
        sage: (a+5)^177
        68210977979428 + 90313850704069*a + 73948093055069*a^2 + O(5^20)
        sage: (a/5+1)^177
        68210977979428/5^177 + 90313850704069/5^177*a + 73948093055069/5^177*a^2 + O(5^-157)

      As of version 3.3, if coefficients of the polynomial are
      non-integral, they are always printed with an explicit power of `p`
      in the denominator. ::

        sage: 5*a + a^2/25
        5*a + 1/5^2*a^2 + O(5^18)

      *print_pos* controls whether to use a balanced representation or
      not. ::

        sage: (a-5)^6
        22864 + 95367431627998*a + 8349*a^2 + O(5^20)
        sage: S.<a> = Qq(125, print_mode='terse', print_pos=False); b = (a-5)^6; b
        22864 - 12627*a + 8349*a^2 + O(5^20)
        sage: (a - 1/5)^6
        -20624/5^6 + 18369/5^5*a + 1353/5^3*a^2 + O(5^14)

      *ram_name* affects how the prime is printed. ::

        sage: T.<a> = Qq(125, print_mode='terse', ram_name='p'); (a - 1/5)^6
        95367431620001/p^6 + 18369/p^5*a + 1353/p^3*a^2 + O(p^14)

      *print_max_terse_terms* controls how many terms of the polynomial
      are shown. ::

        sage: U.<a> = Qq(625, print_mode='terse', print_max_terse_terms=2); (a-1/5)^6
        106251/5^6 + 49994/5^5*a + ... + O(5^14)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` type
      and ``True`` for all other types. ::

        sage: U.<e> = Qq(9, 2, print_mode='terse', show_prec=False); repr(-3*(1+2*e)^4)
        '3 + 9*e'

      *print_sep*, *print_max_ram_terms* and *print_max_unram_terms* have no
      effect.

      Equality again depends on the printing options::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    4. **digits**: This print mode is not available when the residue
       field is not prime.

       It might make sense to have a dictionary for small fields, but
       this isn't implemented.

    5. **bars**: elements are displayed in a similar fashion to
       series, but more compactly. ::

        sage: R.<a> = Qq(125); (a+5)^6
        (4*a^2 + 3*a + 4) + (3*a^2 + 2*a)*5 + (a^2 + a + 1)*5^2 + (3*a + 2)*5^3 + (3*a^2 + a + 3)*5^4 + (2*a^2 + 3*a + 2)*5^5 + O(5^20)
        sage: R.<a> = Qq(125, print_mode='bars', prec=8); repr((a+5)^6)
        '...[2, 3, 2]|[3, 1, 3]|[2, 3]|[1, 1, 1]|[0, 2, 3]|[4, 3, 4]'
        sage: repr((a-5)^6)
        '...[0, 4]|[1, 4]|[2, 0, 2]|[1, 4, 3]|[2, 3, 1]|[4, 4, 3]|[2, 4, 4]|[4, 3, 4]'

      Note that elements with negative valuation are shown with a
      decimal point at valuation 0. ::

        sage: repr((a+1/5)^6)
        '...[3]|[4, 1, 3]|.|[1, 2, 3]|[3, 3]|[0, 0, 3]|[0, 1]|[0, 1]|[1]'
        sage: repr((a+1/5)^2)
        '...[0, 0, 1]|.|[0, 2]|[1]'

      If not enough precision is known, ``'?'`` is used instead. ::

        sage: repr((a+R(1/5,relprec=3))^7)
        '...|.|?|?|?|?|[0, 1, 1]|[0, 2]|[1]'

      Note that it's not possible to read off the precision from the
      representation in this mode. ::

        sage: b = a + 3; repr(b)
        '...[3, 1]'
        sage: c = a + R(3, 4); repr(c)
        '...[3, 1]'
        sage: b.precision_absolute()
        8
        sage: c.precision_absolute()
        4

      *print_pos* controls whether the digits can be negative. ::

        sage: S.<a> = Qq(125, print_mode='bars', print_pos=False); repr((a-5)^6)
        '...[1, -1, 1]|[2, 1, -2]|[2, 0, -2]|[-2, -1, 2]|[0, 0, -1]|[-2]|[-1, -2, -1]'
        sage: repr((a-1/5)^6)
        '...[0, 1, 2]|[-1, 1, 1]|.|[-2, -1, -1]|[2, 2, 1]|[0, 0, -2]|[0, -1]|[0, -1]|[1]'

      *print_max_ram_terms* controls the maximum number of "digits" shown.
      Note that this puts a cap on the relative precision, not the
      absolute precision. ::

        sage: T.<a> = Qq(125, print_max_ram_terms=3, print_pos=False); (a-5)^6
        (-a^2 - 2*a - 1) - 2*5 - a^2*5^2 + ... + O(5^20)
        sage: 5*(a-5)^6 + 50
        (-a^2 - 2*a - 1)*5 - a^2*5^3 + (2*a^2 - a - 2)*5^4 + ... + O(5^21)

      *print_sep* controls the separating character (``'|'`` by default). ::

        sage: U.<a> = Qq(625, print_mode='bars', print_sep=''); b = (a+5)^6; repr(b)
        '...[0, 1][4, 0, 2][3, 2, 2, 3][4, 2, 2, 4][0, 3][1, 1, 3][3, 1, 4, 1]'

      *print_max_unram_terms* controls how many terms are shown in each
      "digit"::

        sage: with local_print_mode(U, {'max_unram_terms': 3}): repr(b)
        '...[0, 1][4,..., 0, 2][3,..., 2, 3][4,..., 2, 4][0, 3][1,..., 1, 3][3,..., 4, 1]'
        sage: with local_print_mode(U, {'max_unram_terms': 2}): repr(b)
        '...[0, 1][4,..., 2][3,..., 3][4,..., 4][0, 3][1,..., 3][3,..., 1]'
        sage: with local_print_mode(U, {'max_unram_terms': 1}): repr(b)
        '...[..., 1][..., 2][..., 3][..., 4][..., 3][..., 3][..., 1]'
        sage: with local_print_mode(U, {'max_unram_terms':0}): repr(b-75*a)
        '...[...][...][...][...][][...][...]'

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'dots'
      (or equivalently ``True``) or 'bigoh'
      The default is ``False`` for the ``'floating-point'`` type
      and ``True`` for all other types. ::

        sage: U.<e> = Qq(9, 2, print_mode='bars', show_prec=True); repr(-3*(1+2*e)^4)
        '...[0, 1]|[1]|[]'

      *ram_name* and *print_max_terse_terms* have no effect.

      Equality depends on printing options::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    EXAMPLES:

    Unlike for ``Qp``, you can't create ``Qq(N)`` when ``N`` is not a prime power.

    However, you can use ``check=False`` to pass in a pair in order to not
    have to factor.  If you do so, you need to use names explicitly
    rather than the ``R.<a>`` syntax. ::

        sage: p = next_prime(2^123)
        sage: k = Qp(p)
        sage: R.<x> = k[]
        sage: K = Qq([(p, 5)], modulus=x^5+x+4, names='a', ram_name='p', print_pos=False, check=False)
        sage: K.0^5
        (-a - 4) + O(p^20)

    In tests on ``sage.math.washington.edu``, the creation of ``K`` as above took an
    average of 1.58ms, while::

        sage: K = Qq(p^5, modulus=x^5+x+4, names='a', ram_name='p', print_pos=False, check=True)

    took an average of 24.5ms.  Of course, with smaller primes these
    savings disappear.

    TESTS:

    Check that :trac:`8162` is resolved::

        sage: R = Qq([(5,3)], names="alpha", check=False); R
        5-adic Unramified Extension Field in alpha defined by x^3 + 3*x + 3
        sage: Qq((5, 3), names="alpha") is R
        True
        sage: Qq(125.factor(), names="alpha") is R
        True

    Check that :trac:`18606` is resolved::

        sage: x = QQ['x'].gen()
        sage: F = Qp(5,20)
        sage: K0 = F.extension(x^2-F(13),names = 'g')
        sage: K1 = F.extension(x^2-13,names = 'g')
        sage: K0 is K1
        True
    """
    if is_Element(q):
        F = Integer(q).factor()
        if len(F) != 1:
            raise ValueError("q must be a prime power")
        q = F
    if isinstance(q, Factorization):
        if len(q) != 1:
            raise ValueError("q must be a factorization of a prime power")
        q = list(q)
    if not isinstance(q, (list, tuple)):
        raise TypeError("q must be an integer, list, tuple or Factorization")
    if len(q) != 2:
        if len(q) != 1:
            raise ValueError("q must have shape [(p,k)]")
        q = q[0]
    if len(q) != 2:
        raise ValueError("q must have shape (p,k)")
    if not isinstance(q, tuple):
        q = tuple(q)

    p,k = q
    if not isinstance(p, Integer):
        p = Integer(p)
    if not isinstance(k, Integer):
        k = Integer(k)

    if check:
        if not p.is_prime() or k <=0:
            raise ValueError("q must be a prime power")

    if prec is not None and not isinstance(prec, Integer):
        prec = Integer(prec)

    base = Qp(p=p, prec=prec, type=type, print_mode=print_mode, names=ram_name, print_pos=print_pos,
              print_sep=print_sep, print_max_terms=print_max_ram_terms, show_prec=show_prec, check=check)

    if k == 1:
        return base

    if isinstance(names, (list, tuple)):
        if len(names) != 1:
            raise ValueError("must provide exactly one generator name")
        names = names[0]
    if names is None:
        raise TypeError("You must specify the name of the generator.")
    if not isinstance(names, str):
        names = str(names)

    if res_name is None:
        res_name = names + '0'

    if modulus is None:
        from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
        modulus = GF((p, k), res_name).modulus().change_ring(ZZ)
    return ExtensionFactory(base=base, modulus=modulus, prec=prec, print_mode=print_mode,
                            names=names, res_name=res_name, ram_name=ram_name, print_pos=print_pos,
                            print_sep=print_sep, print_max_ram_terms=print_max_ram_terms,
                            print_max_unram_terms=print_max_unram_terms,
                            print_max_terse_terms=print_max_terse_terms, show_prec=show_prec, check=check,
                            unram=True, implementation=implementation)

######################################################
# Short constructor names for different types
######################################################

def QpCR(p, prec = None, *args, **kwds):
    r"""
    A shortcut function to create capped relative `p`-adic fields.

    Same functionality as ``Qp``.  See documentation for ``Qp`` for a
    description of the input parameters.

    EXAMPLES::

        sage: QpCR(5, 40)
        5-adic Field with capped relative precision 40
    """
    return Qp(p, prec, 'capped-rel', *args, **kwds)

def QpFP(p, prec = None, *args, **kwds):
    r"""
    A shortcut function to create floating point `p`-adic fields.

    Same functionality as ``Qp``.  See documentation for ``Qp`` for a
    description of the input parameters.

    EXAMPLES::

        sage: QpFP(5, 40)
        5-adic Field with floating precision 40
    """
    return Qp(p, prec, 'floating-point', *args, **kwds)

def QqCR(q, prec = None, *args, **kwds):
    r"""
    A shortcut function to create capped relative unramified `p`-adic
    fields.

    Same functionality as ``Qq``.  See documentation for ``Qq`` for a
    description of the input parameters.

    EXAMPLES::

        sage: R.<a> = QqCR(25, 40); R
        5-adic Unramified Extension Field in a defined by x^2 + 4*x + 2
    """
    return Qq(q, prec, 'capped-rel', *args, **kwds)

def QqFP(q, prec = None, *args, **kwds):
    r"""
    A shortcut function to create floating point unramified `p`-adic
    fields.

    Same functionality as ``Qq``.  See documentation for ``Qq`` for a
    description of the input parameters.

    EXAMPLES::

        sage: R.<a> = QqFP(25, 40); R
        5-adic Unramified Extension Field in a defined by x^2 + 4*x + 2
    """
    return Qq(q, prec, 'floating-point', *args, **kwds)

@experimental(23505)
def QpLC(p, prec = None, *args, **kwds):
    r"""
    A shortcut function to create `p`-adic fields with lattice precision.

    See :func:`ZpLC` for more information about this model of precision.

    EXAMPLES::

        sage: R = QpLC(2)
        sage: R
        2-adic Field with lattice-cap precision
    """
    return Qp(p, prec, 'lattice-cap', *args, **kwds)

@experimental(23505)
def QpLF(p, prec = None, *args, **kwds):
    r"""
    A shortcut function to create `p`-adic fields with lattice precision.

    See :func:`ZpLC` for more information about this model of precision.

    EXAMPLES::

        sage: R = QpLF(2)
        sage: R
        2-adic Field with lattice-float precision
    """
    return Qp(p, prec, 'lattice-float', *args, **kwds)

def QpER(p, prec=None, halt=None, secure=False, *args, **kwds):
    r"""
    A shortcut function to create relaxed `p`-adic fields.

    See :func:`ZpER` for more information about this model of precision.

    EXAMPLES::

        sage: R = QpER(2)
        sage: R
        2-adic Field handled with relaxed arithmetics
    """
    return Qp(p, (prec, halt, secure), 'relaxed', *args, **kwds)

#######################################################################################################
#
#  p-Adic Rings
#  Zp -- base rings
#  Zq -- unramified extension ring of Zp
#  ZpCR, ZpCA, ZpFM, ZpL, ZqCR, ZqCA, ZqFM, ZqL -- shortcuts for precision-type versions of Zp and Zq
#
#######################################################################################################

class Zp_class(UniqueFactory):
    r"""
    A creation function for `p`-adic rings.

    INPUT:

    - ``p`` -- integer: the `p` in `\ZZ_p`

    - ``prec`` -- integer (default: ``20``) the precision cap of the
      ring.  In the lattice capped case, ``prec`` can either be a
      pair (``relative_cap``, ``absolute_cap``) or an integer
      (understood at relative cap).
      In the relaxed case, ``prec`` can be either a
      pair (``default_prec``, ``halting_prec``) or an integer
      (understood at default precision).
      Except for the fixed modulus and floating point cases, individual elements
      keep track of their own precision.  See TYPES and PRECISION
      below.

    - ``type`` -- string (default: ``'capped-rel'``) Valid types are
      ``'capped-rel'``, ``'capped-abs'``, ``'fixed-mod'``,
      ``'floating-point'``, ``'lattice-cap'``, ``'lattice-float'``, ``'relaxed'``
      See TYPES and PRECISION below

    - ``print_mode`` -- string (default: ``None``).  Valid modes are
      ``'series'``, ``'val-unit'``, ``'terse'``, ``'digits'``, and
      ``'bars'``. See PRINTING below

    - ``names`` -- string or tuple (defaults to a string
      representation of `p`).  What to use whenever `p` is printed.

    - ``print_pos`` -- bool (default ``None``) Whether to only use
      positive integers in the representations of elements. See
      PRINTING below.

    - ``print_sep`` -- string (default ``None``) The separator
      character used in the ``'bars'`` mode. See PRINTING below.

    - ``print_alphabet`` -- tuple (default ``None``) The encoding into
      digits for use in the ``'digits'`` mode. See PRINTING below.

    - ``print_max_terms`` -- integer (default ``None``) The maximum
      number of terms shown.  See PRINTING below.

    - ``show_prec`` -- bool (default ``None``) whether to show the precision
      for elements.  See PRINTING below.

    - ``check`` -- bool (default ``True``) whether to check if `p` is
      prime.  Non-prime input may cause seg-faults (but can also be
      useful for base `n` expansions for example)

    - ``label`` -- string (default ``None``) used for lattice precision to
      create parents with different lattices.

    OUTPUT:

    - The corresponding `p`-adic ring.

    TYPES AND PRECISION:

    There are two main types of precision.
    The first is relative precision; it gives the number of known
    `p`-adic digits::

        sage: R = Zp(5, 20, 'capped-rel', 'series'); a = R(675); a
        2*5^2 + 5^4 + O(5^22)
        sage: a.precision_relative()
        20

    The second type of precision is absolute precision, which gives
    the power of `p` that this element is defined modulo::

        sage: a.precision_absolute()
        22

    There are several types of `p`-adic rings, depending on the methods
    used for tracking precision. Namely, we have:

    - capped relative rings (``type='capped-rel'``)

    - capped absolute rings (``type='capped-abs'``)

    - fixed modulus rings (``type='fixed-mod'``)

    - floating point rings (``type='floating-point'``)

    - lattice precision rings (``type='lattice-cap'`` or ``type='lattice-float'``)

    - exact fields with relaxed arithmetics (``type='relaxed'``)

    In the capped relative case, the relative precision of an element
    is restricted to be at most a certain value, specified at the
    creation of the field.  Individual elements also store their own
    precision, so the effect of various arithmetic operations on
    precision is tracked.  When you cast an exact element into a
    capped relative field, it truncates it to the precision cap of the
    field. ::

        sage: R = Zp(5, 5, 'capped-rel', 'series'); a = R(4006); a
        1 + 5 + 2*5^3 + 5^4 + O(5^5)
        sage: b = R(4025); b
        5^2 + 2*5^3 + 5^4 + 5^5 + O(5^7)
        sage: a + b
        1 + 5 + 5^2 + 4*5^3 + 2*5^4 + O(5^5)

    In the capped absolute type, instead of having a cap on the
    relative precision of an element there is instead a cap on the
    absolute precision.  Elements still store their own precisions,
    and as with the capped relative case, exact elements are truncated
    when cast into the ring. ::

        sage: R = Zp(5, 5, 'capped-abs', 'series'); a = R(4005); a
        5 + 2*5^3 + 5^4 + O(5^5)
        sage: b = R(4025); b
        5^2 + 2*5^3 + 5^4 + O(5^5)
        sage: a * b
        5^3 + 2*5^4 + O(5^5)
        sage: (a * b) // 5^3
        1 + 2*5 + O(5^2)

    The fixed modulus type is the leanest of the `p`-adic rings: it is
    basically just a wrapper around `\ZZ / p^n \ZZ` providing a unified
    interface with the rest of the `p`-adics.  This is the type you
    should use if your sole interest is speed.  It does not track
    precision of elements. ::

        sage: R = Zp(5,5,'fixed-mod','series'); a = R(4005); a
        5 + 2*5^3 + 5^4
        sage: a // 5
        1 + 2*5^2 + 5^3

    The floating point case is similar to the fixed modulus type
    in that elements do not trac their own precision.  However, relative
    precision is truncated with each operation rather than absolute precision.

    On the contrary, the lattice type tracks precision using lattices
    and automatic differentiation. It is rather slow but provides sharp
    (often optimal) results regarding precision.
    We refer to the documentation of the function :func:`ZpLC` for a
    small demonstration of the capabilities of this precision model.

    Finally, the model for relaxed `p`-adics is quite different from any of
    the other types. In addition to storing a finite approximation, one
    also stores a method for increasing the precision.
    A quite interesting feature with relaxed `p`-adics is the possibility to
    create (in some cases) self-referent numbers, that are numbers whose
    `n`-th digit is defined by the previous ones.
    We refer to the documentation of the function :func:`ZpL` for a
    small demonstration of the capabilities of this precision model.

    PRINTING:

    There are many different ways to print `p`-adic elements.  The
    way elements of a given ring print is controlled by options
    passed in at the creation of the ring.  There are five basic
    printing modes (series, val-unit, terse, digits and bars), as
    well as various options that either hide some information in
    the print representation or sometimes make print
    representations more compact.  Note that the printing options
    affect whether different `p`-adic fields are considered equal.

    1. **series**: elements are displayed as series in `p`. ::

        sage: R = Zp(5, print_mode='series'); a = R(70700); a
        3*5^2 + 3*5^4 + 2*5^5 + 4*5^6 + O(5^22)
        sage: b = R(-70700); b
        2*5^2 + 4*5^3 + 5^4 + 2*5^5 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + 4*5^20 + 4*5^21 + O(5^22)

      *print_pos* controls whether negatives can be used in the
      coefficients of powers of `p`. ::

        sage: S = Zp(5, print_mode='series', print_pos=False); a = S(70700); a
        -2*5^2 + 5^3 - 2*5^4 - 2*5^5 + 5^7 + O(5^22)
        sage: b = S(-70700); b
        2*5^2 - 5^3 + 2*5^4 + 2*5^5 - 5^7 + O(5^22)

      *print_max_terms* limits the number of terms that appear. ::

        sage: T = Zp(5, print_mode='series', print_max_terms=4); b = R(-70700); b
        2*5^2 + 4*5^3 + 5^4 + 2*5^5 + ... + O(5^22)

      *names* affects how the prime is printed. ::

        sage: U.<p> = Zp(5); p
        p + O(p^21)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` and
      ``'fixed-mod'`` types and ``True`` for all other types. ::

        sage: Zp(5, show_prec=False)(6)
        1 + 5

      *print_sep* and *print_alphabet* have no effect.

      Note that print options affect equality::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    2. **val-unit**: elements are displayed as `p^k u`::

        sage: R = Zp(5, print_mode='val-unit'); a = R(70700); a
        5^2 * 2828 + O(5^22)
        sage: b = R(-707*5); b
        5 * 95367431639918 + O(5^21)

      *print_pos* controls whether to use a balanced representation or
      not. ::

        sage: S = Zp(5, print_mode='val-unit', print_pos=False); b = S(-70700); b
        5^2 * (-2828) + O(5^22)

      *names* affects how the prime is printed. ::

        sage: T = Zp(5, print_mode='val-unit', names='pi'); a = T(70700); a
        pi^2 * 2828 + O(pi^22)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` and
      ``'fixed-mod'`` types and ``True`` for all other types. ::

        sage: Zp(5, print_mode='val-unit', show_prec=False)(30)
        5 * 6

      *print_max_terms*, *print_sep* and *print_alphabet* have no effect.

      Equality again depends on the printing options::

        sage: R == S, R == T, S == T
        (False, False, False)

    3. **terse**: elements are displayed as an integer in base 10::

        sage: R = Zp(5, print_mode='terse'); a = R(70700); a
        70700 + O(5^22)
        sage: b = R(-70700); b
        2384185790944925 + O(5^22)

      *print_pos* controls whether to use a balanced representation or not. ::

        sage: S = Zp(5, print_mode='terse', print_pos=False); b = S(-70700); b
        -70700 + O(5^22)

      *name* affects how the name is printed.  Note that this interacts
      with the choice of shorter string for denominators. ::

        sage: T.<unif> = Zp(5, print_mode='terse'); c = T(-707); c
        95367431639918 + O(unif^20)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` and
      ``'fixed-mod'`` types and ``True`` for all other types. ::

        sage: Zp(5, print_mode='terse', show_prec=False)(30)
        30

      *print_max_terms*, *print_sep* and *print_alphabet* have no effect.

      Equality depends on printing options::

        sage: R == S, R == T, S == T
        (False, False, False)

    4. **digits**: elements are displayed as a string of base `p` digits

      Restriction: you can only use the digits printing mode for small
      primes.  Namely, `p` must be less than the length of the alphabet
      tuple (default alphabet has length 62). ::

        sage: R = Zp(5, print_mode='digits'); a = R(70700); repr(a)
        '...4230300'
        sage: b = R(-70700); repr(b)
        '...4444444444444440214200'

      Note that it's not possible to read off the precision from the
      representation in this mode.

      *print_max_terms* limits the number of digits that are printed. ::

        sage: S = Zp(5, print_max_terms=4); S(-70700)
        2*5^2 + 4*5^3 + 5^4 + 2*5^5 + ... + O(5^22)

      *print_alphabet* controls the symbols used to substitute for digits
      greater than 9.  Defaults to
      ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z')::

        sage: T = Zp(5, print_mode='digits', print_alphabet=('1','2','3','4','5')); repr(T(-70700))
        '...5555555555555551325311'

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'dots'
      (or equivalently ``True``) or 'bigoh'.
      The default is ``False`` for the ``'floating-point'`` and
      ``'fixed-mod'`` types and ``True`` for all other types. ::

        sage: repr(Zp(5, 2, print_mode='digits', show_prec=True)(6))
        '...11'
        sage: repr(Zp(5, 2, print_mode='digits', show_prec='bigoh')(6))
        '11 + O(5^2)'

      *print_pos*, *name* and *print_sep* have no effect.

      Equality depends on printing options::

        sage: R == S, R == T, S == T
        (False, False, False)

    5. **bars**: elements are displayed as a string of base `p` digits
       with separators::

        sage: R = Zp(5, print_mode='bars'); a = R(70700); repr(a)
        '...4|2|3|0|3|0|0'
        sage: b = R(-70700); repr(b)
        '...4|4|4|4|4|4|4|4|4|4|4|4|4|4|4|0|2|1|4|2|0|0'

      Again, note that it's not possible to read off the precision from
      the representation in this mode.

      *print_pos* controls whether the digits can be negative. ::

        sage: S = Zp(5, print_mode='bars',print_pos=False); b = S(-70700); repr(b)
        '...-1|0|2|2|-1|2|0|0'

      *print_max_terms* limits the number of digits that are printed. ::

        sage: T = Zp(5, print_max_terms=4); T(-70700)
        2*5^2 + 4*5^3 + 5^4 + 2*5^5 + ... + O(5^22)

      *print_sep* controls the separation character. ::

        sage: U = Zp(5, print_mode='bars', print_sep=']['); a = U(70700); repr(a)
        '...4][2][3][0][3][0][0'

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'dots'
      (or equivalently ``True``) or 'bigoh'.
      The default is ``False`` for the ``'floating-point'`` and
      ``'fixed-mod'`` types and ``True`` for all other types. ::

        sage: repr(Zp(5, 2, print_mode='bars', show_prec=True)(6))
        '...1|1'
        sage: repr(Zp(5, 2, print_mode='bars', show_prec=False)(6))
        '1|1'

      *name* and *print_alphabet* have no effect.

      Equality depends on printing options::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    EXAMPLES:

    We allow non-prime `p`, but only if ``check = False``.  Note that some
    features will not work. ::

        sage: K = Zp(15, check=False); a = K(999); a
        9 + 6*15 + 4*15^2 + O(15^20)

    We create rings with various parameters::

        sage: Zp(7)
        7-adic Ring with capped relative precision 20
        sage: Zp(9)
        Traceback (most recent call last):
        ...
        ValueError: p must be prime
        sage: Zp(17, 5)
        17-adic Ring with capped relative precision 5
        sage: Zp(17, 5)(-1)
        16 + 16*17 + 16*17^2 + 16*17^3 + 16*17^4 + O(17^5)

    It works even with a fairly huge cap::

        sage: Zp(next_prime(10^50), 100000)
        100000000000000000000000000000000000000000000000151-adic Ring with capped relative precision 100000

    We create each type of ring::

        sage: Zp(7, 20, 'capped-rel')
        7-adic Ring with capped relative precision 20
        sage: Zp(7, 20, 'fixed-mod')
        7-adic Ring of fixed modulus 7^20
        sage: Zp(7, 20, 'capped-abs')
        7-adic Ring with capped absolute precision 20

    We create a capped relative ring with each print mode::

        sage: k = Zp(7, 8, print_mode='series'); k
        7-adic Ring with capped relative precision 8
        sage: k(7*(19))
        5*7 + 2*7^2 + O(7^9)
        sage: k(7*(-19))
        2*7 + 4*7^2 + 6*7^3 + 6*7^4 + 6*7^5 + 6*7^6 + 6*7^7 + 6*7^8 + O(7^9)

    ::

        sage: k = Zp(7, print_mode='val-unit'); k
        7-adic Ring with capped relative precision 20
        sage: k(7*(19))
        7 * 19 + O(7^21)
        sage: k(7*(-19))
        7 * 79792266297611982 + O(7^21)

    ::

        sage: k = Zp(7, print_mode='terse'); k
        7-adic Ring with capped relative precision 20
        sage: k(7*(19))
        133 + O(7^21)
        sage: k(7*(-19))
        558545864083283874 + O(7^21)

    Note that `p`-adic rings are cached (via weak references)::

        sage: a = Zp(7); b = Zp(7)
        sage: a is b
        True

    We create some elements in various rings::

        sage: R = Zp(5); a = R(4); a
        4 + O(5^20)
        sage: S = Zp(5, 10, type = 'capped-abs'); b = S(2); b
        2 + O(5^10)
        sage: a + b
        1 + 5 + O(5^10)
    """
    def create_key(self, p, prec = None, type = 'capped-rel', print_mode = None,
                   names = None, ram_name = None, print_pos = None, print_sep = None, print_alphabet = None,
                   print_max_terms = None, show_prec = None, check = True,
                   label = None):
        r"""
        Creates a key from input parameters for ``Zp``.

        See the documentation for ``Zp`` for more information.

        TESTS::

            sage: Zp.create_key(5,40)
            (5, 40, 'capped-rel', 'series', '5', True, '|', (), -1, 'bigoh', None)
            sage: Zp.create_key(5,40,print_mode='digits')
            (5,
             40,
             'capped-rel',
             'digits',
             '5',
             True,
             '|',
             ('0', '1', '2', '3', '4'),
             -1,
             'dots',
             None)
        """
        if isinstance(names, (int, Integer)):
            # old pickle; names is what used to be halt.
            names = ram_name
            ram_name = print_pos
            print_pos = print_sep
            print_alphabet = print_max_terms
            print_max_terms = check
            check = True
        if label is not None and type not in ['lattice-cap','lattice-float']:
            raise ValueError("label keyword only supported for lattice precision")
        return get_key_base(p, prec, type, print_mode, names, ram_name, print_pos, print_sep, print_alphabet,
                            print_max_terms, show_prec, check,
                            ['capped-rel', 'fixed-mod', 'capped-abs', 'floating-point', 'lattice-cap', 'lattice-float', 'relaxed'],
                            label=label)

    def create_object(self, version, key):
        r"""
        Creates an object using a given key.

        See the documentation for ``Zp`` for more information.

        TESTS::

            sage: Zp.create_object((3,4,2),(5, 41, 'capped-rel', 'series', '5', True, '|', (), -1))
            5-adic Ring with capped relative precision 41
        """
        if (version[0] < 3 or (len(version) > 1 and version[0] == 3 and version[1] < 2) or
            (len(version) > 2 and version[0] == 3 and version[1] == 2 and version[2] < 3)):
            p, prec, type, print_mode, name = key
            print_pos, print_sep, print_alphabet, print_max_terms = None, None, None, None
        elif version[0] < 8:
            p, prec, type, print_mode, name, print_pos, print_sep, print_alphabet, print_max_terms = key
            show_prec = None
            label = None
        else:
            p, prec, type, print_mode, name, print_pos, print_sep, print_alphabet, print_max_terms, show_prec, label = key
        if (version[0] < 4 or (len(version) > 1 and version[0] == 4 and version[1] < 5) or
            (len(version) > 2 and version[0] == 4 and version[1] == 5 and version[2] < 3)):
            # keys changed in order to reduce irrelevant duplications: e.g. two Zps with print_mode 'series'
            # that are identical except for different 'print_alphabet' now return the same object.
            key = get_key_base(p, prec, type, print_mode, name, None, print_pos, print_sep, print_alphabet,
                               print_max_terms, None, False, ['capped-rel', 'fixed-mod', 'capped-abs', 'lattice-cap', 'lattice-float', 'relaxed'])
            try:
                obj = self._cache[version, key]()
                if obj is not None:
                    return obj
            except KeyError:
                pass
            p, prec, type, print_mode, name, print_pos, print_sep, print_alphabet, print_max_terms, show_prec, label = key
        if type == 'capped-rel':
            return pAdicRingCappedRelative(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                     'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name)
        elif type == 'fixed-mod':
            return pAdicRingFixedMod(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                               'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name)
        elif type == 'capped-abs':
            return pAdicRingCappedAbsolute(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                     'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name)
        elif type == 'floating-point':
            return pAdicRingFloatingPoint(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                     'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name)
        elif type == 'relaxed':
            return pAdicRingRelaxed(p, prec, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                           'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name)
        elif type[:8] == 'lattice-':
            subtype = type[8:]
            return pAdicRingLattice(p, prec, subtype, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                                                       'ram_name': name, 'max_ram_terms': print_max_terms, 'show_prec': show_prec}, name, label)
        else:
            raise ValueError("unexpected type")

Zp = Zp_class("Zp")


######################################################
# Zq -- unramified extensions
######################################################

def Zq(q, prec = None, type = 'capped-rel', modulus = None, names=None,
          print_mode=None, ram_name = None, res_name = None, print_pos = None,
       print_sep = None, print_max_ram_terms = None,
       print_max_unram_terms = None, print_max_terse_terms = None, show_prec = None, check = True, implementation = 'FLINT'):
    r"""
    Given a prime power `q = p^n`, return the unique unramified
    extension of `\ZZ_p` of degree `n`.

    INPUT:

    - ``q`` -- integer, list or tuple: the prime power in `\QQ_q`.  Or a
      factorization object, single element list ``[(p, n)]`` where ``p`` is
      a prime and ``n`` a positive integer, or the pair ``(p, n)``.

    - ``prec`` -- integer (default: ``20``) the precision cap of the
      field.  Individual elements keep track of their own precision.
      See TYPES and PRECISION below.

    - ``type`` -- string (default: ``'capped-rel'``) Valid types are
      ``'capped-abs'``, ``'capped-rel'``, ``'fixed-mod'``, and
      ``'floating-point'``.  See TYPES and PRECISION below

    - modulus -- polynomial (default None) A polynomial defining an
      unramified extension of `\ZZ_p`.  See MODULUS below.

    - ``names`` -- string or tuple (``None`` is only allowed when
      `q=p`).  The name of the generator, reducing to a generator of
      the residue field.

    - ``print_mode`` -- string (default: ``None``).  Valid modes are ``'series'``,
      ``'val-unit'``, ``'terse'``, and ``'bars'``. See PRINTING below.

    - ``ram_name`` -- string (defaults to string representation of `p` if
      None).  ``ram_name`` controls how the prime is printed. See PRINTING
      below.

    - ``res_name`` -- string (defaults to ``None``, which corresponds
      to adding a ``'0'`` to the end of the name).  Controls how
      elements of the reside field print.

    - ``print_pos`` -- bool (default ``None``) Whether to only use
      positive integers in the representations of elements. See
      PRINTING below.

    - ``print_sep`` -- string (default ``None``) The separator
      character used in the ``'bars'`` mode. See PRINTING below.

    - ``print_max_ram_terms`` -- integer (default ``None``) The maximum
      number of powers of `p` shown.  See PRINTING below.

    - ``print_max_unram_terms`` -- integer (default ``None``) The
      maximum number of entries shown in a coefficient of `p`.  See
      PRINTING below.

    - ``print_max_terse_terms`` -- integer (default ``None``) The maximum
      number of terms in the polynomial representation of an element
      (using ``'terse'``).  See PRINTING below.

    - ``show_prec`` -- bool (default ``None``) Whether to show the precision
      for elements.  See PRINTING below.

    - ``check`` -- bool (default ``True``) whether to check inputs.

    - ``implementation`` -- string (default ``FLINT``) which
      implementation to use.  ``NTL`` is the other option.

    OUTPUT:

    - The corresponding unramified `p`-adic ring.

    TYPES AND PRECISION:


    There are two types of precision for a `p`-adic element.  The first
    is relative precision (default), which gives the number of known `p`-adic
    digits::

        sage: R.<a> = Zq(25, 20, 'capped-rel', print_mode='series'); b = 25*a; b
        a*5^2 + O(5^22)
        sage: b.precision_relative()
        20

    The second type of precision is absolute precision, which gives
    the power of `p` that this element is defined modulo::

        sage: b.precision_absolute()
        22

    There are many types of `p`-adic rings: capped relative rings
    (``type='capped-rel'``), capped absolute rings
    (``type='capped-abs'``), fixed modulus rings (``type='fixed-mod'``),
    and floating point rings (``type='floating-point'``).

    In the capped relative case, the relative precision of an element
    is restricted to be at most a certain value, specified at the
    creation of the field.  Individual elements also store their own
    precision, so the effect of various arithmetic operations on
    precision is tracked.  When you cast an exact element into a
    capped relative field, it truncates it to the precision cap of the
    field. ::

        sage: R.<a> = Zq(9, 5, 'capped-rel', print_mode='series'); b = (1+2*a)^4; b
        2 + (2*a + 2)*3 + (2*a + 1)*3^2 + O(3^5)
        sage: c = R(3249); c
        3^2 + 3^4 + 3^5 + 3^6 + O(3^7)
        sage: b + c
        2 + (2*a + 2)*3 + (2*a + 2)*3^2 + 3^4 + O(3^5)

    One can invert non-units: the result is in the fraction field. ::

        sage: d = ~(3*b+c); d
        2*3^-1 + (a + 1) + (a + 1)*3 + a*3^3 + O(3^4)
        sage: d.parent()
        3-adic Unramified Extension Field in a defined by x^2 + 2*x + 2

    The capped absolute case is the same as the capped relative case,
    except that the cap is on the absolute precision rather than the
    relative precision. ::

        sage: R.<a> = Zq(9, 5, 'capped-abs', print_mode='series'); b = 3*(1+2*a)^4; b
        2*3 + (2*a + 2)*3^2 + (2*a + 1)*3^3 + O(3^5)
        sage: c = R(3249); c
        3^2 + 3^4 + O(3^5)
        sage: b*c
        2*3^3 + (2*a + 2)*3^4 + O(3^5)
        sage: b*c >> 1
        2*3^2 + (2*a + 2)*3^3 + O(3^4)

    The fixed modulus case is like the capped absolute, except that
    individual elements don't track their precision. ::

        sage: R.<a> = Zq(9, 5, 'fixed-mod', print_mode='series'); b = 3*(1+2*a)^4; b
        2*3 + (2*a + 2)*3^2 + (2*a + 1)*3^3
        sage: c = R(3249); c
        3^2 + 3^4
        sage: b*c
        2*3^3 + (2*a + 2)*3^4
        sage: b*c >> 1
        2*3^2 + (2*a + 2)*3^3

    The floating point case is similar to the fixed modulus type
    in that elements do not trac their own precision.  However, relative
    precision is truncated with each operation rather than absolute precision.

    MODULUS:

    The modulus needs to define an unramified extension of `\ZZ_p`: when it
    is reduced to a polynomial over `\GF{p}` it should be irreducible.

    The modulus can be given in a number of forms.

    1. A **polynomial**.

      The base ring can be `\ZZ`, `\QQ`, `\ZZ_p`, `\GF{p}`, or anything that can
      be converted to `\ZZ_p`. ::

        sage: P.<x> = ZZ[]
        sage: R.<a> = Zq(27, modulus = x^3 + 2*x + 1); R.modulus()
        (1 + O(3^20))*x^3 + O(3^20)*x^2 + (2 + O(3^20))*x + 1 + O(3^20)
        sage: P.<x> = QQ[]
        sage: S.<a> = Zq(27, modulus = x^3 + 2/7*x + 1)
        sage: P.<x> = Zp(3)[]
        sage: T.<a> = Zq(27, modulus = x^3 + 2*x + 1)
        sage: P.<x> = Qp(3)[]
        sage: U.<a> = Zq(27, modulus = x^3 + 2*x + 1)
        sage: P.<x> = GF(3)[]
        sage: V.<a> = Zq(27, modulus = x^3 + 2*x + 1)

      Which form the modulus is given in has no effect on the unramified
      extension produced::

        sage: R == S, R == T, T == U, U == V
        (False, True, True, False)

      unless the modulus is different, or the precision of the modulus
      differs.  In the case of ``V``, the modulus is only given to precision
      ``1``, so the resulting field has a precision cap of ``1``. ::

        sage: V.precision_cap()
        1
        sage: U.precision_cap()
        20
        sage: P.<x> = Zp(3)[]
        sage: modulus = x^3 + (2 + O(3^7))*x + (1 + O(3^10))
        sage: modulus
        (1 + O(3^20))*x^3 + (2 + O(3^7))*x + 1 + O(3^10)
        sage: W.<a> = Zq(27, modulus = modulus); W.precision_cap()
        7

    2. The modulus can also be given as a **symbolic expression**. ::

        sage: x = var('x')
        sage: X.<a> = Zq(27, modulus = x^3 + 2*x + 1); X.modulus()
        (1 + O(3^20))*x^3 + O(3^20)*x^2 + (2 + O(3^20))*x + 1 + O(3^20)
        sage: X == R
        True

      By default, the polynomial chosen is the standard lift of the
      generator chosen for `\GF{q}`. ::

        sage: GF(125, 'a').modulus()
        x^3 + 3*x + 3
        sage: Y.<a> = Zq(125); Y.modulus()
        (1 + O(5^20))*x^3 + O(5^20)*x^2 + (3 + O(5^20))*x + 3 + O(5^20)

      However, you can choose another polynomial if desired (as long as
      the reduction to `\GF{p}[x]` is irreducible). ::

        sage: P.<x> = ZZ[]
        sage: Z.<a> = Zq(125, modulus = x^3 + 3*x^2 + x + 1); Z.modulus()
        (1 + O(5^20))*x^3 + (3 + O(5^20))*x^2 + (1 + O(5^20))*x + 1 + O(5^20)
        sage: Y == Z
        False

    PRINTING:

    There are many different ways to print `p`-adic elements.  The way
    elements of a given field print is controlled by options passed in
    at the creation of the field.  There are four basic printing modes
    (``'series'``, ``'val-unit'``, ``'terse'`` and ``'bars'``; ``'digits'`` is not available), as
    well as various options that either hide some information in the
    print representation or sometimes make print representations more
    compact.  Note that the printing options affect whether different
    `p`-adic fields are considered equal.

    1. **series**: elements are displayed as series in `p`. ::

        sage: R.<a> = Zq(9, 20, 'capped-rel', print_mode='series'); (1+2*a)^4
        2 + (2*a + 2)*3 + (2*a + 1)*3^2 + O(3^20)
        sage: -3*(1+2*a)^4
        3 + a*3^2 + 3^3 + (2*a + 2)*3^4 + (2*a + 2)*3^5 + (2*a + 2)*3^6 + (2*a + 2)*3^7 + (2*a + 2)*3^8 + (2*a + 2)*3^9 + (2*a + 2)*3^10 + (2*a + 2)*3^11 + (2*a + 2)*3^12 + (2*a + 2)*3^13 + (2*a + 2)*3^14 + (2*a + 2)*3^15 + (2*a + 2)*3^16 + (2*a + 2)*3^17 + (2*a + 2)*3^18 + (2*a + 2)*3^19 + (2*a + 2)*3^20 + O(3^21)
        sage: b = ~(3*a+18); b
        (a + 2)*3^-1 + 1 + 2*3 + (a + 1)*3^2 + 3^3 + 2*3^4 + (a + 1)*3^5 + 3^6 + 2*3^7 + (a + 1)*3^8 + 3^9 + 2*3^10 + (a + 1)*3^11 + 3^12 + 2*3^13 + (a + 1)*3^14 + 3^15 + 2*3^16 + (a + 1)*3^17 + 3^18 + O(3^19)
        sage: b.parent() is R.fraction_field()
        True

      *print_pos* controls whether negatives can be used in the
      coefficients of powers of `p`. ::

        sage: S.<b> = Zq(9, print_mode='series', print_pos=False); (1+2*b)^4
        -1 - b*3 - 3^2 + (b + 1)*3^3 + O(3^20)
        sage: -3*(1+2*b)^4
        3 + b*3^2 + 3^3 + (-b - 1)*3^4 + O(3^21)

      *ram_name* controls how the prime is printed. ::

        sage: T.<d> = Zq(9, print_mode='series', ram_name='p'); 3*(1+2*d)^4
        2*p + (2*d + 2)*p^2 + (2*d + 1)*p^3 + O(p^21)

      *print_max_ram_terms* limits the number of powers of `p` that
      appear. ::

        sage: U.<e> = Zq(9, print_mode='series', print_max_ram_terms=4); repr(-3*(1+2*e)^4)
        '3 + e*3^2 + 3^3 + (2*e + 2)*3^4 + ... + O(3^21)'

      *print_max_unram_terms* limits the number of terms that appear in a
      coefficient of a power of `p`. ::

        sage: V.<f> = Zq(128, prec = 8, print_mode='series'); repr((1+f)^9)
        '(f^3 + 1) + (f^5 + f^4 + f^3 + f^2)*2 + (f^6 + f^5 + f^4 + f + 1)*2^2 + (f^5 + f^4 + f^2 + f + 1)*2^3 + (f^6 + f^5 + f^4 + f^3 + f^2 + f + 1)*2^4 + (f^5 + f^4)*2^5 + (f^6 + f^5 + f^4 + f^3 + f + 1)*2^6 + (f + 1)*2^7 + O(2^8)'
        sage: V.<f> = Zq(128, prec = 8, print_mode='series', print_max_unram_terms = 3); repr((1+f)^9)
        '(f^3 + 1) + (f^5 + f^4 + ... + f^2)*2 + (f^6 + f^5 + ... + 1)*2^2 + (f^5 + f^4 + ... + 1)*2^3 + (f^6 + f^5 + ... + 1)*2^4 + (f^5 + f^4)*2^5 + (f^6 + f^5 + ... + 1)*2^6 + (f + 1)*2^7 + O(2^8)'
        sage: V.<f> = Zq(128, prec = 8, print_mode='series', print_max_unram_terms = 2); repr((1+f)^9)
        '(f^3 + 1) + (f^5 + ... + f^2)*2 + (f^6 + ... + 1)*2^2 + (f^5 + ... + 1)*2^3 + (f^6 + ... + 1)*2^4 + (f^5 + f^4)*2^5 + (f^6 + ... + 1)*2^6 + (f + 1)*2^7 + O(2^8)'
        sage: V.<f> = Zq(128, prec = 8, print_mode='series', print_max_unram_terms = 1); repr((1+f)^9)
        '(f^3 + ...) + (f^5 + ...)*2 + (f^6 + ...)*2^2 + (f^5 + ...)*2^3 + (f^6 + ...)*2^4 + (f^5 + ...)*2^5 + (f^6 + ...)*2^6 + (f + ...)*2^7 + O(2^8)'
        sage: V.<f> = Zq(128, prec = 8, print_mode='series', print_max_unram_terms = 0); repr((1+f)^9 - 1 - f^3)
        '(...)*2 + (...)*2^2 + (...)*2^3 + (...)*2^4 + (...)*2^5 + (...)*2^6 + (...)*2^7 + O(2^8)'

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` and
      ``'fixed-mod'`` types and ``True`` for all other types. ::

        sage: U.<e> = Zq(9, 2, show_prec=False); repr(-3*(1+2*e)^4)
        '3 + e*3^2'

      *print_sep* and *print_max_terse_terms* have no effect.

      Note that print options affect equality::

        sage: R == S, R == T, R == U, R == V, S == T, S == U, S == V, T == U, T == V, U == V
        (False, False, False, False, False, False, False, False, False, False)

    2. **val-unit**: elements are displayed as `p^k u`::

        sage: R.<a> = Zq(9, 7, print_mode='val-unit'); b = (1+3*a)^9 - 1; b
        3^3 * (15 + 64*a) + O(3^7)
        sage: ~b
        3^-3 * (41 + a) + O(3)

      *print_pos* controls whether to use a balanced representation or
      not. ::

        sage: S.<a> = Zq(9, 7, print_mode='val-unit', print_pos=False); b = (1+3*a)^9 - 1; b
        3^3 * (15 - 17*a) + O(3^7)
        sage: ~b
        3^-3 * (-40 + a) + O(3)

      *ram_name* affects how the prime is printed. ::

        sage: A.<x> = Zp(next_prime(10^6), print_mode='val-unit')[]
        sage: T.<a> = Zq(next_prime(10^6)^3, 4, print_mode='val-unit', ram_name='p', modulus=x^3+385831*x^2+106556*x+321036); b = (next_prime(10^6)^2*(a^2 + a - 4)^4); b
        p^2 * (87996187118837557387483 + 246348888344392418464080*a + 1353538653775332610349*a^2) + O(p^6)
        sage: b * (a^2 + a - 4)^-4
        p^2 * 1 + O(p^6)

      *print_max_terse_terms* controls how many terms of the polynomial
      appear in the unit part. ::

        sage: U.<a> = Zq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3); b = (17*(a^3-a+14)^6); b
        17 * (12131797 + 12076378*a + 10809706*a^2 + ...) + O(17^7)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` and
      ``'fixed-mod'`` types and ``True`` for all other types. ::

        sage: U.<e> = Zq(9, 2, print_mode='val-unit', show_prec=False); repr(-3*(1+2*e)^4)
        '3 * (1 + 3*e)'

      *print_sep*, *print_max_ram_terms* and *print_max_unram_terms* have no effect.

      Equality again depends on the printing options::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    3. **terse**: elements are displayed as a polynomial of degree less
       than the degree of the extension. ::

        sage: R.<a> = Zq(125, print_mode='terse')
        sage: (a+5)^177
        68210977979428 + 90313850704069*a + 73948093055069*a^2 + O(5^20)
        sage: (a/5+1)^177
        68210977979428/5^177 + 90313850704069/5^177*a + 73948093055069/5^177*a^2 + O(5^-157)

      Note that in this last computation, you get one fewer `p`-adic digit
      than one might expect.  This is because ``R`` is capped absolute, and
      thus 5 is cast in with relative precision 19.

      As of version 3.3, if coefficients of the polynomial are
      non-integral, they are always printed with an explicit power of `p`
      in the denominator. ::

        sage: 5*a + a^2/25
        5*a + 1/5^2*a^2 + O(5^18)

      *print_pos* controls whether to use a balanced representation or
      not. ::

        sage: (a-5)^6
        22864 + 95367431627998*a + 8349*a^2 + O(5^20)
        sage: S.<a> = Zq(125, print_mode='terse', print_pos=False); b = (a-5)^6; b
        22864 - 12627*a + 8349*a^2 + O(5^20)
        sage: (a - 1/5)^6
        -20624/5^6 + 18369/5^5*a + 1353/5^3*a^2 + O(5^14)

      *ram_name* affects how the prime is printed. ::

        sage: T.<a> = Zq(125, print_mode='terse', ram_name='p'); (a - 1/5)^6
        95367431620001/p^6 + 18369/p^5*a + 1353/p^3*a^2 + O(p^14)

      *print_max_terse_terms* controls how many terms of the polynomial
      are shown. ::

        sage: U.<a> = Zq(625, print_mode='terse', print_max_terse_terms=2); (a-1/5)^6
        106251/5^6 + 49994/5^5*a + ... + O(5^14)

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'bigoh'
      (or equivalently ``True``).
      The default is ``False`` for the ``'floating-point'`` and
      ``'fixed-mod'`` types and ``True`` for all other types. ::

        sage: U.<e> = Zq(9, 2, print_mode='terse', show_prec=False); repr(-3*(1+2*e)^4)
        '3 + 9*e'

      *print_sep*, *print_max_ram_terms* and *print_max_unram_terms* have no
      effect.

      Equality again depends on the printing options::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    4. **digits**: This print mode is not available when the residue
       field is not prime.  It might make sense to have a dictionary
       for small fields, but this isn't implemented.

    5. **bars**: elements are displayed in a similar fashion to series,
       but more compactly. ::

        sage: R.<a> = Zq(125); (a+5)^6
        (4*a^2 + 3*a + 4) + (3*a^2 + 2*a)*5 + (a^2 + a + 1)*5^2 + (3*a + 2)*5^3 + (3*a^2 + a + 3)*5^4 + (2*a^2 + 3*a + 2)*5^5 + O(5^20)
        sage: R.<a> = Zq(125, print_mode='bars', prec=8); repr((a+5)^6)
        '...[2, 3, 2]|[3, 1, 3]|[2, 3]|[1, 1, 1]|[0, 2, 3]|[4, 3, 4]'
        sage: repr((a-5)^6)
        '...[0, 4]|[1, 4]|[2, 0, 2]|[1, 4, 3]|[2, 3, 1]|[4, 4, 3]|[2, 4, 4]|[4, 3, 4]'

      Note that it's not possible to read off the precision from the
      representation in this mode. ::

        sage: b = a + 3; repr(b)
        '...[3, 1]'
        sage: c = a + R(3, 4); repr(c)
        '...[3, 1]'
        sage: b.precision_absolute()
        8
        sage: c.precision_absolute()
        4

      *print_pos* controls whether the digits can be negative. ::

        sage: S.<a> = Zq(125, print_mode='bars', print_pos=False); repr((a-5)^6)
        '...[1, -1, 1]|[2, 1, -2]|[2, 0, -2]|[-2, -1, 2]|[0, 0, -1]|[-2]|[-1, -2, -1]'
        sage: repr((a-1/5)^6)
        '...[0, 1, 2]|[-1, 1, 1]|.|[-2, -1, -1]|[2, 2, 1]|[0, 0, -2]|[0, -1]|[0, -1]|[1]'

      *print_max_ram_terms* controls the maximum number of "digits" shown.
      Note that this puts a cap on the relative precision, not the
      absolute precision. ::

        sage: T.<a> = Zq(125, print_max_ram_terms=3, print_pos=False); (a-5)^6
        (-a^2 - 2*a - 1) - 2*5 - a^2*5^2 + ... + O(5^20)
        sage: 5*(a-5)^6 + 50
        (-a^2 - 2*a - 1)*5 - a^2*5^3 + (2*a^2 - a - 2)*5^4 + ... + O(5^21)
        sage: (a-1/5)^6
        5^-6 - a*5^-5 - a*5^-4 + ... + O(5^14)

      *print_sep* controls the separating character (``'|'`` by default). ::

        sage: U.<a> = Zq(625, print_mode='bars', print_sep=''); b = (a+5)^6; repr(b)
        '...[0, 1][4, 0, 2][3, 2, 2, 3][4, 2, 2, 4][0, 3][1, 1, 3][3, 1, 4, 1]'

      *print_max_unram_terms* controls how many terms are shown in each
      ``'digit'``::

        sage: with local_print_mode(U, {'max_unram_terms': 3}): repr(b)
        '...[0, 1][4,..., 0, 2][3,..., 2, 3][4,..., 2, 4][0, 3][1,..., 1, 3][3,..., 4, 1]'
        sage: with local_print_mode(U, {'max_unram_terms': 2}): repr(b)
        '...[0, 1][4,..., 2][3,..., 3][4,..., 4][0, 3][1,..., 3][3,..., 1]'
        sage: with local_print_mode(U, {'max_unram_terms': 1}): repr(b)
        '...[..., 1][..., 2][..., 3][..., 4][..., 3][..., 3][..., 1]'
        sage: with local_print_mode(U, {'max_unram_terms':0}): repr(b-75*a)
        '...[...][...][...][...][][...][...]'

      *show_prec* determines how the precision is printed.
      It can be either 'none' (or equivalently ``False``), 'dots'
      (or equivalently ``True``) or 'bigoh'.
      The default is ``False`` for the ``'floating-point'`` and
      ``'fixed-mod'`` types and ``True`` for all other types. ::

        sage: U.<e> = Zq(9, 2, print_mode='bars', show_prec='bigoh'); repr(-3*(1+2*e)^4)
        '[0, 1]|[1]|[] + O(3^3)'

      *ram_name* and *print_max_terse_terms* have no effect.

      Equality depends on printing options::

        sage: R == S, R == T, R == U, S == T, S == U, T == U
        (False, False, False, False, False, False)

    EXAMPLES:

    Unlike for ``Zp``, you can't create ``Zq(N)`` when ``N`` is not a prime power.

    However, you can use ``check=False`` to pass in a pair in order to not
    have to factor.  If you do so, you need to use names explicitly
    rather than the ``R.<a>`` syntax. ::

        sage: p = next_prime(2^123)
        sage: k = Zp(p)
        sage: R.<x> = k[]
        sage: K = Zq([(p, 5)], modulus=x^5+x+4, names='a', ram_name='p', print_pos=False, check=False)
        sage: K.0^5
        (-a - 4) + O(p^20)

    In tests on sage.math, the creation of ``K`` as above took an average
    of 1.58ms, while::

        sage: K = Zq(p^5, modulus=x^5+x+4, names='a', ram_name='p', print_pos=False, check=True)

    took an average of 24.5ms.  Of course, with smaller primes these
    savings disappear.

    TESTS::

        sage: R = Zq([(5,3)], names="alpha"); R
        5-adic Unramified Extension Ring in alpha defined by x^3 + 3*x + 3
        sage: Zq((5, 3), names="alpha") is R
        True
        sage: Zq(125.factor(), names="alpha") is R
        True

    """
    if check:
        if isinstance(q, Factorization) or isinstance(q, (list, tuple)):
            if not isinstance(q, Factorization) and len(q) == 2:
                F = [(Integer(q[0]), Integer(q[1]))]
            else:
                if len(q) != 1:
                    raise ValueError("q must be a prime power")
                if len(q[0]) != 2:
                    raise ValueError("q must have shape [(p, k)]")
                F = [(Integer(q[0][0]), Integer(q[0][1]))]
            if not F[0][0].is_prime() or F[0][1] <= 0:
                raise ValueError("q must be a prime power")
            q = F[0][0]**F[0][1]
        else:
            q = Integer(q)
            F = q.factor()
            if len(F) != 1:
                raise ValueError("q must be a prime power")
        if prec is not None and not isinstance(prec, Integer):
            prec = Integer(prec)
        if isinstance(names, (list, tuple)):
            names = names[0]
        from sage.structure.element import Expression
        if not (modulus is None or is_Polynomial(modulus) or isinstance(modulus, Expression)):
            raise TypeError("modulus must be a polynomial")
        if names is not None and not isinstance(names, str):
            names = str(names)
            #raise TypeError, "names must be a string"
        q = Integer(q)
        F = q.factor()
        if len(F) != 1:
            raise ValueError("q must be a prime power")
    else:
        F = q
        q = F[0][0]**F[0][1]
    base = Zp(p=F[0][0], prec=prec, type=type, print_mode=print_mode, names=ram_name,
              print_pos=print_pos, print_sep=print_sep, print_max_terms=print_max_ram_terms,
              show_prec=show_prec, check=False)
    if F[0][1] == 1:
        return base
    elif names is None:
        raise TypeError("You must specify the name of the generator.")
    if res_name is None:
        res_name = names + '0'
    if modulus is None:
        from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
        if ram_name is None:
            ram_name = str(F[0][0])
        modulus = GF(q, res_name).modulus().change_ring(ZZ)
    return ExtensionFactory(base=base, modulus=modulus, prec=prec, print_mode=print_mode,
                            names=names, res_name=res_name, ram_name=ram_name, print_pos=print_pos,
                            print_sep=print_sep, print_max_ram_terms=print_max_ram_terms,
                            print_max_unram_terms=print_max_unram_terms,
                            print_max_terse_terms=print_max_terse_terms, show_prec=show_prec, check=check,
                            unram=True, implementation=implementation)

######################################################
# Short constructor names for different types
######################################################

def ZpCR(p, prec = None, *args, **kwds):
    r"""
    A shortcut function to create capped relative `p`-adic rings.

    Same functionality as ``Zp``.  See documentation for ``Zp`` for a
    description of the input parameters.

    EXAMPLES::

        sage: ZpCR(5, 40)
        5-adic Ring with capped relative precision 40
    """
    return Zp(p, prec, 'capped-rel', *args, **kwds)

def ZpCA(p, prec = None, *args, **kwds):
    r"""
    A shortcut function to create capped absolute `p`-adic rings.

    See documentation for :func:`Zp` for a description of the input parameters.

    EXAMPLES::

        sage: ZpCA(5, 40)
        5-adic Ring with capped absolute precision 40
    """
    return Zp(p, prec, 'capped-abs', *args, **kwds)

def ZpFM(p, prec = None, *args, **kwds):
    r"""
    A shortcut function to create fixed modulus `p`-adic rings.

    See documentation for :func:`Zp` for a description of the input parameters.

    EXAMPLES::

        sage: ZpFM(5, 40)
        5-adic Ring of fixed modulus 5^40
    """
    return Zp(p, prec, 'fixed-mod', *args, **kwds)

def ZpFP(p, prec = None, *args, **kwds):
    r"""
    A shortcut function to create floating point `p`-adic rings.

    Same functionality as ``Zp``.  See documentation for ``Zp`` for a
    description of the input parameters.

    EXAMPLES::

        sage: ZpFP(5, 40)
        5-adic Ring with floating precision 40
    """
    return Zp(p, prec, 'floating-point', *args, **kwds)

def ZqCR(q, prec = None, *args, **kwds):
    r"""
    A shortcut function to create capped relative unramified `p`-adic rings.

    Same functionality as ``Zq``.  See documentation for ``Zq`` for a
    description of the input parameters.

    EXAMPLES::

        sage: R.<a> = ZqCR(25, 40); R
        5-adic Unramified Extension Ring in a defined by x^2 + 4*x + 2
    """
    return Zq(q, prec, 'capped-rel', *args, **kwds)

def ZqCA(q, prec = None, *args, **kwds):
    r"""
    A shortcut function to create capped absolute unramified `p`-adic rings.

    See documentation for :func:`Zq` for a description of the input parameters.

    EXAMPLES::

        sage: R.<a> = ZqCA(25, 40); R
        5-adic Unramified Extension Ring in a defined by x^2 + 4*x + 2
    """
    return Zq(q, prec, 'capped-abs', *args, **kwds)

def ZqFM(q, prec = None, *args, **kwds):
    r"""
    A shortcut function to create fixed modulus unramified `p`-adic rings.

    See documentation for :func:`Zq` for a description of the input parameters.

    EXAMPLES::

        sage: R.<a> = ZqFM(25, 40); R
        5-adic Unramified Extension Ring in a defined by x^2 + 4*x + 2
    """
    return Zq(q, prec, 'fixed-mod', *args, **kwds)

def ZqFP(q, prec = None, *args, **kwds):
    r"""
    A shortcut function to create floating point unramified `p`-adic rings.

    Same functionality as ``Zq``.  See documentation for ``Zq`` for a
    description of the input parameters.

    EXAMPLES::

        sage: R.<a> = ZqFP(25, 40); R
        5-adic Unramified Extension Ring in a defined by x^2 + 4*x + 2
    """
    return Zq(q, prec, 'floating-point', *args, **kwds)

@experimental(23505)
def ZpLC(p, prec=None, *args, **kwds):
    r"""
    A shortcut function to create `p`-adic rings with lattice precision
    (precision is encoded by a lattice in a large vector space and tracked
    using automatic differentiation).

    See documentation for :func:`Zp` for a description of the input parameters.

    EXAMPLES:

    Below is a small demo of the features by this model of precision::

        sage: R = ZpLC(3, print_mode='terse')
        sage: R
        3-adic Ring with lattice-cap precision

        sage: x = R(1,10)

    Of course, when we multiply by 3, we gain one digit of absolute
    precision::

        sage: 3*x
        3 + O(3^11)

    The lattice precision machinery sees this even if we decompose
    the computation into several steps::

        sage: y = x+x
        sage: y
        2 + O(3^10)
        sage: x + y
        3 + O(3^11)

    The same works for the multiplication::

        sage: z = x^2
        sage: z
        1 + O(3^10)
        sage: x*z
        1 + O(3^11)

    This can be more surprising when we are working with elements given
    at different precisions::

        sage: R = ZpLC(2, print_mode='terse')
        sage: x = R(1,10)
        sage: y = R(1,5)
        sage: z = x+y; z
        2 + O(2^5)
        sage: t = x-y; t
        O(2^5)
        sage: z+t  # observe that z+t = 2*x
        2 + O(2^11)
        sage: z-t  # observe that z-t = 2*y
        2 + O(2^6)

        sage: x = R(28888,15)
        sage: y = R(204,10)
        sage: z = x/y; z
        242 + O(2^9)
        sage: z*y  # which is x
        28888 + O(2^15)

    The SOMOS sequence is the sequence defined by the recurrence:

    .. MATH::

        u_n = \frac {u_{n-1} u_{n-3} + u_{n-2}^2} {u_{n-4}}

    It is known for its numerical instability.
    On the one hand, one can show that if the initial values are
    invertible in `\ZZ_p` and known at precision `O(p^N)`
    then all the next terms of the SOMOS sequence will be known
    at the same precision as well.
    On the other hand, because of the division, when we unroll
    the recurrence, we loose a lot of precision. Observe::

        sage: R = Zp(2, 30, print_mode='terse')
        sage: a,b,c,d = R(1,15), R(1,15), R(1,15), R(3,15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        4 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        13 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        55 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        21975 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        6639 + O(2^13)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        7186 + O(2^13)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        569 + O(2^13)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        253 + O(2^13)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        4149 + O(2^13)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        2899 + O(2^12)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        3072 + O(2^12)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        349 + O(2^12)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        619 + O(2^12)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        243 + O(2^12)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        3 + O(2^2)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        2 + O(2^2)

    If instead, we use the lattice precision, everything goes well::

        sage: R = ZpLC(2, 30, print_mode='terse')
        sage: a,b,c,d = R(1,15), R(1,15), R(1,15), R(3,15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        4 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        13 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        55 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        21975 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        23023 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        31762 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        16953 + O(2^15)
        sage: a,b,c,d = b,c,d,(b*d+c*c)/a; print(d)
        16637 + O(2^15)

        sage: for _ in range(100):
        ....:     a,b,c,d = b,c,d,(b*d+c*c)/a
        sage: a
        15519 + O(2^15)
        sage: b
        32042 + O(2^15)
        sage: c
        17769 + O(2^15)
        sage: d
        20949 + O(2^15)

    ALGORITHM:

    The precision is global.
    It is encoded by a lattice in a huge vector space whose dimension
    is the number of elements having this parent. Precision is tracked
    using automatic differentiation techniques (see [CRV2014]_ and
    [CRV2018]_).

    Concretely, this precision datum is an instance of the class
    :class:`sage.rings.padic.lattice_precision.PrecisionLattice`.
    It is attached to the parent and is created at the same time
    as the parent.
    (It is actually a bit more subtle because two different parents
    may share the same instance; this happens for instance for a
    `p`-adic ring and its field of fractions.)

    This precision datum is accessible through the method :meth:`precision`::

        sage: R = ZpLC(5, print_mode='terse')
        sage: prec = R.precision()
        sage: prec
        Precision lattice on 0 objects

    This instance knows about all elements of the parent. It is
    automatically updated when a new element (of this parent) is
    created::

        sage: x = R(3513,10)
        sage: prec
        Precision lattice on 1 object
        sage: y = R(176,5)
        sage: prec
        Precision lattice on 2 objects
        sage: z = R.random_element()
        sage: prec
        Precision lattice on 3 objects

    The method :meth:`tracked_elements` provides the list of all
    tracked elements::

        sage: prec.tracked_elements()
        [3513 + O(5^10), 176 + O(5^5), ...]

    Similarly, when a variable is collected by the garbage collector,
    the precision lattice is updated. Note however that the update
    might be delayed. We can force it with the method :meth:`del_elements`::

        sage: z = 0
        sage: prec # random output, could be 2 objects if the garbage collector is fast
        Precision lattice on 3 objects
        sage: prec.del_elements()
        sage: prec
        Precision lattice on 2 objects

    The method :meth:`precision_lattice` returns (a matrix defining)
    the lattice that models the precision. Here we have::

        sage: prec.precision_lattice()
        [9765625       0]
        [      0    3125]

    Observe that `5^10 = 9765625` and `5^5 = 3125`.
    The above matrix then reflects the precision on `x` and `y`.

    Now, observe how the precision lattice changes while performing
    computations::

        sage: x, y = 3*x+2*y, 2*(x-y)
        sage: prec.del_elements()
        sage: prec.precision_lattice()
        [    3125 48825000]
        [       0 48828125]

    The matrix we get is no longer diagonal, meaning that some digits
    of precision are diffused among the two new elements `x` and `y`.
    They nevertheless show up when we compute for instance `x+y`::

        sage: x
        1516 + O(5^5)
        sage: y
        424 + O(5^5)
        sage: x+y
        17565 + O(5^11)

    These diffused digits of precision (which are tracked but
    do not appear on the printing) allow to be always sharp on
    precision.

    NOTE:

    Each elementary operation requires significant manipulations
    on the precision lattice and therefore is costly. Precisely:

    - The creation of a new element has a cost `O(n)` where `n`
      is the number of tracked elements.

    - The destruction of one element has a cost `O(m^2)` where
      `m` is the distance between the destroyed element and
      the last one. Fortunately, it seems that `m` tends to
      be small in general (the dynamics of the list of tracked
      elements is rather close to that of a stack).

    It is nevertheless still possible to manipulate several
    hundred variables (e.g. square matrices of size 5 or
    polynomials of degree 20).

    The class :class:`PrecisionLattice` provides several
    features for introspection, especially concerning timings.
    See :meth:`history` and :meth:`timings` for details.

    .. SEEALSO::

        :func:`ZpLF`
    """
    return Zp(p, prec, 'lattice-cap', *args, **kwds)

@experimental(23505)
def ZpLF(p, prec=None, *args, **kwds):
    r"""
    A shortcut function to create `p`-adic rings where precision
    is encoded by a module in a large vector space.

    See documentation for :func:`Zp` for a description of the input parameters.

    NOTE:

    The precision is tracked using automatic differentiation
    techniques (see [CRV2018]_ and [CRV2014]_).
    Floating point `p`-adic numbers are used for the computation
    of the differential (which is then not exact).

    EXAMPLES::

        sage: R = ZpLF(5, 40)
        sage: R
        5-adic Ring with lattice-float precision

    .. SEEALSO::

        :func:`ZpLC`
    """
    return Zp(p, prec, 'lattice-float', *args, **kwds)

def ZpER(p, prec=None, halt=None, secure=False, *args, **kwds):
    r"""
    A shortcut function to create relaxed `p`-adic rings.

    INPUT:

    - ``prec`` -- an integer (default: ``20``), the default
      precision

    - ``halt`` -- an integer (default: twice ``prec``), the
      halting precision

    - ``secure`` -- a boolean (default: ``False``); if ``False``,
      consider indistinguishable elements at the working precision
      as equal; otherwise, raise an error.

    See documentation for :func:`Zp` for a description of the other
    input parameters.

    A SHORT INTRODUCTION TO RELAXED `p`-ADICS:

    The model for relaxed `p`-adics is quite different from any of the
    other types of `p`-adics. In addition to storing a finite
    approximation, one also stores a method for increasing the
    precision.

    Relaxed `p`-adic rings are created by the constructor :func:`ZpER`::

        sage: R = ZpER(5, print_mode="digits")
        sage: R
        5-adic Ring handled with relaxed arithmetics

    The precision is not capped in `R`::

        sage: R.precision_cap()
        +Infinity

    However, a default precision is fixed. This is the precision
    at which the elements will be printed::

        sage: R.default_prec()
        20

    A default halting precision is also set. It is the default absolute
    precision at which the elements will be compared. By default, it is
    twice the default precision::

        sage: R.halting_prec()
        40

    However, both the default precision and the halting precision can be
    customized at the creation of the parent as follows:

        sage: S = ZpER(5, prec=10, halt=100)
        sage: S.default_prec()
        10
        sage: S.halting_prec()
        100

    One creates elements as usual::

        sage: a = R(17/42)
        sage: a
        ...00244200244200244201

        sage: R.random_element()  # random
        ...21013213133412431402

    Here we notice that 20 digits (that is the default precision) are printed.
    However, the computation model is designed in order to guarantee that more
    digits of `a` will be available on demand.
    This feature is reflected by the fact that, when we ask for the precision
    of `a`, the software answers `+\infty`::

        sage: a.precision_absolute()
        +Infinity

    Asking for more digits is achieved by the methods :meth:`at_precision_absolute`
    and :meth:`at_precision_relative`::

        sage: a.at_precision_absolute(30)
        ...?244200244200244200244200244201

    As a shortcut, one can use the bracket operator::

        sage: a[:30]
        ...?244200244200244200244200244201

    Of course, standard operations are supported::

        sage: b = R(42/17)
        sage: a + b
        ...03232011214322140002
        sage: a - b
        ...42311334324023403400
        sage: a * b
        ...00000000000000000001
        sage: a / b
        ...12442142113021233401
        sage: sqrt(a)
        ...20042333114021142101

    We observe again that only 20 digits are printed but, as before,
    more digits are available on demand::

        sage: sqrt(a)[:30]
        ...?142443342120042333114021142101

    .. RUBRIC:: Equality tests

    Checking equalities between relaxed `p`-adics is a bit subtle and can
    sometimes be puzzling at first glance.

    When the parent is created with ``secure=False`` (which is the
    default), elements are compared at the current precision, or at the
    default halting precision if it is higher::

        sage: a == b
        False

        sage: a == sqrt(a)^2
        True
        sage: a == sqrt(a)^2 + 5^50
        True

    In the above example, the halting precision is `40`; it is the
    reason why a congruence modulo `5^50` is considered as an equality.
    However, if both sides of the equalities have been previously
    computed with more digits, those digits are taken into account.
    Hence comparing two elements at different times can produce
    different results::

        sage: aa = sqrt(a)^2 + 5^50
        sage: a == aa
        True
        sage: a[:60]
        ...?244200244200244200244200244200244200244200244200244200244201
        sage: aa[:60]
        ...?244200244300244200244200244200244200244200244200244200244201
        sage: a == aa
        False

    This annoying situation, where the output of `a == aa` may change
    depending on previous computations, cannot occur when the parent is
    created with ``secure=True``.
    Indeed, in this case, if the equality cannot be decided, an error
    is raised::

        sage: S = ZpER(5, secure=True)
        sage: u = S.random_element()
        sage: uu = u + 5^50
        sage: u == uu
        Traceback (most recent call last):
        ...
        PrecisionError: unable to decide equality; try to bound precision

        sage: u[:60] == uu
        False

    .. RUBRIC:: Self-referent numbers

    A quite interesting feature with relaxed `p`-adics is the possibility to
    create (in some cases) self-referent numbers. Here is an example.
    We first declare a new variable as follows::

        sage: x = R.unknown()
        sage: x
        ...?.0

    We then use the method :meth:`set` to define `x` by writing down an equation
    it satisfies::

        sage: x.set(1 + 5*x^2)
        True

    The variable `x` now contains the unique solution of the equation
    `x = 1 + 5 x^2`::

        sage: x
        ...04222412141121000211

    This works because the `n`-th digit of the right hand size of the
    defining equation only involves the `i`-th digits of `x` with `i < n`
    (this is due to the factor `5`).

    As a comparison, the following does not work::

        sage: y = R.unknown()
        sage: y.set(1 + 3*y^2)
        True
        sage: y
        ...?.0
        sage: y[:20]
        Traceback (most recent call last):
        ...
        RecursionError: definition looks circular

    Self-referent definitions also work with systems of equations::

        sage: u = R.unknown()
        sage: v = R.unknown()
        sage: w = R.unknown()

        sage: u.set(1 + 2*v + 3*w^2 + 5*u*v*w)
        True
        sage: v.set(2 + 4*w + sqrt(1 + 5*u + 10*v + 15*w))
        True
        sage: w.set(3 + 25*(u*v + v*w + u*w))
        True

        sage: u
        ...31203130103131131433
        sage: v
        ...33441043031103114240
        sage: w
        ...30212422041102444403
    """
    return Zp(p, (prec, halt, secure), 'relaxed', *args, **kwds)


#######################################################################################################
#
#  The Extension Factory -- creates extensions of p-adic rings and fields
#
#######################################################################################################

class pAdicExtension_class(UniqueFactory):
    r"""
    A class for creating extensions of `p`-adic rings and fields.

    EXAMPLES::

        sage: R = Zp(5,3)
        sage: S.<x> = ZZ[]
        sage: W.<w> = pAdicExtension(R, x^4-15)
        sage: W
        5-adic Eisenstein Extension Ring in w defined by x^4 - 15
        sage: W.precision_cap()
        12
    """
    def create_key_and_extra_args(self, base, modulus, prec = None, print_mode = None,
                                  names = None, var_name = None, res_name = None,
                                  unram_name = None, ram_name = None, print_pos = None,
                                  print_sep = None, print_alphabet = None, print_max_ram_terms = None,
                                  print_max_unram_terms = None, print_max_terse_terms = None,
                                  show_prec = None, check = True, unram = False, implementation='FLINT'):
        r"""
        Creates a key from input parameters for pAdicExtension.

        See the documentation for ``Qq`` for more information.

        TESTS::

            sage: R = Zp(5,3)
            sage: S.<x> = ZZ[]
            sage: pAdicExtension.create_key_and_extra_args(R, x^4-15,names='w')
            (('e',
              5-adic Ring with capped relative precision 3,
              x^4 - 15,
              ('w', None, None, 'w'),
              12,
              'series',
              True,
              '|',
              (),
              -1,
              -1,
              -1,
              'bigoh',
              'NTL'),
             {'approx_modulus': (1 + O(5^3))*x^4 + O(5^4)*x^3 + O(5^4)*x^2 + O(5^4)*x + 2*5 + 4*5^2 + 4*5^3 + O(5^4)})

            sage: A = Qp(3,5)
            sage: Po.<X> = A[]
            sage: f = Po([3,0,-1])
            sage: K.<a> = A.ext(f)
            sage: -a^2+3
            O(a^12)
            sage: K.defining_polynomial() == f/f.leading_coefficient()
            True

            sage: g = Po([6,3,2])
            sage: H.<b> = A.ext(g)
            sage: 2*b^2+3*b+6
            O(b^12)
            sage: H.defining_polynomial() == g/g.leading_coefficient()
            True
        """
        if print_mode is None:
            print_mode = base.print_mode()
        if print_pos is None:
            print_pos = base._printer._pos()
        if print_sep is None:
            print_sep = base._printer._sep()
        if print_alphabet is None:
            print_alphabet = base._printer._alphabet()
        if print_max_ram_terms is None:
            print_max_ram_terms = base._printer._max_ram_terms()
        if print_max_unram_terms is None:
            print_max_unram_terms = base._printer._max_unram_terms()
        if print_max_terse_terms is None:
            print_max_terse_terms = base._printer._max_terse_terms()
        show_prec = _canonicalize_show_prec(base._prec_type(), print_mode, show_prec)
        from sage.structure.element import Expression
        if check:
            if isinstance(modulus, Expression):
                if len(modulus.variables()) != 1:
                    raise ValueError("symbolic expression must be in only one variable")
                exact_modulus = modulus.polynomial(base.exact_field())
                approx_modulus = modulus.polynomial(base)
            elif is_Polynomial(modulus):
                if modulus.parent().ngens() != 1:
                    raise ValueError("must use univariate polynomial")
                exact_modulus = modulus.change_ring(base.exact_field())
                approx_modulus = modulus.change_ring(base)
            else:
                raise ValueError("modulus must be a polynomial")
            if exact_modulus.degree() <= 1:
                raise NotImplementedError("degree of modulus must be at least 2")
            # need to add more checking here.
            if not unram and not exact_modulus.is_monic():
                exact_modulus = exact_modulus / exact_modulus.leading_coefficient()
                approx_modulus = approx_modulus / approx_modulus.leading_coefficient()
            if names is None:
                if var_name is not None:
                    names = var_name
                else:
                    raise ValueError("must specify name of generator of extension")
            if isinstance(names, tuple):
                names = names[0]
            if not isinstance(names, str):
                names = str(names)
        else:
            exact_modulus = modulus
            approx_modulus = modulus.change_ring(base)

        # We now decide on the extension class: unramified, Eisenstein, two-step or general
        if unram or is_unramified(approx_modulus):
            if unram_name is None:
                unram_name = names
            if res_name is None:
                res_name = unram_name + '0'
            if ram_name is None:
                ram_name = base._printer._uniformizer_name()
            names = (names, res_name, unram_name, ram_name)
            if base.absolute_degree() == 1:
                polytype = 'u'
            else:
                polytype = 'ru'
            if prec is None:
                prec = min([c.precision_absolute() for c in approx_modulus.list()] + [base.precision_cap()])
            elif prec > base.precision_cap():
                raise ValueError("Precision cannot be larger than that of base ring; you may want to call the change method on the base ring.")
            approx_modulus = truncate_to_prec(exact_modulus, base, prec)

        elif is_eisenstein(approx_modulus):
            unram_name = None
            res_name = None
            if ram_name is None:
                ram_name = names
            if base.absolute_degree() == 1:
                unram_name = None
                polytype = 'e'
            else:
                unram_name = base.variable_name()
                polytype = 're'
                implementation = 'Polynomial'
            names = (names, res_name, unram_name, ram_name)
            e = approx_modulus.degree()
            if prec is None:
                prec = min([c.precision_absolute() for c in approx_modulus.list() if not c._is_exact_zero()] + [base.precision_cap()]) * e
            elif prec > base.precision_cap() * e:
                raise ValueError("Precision cannot be larger than that of base ring; you may want to call the change method on the base ring.")
            approx_modulus = truncate_to_prec(exact_modulus, base, (prec/e).ceil() + 1)
        else:
            if unram_name is None:
                unram_name = names + '_u'
            if res_name is None:
                res_name = names + '0'
            if ram_name is None:
                ram_name = names + '_p'
            names = (names, res_name, unram_name, ram_name)
            polytype = 'p'
        if polytype == 'e':
            implementation = "NTL" # for testing - FLINT ramified extensions not implemented yet
        key = (polytype, base, exact_modulus, names, prec, print_mode, print_pos,
               print_sep, tuple(print_alphabet), print_max_ram_terms, print_max_unram_terms,
               print_max_terse_terms, show_prec, implementation)
        return key, {'approx_modulus': approx_modulus}

    def create_object(self, version, key, approx_modulus=None, shift_seed=None):
        r"""
        Creates an object using a given key.

        See the documentation for pAdicExtension for more information.

        TESTS::

            sage: R = Zp(5,3)
            sage: S.<x> = R[]
            sage: pAdicExtension.create_object(version = (6,4,2), key = ('e', R, x^4 - 15, x^4 - 15, ('w', None, None, 'w'), 12, None, 'series', True, '|', (),-1,-1,-1,'NTL'), shift_seed = S(3 + O(5^3)))
            5-adic Eisenstein Extension Ring in w defined by x^4 - 15
        """
        polytype = key[0]
        if version[0] < 6 or version[0] == 6 and version[1] < 1:
            key = list(key)
            key.append('NTL')
        if version[0] < 8:
            (polytype, base, premodulus, approx_modulus, names, prec, halt, print_mode, print_pos, print_sep,
             print_alphabet, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms, implementation) = key
            from sage.structure.element import Expression
            if isinstance(premodulus, Expression):
                exact_modulus = premodulus.polynomial(base.exact_field())
            elif is_Polynomial(premodulus):
                exact_modulus = premodulus.change_ring(base.exact_field())
            show_prec = None
        else:
            (polytype, base, exact_modulus, names, prec, print_mode, print_pos,
             print_sep, print_alphabet, print_max_ram_terms, print_max_unram_terms,
             print_max_terse_terms, show_prec, implementation) = key
            if polytype in ('e', 're'):
                unif = exact_modulus.base_ring()(base.uniformizer())
                shift_seed = (-exact_modulus[:exact_modulus.degree()] / unif).change_ring(base)
            if not krasner_check(exact_modulus, prec):
                raise ValueError("polynomial does not determine a unique extension.  Please specify more precision or use parameter check=False.")

        if show_prec is None:
            show_prec = base._printer._show_prec()
        if polytype == 'p':
            raise NotImplementedError("Extensions by general polynomials not yet supported.  Please use an unramified or Eisenstein polynomial.")
        T = ext_table[polytype, type(base.ground_ring_of_tower()).__base__]
        return T(exact_modulus, approx_modulus, prec,
                 {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet,
                  'max_ram_terms': print_max_ram_terms, 'max_unram_terms': print_max_unram_terms, 'max_terse_terms': print_max_terse_terms, 'show_prec': show_prec},
                 shift_seed, names, implementation)

ExtensionFactory = pAdicExtension = pAdicExtension_class("pAdicExtension")

######################################################
# Helper functions for the Extension Factory
######################################################

def split(poly, prec):
    r"""
    Given a polynomial ``poly`` and a desired precision ``prec``, computes
    ``upoly`` and epoly so that the extension defined by ``poly`` is isomorphic
    to the extension defined by first taking an extension by the unramified
    polynomial ``upoly``, and then an extension by the Eisenstein polynomial
    ``epoly``.

    We need better `p`-adic factoring in Sage before this function can be
    implemented.

    EXAMPLES::

        sage: k = Qp(13)
        sage: x = polygen(k)
        sage: f = x^2+1
        sage: sage.rings.padics.factory.split(f, 10)
        Traceback (most recent call last):
        ...
        NotImplementedError: Extensions by general polynomials not yet supported. Please use an unramified or Eisenstein polynomial.

    TESTS:

    This checks that :trac:`6186` is still fixed::

        sage: k = Qp(13)
        sage: x = polygen(k)
        sage: f = x^2+1
        sage: L.<a> = k.extension(f)
        Traceback (most recent call last):
        ...
        NotImplementedError: Extensions by general polynomials not yet supported. Please use an unramified or Eisenstein polynomial.

    """
    raise NotImplementedError("Extensions by general polynomials not yet supported.  Please use an unramified or Eisenstein polynomial.")

def truncate_to_prec(poly, R, absprec):
    r"""
    Truncates the unused precision off of a polynomial.

    EXAMPLES::

        sage: R = Zp(5)
        sage: S.<x> = R[]
        sage: from sage.rings.padics.factory import truncate_to_prec
        sage: f = x^4 + (3+O(5^6))*x^3 + O(5^4)
        sage: truncate_to_prec(f, R, 5)
        (1 + O(5^5))*x^4 + (3 + O(5^5))*x^3 + O(5^5)*x^2 + O(5^5)*x + O(5^4)
    """
    return R[poly.variable_name()]([R(a, absprec=absprec) for a in poly.list()]) # Is this quite right?  We don't want flat necessarily...

def krasner_check(poly, prec):
    r"""
    Returns True iff poly determines a unique isomorphism class of
    extensions at precision prec.

    Currently just returns True (thus allowing extensions that are not
    defined to high enough precision in order to specify them up to
    isomorphism).  This will change in the future.

    EXAMPLES::

        sage: from sage.rings.padics.factory import krasner_check
        sage: krasner_check(1,2) #this is a stupid example.
        True
    """
    return True #This needs to be implemented

def is_eisenstein(poly):
    r"""
    Returns True iff this monic polynomial is Eisenstein.

    A polynomial is Eisenstein if it is monic, the constant term has
    valuation 1 and all other terms have positive valuation.

    EXAMPLES::

        sage: R = Zp(5)
        sage: S.<x> = R[]
        sage: from sage.rings.padics.factory import is_eisenstein
        sage: f = x^4 - 75*x + 15
        sage: is_eisenstein(f)
        True
        sage: g = x^4 + 75
        sage: is_eisenstein(g)
        False
        sage: h = x^7 + 27*x -15
        sage: is_eisenstein(h)
        False
    """
    if poly[0].valuation() != 1:
        return False
    if reduce(lambda a, b: a or b, [(c.valuation() < 1) for c in poly.list()[1:poly.degree()]]):
        return False
    return True

def is_unramified(poly):
    r"""
    Returns true iff this monic polynomial is unramified.

    A polynomial is unramified if its reduction modulo the maximal
    ideal is irreducible.

    EXAMPLES::

        sage: R = Zp(5)
        sage: S.<x> = R[]
        sage: from sage.rings.padics.factory import is_unramified
        sage: f = x^4 + 14*x + 9
        sage: is_unramified(f)
        True
        sage: g = x^6 + 17*x + 6
        sage: is_unramified(g)
        False
    """
    if poly[0].valuation() > 0:
        return False
    if reduce(lambda a, b: a or b, [(c.valuation() < 0) for c in poly.list()[1:poly.degree()]]):
        return False
    F = poly.parent().change_ring(poly.base_ring().residue_class_field())(poly).factor()
    if len(F) != 1 or F[0][1] != 1:
        return False
    return True
