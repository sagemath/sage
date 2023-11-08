"""
Tests for deprecations of imports in global namespace from :trac:`28064`::

    sage: igusa_clebsch_invariants
    doctest:warning...:
    DeprecationWarning:
    Importing igusa_clebsch_invariants from here is deprecated;
    please use "from sage.schemes.hyperelliptic_curves.invariants import igusa_clebsch_invariants" instead.
    See https://github.com/sagemath/sage/issues/28064 for details.
    ...

    sage: absolute_igusa_invariants_kohel
    doctest:warning...:
    DeprecationWarning:
    Importing absolute_igusa_invariants_kohel from here is deprecated;
    please use "from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_kohel" instead.
    See https://github.com/sagemath/sage/issues/28064 for details.
    ...

    sage: absolute_igusa_invariants_wamelen
    doctest:warning...:
    DeprecationWarning:
    Importing absolute_igusa_invariants_wamelen from here is deprecated;
    please use "from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_wamelen" instead.
    See https://github.com/sagemath/sage/issues/28064 for details.
    ...

    sage: clebsch_invariants
    doctest:warning...:
    DeprecationWarning:
    Importing clebsch_invariants from here is deprecated;
    please use "from sage.schemes.hyperelliptic_curves.invariants import clebsch_invariants" instead.
    See https://github.com/sagemath/sage/issues/28064 for details.
    ...
"""
from sage.misc.lazy_import import lazy_import

from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.schemes.hyperelliptic_curves.kummer_surface import KummerSurface
lazy_import('sage.schemes.hyperelliptic_curves.invariants',
            ['igusa_clebsch_invariants', 'absolute_igusa_invariants_kohel',
             'absolute_igusa_invariants_wamelen', 'clebsch_invariants'],
            deprecation=28064)
from sage.schemes.hyperelliptic_curves.mestre import (Mestre_conic, HyperellipticCurve_from_invariants)
from sage.schemes.hyperelliptic_curves import monsky_washnitzer
