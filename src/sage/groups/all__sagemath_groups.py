# sage_setup: distribution = sagemath-groups
from sage.groups.all__sagemath_modules import *

try:
    from sage.groups.all__sagemath_pari import *
except ImportError:
    pass

from sage.groups.all__sagemath_gap import *

from sage.misc.lazy_import import lazy_import

lazy_import('sage.groups.class_function', 'ClassFunction')

lazy_import('sage.groups.conjugacy_classes', ['ConjugacyClass', 'ConjugacyClassGAP'])

lazy_import('sage.groups.free_group', 'FreeGroup')
lazy_import('sage.groups.braid', 'BraidGroup')
lazy_import('sage.groups.cubic_braid', 'CubicBraidGroup')
lazy_import('sage.groups.cubic_braid', 'AssionGroupU')
lazy_import('sage.groups.cubic_braid', 'AssionGroupS')

lazy_import('sage.groups.artin', 'ArtinGroup')
lazy_import('sage.groups.raag', 'RightAngledArtinGroup')

lazy_import('sage.groups.semimonomial_transformations.semimonomial_transformation_group',
            'SemimonomialTransformationGroup')

lazy_import('sage.groups.group_exp', ['GroupExp', 'GroupExp_Class', 'GroupExpElement'])

lazy_import('sage.groups.group_semidirect_product', [
            'GroupSemidirectProduct', 'GroupSemidirectProductElement'])
del lazy_import
