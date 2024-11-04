from sage.groups.all__sagemath_categories import *

from sage.misc.lazy_import import lazy_import

from sage.groups.pari_group import PariGroup

from sage.groups.matrix_gps.all import *
from sage.groups.abelian_gps.all import *

from sage.groups.perm_gps.all import *

lazy_import('sage.groups.class_function', 'ClassFunction')

from sage.groups.additive_abelian.all import *

lazy_import('sage.groups.conjugacy_classes', ['ConjugacyClass', 'ConjugacyClassGAP'])

lazy_import('sage.groups.free_group', 'FreeGroup')
lazy_import('sage.groups.braid', 'BraidGroup')
lazy_import('sage.groups.cubic_braid', 'CubicBraidGroup')
lazy_import('sage.groups.cubic_braid', 'AssionGroupU')
lazy_import('sage.groups.cubic_braid', 'AssionGroupS')

lazy_import('sage.groups.affine_gps.affine_group', 'AffineGroup')
lazy_import('sage.groups.affine_gps.euclidean_group', 'EuclideanGroup')

lazy_import('sage.groups.artin', 'ArtinGroup')
lazy_import('sage.groups.raag', 'RightAngledArtinGroup')

lazy_import('sage.groups.semimonomial_transformations.semimonomial_transformation_group', 'SemimonomialTransformationGroup')

lazy_import('sage.groups.group_exp', 'GroupExp')
lazy_import('sage.groups.group_exp', ['GroupExp_Class', 'GroupExpElement'],
            deprecation=38238)
lazy_import('sage.groups.group_semidirect_product', 'GroupSemidirectProduct')
lazy_import('sage.groups.group_semidirect_product', 'GroupSemidirectProductElement',
            deprecation=38238)

del lazy_import
