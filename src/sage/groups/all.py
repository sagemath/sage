from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

from .pari_group import PariGroup

from .matrix_gps.all import *
from .abelian_gps.all import *

from .perm_gps.all import *

from .generic import *

lazy_import('sage.groups.class_function', 'ClassFunction')

from .additive_abelian.all import *

lazy_import('sage.groups.conjugacy_classes', ['ConjugacyClass', 'ConjugacyClassGAP'])

lazy_import('sage.groups.free_group', 'FreeGroup')
lazy_import('sage.groups.braid', 'BraidGroup')

lazy_import('sage.groups.affine_gps.affine_group', 'AffineGroup')
lazy_import('sage.groups.affine_gps.euclidean_group', 'EuclideanGroup')

lazy_import('sage.groups.raag', 'RightAngledArtinGroup')

lazy_import('sage.groups', 'groups_catalog', 'groups')

lazy_import('sage.groups.semimonomial_transformations.semimonomial_transformation_group', 'SemimonomialTransformationGroup')

lazy_import('sage.groups.group_exp', ['GroupExp', 'GroupExp_Class', 'GroupExpElement'])

lazy_import('sage.groups.group_semidirect_product', ['GroupSemidirectProduct', 'GroupSemidirectProductElement'])
