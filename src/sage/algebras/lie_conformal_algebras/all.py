# sage_setup: distribution = sagemath-combinat
# ******************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************


from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.lie_conformal_algebras.lie_conformal_algebra',
            'LieConformalAlgebra')
lazy_import('sage.algebras.lie_conformal_algebras',
            'examples', 'lie_conformal_algebras')
del lazy_import
