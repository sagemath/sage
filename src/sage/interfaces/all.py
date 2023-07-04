# interfaces to other interpreters

from .sage0 import sage0, sage0_version, Sage
from .gap import gap, gap_reset_workspace, Gap
from .gp import gp, gp_version, Gp
# import problems
#from maxima_lib import maxima_lib
from .maxima import maxima, Maxima
from .singular import singular, singular_version, Singular

from .magma import magma, Magma
from .polymake import polymake

from sage.misc.lazy_import import lazy_import

lazy_import('sage.interfaces.axiom', ['Axiom', 'axiom'])
lazy_import('sage.interfaces.ecm', ['ECM', 'ecm'])
lazy_import('sage.interfaces.four_ti_2', 'four_ti_2')
lazy_import('sage.interfaces.fricas', ['FriCAS', 'fricas'])
lazy_import('sage.interfaces.frobby', 'frobby')
lazy_import('sage.interfaces.gap3', ['gap3', 'gap3_version', 'Gap3'])
lazy_import('sage.interfaces.genus2reduction', ['genus2reduction', 'Genus2reduction'])
lazy_import('sage.interfaces.gfan', ['gfan', 'Gfan'])
lazy_import('sage.interfaces.giac', ['giac', 'Giac'])
lazy_import('sage.interfaces.gnuplot', 'gnuplot')
lazy_import('sage.interfaces.kash', ['kash', 'kash_version', 'Kash'])
lazy_import('sage.interfaces.lie', ['lie', 'LiE'])
lazy_import('sage.interfaces.lisp', ['lisp', 'Lisp'])
lazy_import('sage.interfaces.macaulay2', ['macaulay2', 'Macaulay2'])
lazy_import('sage.interfaces.magma_free', 'magma_free')
lazy_import('sage.interfaces.maple', ['maple', 'Maple'])
lazy_import('sage.interfaces.mathematica', ['mathematica', 'Mathematica'])
lazy_import('sage.interfaces.mathics', ['mathics', 'Mathics'])
lazy_import('sage.interfaces.matlab', ['matlab', 'matlab_version', 'Matlab'])
lazy_import('sage.interfaces.mupad', ['mupad', 'Mupad'])  # NOT functional yet
lazy_import('sage.interfaces.mwrank', ['mwrank', 'Mwrank'])
lazy_import('sage.interfaces.octave', ['octave', 'Octave'])
lazy_import('sage.interfaces.povray', 'povray')
lazy_import('sage.interfaces.psage', 'PSage')
lazy_import('sage.interfaces.qepcad', ['qepcad', 'qepcad_version', 'qepcad_formula'])
lazy_import('sage.interfaces.r', ['r', 'R', 'r_version'])
lazy_import('sage.interfaces.read_data', 'read_data')
lazy_import('sage.interfaces.scilab', 'scilab')
lazy_import('sage.interfaces.tachyon', 'tachyon_rt')

# The following variable is used by sage-shell-mode in emacs:
interfaces = ['gap', 'gap3', 'giac', 'gp', 'mathematica', 'gnuplot',
              'kash', 'magma', 'macaulay2', 'maple', 'maxima',
              'mathematica', 'mwrank', 'octave', 'r', 'singular',
              'sage0', 'sage']
