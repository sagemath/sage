# interfaces to other interpreters

from sage.misc.lazy_import import lazy_import

from .frobby import frobby
from .four_ti_2 import four_ti_2
from .axiom import Axiom, axiom
from .fricas import FriCAS, fricas

from .gap import gap, gap_reset_workspace, Gap
from .gap3 import gap3, gap3_version, Gap3
lazy_import('sage.interfaces.genus2reduction', ['genus2reduction', 'Genus2reduction'])
from .gfan import gfan,  Gfan
from .giac import giac, Giac
from .gp import gp, gp_version, Gp
from .gnuplot import gnuplot
from .kash import  kash, kash_version, Kash
from .lisp import lisp, Lisp
from .magma import magma, Magma
from .magma_free import magma_free
from .macaulay2 import macaulay2, Macaulay2
from .maple import maple, Maple
from .maxima import maxima, Maxima
# import problems
#from maxima_lib import maxima_lib
from .mathematica import mathematica, Mathematica
from .mathics import mathics, Mathics
from .matlab import matlab, matlab_version, Matlab
from .mupad import mupad, Mupad  # NOT functional yet
from .mwrank import mwrank, Mwrank
from .octave import octave, Octave
from .polymake import polymake
from .qepcad import qepcad, qepcad_version, qepcad_formula
from .qsieve import qsieve
from .singular import singular, singular_version, Singular
from .sage0 import sage0 as sage0, sage0_version, Sage
from .scilab import scilab
from .tachyon import tachyon_rt
from .psage import PSage
from .ecm import ECM, ecm
from .povray import povray
from .lie import lie, LiE
from .r import r, R, r_version
from .read_data import read_data

interfaces = ['gap', 'gap3', 'giac', 'gp', 'mathematica', 'gnuplot', \
              'kash', 'magma', 'macaulay2', 'maple', 'maxima', \
              'mathematica', 'mwrank', 'octave', 'r', \
              'singular', 'sage0', 'sage']


from sage.repl.rich_output.display_manager import get_display_manager
if get_display_manager().is_in_terminal():
    from .axiom import axiom_console
    from .fricas import fricas_console
    from .gap import gap_console
    from .gap3 import gap3_console
    from .giac import giac_console
    from .gp import gp_console
    from .gnuplot import gnuplot_console
    from .kash import  kash_console
    from .lisp import lisp_console
    from .magma import magma_console
    from .macaulay2 import macaulay2_console
    from .maple import maple_console
    from .maxima_abstract import maxima_console
    from .mathematica import mathematica_console
    from .mathics import mathics_console
    from .matlab import matlab_console
    from .mupad import mupad_console
    from .mwrank import mwrank_console
    from .octave import octave_console
    from .qepcad import qepcad_console
    from .singular import singular_console
    from .sage0 import sage0_console
    from .lie import lie_console
    from .r import r_console
