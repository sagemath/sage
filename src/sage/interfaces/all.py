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
lazy_import('sage.interfaces.qsieve', 'qsieve')
lazy_import('sage.interfaces.r', ['r', 'R', 'r_version'])
lazy_import('sage.interfaces.read_data', 'read_data')
lazy_import('sage.interfaces.scilab', 'scilab')
lazy_import('sage.interfaces.tachyon', 'tachyon_rt')

try:
    from sage.repl.rich_output.display_manager import get_display_manager as _get_display_manager
except ImportError:
    pass
else:
    if _get_display_manager().is_in_terminal():
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
