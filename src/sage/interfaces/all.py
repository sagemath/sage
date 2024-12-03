# interfaces to other interpreters

# import problems
try:
    # from maxima_lib import maxima_lib
    from sage.interfaces.maxima import maxima, Maxima
except ImportError:
    pass

from sage.misc.lazy_import import lazy_import

lazy_import('sage.interfaces.sage0', ['sage0', 'sage0_version', 'Sage'])
lazy_import('sage.interfaces.axiom', ['Axiom', 'axiom'])
lazy_import('sage.interfaces.ecm', ['ECM', 'ecm'])
lazy_import('sage.interfaces.four_ti_2', 'four_ti_2')
lazy_import('sage.interfaces.fricas', ['FriCAS', 'fricas'])
lazy_import('sage.interfaces.frobby', 'frobby')
lazy_import('sage.interfaces.gap', ['gap', 'gap_reset_workspace', 'Gap'])
lazy_import('sage.interfaces.gap3', ['gap3', 'gap3_version', 'Gap3'])
lazy_import('sage.interfaces.genus2reduction', ['genus2reduction', 'Genus2reduction'])
lazy_import('sage.interfaces.gfan', ['gfan', 'Gfan'])
lazy_import('sage.interfaces.giac', ['giac', 'Giac'])
lazy_import('sage.interfaces.gnuplot', 'gnuplot')
lazy_import('sage.interfaces.gp', ['gp', 'gp_version', 'Gp'])
lazy_import('sage.interfaces.kash', ['kash', 'kash_version', 'Kash'])
lazy_import('sage.interfaces.lie', ['lie', 'LiE'])
lazy_import('sage.interfaces.lisp', ['lisp', 'Lisp'])
lazy_import('sage.interfaces.macaulay2', ['macaulay2', 'Macaulay2'])
lazy_import('sage.interfaces.magma', ['magma', 'Magma'])
lazy_import('sage.interfaces.magma_free', 'magma_free')
lazy_import('sage.interfaces.maple', ['maple', 'Maple'])
lazy_import('sage.interfaces.mathematica', ['mathematica', 'Mathematica'])
lazy_import('sage.interfaces.mathics', ['mathics', 'Mathics'])
lazy_import('sage.interfaces.matlab', ['matlab', 'matlab_version', 'Matlab'])
lazy_import('sage.interfaces.mupad', ['mupad', 'Mupad'])  # NOT functional yet
lazy_import('sage.interfaces.mwrank', ['mwrank', 'Mwrank'])
lazy_import('sage.interfaces.octave', ['octave', 'Octave'])
lazy_import('sage.interfaces.polymake', 'polymake')
lazy_import('sage.interfaces.povray', 'povray')
lazy_import('sage.interfaces.psage', 'PSage')
lazy_import('sage.interfaces.qepcad', ['qepcad', 'qepcad_version', 'qepcad_formula'])
lazy_import('sage.interfaces.r', ['r', 'R', 'r_version'])
lazy_import('sage.interfaces.read_data', 'read_data')
lazy_import('sage.interfaces.scilab', 'scilab')
lazy_import('sage.interfaces.singular', ['singular', 'singular_version', 'Singular'])
lazy_import('sage.interfaces.tachyon', 'tachyon_rt')

# The following variable is used by sage-shell-mode in emacs:
interfaces = ['gap', 'gap3', 'giac', 'gp', 'mathematica', 'gnuplot',
              'kash', 'magma', 'macaulay2', 'maple', 'maxima',
              'mathematica', 'mwrank', 'octave', 'r', 'singular',
              'sage0', 'sage']
del lazy_import
