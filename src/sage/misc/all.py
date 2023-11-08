#from .all__sagemath_objects import *
from .all__sagemath_environment import *
from .all__sagemath_modules import *
from .all__sagemath_repl import *

from .remote_file import get_remote_file

lazy_import('sage.misc.dist', 'install_scripts', deprecation=34259)

from sage.misc.classgraph import class_graph

lazy_import('sage.repl.interpreter', 'logstr', deprecation=34259)
