#from .all__sagemath_objects import *
from .all__sagemath_environment import *
from .all__sagemath_modules import *
from .all__sagemath_repl import *

from .remote_file import get_remote_file

from .dist import install_scripts

from .classgraph import class_graph


class logstr(str):
    def __repr__(self):
        return self

    def _latex_(self):
        # return "\\begin{verbatim}%s\\end{verbatim}"%self
        if '#' not in self:
            delim = '#'
        elif '@' not in self:
            delim = '@'
        elif '~' not in self:
            delim = '~'
        return r"""\verb%s%s%s""" % (delim, self.replace('\n\n', '\n').replace('\n', '; '), delim)
