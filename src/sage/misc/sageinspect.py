# sage_setup: distribution = sagemath-objects
r"""
Inspect Python, Sage, and Cython objects

This module extends parts of Python's inspect module to Cython objects.

EXAMPLES::

    sage: from sage.misc.sageinspect import *

Test introspection of modules defined in Python and Cython files:

Cython modules::

    sage: sage_getfile(sage.rings.rational)
    '.../rational.pyx'
    sage: sage_getdoc(sage.rings.rational).lstrip()
    'Rational Numbers...'
    sage: sage_getsource(sage.rings.rational)
    '# distutils: ...Rational Numbers...'

Python modules::

    sage: sage_getfile(sage.misc.sageinspect)
    '.../sageinspect.py'
    sage: print(sage_getdoc(sage.misc.sageinspect).lstrip()[:40])
    Inspect Python, Sage, and Cython objects
    sage: sage_getsource(sage.misc.sageinspect).lstrip()[51:-1]
    'Inspect Python, Sage, and Cython objects...'

Test introspection of classes defined in Python and Cython files:

Cython classes::

    sage: sage_getfile(sage.rings.rational.Rational)
    '.../rational.pyx'
    sage: sage_getdoc(sage.rings.rational.Rational).lstrip()
    'A rational number...'
    sage: sage_getsource(sage.rings.rational.Rational)
    'cdef class Rational...'

Python classes::

    sage: sage_getfile(BlockFinder)
    '.../sage/misc/sageinspect.py'
    sage: sage_getdoc(BlockFinder).lstrip()[:50]                                        # needs sphinx
    'Provide a "tokeneater()" method to detect the end '
    sage: sage_getsource(BlockFinder)
    'class BlockFinder:...'

Test introspection of functions defined in Python and Cython files:

Cython functions::

    sage: sage_getdef(sage.rings.rational.make_rational, obj_name='mr')
    'mr(s)'
    sage: sage_getfile(sage.rings.rational.make_rational)
    '.../rational.pyx'
    sage: sage_getdoc(sage.rings.rational.make_rational).lstrip()
    'Make a rational number ...'
    sage: sage_getsource(sage.rings.rational.make_rational)
    '@cython.binding(True)\ndef make_rational(s):...'

Python functions::

    sage: sage_getdef(sage.misc.sageinspect.sage_getfile, obj_name='sage_getfile')
    'sage_getfile(obj)'
    sage: sage_getfile(sage.misc.sageinspect.sage_getfile)
    '.../sageinspect.py'
    sage: sage_getdoc(sage.misc.sageinspect.sage_getfile).lstrip()
    'Get the full file name associated to "obj" as a string...'
    sage: sage_getsource(sage.misc.sageinspect.sage_getfile)[4:]
    'sage_getfile(obj):...'

Unfortunately, no argspec is extractable from builtins. Hence, we use a
generic argspec::

    sage: sage_getdef(''.find, 'find')
    'find(*args, **kwds)'
    sage: sage_getdef(str.find, 'find')
    'find(*args, **kwds)'

By :issue:`9976` and :issue:`14017`, introspection also works for interactively
defined Cython code, and with rather tricky argument lines::

    sage: # needs sage.misc.cython
    sage: cython('def foo(unsigned int x=1, a=\')"\', b={not (2+1==3):\'bar\'}, *args, **kwds): return')
    sage: print(sage_getsource(foo))
    def foo(unsigned int x=1, a=')"', b={not (2+1==3):'bar'}, *args, **kwds): return
    sage: sage_getargspec(foo)
    FullArgSpec(args=['x', 'a', 'b'], varargs='args', varkw='kwds', defaults=(1, ')"', {False: 'bar'}), kwonlyargs=[], kwonlydefaults=None, annotations={})

AUTHORS:

- Originally taken from Fernando Perez's IPython
- William Stein: extensive modifications
- William Stein: in :func:`_sage_getargspec_cython`, a modified version of
  ``inspect.getargspec`` from the Python Standard Library, which was taken from
  IPython for use in Sage
- Nick Alexander: extensions, testing
- Simon King: some extension for Cython, generalisation of SageArgSpecVisitor
- Simon King: in :func:`sage_getsourcelines`, if a class has no docstring then let the
  class definition be found starting from the ``__init__`` method.
- Simon King: in :func:`sage_getsourcelines`, get source lines for dynamic classes
- Simon King: in :func:`_sage_getargspec_cython`, return an ``ArgSpec``, fix some bugs
- Simon King (2011-09): added :func:`_sage_getsourcelines_name_with_dot`
- Simon King (2013-02): in :func:`_sage_getargspec_cython`, recognise varargs and
  default values in cython code, and return an ``ArgSpec``
"""

import ast
import inspect
import functools
import os
import tokenize
import re
from inspect import Signature, Parameter

try:
    import importlib.machinery as import_machinery
except ImportError:
    pass


def is_function_or_cython_function(obj):
    """
    Check whether something is a function.

    This is a variant of :func:`inspect.isfunction`:
    We assume that anything which has a genuine ``__code__``
    attribute (not using ``__getattr__`` overrides) is a function.
    This is meant to support Cython functions.

    Think twice before using this function (or any function from the
    :mod:`inspect` or :mod:`sage.misc.sageinspect` modules).  Most uses of
    :func:`inspect.isfunction` in ordinary library code can be replaced by
    :func:`callable`.

    EXAMPLES::

        sage: from sage.misc.sageinspect import is_function_or_cython_function
        sage: def f(): pass
        sage: is_function_or_cython_function(f)
        True
        sage: is_function_or_cython_function(lambda x:x)
        True
        sage: from sage.categories.coercion_methods import _mul_parent
        sage: is_function_or_cython_function(_mul_parent)
        True
        sage: is_function_or_cython_function(Integer.digits)     # unbound method
        False
        sage: is_function_or_cython_function(Integer(1).digits)  # bound method
        False

    TESTS:

    Verify that ipywidgets can correctly determine signatures of Cython
    functions::

        sage: from ipywidgets.widgets.interaction import signature
        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandelbrot_plot   # needs sage.symbolic
        sage: signature(fast_mandelbrot_plot)  # random                                 # needs sage.symbolic
        <IPython.utils._signatures.Signature object at 0x7f3ec8274e10>
    """
    # We use type(obj) instead of just obj to avoid __getattr__().
    # Some types, like methods, will return the __code__ of the
    # underlying function in __getattr__() but we don't want to
    # detect those as functions.
    return hasattr(type(obj), "__code__")


def isclassinstance(obj):
    r"""
    Check if argument is instance of non built-in class.

    INPUT:

    - ``obj`` -- object

    EXAMPLES::

        sage: from sage.misc.sageinspect import isclassinstance
        sage: isclassinstance(int)
        False
        sage: class myclass: pass
        sage: isclassinstance(myclass())
        True
        sage: isclassinstance(myclass)
        False
        sage: class mymetaclass(type): pass
        sage: class myclass2(metaclass=mymetaclass): pass
        sage: isclassinstance(myclass2)
        False
    """
    builtin_mods = {'__builtin__', 'builtins', 'exceptions'}

    return (not inspect.isclass(obj) and
            hasattr(obj, '__class__') and
            hasattr(obj.__class__, '__module__') and
            obj.__class__.__module__ not in builtin_mods and
            # Starting with Cython 3, Cython's builtin types have __module__ set
            # to the shared module names like _cython_3_0_0.
            not (isinstance(obj.__class__.__module__, str) and
                 obj.__class__.__module__.startswith('_cython_')))


# Parse strings of form "File: sage/rings/rational.pyx (starting at line 1080)"
# "\ " protects a space in re.VERBOSE mode.
__embedded_position_re = re.compile(r'''
^                                           # anchor to the beginning of the line
File:\ (?P<FILENAME>.*?)                    # match File: then filename
\ \(starting\ at\ line\ (?P<LINENO>\d+)\)   # match line number
\n?                                         # if there is a newline, eat it
(?P<ORIGINAL>.*)                            # the original docstring is the end
\Z                                          # anchor to the end of the string
''', re.MULTILINE | re.DOTALL | re.VERBOSE)

# Parse Python identifiers
__identifier_re = re.compile(r"^[^\d\W]\w*")


def _extract_embedded_position(docstring):
    r"""
    If docstring has a Cython embedded position, return a tuple
    (original_docstring, filename, line).  If not, return None.

    INPUT:

    - ``docstring`` -- string

    EXAMPLES::

       sage: from sage.misc.sageinspect import _extract_embedded_position
       sage: import inspect
       sage: _extract_embedded_position(inspect.getdoc(var))[1][-21:]                   # needs sage.symbolic
       'sage/calculus/var.pyx'

    TESTS:

    The following has been fixed in :issue:`13916`::

        sage: cython('''cpdef test_funct(x, y): return''')                               # needs sage.misc.cython
        sage: func_doc = inspect.getdoc(test_funct)                                     # needs sage.misc.cython
        sage: with open(_extract_embedded_position(func_doc)[1]) as f:                  # needs sage.misc.cython
        ....:     print(f.read())
        cpdef test_funct(x, y): return

    Ensure that the embedded filename of the compiled function is
    correct.  In particular it should be relative to ``spyx_tmp()`` in
    order for certain documentation functions to work properly.  See
    :issue:`24097`::

        sage: from sage.env import DOT_SAGE
        sage: from sage.misc.sage_ostools import restore_cwd
        sage: with restore_cwd(DOT_SAGE):                                               # needs sage.misc.cython
        ....:     cython('''cpdef test_funct(x, y): return''')
        sage: func_doc = inspect.getdoc(test_funct)                                     # needs sage.misc.cython
        sage: with open(_extract_embedded_position(func_doc)[1]) as f:                  # needs sage.misc.cython
        ....:     print(f.read())
        cpdef test_funct(x, y): return
    """
    try:
        res = __embedded_position_re.search(docstring)
    except TypeError:
        return None

    if res is None:
        return None

    raw_filename = res.group('FILENAME')
    filename = raw_filename

    if not os.path.isabs(filename):
        # Try some common path prefixes for Cython modules built by/for Sage
        # 1) Module in the sage src tree
        # 2) Module compiled by Sage's inline cython() compiler
        from sage.misc.temporary_file import spyx_tmp
        if raw_filename.startswith('sage/'):
            import sage
            try_filenames = [os.path.join(directory, raw_filename[5:])
                             for directory in sage.__path__]
        else:
            try_filenames = []
        try_filenames.append(
            os.path.join(spyx_tmp(), '_'.join(raw_filename.split('_')[:-1]),
                         raw_filename))
        for try_filename in try_filenames:
            if os.path.exists(try_filename):
                filename = try_filename
                break
        # Otherwise we keep the relative path and just hope it's relative to
        # the cwd; otherwise there's no way to be sure.

    lineno = int(res.group('LINENO'))
    original = res.group('ORIGINAL')
    return (original, filename, lineno)


def _extract_embedded_signature(docstring, name):
    r"""
    If docstring starts with the embedded of a method called ``name``, return
    a tuple (original_docstring, argspec).  If not, return (docstring, None).

    See :issue:`17814`.

    INPUT:

    - ``docstring`` -- string

    EXAMPLES::

        sage: from sage.misc.sageinspect import _extract_embedded_signature
        sage: from sage.misc.nested_class import MainClass
        sage: print(_extract_embedded_signature(MainClass.NestedClass.NestedSubClass.dummy.__doc__, 'dummy')[0])
        File: ...sage/misc/nested_class.pyx (starting at line ...)
        ...
        sage: _extract_embedded_signature(MainClass.NestedClass.NestedSubClass.dummy.__doc__, 'dummy')[1]
        FullArgSpec(args=['self', 'x', 'r'], varargs='args', varkw='kwds', defaults=((1, 2, 3.4),), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _extract_embedded_signature(range.__call__.__doc__, '__call__')
        ('Call self as a function.', None)

    """
    # If there is an embedded signature, it is in the first line
    L = docstring.split(os.linesep, 1)
    firstline = L[0]
    # It is possible that the signature is of the form ClassName.method_name,
    # and thus we need to do the following:
    if name not in firstline:
        return docstring, None
    signature = firstline.split(name, 1)[-1]
    if signature.startswith("(") and signature.endswith(")"):
        docstring = L[1] if len(L) > 1 else ''  # Remove first line, keep the rest
        def_string = "def " + name + signature + ": pass"
        try:
            return docstring, inspect.FullArgSpec(*_sage_getargspec_cython(def_string))
        except SyntaxError:
            docstring = os.linesep.join(L)
    return docstring, None


class BlockFinder:
    """
    Provide a :meth:`tokeneater` method to detect the end of a code block.

    This is the Python library's :class:`inspect.BlockFinder` modified
    to recognize Cython definitions.
    """
    def __init__(self):
        self.indent = 0
        self.islambda = False
        self.started = False
        self.passline = False
        self.last = 1

    def tokeneater(self, type, token, srow_scol, erow_ecol, line):
        srow, scol = srow_scol
        erow, ecol = erow_ecol
        if not self.started:
            # look for the first "(cp)def", "class" or "lambda"
            if token in ("def", "cpdef", "class", "lambda"):
                if token == "lambda":
                    self.islambda = True
                self.started = True
            self.passline = True    # skip to the end of the line
        elif type == tokenize.NEWLINE:
            self.passline = False   # stop skipping when a NEWLINE is seen
            self.last = srow
            if self.islambda:       # lambdas always end at the first NEWLINE
                raise inspect.EndOfBlock
        elif self.passline:
            pass
        elif type == tokenize.INDENT:
            self.indent = self.indent + 1
            self.passline = True
        elif type == tokenize.DEDENT:
            self.indent = self.indent - 1
            # the end of matching indent/dedent pairs end a block
            # (note that this only works for "def"/"class" blocks,
            #  not e.g. for "if: else:" or "try: finally:" blocks)
            if self.indent <= 0:
                raise inspect.EndOfBlock
        elif self.indent == 0 and type not in (tokenize.COMMENT, tokenize.NL):
            # any other token on the same indentation level end the previous
            # block as well, except the pseudo-tokens COMMENT and NL.
            raise inspect.EndOfBlock


def _getblock(lines):
    """
    Extract the block of code at the top of the given list of lines.

    This is the Python library's :func:`inspect.getblock`, except that
    it uses an instance of our custom :class:`BlockFinder`.
    """
    blockfinder = BlockFinder()
    iter_lines = iter(lines)
    tokenizer = tokenize.tokenize

    def readline():
        return next(iter_lines).encode('utf-8')
    try:
        for tok in tokenizer(readline):
            blockfinder.tokeneater(*tok)
    except (inspect.EndOfBlock, IndentationError):
        pass
    return lines[:blockfinder.last]


def _extract_source(lines, lineno):
    r"""
    Given a list of lines or a multiline string and a starting lineno,
    _extract_source returns [source_lines].  [source_lines] is the smallest
    indentation block starting at lineno.

    INPUT:

    - ``lines`` -- string or list of strings
    - ``lineno`` -- positive integer

    EXAMPLES::

        sage: from sage.misc.sageinspect import _extract_source
        sage: s2 = "#hello\n\n  class f():\n    pass\n\n#goodbye"
        sage: _extract_source(s2, 3)
        ['  class f():\n', '    pass\n']
    """
    if lineno < 1:
        raise ValueError("Line numbering starts at 1! (tried to extract line {})".format(lineno))
    lineno -= 1

    if isinstance(lines, str):
        lines = lines.splitlines(True)  # true keeps the '\n'
    if len(lines):
        # Fixes an issue with getblock
        lines[-1] += '\n'

    return _getblock(lines[lineno:])


class SageArgSpecVisitor(ast.NodeVisitor):
    """
    A simple visitor class that walks an abstract-syntax tree (AST)
    for a Python function's argspec.  It returns the contents of nodes
    representing the basic Python types: None, booleans, numbers,
    strings, lists, tuples, and dictionaries.  We use this class in
    :func:`_sage_getargspec_from_ast` to extract an argspec from a
    function's or method's source code.

    EXAMPLES::

        sage: import ast, sage.misc.sageinspect as sms
        sage: visitor = sms.SageArgSpecVisitor()
        sage: visitor.visit(ast.parse('[1,2,3]').body[0].value)
        [1, 2, 3]
        sage: v = visitor.visit(ast.parse("{'a':('e',2,[None,({False:True},'pi')]), 37.0:'temp'}").body[0].value)
        sage: sorted(v.items(), key=lambda x: str(x[0]))
        [(37.0, 'temp'), ('a', ('e', 2, [None, ({False: True}, 'pi')]))]
        sage: v = ast.parse("jc = ['veni', 'vidi', 'vici']").body[0]; v
        <...ast.Assign object at ...>
        sage: attrs = [x for x in dir(v) if not x.startswith('__')]
        sage: '_attributes' in attrs and '_fields' in attrs and 'col_offset' in attrs
        True
        sage: visitor.visit(v.targets[0])
        'jc'
        sage: visitor.visit(v.value)
        ['veni', 'vidi', 'vici']
    """
    def visit_Name(self, node):
        """
        Visit a Python AST :class:`ast.Name` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: ``None``, ``True``, ``False``, or the ``node``'s name as a string

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Name(ast.parse(x).body[0].value)
            sage: [vis(n) for n in ['foo', 'bar']]
            ['foo', 'bar']
            sage: [type(vis(n)) for n in ['foo', 'bar']]
            [<class 'str'>, <class 'str'>]
        """
        return node.id

    def visit_NameConstant(self, node):
        """
        Visit a Python AST :class:`ast.NameConstant` node.

        This is an optimization added in Python 3.4 for the special cases
        of True, False, and None.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: ``None``, ``True``, ``False``

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_NameConstant(ast.parse(x).body[0].value)
            sage: [vis(n) for n in ['True', 'False', 'None']]
            [True, False, None]
            sage: [type(vis(n)) for n in ['True', 'False', 'None']]
            [<class 'bool'>, <class 'bool'>, <class 'NoneType'>]
        """
        return node.value

    def visit_arg(self, node):
        r"""
        Visit a Python AST :class:`ast.arg` node.

        This node type is only on Python 3, where function arguments are
        more complex than just an identifier (e.g. they may also include
        annotations).

        For now we simply return the argument identifier as a string.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the argument name

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: s = "def f(a, b=2, c={'a': [4, 5.5, False]}, d=(None, True)):\n    return"
            sage: visitor = sms.SageArgSpecVisitor()
            sage: args = ast.parse(s).body[0].args.args
            sage: [visitor.visit_arg(n) for n in args]
            ['a', 'b', 'c', 'd']
        """
        return node.arg

    def visit_Num(self, node):
        """
        Visit a Python AST :class:`ast.Num` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the number the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Num(ast.parse(x).body[0].value)
            sage: [vis(n) for n in ['123', '0.0']]
            [123, 0.0]

        .. NOTE::

            On Python 3 negative numbers are parsed first, for some reason, as
            a UnaryOp node.
        """
        return node.value

    def visit_Str(self, node):
        r"""
        Visit a Python AST :class:`ast.Str` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the string the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Str(ast.parse(x).body[0].value)
            sage: [vis(s) for s in ['"abstract"', "'syntax'", r'''r"tr\ee"''']]
            ['abstract', 'syntax', 'tr\\ee']
        """
        return node.value

    def visit_List(self, node):
        """
        Visit a Python AST :class:`ast.List` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the list the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_List(ast.parse(x).body[0].value)
            sage: [vis(l) for l in ['[]', "['s', 't', 'u']", '[[e], [], [pi]]']]
            [[], ['s', 't', 'u'], [['e'], [], ['pi']]]
        """
        return [self.visit(n) for n in node.elts]

    def visit_Tuple(self, node):
        """
        Visit a Python AST :class:`ast.Tuple` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the tuple the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Tuple(ast.parse(x).body[0].value)
            sage: [vis(t) for t in ['()', '(x,y)', '("Au", "Al", "Cu")']]
            [(), ('x', 'y'), ('Au', 'Al', 'Cu')]
        """
        return tuple(self.visit(n) for n in node.elts)

    def visit_Dict(self, node):
        """
        Visit a Python AST :class:`ast.Dict` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the dictionary the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Dict(ast.parse(x).body[0].value)
            sage: v = [vis(d) for d in ['{}', "{1:one, 'two':2, other:bother}"]]
            sage: [sorted(d.items(), key=lambda x: str(x[0])) for d in v]
            [[], [(1, 'one'), ('other', 'bother'), ('two', 2)]]
        """
        d = {}
        for k, v in zip(node.keys, node.values):
            d[self.visit(k)] = self.visit(v)
        return d

    def visit_BoolOp(self, node):
        """
        Visit a Python AST :class:`ast.BoolOp` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the result that ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit(ast.parse(x).body[0].value)
            sage: [vis(d) for d in ['True and 1', 'False or 3 or None', '3 and 4']] #indirect doctest
            [1, 3, 4]
        """
        op = node.op.__class__.__name__
        L = list(node.values)
        out = self.visit(L.pop(0))
        if op == 'And':
            while L:
                next = self.visit(L.pop(0))
                out = out and next
            return out
        if op == 'Or':
            while L:
                next = self.visit(L.pop(0))
                out = out or next
            return out

    def visit_Compare(self, node):
        """
        Visit a Python AST :class:`ast.Compare` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the result that ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Compare(ast.parse(x).body[0].value)
            sage: [vis(d) for d in ['1<2==2!=3', '1==1>2', '1<2>1', '1<3<2<4']]
            [True, False, True, False]
        """
        left = self.visit(node.left)
        ops = list(node.ops)
        comparators = list(node.comparators)  # the things to be compared with.
        while ops:
            op = ops.pop(0).__class__.__name__
            right = self.visit(comparators.pop(0))
            if op == 'Lt':
                if not left < right:
                    return False
            elif op == 'LtE':
                if not left <= right:
                    return False
            elif op == 'Gt':
                if not left > right:
                    return False
            elif op == 'GtE':
                if not left >= right:
                    return False
            elif op == 'Eq':
                if not left == right:
                    return False
            elif op == 'NotEq':
                if not left != right:
                    return False
            left = right
        return True

    def visit_BinOp(self, node):
        """
        Visit a Python AST :class:`ast.BinOp` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the result that ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit(ast.parse(x).body[0].value)
            sage: [vis(d) for d in ['(3+(2*4))', '7|8', '5^3', '7/3', '7//3', '3<<4']] #indirect doctest
            [11, 15, 6, 2.3333333333333335, 2, 48]
        """
        op = node.op.__class__.__name__
        if op == 'Add':
            return self.visit(node.left) + self.visit(node.right)
        if op == 'Mult':
            return self.visit(node.left) * self.visit(node.right)
        if op == 'BitAnd':
            return self.visit(node.left) & self.visit(node.right)
        if op == 'BitOr':
            return self.visit(node.left) | self.visit(node.right)
        if op == 'BitXor':
            return self.visit(node.left) ^ self.visit(node.right)
        if op == 'Div':
            return self.visit(node.left) / self.visit(node.right)
        if op == 'Eq':
            return self.visit(node.left) == self.visit(node.right)
        if op == 'FloorDiv':
            return self.visit(node.left) // self.visit(node.right)
        if op == 'NotEq':
            return self.visit(node.left) != self.visit(node.right)
        if op == 'NotIn':
            return self.visit(node.left) not in self.visit(node.right)
        if op == 'Pow':
            return self.visit(node.left) ** self.visit(node.right)
        if op == 'RShift':
            return self.visit(node.left) >> self.visit(node.right)
        if op == 'LShift':
            return self.visit(node.left) << self.visit(node.right)
        if op == 'Sub':
            return self.visit(node.left) - self.visit(node.right)
        if op == 'Gt':
            return self.visit(node.left) > self.visit(node.right)
        if op == 'GtE':
            return self.visit(node.left) >= self.visit(node.right)
        if op == 'In':
            return self.visit(node.left) in self.visit(node.right)
        if op == 'Is':
            return self.visit(node.left) is self.visit(node.right)
        if op == 'IsNot':
            return self.visit(node.left) is not self.visit(node.right)
        if op == 'Lt':
            return self.visit(node.left) < self.visit(node.right)
        if op == 'LtE':
            return self.visit(node.left) <= self.visit(node.right)
        if op == 'Mod':
            return self.visit(node.left) % self.visit(node.right)

    def visit_UnaryOp(self, node):
        """
        Visit a Python AST :class:`ast.BinOp` node.

        INPUT:

        - ``node`` -- the node instance to visit

        OUTPUT: the result that ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_UnaryOp(ast.parse(x).body[0].value)
            sage: [vis(d) for d in ['+(3*2)', '-(3*2)']]
            [6, -6]
        """
        op = node.op.__class__.__name__
        if op == 'Not':
            return not self.visit(node.operand)
        if op == 'UAdd':
            return self.visit(node.operand)
        if op == 'USub':
            return -self.visit(node.operand)


def _grep_first_pair_of_parentheses(s):
    r"""
    Return the first matching pair of parentheses in a code string.

    INPUT:

    - ``s`` -- string

    OUTPUT:

    A substring of the input, namely the part between the first
    (outmost) matching pair of parentheses (including the
    parentheses).

    Parentheses between single or double quotation marks do not
    count. If no matching pair of parentheses can be found, a
    :exc:`SyntaxError` is raised.

    EXAMPLES::

        sage: from sage.misc.sageinspect import _grep_first_pair_of_parentheses
        sage: code = 'def foo(a="\'):", b=4):\n    return'
        sage: _grep_first_pair_of_parentheses(code)
        '(a="\'):", b=4)'
        sage: code = 'def foo(a="%s):", \'b=4):\n    return'%("'")
        sage: _grep_first_pair_of_parentheses(code)
        Traceback (most recent call last):
        ...
        SyntaxError: The given string does not contain balanced parentheses
    """
    out = []
    single_quote = False
    double_quote = False
    escaped = False
    level = 0
    for c in s:
        if level > 0:
            out.append(c)
        if c == '(' and not single_quote and not double_quote and not escaped:
            level += 1
        elif c == '"' and not single_quote and not escaped:
            double_quote = not double_quote
        elif c == "'" and not double_quote and not escaped:
            single_quote = not single_quote
        elif c == ')' and not single_quote and not double_quote and not escaped:
            if level == 1:
                return '(' + ''.join(out)
            level -= 1
        elif c == "\\" and (single_quote or double_quote):
            escaped = not escaped
        else:
            escaped = False
    raise SyntaxError("The given string does not contain balanced parentheses")


def _split_syntactical_unit(s):
    """
    Split off a sub-expression from the start of a given string.

    INPUT:

    - ``s`` -- string

    OUTPUT:

    A pair ``unit, s2``, such that ``unit`` is the string representation of a
    string (single or double quoted) or of a sub-expression surrounded by
    brackets (round, square or curly brackets), or of an identifier, or a
    single character, if none of the above is available. The given string ``s``
    is obtained by appending some whitespace followed by ``s2`` to ``unit``.

    Blank space between the units is removed.

    EXAMPLES::

        sage: from sage.misc.sageinspect import _split_syntactical_unit
        sage: s = "(Hel) lo_1=[)\"!\" ] '''? {world} '''?"
        sage: while s:
        ....:     u, s = _split_syntactical_unit(s)
        ....:     print(u)
        (Hel)
        lo_1
        =
        [)"!"]
        '''? {world} '''
        ?

    If the string ends before the unit is completed (mispatching parentheses
    or missing quotation mark), then a syntax error is raised::

        sage: s = "'''({SAGE}]"
        sage: _split_syntactical_unit(s)
        Traceback (most recent call last):
        ...
        SyntaxError: EOF while scanning string literal
        sage: s = "({SAGE}]"
        sage: _split_syntactical_unit(s)
        Traceback (most recent call last):
        ...
        SyntaxError: Syntactical group starting with '(' did not end with ')'

    Numbers are not recognised::

        sage: _split_syntactical_unit('123')
        ('1', '23')

    TESTS:

    The following was fixed in :issue:`16309`::

        sage: _split_syntactical_unit('()): pass')
        ('()', '): pass')
    """
    s = s.strip()
    if not s:
        return s

    # Split a given string at the next unescaped quotation mark
    def split_string(s, quot):
        escaped = False
        l = len(quot)
        for i in range(len(s)):
            if s[i] == '\\':
                escaped = not escaped
                continue
            if not escaped and s[i:i + l] == quot:
                return s[:i], s[i + l:]
            escaped = False
        raise SyntaxError("EOF while scanning string literal")
    # 1. s is a triple-quoted string
    if s.startswith('"""'):
        a, b = split_string(s[3:], '"""')
        return '"""' + a + '"""', b.strip()
    if s.startswith('r"""'):
        a, b = split_string(s[4:], '"""')
        return 'r"""'+a+'"""', b.strip()
    if s.startswith("'''"):
        a, b = split_string(s[3:], "'''")
        return "'''"+a+"'''", b.strip()
    if s.startswith("r'''"):
        a, b = split_string(s[4:], "'''")
        return "r'''"+a+"'''", b.strip()

    # 2. s is a single-quoted string
    if s.startswith('"'):
        a, b = split_string(s[1:], '"')
        return '"'+a+'"', b.strip()
    if s.startswith("'"):
        a, b = split_string(s[1:], "'")
        return "'"+a+"'", b.strip()
    if s.startswith('r"'):
        a, b = split_string(s[2:], '"')
        return 'r"'+a+'"', b.strip()
    if s.startswith("r'"):
        a, b = split_string(s[2:], "'")
        return "r'"+a+"'", b.strip()

    # 3. s is not a string
    start = s[0]
    out = [start]
    if start == '(':
        stop = ')'
    elif start == '[':
        stop = ']'
    elif start == '{':
        stop = '}'
    elif start == '\\':
        # note that python would raise a syntax error
        # if the line contains anything but whitespace
        # after the backslash. But we assume here that
        # the input is syntactically correct.
        return _split_syntactical_unit(s[1:])
    elif start == '#':
        linebreak = s.index(os.linesep)
        if linebreak == -1:
            return '', ''
        return '', s[linebreak:].strip()
    else:
        M = __identifier_re.search(s)
        if M is None:
            return s[0], s[1:].strip()
        return M.group(), s[M.end():].strip()

    s = s[1:]
    while s:
        tmp_group, s = _split_syntactical_unit(s)
        out.append(tmp_group)
        s = s.strip()
        if tmp_group == stop:
            return ''.join(out), s
        elif s.startswith(stop):
            out.append(stop)
            return ''.join(out), s[1:].strip()
    raise SyntaxError("Syntactical group starting with %s did not end with %s" % (repr(start), repr(stop)))


def _sage_getargspec_from_ast(source):
    r"""
    Return an argspec for a Python function or method by compiling its
    source to an abstract-syntax tree (AST) and walking its ``args``
    subtrees with :class:`SageArgSpecVisitor`.  We use this in
    :func:`_sage_getargspec_cython`.

    INPUT:

    - ``source`` -- string; the function's (or method's) source code
      definition. The function's body is ignored.

    OUTPUT: an instance of :obj:`inspect.ArgSpec`, i.e., a named tuple

    EXAMPLES::

        sage: import inspect, sage.misc.sageinspect as sms
        sage: from_ast = sms._sage_getargspec_from_ast
        sage: s = "def f(a, b=2, c={'a': [4, 5.5, False]}, d=(None, True)):\n    return"
        sage: from_ast(s)
        FullArgSpec(args=['a', 'b', 'c', 'd'], varargs=None, varkw=None, defaults=(2, {'a': [4, 5.5, False]}, (None, True)), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: context = {}
        sage: exec(compile(s, '<string>', 'single'), context)
        sage: inspect.getfullargspec(context['f'])
        FullArgSpec(args=['a', 'b', 'c', 'd'], varargs=None, varkw=None, defaults=(2, {'a': [4, 5.5, False]}, (None, True)), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: from_ast(s) == inspect.getfullargspec(context['f'])
        True
        sage: set(from_ast(sms.sage_getsource(x)) == inspect.getfullargspec(x) for x in [factor, identity_matrix, Graph.__init__])                              # needs sage.graphs sage.modules
        {True}
    """
    ast_args = ast.parse(source.lstrip()).body[0].args

    visitor = SageArgSpecVisitor()
    args = [visitor.visit(a) for a in ast_args.args]
    defaults = [visitor.visit(d) for d in ast_args.defaults]

    # vararg and kwarg may be None
    vararg = getattr(ast_args.vararg, 'arg', None)
    kwarg = getattr(ast_args.kwarg, 'arg', None)

    return inspect.FullArgSpec(args, vararg, kwarg,
                               tuple(defaults) if defaults else None,
                               kwonlyargs=[], kwonlydefaults=None, annotations={})


def _sage_getargspec_cython(source):
    r"""
    inspect.getargspec from source code.  That is, get the names and
    default values of a function's arguments.

    INPUT:

    - ``source`` -- string; the function's (or method's) source code
      definition.  The function's body is ignored. The definition may
      contain type definitions for the function arguments.

    OUTPUT: an instance of :class:`inspect.FullArgSpec`, i.e., a named tuple

    EXAMPLES::

        sage: from sage.misc.sageinspect import _sage_getargspec_cython as sgc
        sage: sgc("cpdef double abc(self, Element x=None, Parent base=0):")
        FullArgSpec(args=['self', 'x', 'base'], varargs=None, varkw=None, defaults=(None, 0), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sgc("def __init__(self, x=None, unsigned int base=0):")
        FullArgSpec(args=['self', 'x', 'base'], varargs=None, varkw=None, defaults=(None, 0), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sgc('def o(p, r={}, *q, **s) except? -1:')
        FullArgSpec(args=['p', 'r'], varargs='q', varkw='s', defaults=({},), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sgc('cpdef how(r=(None, "u:doing?")):')
        FullArgSpec(args=['r'], varargs=None, varkw=None, defaults=((None, 'u:doing?'),), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sgc('def _(x="):"):')
        FullArgSpec(args=['x'], varargs=None, varkw=None, defaults=('):',), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sgc('def f(z = {(1, 2, 3): True}):\n    return z')
        FullArgSpec(args=['z'], varargs=None, varkw=None, defaults=({(1, 2, 3): True},), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sgc('def f(double x, z = {(1, 2, 3): True}):\n    return z')
        FullArgSpec(args=['x', 'z'], varargs=None, varkw=None, defaults=({(1, 2, 3): True},), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sgc('def f(*args): pass')
        FullArgSpec(args=[], varargs='args', varkw=None, defaults=None, kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sgc('def f(**args): pass')
        FullArgSpec(args=[], varargs=None, varkw='args', defaults=None, kwonlyargs=[], kwonlydefaults=None, annotations={})

    Some malformed input is detected::

        sage: sgc('def f(x, y')
        Traceback (most recent call last):
        ...
        SyntaxError: Unexpected EOF while parsing argument list
        sage: sgc('def f(*x = 5, z = {(1,2,3): True}): pass')
        Traceback (most recent call last):
        ...
        SyntaxError: invalid ...
        sage: sgc('def f(int *x = 5, z = {(1,2,3): True}): pass')
        Traceback (most recent call last):
        ...
        SyntaxError: Pointer types not allowed in def or cpdef functions
        sage: sgc('def f(x = , z = {(1,2,3): True}): pass')
        Traceback (most recent call last):
        ...
        SyntaxError: Definition of a default argument expected
        sage: sgc('def f(int x = 5, , z = {(1,2,3): True}): pass')
        Traceback (most recent call last):
        ...
        SyntaxError: invalid ...

    TESTS:

    Some input that is malformed in Python 2 but well formed in Cython or
    Python 3 is correctly parsed::

        sage: def dummy_python(self, *args, x=1): pass
        sage: sgc("def dummy_python(self, *args, x=1): pass")
        FullArgSpec(args=['self', 'x'], varargs='args', varkw=None, defaults=(1,),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: cython("def dummy_cython(self, *args, x=1): pass")                        # needs sage.misc.cython
        sage: sgc("def dummy_cython(self, *args, x=1): pass")
        FullArgSpec(args=['self', 'x'], varargs='args', varkw=None, defaults=(1,),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})

    In some examples above, a syntax error was raised when a type
    definition contains a pointer. An exception is made for ``char*``,
    since C strings are acceptable input in public Cython functions::

        sage: sgc('def f(char *x = "a string", z = {(1,2,3): True}): pass')
        FullArgSpec(args=['x', 'z'], varargs=None, varkw=None,
                    defaults=('a string', {(1, 2, 3): True}),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})
    """
    defpos = source.find('def ')
    assert defpos > -1, "The given source does not contain 'def'"
    s = source[defpos:].strip()
    while s:
        if s.startswith('('):
            break
        _, s = _split_syntactical_unit(s)
    s = s[1:].strip()
    if not s:
        raise SyntaxError("Function definition must contain an argument list")

    # We remove the type declarations, build a dummy Python function, and
    # then call _get_argspec_from_ast. This should be
    # better than creating a complete parser for Cython syntax,
    # even though _split_syntactical_unit does part of the parsing work anyway.

    cy_units = []
    while not s.startswith(')'):
        if not s:
            raise SyntaxError("Unexpected EOF while parsing argument list")
        u, s = _split_syntactical_unit(s)
        cy_units.append(u)

    py_units = []
    name = None
    i = 0
    l = len(cy_units)
    expect_default = False
    nb_stars = 0
    varargs = None
    keywords = None
    while (i < l):
        unit = cy_units[i]
        if expect_default:
            if unit in ('=', '*', ','):
                raise SyntaxError("Definition of a default argument expected")
            while unit != ',':
                py_units.append(unit)
                i += 1
                if i == l:
                    break
                unit = cy_units[i]
            expect_default = False
            name = None
            if nb_stars:
                raise SyntaxError("The %s argument has no default" % ('varargs' if nb_stars == 1 else 'keywords'))
            continue
        i += 1
        if unit == '*':
            if name:
                if name != 'char':
                    raise SyntaxError("Pointer types not allowed in def or cpdef functions")
                else:
                    continue
            else:
                nb_stars += 1
            continue
        elif unit == ',':
            if expect_default:
                raise SyntaxError("Unexpected EOF while parsing argument list")
            name = None
            if nb_stars:
                nb_stars = 0
                continue
        elif unit == '=':
            expect_default = True
            name = None
            if nb_stars:
                raise SyntaxError("The %s argument has no default" % ('varargs' if nb_stars == 1 else 'keywords'))
        else:
            name = unit
        if name is not None:
            # Is "name" part of a type definition?
            # If it is the last identifier before '=' or ',',
            # then it *is* a variable name,
            if i == l or cy_units[i] in ('=', ','):
                if nb_stars == 0:
                    py_units.append(name)
                elif nb_stars == 1:
                    if varargs is None:
                        varargs = name
                        # skip the "=" or ",", since varargs
                        # is treated separately
                        i += 1
                        name = None
                        nb_stars = 0
                    else:
                        raise SyntaxError("varargs cannot be defined twice")
                elif nb_stars == 2:
                    if keywords is None:
                        keywords = name
                        # skip the "=" or ",", since varargs
                        # is treated separately
                        i += 1
                        name = None
                        nb_stars = 0
                    else:
                        raise SyntaxError("varargs cannot be defined twice")
                else:
                    raise SyntaxError("variable declaration comprises at most two '*'")
        else:
            py_units.append(unit)
    if varargs is None:
        varargs = ''
    elif not py_units or py_units[-1] == ',':
        varargs = '*' + varargs
    else:
        varargs = ',*' + varargs
    if keywords is None:
        keywords = ''
    elif varargs or (py_units and py_units[-1] != ','):
        keywords = ',**' + keywords
    else:
        keywords = '**' + keywords
    return _sage_getargspec_from_ast('def dummy(' + ''.join(py_units) +
                                     varargs + keywords + '): pass')


def sage_getfile(obj):
    r"""
    Get the full file name associated to ``obj`` as a string.

    INPUT:

    - ``obj`` -- a Sage object, module, etc.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getfile
        sage: sage_getfile(sage.rings.rational)
        '...sage/rings/rational.pyx'
        sage: from sage.algebras.steenrod.steenrod_algebra import Sq                    # needs sage.combinat sage.modules
        sage: sage_getfile(Sq)                                                          # needs sage.combinat sage.modules
        '...sage/algebras/steenrod/steenrod_algebra.py'
        sage: sage_getfile(x)                                                           # needs sage.symbolic
        '...sage/symbolic/expression.pyx'

    The following tests against some bugs fixed in :issue:`9976`::

        sage: obj = sage.combinat.partition_algebra.SetPartitionsAk                     # needs sage.combinat sage.modules
        sage: sage_getfile(obj)                                                         # needs sage.combinat sage.modules
        '...sage/combinat/partition_algebra.py'

    And here is another bug, fixed in :issue:`11298`::

        sage: P.<x,y> = QQ[]
        sage: sage_getfile(P)                                                           # needs sage.libs.singular
        '...sage/rings/polynomial/multi_polynomial_libsingular...'

    A problem fixed in :issue:`16309`::

        sage: cython(                                                                   # needs sage.misc.cython
        ....: '''
        ....: class Bar: pass
        ....: cdef class Foo: pass
        ....: ''')
        sage: sage_getfile(Bar)                                                         # needs sage.misc.cython
        '...pyx'
        sage: sage_getfile(Foo)                                                         # needs sage.misc.cython
        '...pyx'

    By :issue:`18249`, we return an empty string for Python builtins. In that
    way, there is no error when the user types, for example, ``range?``::

        sage: sage_getfile(range)
        ''
    """
    # We try to extract from docstrings, but not using Python's inspect
    # because _sage_getdoc_unformatted is more robust.
    d = _sage_getdoc_unformatted(obj)
    pos = _extract_embedded_position(d)
    if pos is not None:
        (_, filename, _) = pos
        return filename

    # The instance case
    if isclassinstance(obj):
        if isinstance(obj, functools.partial):
            return sage_getfile(obj.func)
        return sage_getfile(obj.__class__)  # inspect.getabsfile(obj.__class__)

    # No go? fall back to inspect.
    try:
        sourcefile = inspect.getabsfile(obj)
    except TypeError:  # this happens for Python builtins
        return ''
    for suffix in import_machinery.EXTENSION_SUFFIXES:
        if sourcefile.endswith(suffix):
            return sourcefile[:-len(suffix)]+os.path.extsep+'pyx'
    return sourcefile


def sage_getfile_relative(obj):
    r"""
    Get the file name associated to ``obj`` as a string.

    This is the same as :func:`sage_getfile`, but
    if the source file is part of the ``sage.*`` namespace, it
    makes the file name relative so that it starts with ``sage/``.

    INPUT:

    - ``obj`` -- a Sage object, module, etc.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getfile_relative
        sage: sage_getfile_relative(sage.rings.rational)
        'sage/rings/rational.pyx'
        sage: from sage.algebras.steenrod.steenrod_algebra import Sq                    # needs sage.combinat sage.modules
        sage: sage_getfile_relative(Sq)                                                 # needs sage.combinat sage.modules
        'sage/algebras/steenrod/steenrod_algebra.py'
        sage: sage_getfile_relative(x)                                                  # needs sage.symbolic
        'sage/symbolic/expression.pyx'
        sage: sage_getfile_relative(range)
        ''
    """
    filename = sage_getfile(obj)
    if not filename:
        return filename

    from os.path import relpath, normpath, commonprefix

    def directories():
        try:
            from sage.env import SAGE_SRC
        except ImportError:
            pass
        else:
            if SAGE_SRC:
                yield normpath(os.path.join(SAGE_SRC, 'sage'))
        import sage
        yield from sage.__path__

    for directory in directories():
        if commonprefix([filename, directory]) == directory:
            return os.path.join('sage', relpath(filename, directory))

    return filename


def sage_getargspec(obj):
    r"""
    Return the names and default values of a function's arguments.

    INPUT:

    - ``obj`` -- any callable object

    OUTPUT:

    A named tuple :class:`FullArgSpec` is returned, as specified by the
    Python library function :func:`inspect.getfullargspec`.

    NOTE:

    If the object has a method ``_sage_argspec_``, then the output of
    that method is transformed into a named tuple and then returned.

    If a class instance has a method ``_sage_src_``, then its output
    is studied to determine the argspec. This is because currently
    the :class:`~sage.misc.cachefunc.CachedMethod` decorator has
    no ``_sage_argspec_`` method.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getargspec
        sage: def f(x, y, z=1, t=2, *args, **keywords):
        ....:     pass
        sage: sage_getargspec(f)
        FullArgSpec(args=['x', 'y', 'z', 't'], varargs='args', varkw='keywords',
                    defaults=(1, 2), kwonlyargs=[], kwonlydefaults=None, annotations={})

    We now run :func:`sage_getargspec` on some functions from the Sage library::

        sage: sage_getargspec(identity_matrix)                                          # needs sage.modules
        FullArgSpec(args=['ring', 'n', 'sparse'], varargs=None, varkw=None,
                    defaults=(0, False), kwonlyargs=[], kwonlydefaults=None,
                    annotations={})
        sage: sage_getargspec(factor)
        FullArgSpec(args=['n', 'proof', 'int_', 'algorithm', 'verbose'],
                    varargs=None, varkw='kwds', defaults=(None, False, 'pari', 0),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})

    In the case of a class or a class instance, the :class:`FullArgSpec` of the
    ``__new__``, ``__init__`` or ``__call__`` method is returned::

        sage: P.<x,y> = QQ[]
        sage: sage_getargspec(P)                                                        # needs sage.libs.singular
        FullArgSpec(args=['base_ring', 'n', 'names', 'order'],
                    varargs=None, varkw=None, defaults=('degrevlex',),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sage_getargspec(P.__class__)                                              # needs sage.libs.singular
        FullArgSpec(args=['self', 'x'], varargs='args', varkw='kwds', defaults=(0,),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})

    The following tests against various bugs that were fixed in
    :issue:`9976`::

        sage: from sage.rings.polynomial.real_roots import bernstein_polynomial_factory_ratlist     # needs sage.modules
        sage: sage_getargspec(bernstein_polynomial_factory_ratlist.coeffs_bitsize)                  # needs sage.modules
        FullArgSpec(args=['self'], varargs=None, varkw=None, defaults=None,
                    kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: from sage.rings.polynomial.pbori.pbori import BooleanMonomialMonoid       # needs sage.rings.polynomial.pbori
        sage: sage_getargspec(BooleanMonomialMonoid.gen)                                # needs sage.rings.polynomial.pbori
        FullArgSpec(args=['self', 'i'], varargs=None, varkw=None, defaults=(0,),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: I = P*[x,y]
        sage: sage_getargspec(I.groebner_basis)                                         # needs sage.libs.singular
        FullArgSpec(args=['self', 'algorithm', 'deg_bound', 'mult_bound', 'prot'],
                    varargs='args', varkw='kwds', defaults=('', None, None, False),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: cython("cpdef int foo(x,y) except -1: return 1")                          # needs sage.misc.cython
        sage: sage_getargspec(foo)                                                      # needs sage.misc.cython
        FullArgSpec(args=['x', 'y'], varargs=None, varkw=None, defaults=None,
                    kwonlyargs=[], kwonlydefaults=None, annotations={})

    If a :func:`functools.partial` instance is involved, we see no other meaningful solution
    than to return the argspec of the underlying function::

        sage: def f(a, b, c, d=1):
        ....:     return a + b + c + d
        sage: import functools
        sage: f1 = functools.partial(f, 1, c=2)
        sage: sage_getargspec(f1)
        FullArgSpec(args=['a', 'b', 'c', 'd'], varargs=None, varkw=None, defaults=(1,),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})

    TESTS:

    By :issue:`9976`, rather complicated cases work. In the
    following example, we dynamically create an extension class
    that returns some source code, and the example shows that
    the source code is taken for granted, i.e., the argspec of
    an instance of that class does not coincide with the argspec
    of its call method. That behaviour is intended, since a
    decorated method appears to have the generic signature
    ``*args, **kwds``, but in fact it is only supposed to be called
    with the arguments requested by the underlying undecorated
    method. We saw an easy example above, namely ``I.groebner_basis``.
    Here is a more difficult one::

        sage: # needs sage.misc.cython
        sage: cython_code = [
        ....: 'cdef class MyClass:',
        ....: '    def _sage_src_(self):',
        ....: '        return "def foo(x, a=\\\')\\\"\\\', b={(2+1):\\\'bar\\\', not 1:3, 3<<4:5}): return\\n"',
        ....: '    def __call__(self, m, n): return "something"']
        sage: cython('\n'.join(cython_code))
        sage: O = MyClass()
        sage: print(sage.misc.sageinspect.sage_getsource(O))
        def foo(x, a=')"', b={(2+1):'bar', not 1:3, 3<<4:5}): return
        sage: spec = sage.misc.sageinspect.sage_getargspec(O)
        sage: spec.args, spec.varargs, spec.varkw
        (['x', 'a', 'b'], None, None)
        sage: spec.defaults[0]
        ')"'
        sage: sorted(spec.defaults[1].items(), key=lambda x: str(x))
        [(3, 'bar'), (48, 5), (False, 3)]
        sage: sage.misc.sageinspect.sage_getargspec(O.__call__)
        FullArgSpec(args=['self', 'm', 'n'], varargs=None, varkw=None, defaults=None,
                    kwonlyargs=[], kwonlydefaults=None, annotations={})

    ::

        sage: cython('def foo(x, a=\'\\\')"\', b={not (2+1==3):\'bar\'}): return')      # needs sage.misc.cython
        sage: print(sage.misc.sageinspect.sage_getsource(foo))                          # needs sage.misc.cython
        def foo(x, a='\')"', b={not (2+1==3):'bar'}): return
        <BLANKLINE>
        sage: sage.misc.sageinspect.sage_getargspec(foo)                                # needs sage.misc.cython
        FullArgSpec(args=['x', 'a', 'b'], varargs=None, varkw=None,
                    defaults=('\')"', {False: 'bar'}),
                    kwonlyargs=[], kwonlydefaults=None, annotations={})

    The following produced a syntax error before the patch at :issue:`11913`,
    see also :issue:`26906`::

        sage: sage.misc.sageinspect.sage_getargspec(r.lm)       # optional - rpy2
        FullArgSpec(args=['self'], varargs='args', varkw='kwds', defaults=None,
                    kwonlyargs=[], kwonlydefaults=None, annotations={})

    The following was fixed in :issue:`16309`::

        sage: # needs sage.misc.cython
        sage: cython(
        ....: '''
        ....: class Foo:
        ....:     @staticmethod
        ....:     def join(categories, bint as_list=False, tuple ignore_axioms=(), tuple axioms=()): pass
        ....: cdef class Bar:
        ....:     @staticmethod
        ....:     def join(categories, bint as_list=False, tuple ignore_axioms=(), tuple axioms=()): pass
        ....:     cpdef meet(categories, bint as_list=False, tuple ignore_axioms=(), tuple axioms=()): pass
        ....: ''')
        sage: sage_getargspec(Foo.join)
        FullArgSpec(args=['categories', 'as_list', 'ignore_axioms', 'axioms'], varargs=None, varkw=None,
                    defaults=(False, (), ()), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sage_getargspec(Bar.join)
        FullArgSpec(args=['categories', 'as_list', 'ignore_axioms', 'axioms'], varargs=None, varkw=None,
                    defaults=(False, (), ()), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: sage_getargspec(Bar.meet)
        FullArgSpec(args=['categories', 'as_list', 'ignore_axioms', 'axioms'], varargs=None, varkw=None,
                    defaults=(False, (), ()), kwonlyargs=[], kwonlydefaults=None, annotations={})

    Test that :issue:`17009` is fixed::

        sage: sage_getargspec(gap)                                                      # needs sage.libs.gap
        FullArgSpec(args=['self', 'x', 'name'], varargs=None, varkw=None,
                    defaults=(None,), kwonlyargs=[], kwonlydefaults=None, annotations={})

    By :issue:`17814`, the following gives the correct answer (previously, the
    defaults would have been found ``None``)::

        sage: from sage.misc.nested_class import MainClass
        sage: sage_getargspec(MainClass.NestedClass.NestedSubClass.dummy)
        FullArgSpec(args=['self', 'x', 'r'], varargs='args', varkw='kwds',
                    defaults=((1, 2, 3.4),), kwonlyargs=[], kwonlydefaults=None, annotations={})

    In :issue:`18249` was decided to return a generic signature for Python
    builtin functions, rather than to raise an error (which is what Python's
    inspect module does)::

        sage: import inspect
        sage: sage_getargspec(range)
        FullArgSpec(args=[], varargs='args', varkw='kwds', defaults=None, kwonlyargs=[], kwonlydefaults=None, annotations={})

    Test that :issue:`28524` is fixed::

        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: shell.run_cell(
        ....:     'class Foo:\n'
        ....:     '    def __call__(self):\n'
        ....:     '        return None\n'
        ....:     '    def __module__(self):\n'
        ....:     '        return "sage.misc.sageinspect"\n'
        ....:     '    def _sage_src_(self):\n'
        ....:     '        return "the source code string"')
        sage: shell.run_cell('f = Foo()')
        sage: shell.run_cell('f??')
        ...the source code string...
    """
    from sage.misc.lazy_attribute import lazy_attribute
    from sage.misc.abstract_method import AbstractMethod
    if inspect.isclass(obj):
        return sage_getargspec(obj.__call__)
    if isinstance(obj, (lazy_attribute, AbstractMethod)):
        source = sage_getsource(obj)
        return inspect.FullArgSpec(*_sage_getargspec_cython(source))
    if not callable(obj):
        raise TypeError("obj is not a code object")
    try:
        return inspect.FullArgSpec(*obj._sage_argspec_())
    except (AttributeError, TypeError):
        pass
    # If we are lucky, the function signature is embedded in the docstring.
    docstring = _sage_getdoc_unformatted(obj)
    try:
        name = obj.__name__
    except AttributeError:
        name = type(obj).__name__
    argspec = _extract_embedded_signature(docstring, name)[1]
    if argspec is not None:
        return argspec
    if hasattr(obj, '__code__'):
        # Note that this may give a wrong result for the constants!
        try:
            args, varargs, varkw = inspect.getargs(obj.__code__)
            return inspect.FullArgSpec(args, varargs, varkw, obj.__defaults__,
                                       kwonlyargs=[], kwonlydefaults=None, annotations={})
        except (TypeError, AttributeError):
            pass
    if isclassinstance(obj):
        if hasattr(obj, '_sage_src_'):  # it may be a decorator!
            source = sage_getsource(obj)
            try:
                # we try to find the definition and parse it by
                # _sage_getargspec_ast
                proxy = 'def dummy' + _grep_first_pair_of_parentheses(source) \
                        + ':\n    return'
                return _sage_getargspec_from_ast(proxy)
            except SyntaxError:
                # To fix trac #10860. See #11913 for more information.
                # See also #26906 and #28524.
                pass
        if isinstance(obj, functools.partial):
            base_spec = sage_getargspec(obj.func)
            return base_spec
        return sage_getargspec(obj.__class__.__call__)
    elif (hasattr(obj, '__objclass__') and hasattr(obj, '__name__') and
          obj.__name__ == 'next'):
        # Handle sage.rings.ring.FiniteFieldIterator.next and similar
        # slot wrappers.  This is mainly to suppress Sphinx warnings.
        return ['self'], None, None, None
    else:
        # We try to get the argspec by reading the source, which may be
        # expensive, but should only be needed for functions defined outside
        # of the Sage library (since otherwise the signature should be
        # embedded in the docstring)
        try:
            source = sage_getsource(obj)
        except TypeError:  # happens for Python builtins
            source = ''
        if source:
            return inspect.FullArgSpec(*_sage_getargspec_cython(source))
        else:
            func_obj = obj

    # Otherwise we're (hopefully!) plain Python, so use inspect
    try:
        args, varargs, varkw = inspect.getargs(func_obj.__code__)
    except AttributeError:
        try:
            args, varargs, varkw = inspect.getargs(func_obj)
        except TypeError:  # arg is not a code object
            # The above "hopefully" was wishful thinking:
            try:
                return inspect.FullArgSpec(*_sage_getargspec_cython(sage_getsource(obj)))
            except TypeError:  # This happens for Python builtins
                # The best we can do is to return a generic argspec
                args = []
                varargs = 'args'
                varkw = 'kwds'
    try:
        defaults = func_obj.__defaults__
    except AttributeError:
        defaults = None
    return inspect.FullArgSpec(args, varargs, varkw, defaults,
                               kwonlyargs=[], kwonlydefaults=None, annotations={})


def _fullargspec_to_signature(fullargspec):
    """
    Converts a :class:`FullArgSpec` instance to a :class:`Signature` instance by best effort.

    EXAMPLES::

        sage: from sage.misc.sageinspect import _fullargspec_to_signature
        sage: from inspect import FullArgSpec
        sage: fullargspec = FullArgSpec(args=['self', 'x', 'base'], varargs=None, varkw=None, defaults=(None, 0), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (self, x=None, base=0)>

    TESTS::

        sage: fullargspec = FullArgSpec(args=['p', 'r'], varargs='q', varkw='s', defaults=({},), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (p, r={}, *q, **s)>
        sage: fullargspec = FullArgSpec(args=['r'], varargs=None, varkw=None, defaults=((None, 'u:doing?'),), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (r=(None, 'u:doing?'))>
        sage: fullargspec = FullArgSpec(args=['x'], varargs=None, varkw=None, defaults=('):',), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (x='):')>
        sage: fullargspec = FullArgSpec(args=['z'], varargs=None, varkw=None, defaults=({(1, 2, 3): True},), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (z={(1, 2, 3): True})>
        sage: fullargspec = FullArgSpec(args=['x', 'z'], varargs=None, varkw=None, defaults=({(1, 2, 3): True},), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (x, z={(1, 2, 3): True})>
        sage: fullargspec = FullArgSpec(args=[], varargs='args', varkw=None, defaults=None, kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (*args)>
        sage: fullargspec = FullArgSpec(args=[], varargs=None, varkw='args', defaults=None, kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (**args)>
        sage: fullargspec = FullArgSpec(args=['self', 'x'], varargs='args', varkw=None, defaults=(1,), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (self, x=1, *args)>
        sage: fullargspec = FullArgSpec(args=['self', 'x'], varargs='args', varkw=None, defaults=(1,), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (self, x=1, *args)>
        sage: fullargspec = FullArgSpec(args=['x', 'z'], varargs=None, varkw=None, defaults=('a string', {(1, 2, 3): True}), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _fullargspec_to_signature(fullargspec)
        <Signature (x='a string', z={(1, 2, 3): True})>
        sage: _fullargspec_to_signature(FullArgSpec(args=['a'], varargs=None, varkw=None, defaults=None, kwonlyargs=['b', 'c'], kwonlydefaults={'b': 1}, annotations={}))
        <Signature (a, *, b=1, c)>
        sage: _fullargspec_to_signature(FullArgSpec(args=['a', 'b'], varargs=None, varkw=None, defaults=(1,), kwonlyargs=['c'], kwonlydefaults=None, annotations={}))
        <Signature (a, b=1, *, c)>
    """
    parameters = []
    defaults_start = len(fullargspec.args) - len(fullargspec.defaults) if fullargspec.defaults else None

    for i, arg in enumerate(fullargspec.args):
        default = fullargspec.defaults[i - defaults_start] if defaults_start is not None and i >= defaults_start else Parameter.empty
        param = Parameter(arg, Parameter.POSITIONAL_OR_KEYWORD, default=default)
        parameters.append(param)

    if fullargspec.varargs:
        param = Parameter(fullargspec.varargs, Parameter.VAR_POSITIONAL)
        parameters.append(param)

    if fullargspec.varkw:
        param = Parameter(fullargspec.varkw, Parameter.VAR_KEYWORD)
        parameters.append(param)

    for arg in fullargspec.kwonlyargs:
        param = Parameter(arg, Parameter.KEYWORD_ONLY, default=Parameter.empty if fullargspec.kwonlydefaults is None else
                          fullargspec.kwonlydefaults.get(arg, Parameter.empty))
        parameters.append(param)

    return Signature(parameters)


def sage_signature(obj):
    r"""
    Return the names and default values of a function's arguments.

    INPUT:

    - ``obj`` -- any callable object

    OUTPUT:

    A :class:`Signature` is returned, as specified by the
    Python library function :func:`inspect.signature`.


    .. NOTE::

        Currently the type information is not returned, because the output
        is converted from the return value of :func:`sage_getargspec`.
        This should be changed in the future.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_signature
        sage: def f(x, y, z=1, t=2, *args, **keywords):
        ....:     pass
        sage: sage_signature(f)
        <Signature (x, y, z=1, t=2, *args, **keywords)>

    We now run :func:`sage_signature` on some functions from the Sage library::

        sage: sage_signature(identity_matrix)                                          # needs sage.modules
        <Signature (ring, n=0, sparse=False)>
        sage: sage_signature(factor)
        <Signature (n, proof=None, int_=False, algorithm='pari', verbose=0, **kwds)>

    In the case of a class or a class instance, the `Signature` of the
    `__new__`, `__init__` or `__call__` method is returned::

        sage: P.<x,y> = QQ[]
        sage: sage_signature(P)                                                        # needs sage.libs.singular
        <Signature (base_ring, n, names, order='degrevlex')>
        sage: sage_signature(P.__class__)                                              # needs sage.libs.singular
        <Signature (self, x=0, *args, **kwds)>

    The following tests against various bugs that were fixed in
    :issue:`9976`::

        sage: from sage.rings.polynomial.real_roots import bernstein_polynomial_factory_ratlist     # needs sage.modules
        sage: sage_signature(bernstein_polynomial_factory_ratlist.coeffs_bitsize)                  # needs sage.modules
        <Signature (self)>
        sage: from sage.rings.polynomial.pbori.pbori import BooleanMonomialMonoid       # needs sage.rings.polynomial.pbori
        sage: sage_signature(BooleanMonomialMonoid.gen)                                # needs sage.rings.polynomial.pbori
        <Signature (self, i=0)>
        sage: I = P*[x,y]
        sage: sage_signature(I.groebner_basis)                                         # needs sage.libs.singular
        <Signature (self, algorithm='', deg_bound=None, mult_bound=None, prot=False, *args, **kwds)>
        sage: cython("cpdef int foo(x,y) except -1: return 1")                          # needs sage.misc.cython
        sage: sage_signature(foo)                                                      # needs sage.misc.cython
        <Signature (x, y)>

    If a `functools.partial` instance is involved, we see no other meaningful solution
    than to return the signature of the underlying function::

        sage: def f(a, b, c, d=1):
        ....:     return a + b + c + d
        sage: import functools
        sage: f1 = functools.partial(f, 1, c=2)
        sage: sage_signature(f1)
        <Signature (a, b, c, d=1)>
    """
    return _fullargspec_to_signature(sage_getargspec(obj))


def formatannotation(annotation, base_module=None):
    """
    This is taken from Python 3.7's inspect.py; the only change is to
    add documentation.

    INPUT:

    - ``annotation`` -- annotation for a function
    - ``base_module`` -- (default: ``None``)

    This is only relevant with Python 3, so the doctests are marked
    accordingly.

    EXAMPLES::

        sage: from sage.misc.sageinspect import formatannotation
        sage: import inspect
        sage: def foo(a, *, b:int, **kwargs):
        ....:     pass
        sage: s = inspect.signature(foo)

        sage: a = s.parameters['a'].annotation
        sage: a
        <class 'inspect._empty'>
        sage: formatannotation(a)
        'inspect._empty'

        sage: b = s.parameters['b'].annotation
        sage: b
        <class 'int'>
        sage: formatannotation(b)
        'int'
    """
    if getattr(annotation, '__module__', None) == 'typing':
        return repr(annotation).replace('typing.', '')
    if isinstance(annotation, type):
        if annotation.__module__ in ('builtins', base_module):
            return annotation.__qualname__
        return annotation.__module__ + '.' + annotation.__qualname__
    return repr(annotation)


_formatannotation = formatannotation


def sage_formatargspec(args, varargs=None, varkw=None, defaults=None,
                       kwonlyargs=(), kwonlydefaults=None, annotations={},
                       formatarg=str,
                       formatvarargs=None,
                       formatvarkw=None,
                       formatvalue=None,
                       formatreturns=None,
                       formatannotation=None):
    """
    Format an argument spec from the values returned by getfullargspec.

    The first seven arguments are (args, varargs, varkw, defaults,
    kwonlyargs, kwonlydefaults, annotations).  The other five arguments
    are the corresponding optional formatting functions that are called to
    turn names and values into strings.  The last argument is an optional
    function to format the sequence of arguments.

    This is taken from Python 3.7's inspect.py, where it is
    deprecated. The only change, aside from documentation (this
    paragraph and the next, plus doctests), is to remove the
    deprecation warning.

    Sage uses this function to format arguments, as obtained by
    :func:`sage_getargspec`. Since :func:`sage_getargspec` works for
    Cython functions while Python's inspect module does not, it makes
    sense to keep this function for formatting instances of
    ``inspect.FullArgSpec``.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_formatargspec
        sage: args = ['a', 'b', 'c']
        sage: defaults = [3]
        sage: sage_formatargspec(args, defaults=defaults)
        '(a, b, c=3)'
    """
    if formatvarargs is None:
        formatvarargs = lambda name: '*' + name
    if formatvarkw is None:
        formatvarkw = lambda name: '**' + name
    if formatvalue is None:
        formatvalue = lambda value: '=' + repr(value)
    if formatreturns is None:
        formatreturns = lambda text: ' -> ' + text
    if formatannotation is None:
        formatannotation = _formatannotation

    def formatargandannotation(arg):
        result = formatarg(arg)
        if arg in annotations:
            result += ': ' + formatannotation(annotations[arg])
        return result
    specs = []
    if defaults:
        firstdefault = len(args) - len(defaults)
    for i, arg in enumerate(args):
        spec = formatargandannotation(arg)
        if defaults and i >= firstdefault:
            spec = spec + formatvalue(defaults[i - firstdefault])
        specs.append(spec)
    if varargs is not None:
        specs.append(formatvarargs(formatargandannotation(varargs)))
    else:
        if kwonlyargs:
            specs.append('*')
    if kwonlyargs:
        for kwonlyarg in kwonlyargs:
            spec = formatargandannotation(kwonlyarg)
            if kwonlydefaults and kwonlyarg in kwonlydefaults:
                spec += formatvalue(kwonlydefaults[kwonlyarg])
            specs.append(spec)
    if varkw is not None:
        specs.append(formatvarkw(formatargandannotation(varkw)))
    result = '(' + ', '.join(specs) + ')'
    if 'return' in annotations:
        result += formatreturns(formatannotation(annotations['return']))
    return result


def sage_getdef(obj, obj_name=''):
    r"""
    Return the definition header for any callable object.

    INPUT:

    - ``obj`` -- function
    - ``obj_name`` -- string (default: ``''``); prepended to the output

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getdef
        sage: sage_getdef(identity_matrix)                                              # needs sage.modules
        '(ring, n=0, sparse=False)'
        sage: sage_getdef(identity_matrix, 'identity_matrix')                           # needs sage.modules
        'identity_matrix(ring, n=0, sparse=False)'

    Check that :issue:`6848` has been fixed::

        sage: sage_getdef(RDF.random_element)
        '(min=-1, max=1)'

    If an exception is generated, None is returned instead and the
    exception is suppressed.
    """
    try:
        spec = sage_getargspec(obj)
        s = str(sage_formatargspec(*spec))
        s = s.strip('(').strip(')').strip()
        if s[:4] == 'self':
            s = s[4:]
        s = s.lstrip(',').strip()
        # for use with typesetting the definition with the notebook:
        # sometimes s contains "*args" or "**keywds", and the
        # asterisks confuse ReST/sphinx/docutils, so escape them:
        # change * to \*, and change ** to \**.
        return obj_name + '(' + s + ')'
    except (AttributeError, TypeError, ValueError):
        return '%s( [noargspec] )' % obj_name


def _sage_getdoc_unformatted(obj):
    r"""
    Return the unformatted docstring associated to ``obj`` as a
    string.

    If ``obj`` is a Cython object with an embedded position in its
    docstring, the embedded position is **not** stripped.

    INPUT:

    - ``obj`` -- a function, module, etc.: something with a docstring

    EXAMPLES::

        sage: from sage.misc.sageinspect import _sage_getdoc_unformatted
        sage: print(_sage_getdoc_unformatted(sage.rings.integer.Integer))
        Integer(x=None, base=0)
        File: ...sage/rings/integer.pyx (starting at line ...)
        <BLANKLINE>
            The :class:`Integer` class represents arbitrary precision
            integers. It derives from the :class:`Element` class, so
            integers can be used as ring elements anywhere in Sage.
        ...

    TESTS:

    Test that we suppress useless built-in output (:issue:`3342`)::

        sage: from sage.misc.sageinspect import _sage_getdoc_unformatted
        sage: _sage_getdoc_unformatted(isinstance.__class__)
        ''

    """
    if obj is None:
        return ''
    r = obj.__doc__

    # Check if the __doc__ attribute was actually a string, and
    # not a 'getset_descriptor' or similar.
    if isinstance(r, str):
        return r
    else:
        # Not a string of any kind
        return ''


def sage_getdoc_original(obj):
    r"""
    Return the unformatted docstring associated to ``obj`` as a
    string.

    If ``obj`` is a Cython object with an embedded position or signature in
    its docstring, the embedded information is stripped. If the stripped
    docstring is empty, then the stripped docstring of ``obj.__init__`` is
    returned instead.

    Feed the results from this into the function
    :func:`sage.misc.sagedoc.format` for printing to the screen.

    INPUT:

    - ``obj`` -- a function, module, etc.: something with a docstring

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getdoc_original

    Here is a class that has its own docstring::

        sage: print(sage_getdoc_original(sage.rings.integer.Integer))
        <BLANKLINE>
            The :class:`Integer` class represents arbitrary precision
            integers. It derives from the :class:`Element` class, so
            integers can be used as ring elements anywhere in Sage.
        ...

    If the class does not have a docstring, the docstring of the
    ``__init__`` method is used, but not the ``__init__`` method
    of the base class (this was fixed in :issue:`24936`)::

        sage: from sage.categories.category import Category
        sage: class A(Category):
        ....:     def __init__(self):
        ....:         '''The __init__ docstring'''
        sage: sage_getdoc_original(A)
        'The __init__ docstring'
        sage: class B(Category):
        ....:     pass
        sage: sage_getdoc_original(B)
        ''

    Old-style classes are supported::

        sage: class OldStyleClass:
        ....:     def __init__(self):
        ....:         '''The __init__ docstring'''
        ....:         pass
        sage: print(sage_getdoc_original(OldStyleClass))
        The __init__ docstring

    When there is no ``__init__`` method, we just get an empty string::

        sage: class OldStyleClass:
        ....:     pass
        sage: sage_getdoc_original(OldStyleClass)
        ''

    If an instance of a class does not have its own docstring, the docstring
    of its class results::

        sage: sage_getdoc_original(sage.plot.colors.aliceblue) == sage_getdoc_original(sage.plot.colors.Color)          # needs sage.plot
        True
    """
    # typ is the type corresponding to obj, which is obj itself if
    # that was a type or old-style class
    if isinstance(obj, type):
        typ = obj
    else:
        typ = type(obj)

    s, argspec = _extract_embedded_signature(_sage_getdoc_unformatted(obj), typ.__name__)
    if s:
        pos = _extract_embedded_position(s)
        if pos is not None:
            s = pos[0]
    if not s:
        # The docstring of obj is empty. To get something, we want to use
        # the documentation of the __init__ method, but only if it belongs
        # to (the type of) obj.
        init = typ.__dict__.get("__init__")
        if init:
            return sage_getdoc_original(init)
    return s


def sage_getdoc(obj, obj_name='', embedded=False):
    r"""
    Return the docstring associated to ``obj`` as a string.

    If ``obj`` is a Cython object with an embedded position in its
    docstring, the embedded position is stripped.

    The optional boolean argument ``embedded`` controls the
    string formatting. It is False by default.

    INPUT:

    - ``obj`` -- a function, module, etc.: something with a docstring

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getdoc
        sage: sage_getdoc(identity_matrix)[87:124]                                      # needs sage.modules
        '...the n x n identity matrix...'
        sage: def f(a, b, c, d=1): return a+b+c+d
        ...
        sage: import functools
        sage: f1 = functools.partial(f, 1,c=2)
        sage: f.__doc__ = "original documentation"
        sage: f1.__doc__ = "specialised documentation"
        sage: sage_getdoc(f)
        'original documentation\n'
        sage: sage_getdoc(f1)
        'specialised documentation\n'
    """
    import sage.misc.sagedoc
    if obj is None:
        return ''
    r = sage_getdoc_original(obj)
    s = sage.misc.sagedoc.format(r, embedded=embedded)
    f = sage_getfile(obj)
    if f and os.path.exists(f):
        from sage.doctest.control import skipfile
        skip = skipfile(f)
        if isinstance(skip, str):
            warn = """WARNING: the enclosing module is marked '{}',
so doctests may not pass.""".format(skip)
            s = warn + "\n\n" + s

    # Fix object naming
    if obj_name != '':
        i = obj_name.find('.')
        if i != -1:
            obj_name = obj_name[:i]
        s = s.replace('self.', '%s.' % obj_name)

    return s


def sage_getsource(obj):
    r"""
    Return the source code associated to obj as a string, or None.

    INPUT:

    - ``obj`` -- function, etc.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getsource
        sage: sage_getsource(identity_matrix)[19:60]                                    # needs sage.modules
        'identity_matrix(ring, n=0, sparse=False):'
        sage: sage_getsource(identity_matrix)[19:60]                                    # needs sage.modules
        'identity_matrix(ring, n=0, sparse=False):'
    """
    # First we should check if the object has a _sage_src_
    # method.  If it does, we just return the output from
    # that.  This is useful for getting pexpect interface
    # elements to behave similar to regular Python objects
    # with respect to introspection.
    try:
        return obj._sage_src_()
    except (AttributeError, TypeError):
        pass

    t = sage_getsourcelines(obj)
    if not t:
        return None
    (source_lines, lineno) = t
    return ''.join(source_lines)


def _sage_getsourcelines_name_with_dot(obj):
    r"""
    Get the source lines of an object whose name
    contains a dot and whose source lines can not
    be obtained by different methods.

    EXAMPLES::

        sage: C = Rings()
        sage: from sage.misc.sageinspect import sage_getsource
        sage: print(sage_getsource(C.parent_class))  #indirect doctest
        class ParentMethods:
        ...
                Return the Lie bracket `[x, y] = x y - y x` of `x` and `y`.
        ...

    TESTS:

    The following was fixed in :issue:`16309`::

        sage: # needs sage.misc.cython
        sage: cython(
        ....: '''
        ....: class A:
        ....:     def __init__(self):
        ....:         "some init doc"
        ....:         pass
        ....: class B:
        ....:     "some class doc"
        ....:     class A(A):
        ....:         pass
        ....: ''')
        sage: B.A.__name__
        'A'
        sage: B.A.__qualname__
        'B.A'
        sage: sage_getsource(B.A)
        '    class A(A):\n        pass\n\n'

    Note that for this example to work, it is essential that the class ``B``
    has a docstring. Otherwise, the code of ``B`` could not be found (Cython
    inserts embedding information into the docstring) and thus the code of
    ``B.A`` couldn't be found either.
    """
    # First, split the name:
    if '.' in obj.__name__:
        splitted_name = obj.__name__.split('.')
    elif hasattr(obj, '__qualname__'):
        splitted_name = obj.__qualname__.split('.')
    else:
        splitted_name = obj.__name__
    path = obj.__module__.split('.')+splitted_name[:-1]
    name = splitted_name[-1]
    try:
        M = __import__(path.pop(0))
    except ImportError:
        try:
            B = obj.__base__
            if B is None:
                raise AttributeError
        except AttributeError:
            raise OSError("could not get source code")
        return sage_getsourcelines(B)
    # M should just be the top-most module.
    # Hence, normally it is just 'sage'
    try:
        while path:
            M = getattr(M, path.pop(0))
    except AttributeError:
        try:
            B = obj.__base__
            if B is None:
                raise AttributeError
        except AttributeError:
            raise OSError("could not get source code")
        return sage_getsourcelines(B)

    lines, base_lineno = sage_getsourcelines(M)
    # the rest of the function is copied from
    # inspect.findsource
    if not lines:
        raise OSError('could not get source code')

    if inspect.ismodule(obj):
        return lines, base_lineno

    if inspect.isclass(obj):
        pat = re.compile(r'^(\s*)class\s*' + name + r'\b')
        # make some effort to find the best matching class definition:
        # use the one with the least indentation, which is the one
        # that's most probably not inside a function definition.
        candidates = []
        for i in range(len(lines)):
            match = pat.match(lines[i])
            if match:
                # if it's at toplevel, it's already the best one
                if lines[i][0] == 'c':
                    return inspect.getblock(lines[i:]), i+base_lineno
                # else add whitespace to candidate list
                candidates.append((match.group(1), i))
        if candidates:
            # this will sort by whitespace, and by line number,
            # less whitespace first
            candidates.sort()
            return inspect.getblock(lines[candidates[0][1]:]), candidates[0][1]+base_lineno
        else:
            raise OSError('could not find class definition')

    if inspect.ismethod(obj):
        obj = obj.__func__
    if is_function_or_cython_function(obj):
        obj = obj.__code__
    if inspect.istraceback(obj):
        obj = obj.tb_frame
    if inspect.isframe(obj):
        obj = obj.f_code
    if inspect.iscode(obj):
        if not hasattr(obj, 'co_firstlineno'):
            raise OSError('could not find function definition')
        pat = re.compile(r'^(\s*def\s)|(.*(?<!\w)lambda(:|\s))|^(\s*@)')
        pmatch = pat.match
        # fperez - fix: sometimes, co_firstlineno can give a number larger than
        # the length of lines, which causes an error.  Safeguard against that.
        lnum = min(obj.co_firstlineno, len(lines)) - 1
        while lnum > 0:
            if pmatch(lines[lnum]):
                break
            lnum -= 1

        return inspect.getblock(lines[lnum:]), lnum+base_lineno
    raise OSError('could not find code object')


def sage_getsourcelines(obj):
    r"""
    Return a pair ([source_lines], starting line number) of the source
    code associated to obj, or None.

    INPUT:

    - ``obj`` -- function, etc.

    OUTPUT:

    (source_lines, lineno) or None: ``source_lines`` is a list of
    strings, and ``lineno`` is an integer.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getsourcelines

        sage: # needs sage.modules
        sage: sage_getsourcelines(matrix)[1]
        21
        sage: sage_getsourcelines(matrix)[0][0]
        'def matrix(*args, **kwds):\n'

    Some classes customize this using a ``_sage_src_lines_`` method,
    which gives the source lines of a class instance, but not the class
    itself. We demonstrate this for :class:`CachedFunction`::

        sage: # needs sage.combinat
        sage: cachedfib = cached_function(fibonacci)
        sage: sage_getsourcelines(cachedfib)[0][0]
        "def fibonacci(n, algorithm='pari') -> Integer:\n"
        sage: sage_getsourcelines(type(cachedfib))[0][0]
        'cdef class CachedFunction():\n'

    TESTS::

        sage: # needs sage.misc.cython
        sage: cython('''cpdef test_funct(x, y): return''')
        sage: sage_getsourcelines(test_funct)
        (['cpdef test_funct(x, y): return\n'], 1)

    The following tests that an instance of ``functools.partial`` is correctly
    dealt with (see :issue:`9976`)::

        sage: from sage.tests.functools_partial_src import test_func
        sage: sage_getsourcelines(test_func)
        (['def base(x):\n',
        ...
        '    return x\n'], 8)

    Here are some cases that were covered in :issue:`11298`;
    note that line numbers may easily change, and therefore we do
    not test them::

        sage: P.<x,y> = QQ[]
        sage: I = P*[x,y]
        sage: sage_getsourcelines(P)                                                    # needs sage.libs.singular
        (['cdef class MPolynomialRing_libsingular(MPolynomialRing_base):\n',
          '\n',
          '    def __cinit__(self):\n',
        ...)
        sage: sage_getsourcelines(I)                                                    # needs sage.libs.singular
        ([...'class MPolynomialIdeal(MPolynomialIdeal_singular_repr,\n',
        ...)
        sage: x = var('x')                                                              # needs sage.symbolic
        sage: lines, lineno = sage_getsourcelines(x); lines[0:5]                        # needs sage.symbolic
        ['cdef class Expression(...):\n',
         '\n',
         '    cdef GEx _gobj\n',
         '\n',
         '    cpdef object pyobject(self):\n']
        sage: lines[-1]    # last line                                                  # needs sage.symbolic
        '        return S\n'

    We show some enhancements provided by :issue:`11768`. First, we
    use a dummy parent class that has defined an element class by a
    nested class definition::

        sage: from sage.misc.test_nested_class import TestNestedParent
        sage: from sage.misc.sageinspect import sage_getsource
        sage: P = TestNestedParent()
        sage: E = P.element_class
        sage: E.__bases__
        (<class 'sage.misc.test_nested_class.TestNestedParent.Element'>,
         <class 'sage.categories.sets_cat.Sets.element_class'>)
        sage: print(sage_getsource(E))
            class Element:
                "This is a dummy element class"
                pass
        sage: print(sage_getsource(P))
        class TestNestedParent(UniqueRepresentation, Parent):
            ...
            class Element:
                "This is a dummy element class"
                pass

    Here is another example that relies on a nested class definition
    in the background::

        sage: C = AdditiveMagmas()
        sage: HC = C.Homsets()
        sage: sage_getsourcelines(HC)
        (['    class Homsets(HomsetsCategory):\n', ...], ...)

    Testing against a bug that has occurred during work on :issue:`11768`::

        sage: P.<x,y> = QQ[]
        sage: I = P*[x,y]
        sage: sage_getsourcelines(I)
        ([...'class MPolynomialIdeal(MPolynomialIdeal_singular_repr,\n',
          '                       MPolynomialIdeal_macaulay2_repr,\n',
          '                       MPolynomialIdeal_magma_repr,\n',
          '                       Ideal_generic):\n',
          '    def __init__(self, ring, gens, coerce=True):\n',
          ...)
    """
    # First try the method _sage_src_lines_(), which is meant to give
    # the source lines of an object (not of its type!).
    try:
        sage_src_lines = obj._sage_src_lines_
    except AttributeError:
        pass
    else:
        try:
            return sage_src_lines()
        except (NotImplementedError, TypeError):
            # NotImplementedError can be raised by _sage_src_lines_()
            # to indicate that it didn't find the source lines.
            #
            # TypeError can happen when obj is a type and
            # obj._sage_src_lines_ is an unbound method. In this case,
            # we don't want to use _sage_src_lines_(), we just want to
            # get the source of the type itself.
            pass

    # Check if we deal with an instance
    if isclassinstance(obj):
        if isinstance(obj, functools.partial):
            return sage_getsourcelines(obj.func)
        else:
            return sage_getsourcelines(obj.__class__)

    # First, we deal with nested classes. Their name contains a dot, and we
    # have a special function for that purpose.
    # This is the case for ParentMethods of categories, for example.
    if (inspect.isclass(obj) and
            ('.' in obj.__name__ or '.' in getattr(obj, '__qualname__', ''))):
        return _sage_getsourcelines_name_with_dot(obj)

    # Next, we try _sage_getdoc_unformatted()
    d = _sage_getdoc_unformatted(obj)
    pos = _extract_embedded_position(d)
    if pos is None:
        try:
            return inspect.getsourcelines(obj)
        except (OSError, TypeError) as err:
            try:
                objinit = obj.__init__
            except AttributeError:
                pass
            else:
                d = _sage_getdoc_unformatted(objinit)
                pos = _extract_embedded_position(d)
                if pos is None:
                    if inspect.isclass(obj):
                        try:
                            B = obj.__base__
                        except AttributeError:
                            B = None
                        if B is not None and B is not obj:
                            return sage_getsourcelines(B)
                    if obj.__class__ != type:
                        return sage_getsourcelines(obj.__class__)
                    raise err

    (orig, filename, lineno) = pos
    try:
        with open(filename) as f:
            source_lines = f.readlines()
    except OSError:
        try:
            from sage.misc.temporary_file import spyx_tmp
            raw_name = filename.split('/')[-1]
            newname = os.path.join(spyx_tmp(), '_'.join(raw_name.split('_')[:-1]), raw_name)
            with open(newname) as f:
                source_lines = f.readlines()
        except OSError:
            return None

    # It is possible that the source lines belong to the __init__ method,
    # rather than to the class. So, we try to look back and find the class
    # definition.
    first_line = source_lines[lineno-1]
    leading_blanks = len(first_line)-len(first_line.lstrip())
    if first_line.lstrip().startswith('def ') and "__init__" in first_line and obj.__name__ != '__init__':
        ignore = False
        double_quote = None
        for lnb in range(lineno, 0, -1):
            new_first_line = source_lines[lnb-1]
            nfl_strip = new_first_line.lstrip()
            if nfl_strip.startswith('"""'):
                if double_quote is None:
                    double_quote = True
                if double_quote:
                    ignore = not ignore
            elif nfl_strip.startswith("'''"):
                if double_quote is None:
                    double_quote = False
                if double_quote is False:
                    ignore = not ignore
            if ignore:
                continue
            if len(new_first_line)-len(nfl_strip) < leading_blanks and nfl_strip:
                # We are not inside a doc string. So, if the indentation
                # is less than the indentation of the __init__ method
                # then we must be at the class definition!
                lineno = lnb
                break
    return _extract_source(source_lines, lineno), lineno


def sage_getvariablename(self, omit_underscore_names=True):
    """
    Attempt to get the name of a Sage object.

    INPUT:

    - ``self`` -- any object

    - ``omit_underscore_names`` -- boolean (default: ``True``)

    OUTPUT:

    If the user has assigned an object ``obj`` to a variable name,
    then return that variable name.  If several variables point to
    ``obj``, return a sorted list of those names.  If
    ``omit_underscore_names`` is ``True`` (the default) then omit names
    starting with an underscore "_".

    EXAMPLES::

        sage: # needs sage.modules
        sage: from sage.misc.sageinspect import sage_getvariablename
        sage: A = random_matrix(ZZ, 100)
        sage: sage_getvariablename(A)
        'A'
        sage: B = A
        sage: sage_getvariablename(A)
        ['A', 'B']

    If an object is not assigned to a variable, an empty list is returned::

        sage: sage_getvariablename(random_matrix(ZZ, 60))                               # needs sage.modules
        []
    """
    # This is a modified version of code taken from
    # https://web.archive.org/web/20100416095847/http://pythonic.pocoo.org/2009/5/30/finding-objects-names
    # written by Georg Brandl.
    result = []
    for frame in inspect.stack():
        for name, obj in frame[0].f_globals.items():
            if obj is self:
                result.append(name)
    if len(result) == 1:
        return result[0]
    else:
        return sorted(result)


__internal_teststring = '''
import os                                  # 1
# preceding comment not include            # 2
def test1(a, b=2):                         # 3
    if a:                                  # 4
        return 1                           # 5
    return b                               # 6
# intervening comment not included         # 7
class test2():                             # 8
    pass                                   # 9
    # indented comment not included        # 10
# trailing comment not included            # 11
def test3(b,                               # 12
          a=2):                            # 13
    pass # EOF                             # 14'''


def __internal_tests():
    r"""
    Test internals of the sageinspect module.

    EXAMPLES::

        sage: from sage.misc.sageinspect import *
        sage: from sage.misc.sageinspect import _extract_source, _extract_embedded_position, _sage_getargspec_cython, __internal_teststring

    If docstring is None, nothing bad happens::

        sage: sage_getdoc(None)
        ''

        sage: import sage.all__sagemath_objects
        sage: sage_getsource(sage.all__sagemath_objects)
        '...all...'

    A cython function with default arguments (one of which is a string)::

        sage: sage_getdef(sage.rings.integer.Integer.factor, obj_name='factor')
        "factor(algorithm='pari', proof=None, limit=None, int_=False, verbose=0)"

    This used to be problematic, but was fixed in :issue:`10094`::

        sage: sage_getsource(sage.rings.integer.Integer.__init__)
        '    def __init__(self, x=None, base=0):\n...'
        sage: sage_getdef(sage.rings.integer.Integer.__init__, obj_name='__init__')
        '__init__(x=None, base=0)'

    Test _extract_source with some likely configurations, including no trailing
    newline at the end of the file::

        sage: s = __internal_teststring.strip()
        sage: es = lambda ls, l: ''.join(_extract_source(ls, l)).rstrip()

        sage: print(es(s, 3))
        def test1(a, b=2):                         # 3
            if a:                                  # 4
                return 1                           # 5
            return b                               # 6

        sage: print(es(s, 8))
        class test2():                             # 8
            pass                                   # 9

        sage: print(es(s, 12))
        def test3(b,                               # 12
                  a=2):                            # 13
            pass # EOF                             # 14

    Test _sage_getargspec_cython with multiple default arguments and a type::

        sage: _sage_getargspec_cython("def init(self, x=None, base=0):")
        FullArgSpec(args=['self', 'x', 'base'], varargs=None, varkw=None, defaults=(None, 0), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _sage_getargspec_cython("def __init__(self, x=None, base=0):")
        FullArgSpec(args=['self', 'x', 'base'], varargs=None, varkw=None, defaults=(None, 0), kwonlyargs=[], kwonlydefaults=None, annotations={})
        sage: _sage_getargspec_cython("def __init__(self, x=None, unsigned int base=0, **keys):")
        FullArgSpec(args=['self', 'x', 'base'], varargs=None, varkw='keys', defaults=(None, 0), kwonlyargs=[], kwonlydefaults=None, annotations={})

    Test _extract_embedded_position:

    We cannot test the filename since it depends on the installation location.

    Make sure things work with no trailing newline::

        sage: _extract_embedded_position('File: sage/rings/rational.pyx (starting at line 1080)')
        ('', '.../rational.pyx', 1080)

    And with a trailing newline::

        sage: s = 'File: sage/rings/rational.pyx (starting at line 1080)\n'
        sage: _extract_embedded_position(s)
        ('', '.../rational.pyx', 1080)

    And with an original docstring::

        sage: s = 'File: sage/rings/rational.pyx (starting at line 1080)\noriginal'
        sage: _extract_embedded_position(s)
        ('original', '.../rational.pyx', 1080)

    And with a complicated original docstring::

        sage: s = 'File: sage/rings/rational.pyx (starting at line 1080)\n\n\noriginal test\noriginal'
        sage: _extract_embedded_position(s)
        ('\n\noriginal test\noriginal', ..., 1080)

        sage: s = 'no embedded position'
        sage: _extract_embedded_position(s) is None
        True
    """
