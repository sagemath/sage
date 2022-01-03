r"""
Some tools for developers

AUTHORS:

- Nicolas M. Thiery: initial version

- Vincent Delecroix (2012 and 2013): improve import_statements
"""
# ****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import os
import re
import sys

from collections import defaultdict


def runsnake(command):
    """
    Graphical profiling with ``runsnake``

    INPUT:

    - ``command`` -- the command to be run as a string.

    EXAMPLES::

        sage: runsnake("list(SymmetricGroup(3))")        # optional - runsnake

    ``command`` is first preparsed (see :func:`preparse`)::

        sage: runsnake('for x in range(1,4): print(x^2)') # optional - runsnake
        1
        4
        9

    :func:`runsnake` requires the program ``runsnake``. Due to non
    trivial dependencies (python-wxgtk, ...), installing it within the
    Sage distribution is unpractical. Hence, we recommend installing
    it with the system wide Python. On Ubuntu 10.10, this can be done
    with::

        > sudo apt-get install python-profiler python-wxgtk2.8 python-setuptools
        > sudo easy_install RunSnakeRun

    See the ``runsnake`` website for instructions for other platforms.

    :func:`runsnake` further assumes that the system wide Python is
    installed in ``/usr/bin/python``.

    .. SEEALSO::

        - `The runsnake website <http://www.vrplumber.com/programming/runsnakerun/>`_
        - ``%prun``
        - :class:`Profiler`

    """
    import cProfile
    from sage.misc.temporary_file import tmp_filename
    from sage.misc.misc import get_main_globals
    from sage.repl.preparse import preparse
    tmpfile = tmp_filename()
    cProfile.runctx(preparse(command.lstrip().rstrip()), get_main_globals(),
                    locals(), filename=tmpfile)
    os.system("/usr/bin/python -E `which runsnake` %s &" % tmpfile)


def import_statement_string(module, names, lazy):
    r"""
    Return a (lazy) import statement for ``names`` from ``module``.

    INPUT:

    - ``module`` -- the name of a module

    - ``names`` -- a list of 2-tuples containing names and alias to
      import

    - ``lazy`` -- a boolean: whether to return a lazy import statement

    EXAMPLES::

        sage: import sage.misc.dev_tools as dt
        sage: modname = 'sage.misc.dev_tools'
        sage: names_and_aliases = [('import_statement_string', 'iss')]
        sage: dt.import_statement_string(modname, names_and_aliases, False)
        'from sage.misc.dev_tools import import_statement_string as iss'
        sage: dt.import_statement_string(modname, names_and_aliases, True)
        "lazy_import('sage.misc.dev_tools', 'import_statement_string', 'iss')"
        sage: dt.import_statement_string(modname, [('a','b'),('c','c'),('d','e')], False)
        'from sage.misc.dev_tools import a as b, c, d as e'
        sage: dt.import_statement_string(modname, [(None,None)], False)
        'import sage.misc.dev_tools'
    """
    if lazy:
        if len(names) == 1:
            name, alias = names[0]
            if name == alias:
                if name is None:
                    raise ValueError("cannot lazy import modules")
                return "lazy_import('%s', '%s')" % (module, name)
            else:
                return "lazy_import('%s', '%s', '%s')" % (module, name, alias)
        obj_names = "[" + ", ".join("'" + name[0] + "'" for name in names) + "]"
        obj_aliases = "[" + ", ".join("'" + name[1] + "'" for name in names) + "]"
        return "lazy_import('%s', %s, %s)" % (module, obj_names, obj_aliases)
    else:
        import_module = False
        name_list = []
        for name, alias in names:
            if name == alias:
                if name is None:
                    import_module = True
                    continue
                name_list.append(name)
            else:
                name_list.append("%s as %s" % (name, alias))
        res = []
        if import_module:
            res.append("import %s" % module)
        if name_list:
            res.append("from %s import %s" % (module, ', '.join(name_list)))
        return "\n".join(res)


def load_submodules(module=None, exclude_pattern=None):
    r"""
    Load all submodules of a given modules.

    This method is intended to be used by developers and especially the one
    who uses :func:`import_statements`. By default it load the sage library and
    it takes around a minute.

    INPUT:

    - ``module`` - an optional module

    - ``exclude_pattern`` - an optional regular expression pattern of module
      names that have to be excluded.

    EXAMPLES::

        sage: sage.misc.dev_tools.load_submodules(sage.combinat)
        load sage.combinat.algebraic_combinatorics... succeeded
        ...
        load sage.combinat.words.suffix_trees... succeeded

    Calling a second time has no effect (since the function does not import
    modules already imported)::

        sage: sage.misc.dev_tools.load_submodules(sage.combinat)

    The second argument allows to exclude a pattern::

        sage: sage.misc.dev_tools.load_submodules(sage.geometry, "database$|lattice")
        load sage.geometry.cone_catalog... succeeded
        load sage.geometry.fan_isomorphism... succeeded
        ...
        load sage.geometry.riemannian_manifolds.surface3d_generators... succeeded

        sage: sage.misc.dev_tools.load_submodules(sage.geometry)
        load sage.geometry.polyhedron.lattice_euclidean_group_element... succeeded
        load sage.geometry.polyhedron.palp_database... succeeded
        load sage.geometry.polyhedron.ppl_lattice_polygon... succeeded
    """
    import pkgutil

    if module is None:
        import sage
        module = sage
        exclude_pattern = r"^sage\.libs|^sage\.tests|tests$|^sage\.all_|all$|sage\.interacts$|^sage\.misc\.benchmark$"

    if exclude_pattern:
        exclude = re.compile(exclude_pattern)
    else:
        exclude = None

    for importer, module_name, ispkg in pkgutil.walk_packages(module.__path__, module.__name__ + '.'):
        if ispkg or module_name in sys.modules:
            continue

        # we exclude several sage components because loading them is much of a
        # problem...
        if exclude and exclude.search(module_name):
            continue

        try:
            sys.stdout.write("load %s..." % module_name)
            sys.stdout.flush()
            loader = importer.find_module(module_name)
            loader.load_module(module_name)
            sys.stdout.write(" succeeded\n")
        except (ValueError, AttributeError, TypeError, ImportError):
            # we might get error because of cython code that has been
            # compiled but with source removed
            sys.stdout.write("failed\n")


def find_objects_from_name(name, module_name=None):
    r"""
    Return the list of objects from ``module_name`` whose name is ``name``.

    If ``module_name`` is ``None``, the function runs through all
    loaded modules and returns the list of objects whose name matches ``name``.

    If ``module_name`` is not ``None``, then search only in submodules of
    ``module_name``.

    In order to search through more modules you might use the function
    :func:`load_submodules`.

    EXAMPLES::

        sage: import sage.misc.dev_tools as dt
        sage: dt.find_objects_from_name('FareySymbol')
        [<class 'sage.modular.arithgroup.farey_symbol.Farey'>]

        sage: import sympy
        sage: dt.find_objects_from_name('RR')
        [Real Field with 53 bits of precision, RR]
        sage: dt.find_objects_from_name('RR', 'sage')
        [Real Field with 53 bits of precision]
        sage: dt.find_objects_from_name('RR', 'sympy')
        [RR]

    Examples that do not belong to the global namespace but in a loaded module::

        sage: 'find_objects_from_name' in globals()
        False
        sage: objs = dt.find_objects_from_name('find_objects_from_name')
        sage: len(objs)
        1
        sage: dt.find_objects_from_name is dt.find_objects_from_name
        True

    .. NOTE::

        It might be a good idea to move this function into
        :mod:`sage.misc.sageinspect`.
    """

    obj = []
    for smodule_name, smodule in sys.modules.items():
        if module_name and not smodule_name.startswith(module_name):
            continue
        if hasattr(smodule, '__dict__') and name in smodule.__dict__:
            u = smodule.__dict__[name]
            if all(v is not u for v in obj):
                obj.append(u)

    return obj


def find_object_modules(obj):
    r"""
    Return a dictionary whose keys are the names of the modules where ``obj``
    appear and the value at a given module name is the list of names that
    ``obj`` have in that module.

    It is very unlikely that the output dictionary has several keys except when
    ``obj`` is an instance of a class.

    EXAMPLES::

        sage: from sage.misc.dev_tools import find_object_modules
        sage: find_object_modules(RR)
        {'sage.rings.real_mpfr': ['RR']}
        sage: find_object_modules(ZZ)
        {'sage.rings.integer_ring': ['Z', 'ZZ']}

    .. NOTE::

        It might be a good idea to move this function in
        :mod:`sage.misc.sageinspect`.
    """
    from sage.misc import sageinspect

    # see if the object is defined in its own module
    # might be wrong for class instances as the instanciation might appear
    # outside of the module !!
    module_name = None
    if sageinspect.isclassinstance(obj):
        module_name = obj.__class__.__module__
    elif hasattr(obj, '__module__') and obj.__module__:
        module_name = obj.__module__

    if module_name:
        if module_name not in sys.modules:
            raise ValueError("This should not happen!")
        d = sys.modules[module_name].__dict__
        matching = sorted(key for key in d if d[key] is obj)
        if matching:
            return {module_name: matching}

    # otherwise, we parse all (already loaded) modules and hope to find
    # something
    module_to_obj = {}
    for module_name, module in sys.modules.items():
        if module_name != '__main__' and hasattr(module, '__dict__'):
            d = module.__dict__
            names = [key for key in d if d[key] is obj]
            if names:
                module_to_obj[module_name] = names

    # if the object is an instance, we try to guess where it is defined
    if sageinspect.isclassinstance(obj):
        dec_pattern = re.compile(r"^(\w[\w0-9\_]*)\s*=", re.MULTILINE)
        module_to_obj2 = {}
        for module_name, obj_names in module_to_obj.items():
            module_to_obj2[module_name] = []
            try:
                src = sageinspect.sage_getsource(sys.modules[module_name])
            except TypeError:
                pass
            else:
                m = dec_pattern.search(src)
                while m:
                    if m.group(1) in obj_names:
                        module_to_obj2[module_name].append(m.group(1))
                    m = dec_pattern.search(src, m.end())
            if not module_to_obj2[module_name]:
                del module_to_obj2[module_name]

        if module_to_obj2:
            return module_to_obj2

    return module_to_obj


def import_statements(*objects, **kwds):
    r"""
    Print import statements for the given objects.

    INPUT:

    - ``*objects`` -- a sequence of objects or names.

    - ``lazy`` -- a boolean (default: ``False``)
      Whether to print a lazy import statement.

    - ``verbose`` -- a boolean (default: ``True``)
      Whether to print information in case of ambiguity.

    - ``answer_as_str`` -- a boolean (default: ``False``)
      If ``True`` return a string instead of printing the statement.

    EXAMPLES::

        sage: import_statements(WeylGroup, lazy_attribute)
        from sage.combinat.root_system.weyl_group import WeylGroup
        from sage.misc.lazy_attribute import lazy_attribute

        sage: import_statements(IntegerRing)
        from sage.rings.integer_ring import IntegerRing

    If ``lazy`` is True, then :func:`lazy_import` statements are
    displayed instead::

        sage: import_statements(WeylGroup, lazy_attribute, lazy=True)
        from sage.misc.lazy_import import lazy_import
        lazy_import('sage.combinat.root_system.weyl_group', 'WeylGroup')
        lazy_import('sage.misc.lazy_attribute', 'lazy_attribute')

    In principle, the function should also work on object which are instances.
    In case of ambiguity, one or two warning lines are printed::

        sage: import_statements(RDF)
        from sage.rings.real_double import RDF

        sage: import_statements(ZZ)
        # ** Warning **: several names for that object: Z, ZZ
        from sage.rings.integer_ring import Z

        sage: import_statements(euler_phi)
        from sage.arith.misc import euler_phi

        sage: import_statements(x)
        from sage.calculus.predefined import x

    If you don't like the warning you can disable them with the option ``verbose``::

        sage: import_statements(ZZ, verbose=False)
        from sage.rings.integer_ring import Z

        sage: import_statements(x, verbose=False)
        from sage.calculus.predefined import x

    If the object has several names, an other way to get the import
    statement you expect is to use a string instead of the object::

        sage: import_statements(matrix)
        # ** Warning **: several names for that object: Matrix, matrix
        from sage.matrix.constructor import Matrix

        sage: import_statements('cached_function')
        from sage.misc.cachefunc import cached_function
        sage: import_statements('Z')
        # **Warning**: distinct objects with name 'Z' in:
        #   - sage.calculus.predefined
        #   - sage.rings.integer_ring
        from sage.rings.integer_ring import Z

    Specifying a string is also useful for objects that are not
    imported in the Sage interpreter namespace by default. In this
    case, an object with that name is looked up in all the modules
    that have been imported in this session::

        sage: import_statement_string
        Traceback (most recent call last):
        ...
        NameError: name 'import_statement_string' is not defined

        sage: import_statements("import_statement_string")
        from sage.misc.dev_tools import import_statement_string

    Sometimes objects are imported as an alias (from XXX import YYY as ZZZ) or
    are affected (XXX = YYY) and the function might detect it::

        sage: import_statements('FareySymbol')
        from sage.modular.arithgroup.farey_symbol import Farey as FareySymbol

        sage: import_statements('power')
        from sage.arith.power import generic_power as power

    In order to be able to detect functions that belong to a non-loaded module,
    you might call the helper :func:`load_submodules` as in the following::

        sage: import_statements('HeckeMonoid')
        Traceback (most recent call last):
        ...
        LookupError: no object named 'HeckeMonoid'
        sage: from sage.misc.dev_tools import load_submodules
        sage: load_submodules(sage.monoids)
        load sage.monoids.automatic_semigroup... succeeded
        load sage.monoids.hecke_monoid... succeeded
        load sage.monoids.indexed_free_monoid... succeeded
        sage: import_statements('HeckeMonoid')
        from sage.monoids.hecke_monoid import HeckeMonoid

    We test different objects which have no appropriate answer::

        sage: import_statements('my_tailor_is_rich')
        Traceback (most recent call last):
        ...
        LookupError: no object named 'my_tailor_is_rich'
        sage: import_statements(5)
        Traceback (most recent call last):
        ...
        ValueError: no import statement found for '5'.

    We test that it behaves well with lazy imported objects (:trac:`14767`)::

        sage: import_statements(NN)
        from sage.rings.semirings.non_negative_integer_semiring import NN
        sage: import_statements('NN')
        from sage.rings.semirings.non_negative_integer_semiring import NN

    Deprecated lazy imports are ignored (see :trac:`17458`)::

        sage: lazy_import('sage.all', 'RR', 'deprecated_RR', namespace=sage.__dict__, deprecation=17458)
        sage: import_statements('deprecated_RR')
        Traceback (most recent call last):
        ...
        LookupError: object named 'deprecated_RR' is deprecated (see trac ticket 17458)
        sage: lazy_import('sage.all', 'RR', namespace=sage.__dict__, deprecation=17458)
        sage: import_statements('RR')
        from sage.rings.real_mpfr import RR

    The following were fixed with :trac:`15351`::

        sage: import_statements('Rationals')
        from sage.rings.rational_field import RationalField as Rationals
        sage: import_statements(sage.combinat.partition_algebra.SetPartitionsAk)
        from sage.combinat.partition_algebra import SetPartitionsAk
        sage: import_statements(CIF)
        from sage.rings.cif import CIF
        sage: import_statements(NaN)
        from sage.symbolic.constants import NaN
        sage: import_statements(pi)
        from sage.symbolic.constants import pi
        sage: import_statements('SAGE_ENV')
        from sage.env import SAGE_ENV
        sage: import_statements('graph_decompositions')
        import sage.graphs.graph_decompositions

    Check that a name from the global namespace is properly found (see
    :trac:`23779`)::

        sage: import_statements('log')
        from sage.misc.functional import log

    .. NOTE::

        The programmers try to made this function as smart as possible.
        Nevertheless it is far from being perfect (for example it does not
        detect deprecated stuff). So, if you use it, double check the answer and
        report weird behaviors.
    """
    import inspect
    from sage.misc.lazy_import import LazyImport

    answer = defaultdict(list)
    module_name = None
    # a dictionary module -> [(name1,alias1), (name2,alias2) ...]
    # where "nameX" is an object in "module" that has to be
    # imported with the alias "aliasX"

    lazy = kwds.pop("lazy", False)
    verbose = kwds.pop("verbose", True)
    answer_as_str = kwds.pop("answer_as_str", False)

    if kwds:
        raise TypeError("Unexpected '{}' argument".format(next(iter(kwds))))

    for obj in objects:
        name = None    # the name of the object

        # 1. if obj is a string, we look for an object that has that name
        if isinstance(obj, str):
            from sage.all import sage_globals
            G = sage_globals()
            name = obj
            if name in G:
                # 1.a. object in the sage namespace
                obj = [G[name]]
            else:
                # 1.b. object inside a submodule of sage
                obj = find_objects_from_name(name, 'sage')
                if not obj:
                    # 1.c. object from something already imported
                    obj = find_objects_from_name(name)

            # remove lazy imported objects from list obj
            i = 0
            deprecation = None
            while i < len(obj):
                if isinstance(obj[i], LazyImport):
                    tmp = obj.pop(i)
                    # Ignore deprecated lazy imports
                    tmp_deprecation = tmp._get_deprecation_ticket()
                    if tmp_deprecation:
                        deprecation = tmp_deprecation
                    else:
                        tmp = tmp._get_object()
                        if all(u is not tmp for u in obj):
                            obj.append(tmp)
                else:
                    i += 1

            if verbose and len(obj) > 1:
                modules = set()
                for o in obj:
                    modules.update(find_object_modules(o))
                print("# **Warning**: distinct objects with name '{}' "
                      "in:".format(name))
                for mod in sorted(modules):
                    print("#   - {}".format(mod))

            # choose a random object among the potentially enormous list of
            # objects we get from "name"
            try:
                obj = obj[0]
            except IndexError:
                if deprecation:
                    raise LookupError(
                        "object named {!r} is deprecated (see trac ticket "
                        "{})".format(name, deprecation))
                else:
                    raise LookupError("no object named {!r}".format(name))

        # 1'. if obj is a LazyImport we recover the real object
        if isinstance(obj, LazyImport):
            obj = obj._get_object()

        # 2. Find out in which modules obj lives
        # and update answer with a couple of strings "(name,alias)" where "name" is
        # the name of the object in the module and "alias" is the name of the
        # object

        # easy case: the object is itself a module
        if inspect.ismodule(obj):
            module_name = obj.__name__
            answer[module_name].append((None, None))
            continue

        modules = find_object_modules(obj)
        if '__main__' in modules:
            del modules['__main__']
        if '__mp_main__' in modules:
            del modules['__mp_main__']

        if not modules:
            raise ValueError("no import statement found for '{}'.".format(obj))

        if name is None:
            # if the object is available under both ascii and unicode names,
            # prefer the ascii version.
            def is_ascii(s):
                """
                Equivalent of `str.isascii` in Python >= 3.7
                """
                return all(ord(c) < 128 for c in s)
            if any(is_ascii(s)
                   for (module_name, obj_names) in modules.items()
                   for s in obj_names):
                for module_name, obj_names in list(modules.items()):
                    if any(not is_ascii(s) for s in obj_names):
                        obj_names = [name for name in obj_names if is_ascii(name)]
                        if not obj_names:
                            del modules[module_name]
                        else:
                            modules[module_name] = obj_names

        if len(modules) == 1:  # the module is well defined
            (module_name, obj_names), = modules.items()
            if name is None:
                if verbose and len(obj_names) > 1:
                    print("# ** Warning **: several names for that object: "
                          "{}".format(', '.join(sorted(obj_names))))
                name = alias = obj_names[0]
            elif name in modules[module_name]:
                alias = name
            else:
                alias = name
                name = obj_names[0]

            answer[module_name].append((name, alias))
            continue

        # here modules contain several answers and we first try to see if there
        # is a best one (i.e. the object "obj" is contained in the module and
        # has name "name")
        if name is not None:
            good_modules = []
            for mod in modules:
                if name in modules[mod]:
                    good_modules.append(mod)

            if len(good_modules) == 1:
                answer[good_modules[0]].append((name, name))
                continue

        # if the object is a class instance, it is likely that it is defined in
        # some XYZ.all module
        from .sageinspect import isclassinstance
        if isclassinstance(obj):
            module_name = type(obj).__module__
            i = module_name.rfind('.')
            all_module_name = module_name[:i] + '.all'
            if all_module_name in modules:
                module_name = all_module_name
                modules[module_name][0]
            else:
                module_name = None

        if module_name is None:
            # here, either "obj" is a class instance but there is no natural
            # candidate for its module or "obj" is not a class instance.
            all_re = re.compile(r'.+\.all(?:_\w+)?$')
            not_all_modules = [mod for mod in modules
                               if not all_re.match(mod)]
            if not not_all_modules:
                print("# ** Warning **: the object {} is only defined in "
                      ".all modules".format(obj))
                module_name = next(iter(modules))
            else:
                if len(not_all_modules) > 1:
                    print("# ** Warning **: several modules for the object "
                          "{}: {}".format(obj, ', '.join(sorted(modules))))
                module_name = not_all_modules[0]

        # 3. Now that we found the module, we fix the problem of the alias
        if name is None:
            alias = name = modules[module_name][0]
        else:
            alias = name
            name = modules[module_name][0]

        answer[module_name].append((name, alias))

    res = []

    if lazy:
        res.append("from sage.misc.lazy_import import lazy_import")

    for module_name in sorted(answer):
        res.append(import_statement_string(module_name, answer[module_name],
                                           lazy))

    if answer_as_str:
        return '\n'.join(res)
    else:
        print('\n'.join(res))
