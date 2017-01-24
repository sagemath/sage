"""
Read and parse the file pari.desc

Run tests from the ``SAGE_SRC`` directory::

    sage: from sage.env import SAGE_SRC
    sage: os.chdir(SAGE_SRC)
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os, re

from sage_setup.autogen.pari.args import pari_arg_types
from sage_setup.autogen.pari.ret import pari_ret_types


def sage_src_pari():
    """
    Return the directory in the Sage source tree where the interface to
    the PARI library is and where the auto-generated files should be
    stored.

    EXAMPLES::

        sage: from sage_setup.autogen.pari.parser import sage_src_pari
        sage: sage_src_pari()
        'sage/libs/cypari2'
    """
    return os.path.join('sage', 'libs', 'cypari2')


def pari_share():
    r"""
    Return the directory where the PARI data files are stored.

    EXAMPLES::

        sage: from sage_setup.autogen.pari.parser import pari_share
        sage: pari_share()
        '.../share/pari'
    """
    from subprocess import Popen, PIPE
    gp = Popen(["gp", "-f", "-q"], stdin=PIPE, stdout=PIPE)
    out = gp.communicate(b"print(default(datadir))")[0]
    datadir = out.strip()
    if not os.path.isdir(datadir):
        raise EnvironmentError("PARI data directory {!r} does not exist".format(datadir))
    return datadir


paren_re = re.compile(r"[(](.*)[)]")
argname_re = re.compile(r"[ {]*([A-Za-z_][A-Za-z0-9_]*)")

def read_pari_desc():
    """
    Read and parse the file ``pari.desc``.

    The output is a dictionary where the keys are GP function names
    and the corresponding values are dictionaries containing the
    ``(key, value)`` pairs from ``pari.desc``.

    EXAMPLES::

        sage: from sage_setup.autogen.pari.parser import read_pari_desc
        sage: D = read_pari_desc()
        sage: D["cos"]
        {'class': 'basic',
         'cname': 'gcos',
         'doc': 'cosine of $x$.',
         'function': 'cos',
         'help': 'cos(x): cosine of x.',
         'prototype': 'Gp',
         'section': 'transcendental'}
    """
    with open(os.path.join(pari_share(), b'pari.desc')) as f:
        lines = f.readlines()

    n = 0
    N = len(lines)

    functions = {}
    while n < N:
        fun = {}
        while True:
            L = lines[n]; n += 1
            if L == "\n":
                break
            # As long as the next lines start with a space, append them
            while lines[n].startswith(" "):
                L += (lines[n])[1:]; n += 1
            key, value = L.split(":", 1)
            # Change key to an allowed identifier name
            key = key.lower().replace("-", "")
            fun[key] = value.strip()

        name = fun["function"]
        functions[name] = fun

    return functions


decl_re = re.compile(" ([A-Za-z][A-Za-z0-9_]*)[(]")

def read_decl():
    """
    Read the files ``paridecl.pxd`` and ``declinl.pxi`` and return a set
    of all declared PARI library functions.

    We do a simple regexp search, so there might be false positives.
    The main use is to skip undeclared functions.

    EXAMPLES::

        sage: from sage_setup.autogen.pari.parser import read_decl
        sage: read_decl()
        {'ABC_to_bnr', ..., 'zx_to_zv'}
    """
    s = set()
    with open(os.path.join(sage_src_pari(), "paridecl.pxd")) as f:
        s.update(decl_re.findall(f.read()))
    with open(os.path.join(sage_src_pari(), "declinl.pxi")) as f:
        s.update(decl_re.findall(f.read()))
    return s


def parse_prototype(proto, help, initial_args=[]):
    """
    Parse arguments and return type of a PARI function.

    INPUT:

    - ``proto`` -- a PARI prototype like ``"GD0,L,DGDGDG"``

    - ``help`` -- a PARI help string like
      ``"qfbred(x,{flag=0},{d},{isd},{sd})"``

    - ``initial_args`` -- other arguments to this function which come
      before the PARI arguments, for example a ``self`` argument.

    OUTPUT: a tuple ``(args, ret)`` where

    - ``args`` is a list consisting of ``initial_args`` followed by
      :class:`PariArgument` instances with all arguments of this
      function.

    - ``ret`` is a :class:`PariReturn` instance with the return type of
      this function.

    EXAMPLES::

        sage: from sage_setup.autogen.pari.parser import parse_prototype
        sage: proto = 'GD0,L,DGDGDG'
        sage: help = 'qfbred(x,{flag=0},{d},{isd},{sd})'
        sage: parse_prototype(proto, help)
        ([GEN x, long flag=0, GEN d=NULL, GEN isd=NULL, GEN sd=NULL], GEN)
        sage: parse_prototype("lp", "foo()", ["TEST"])
        (['TEST', prec precision=0], long)
    """
    # Use the help string just for the argument names.
    # "names" should be an iterator over the argument names.
    m = paren_re.search(help)
    if m is None:
        names = iter([])
    else:
        s = m.groups()[0]
        matches = [argname_re.match(x) for x in s.split(",")]
        names = (m.groups()[0] for m in matches if m is not None)

    # First, handle the return type
    try:
        c = proto[0]
        t = pari_ret_types[c]
        n = 1  # index in proto
    except (IndexError, KeyError):
        t = pari_ret_types[""]
        n = 0  # index in proto
    ret = t()

    # Go over the prototype characters and build up the arguments
    args = list(initial_args)
    default = None
    while n < len(proto):
        c = proto[n]; n += 1

        # Parse default value
        if c == "D":
            default = ""
            if proto[n] not in pari_arg_types:
                while True:
                    c = proto[n]; n += 1
                    if c == ",":
                        break
                    default += c
            c = proto[n]; n += 1
        else:
            default = None

        try:
            t = pari_arg_types[c]
            if t is None:
                raise NotImplementedError('unsupported prototype character %r' % c)
        except KeyError:
            if c == ",":
                continue  # Just skip additional commas
            else:
                raise ValueError('unknown prototype character %r' % c)

        arg = t(names, default, index=len(args))
        args.append(arg)

    return (args, ret)
