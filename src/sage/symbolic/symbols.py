r"""
Symbol table
"""
# ****************************************************************************
#       Copyright (C) 2009 Mike Hansen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

symbol_table = {'functions': {}}


def register_symbol(obj, conversions, nargs=None):
    """
    Add an object to the symbol table, along with how to convert it to
    other systems such as Maxima, Mathematica, etc.  This table is used
    to convert *from* other systems back to Sage.

    INPUT:

    - ``obj`` -- a symbolic object or function

    - ``conversions`` -- dictionary of conversions, where the keys
      are the names of interfaces (e.g., ``'maxima'``), and the values
      are the string representation of ``obj`` in that system

    - ``nargs`` -- (optional) number of arguments; for most functions,
      this can be deduced automatically

    EXAMPLES::

        sage: from sage.symbolic.symbols import register_symbol as rs
        sage: rs(SR(5), {'maxima': 'five'})                                             # needs sage.symbolic
        sage: SR(maxima_calculus('five'))                                               # needs sage.symbolic
        5
    """
    conversions = dict(conversions)
    try:
        conversions['sage'] = obj.name()
    except AttributeError:
        pass
    if nargs is None:
        try:
            nargs = obj.number_of_arguments()
        except AttributeError:
            nargs = -1  # meaning unknown number of arguments
    for system, name in conversions.items():
        system_table = symbol_table.get(system, None)
        if system_table is None:
            symbol_table[system] = system_table = {}
        system_table[(name, nargs)] = obj
