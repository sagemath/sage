r"""
Unicode Characters

This module provides Python identifiers for Unicode characters corresponding
to various mathematical symbols.

The syntax is ``unicode_XXX``, where ``XXX`` is the LaTeX name of the symbol,
stripped of any backslash or any curly brace (e.g. ``unicode_mathbbR`` for
``ℝ``).

EXAMPLES:

Operators::

    sage: from sage.typeset.unicode_characters import unicode_otimes
    sage: unicode_otimes
    '⊗'
    sage: hex(ord(_))
    '0x2297'
    sage: from sage.typeset.unicode_characters import unicode_bigotimes
    sage: unicode_bigotimes
    '⨂'
    sage: hex(ord(_))
    '0x2a02'
    sage: from sage.typeset.unicode_characters import unicode_wedge
    sage: unicode_wedge
    '∧'
    sage: hex(ord(_))
    '0x2227'
    sage: from sage.typeset.unicode_characters import unicode_bigwedge
    sage: unicode_bigwedge
    '⋀'
    sage: hex(ord(_))
    '0x22c0'
    sage: from sage.typeset.unicode_characters import unicode_partial
    sage: unicode_partial
    '∂'
    sage: hex(ord(_))
    '0x2202'

Arrows::

    sage: from sage.typeset.unicode_characters import unicode_to
    sage: unicode_to
    '→'
    sage: hex(ord(_))
    '0x2192'
    sage: from sage.typeset.unicode_characters import unicode_mapsto
    sage: unicode_mapsto
    '↦'
    sage: hex(ord(_))
    '0x21a6'

Letters::

    sage: from sage.typeset.unicode_characters import unicode_mathbbR
    sage: unicode_mathbbR
    'ℝ'
    sage: hex(ord(_))
    '0x211d'
    sage: from sage.typeset.unicode_characters import unicode_mathbbC
    sage: unicode_mathbbC
    'ℂ'
    sage: hex(ord(_))
    '0x2102'

"""

# ****************************************************************************
#       Copyright (C) 2021 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# Operators:

unicode_otimes = '\u2297'          # '⊗'

unicode_bigotimes = '\u2A02'       # '⨂'

unicode_wedge = '\u2227'           # '∧'

unicode_bigwedge = '\u22C0'        # '⋀'

unicode_partial = '\u2202'         # '∂'

# Arrows:

unicode_to = '\u2192'              # '→'

unicode_mapsto = '\u21A6'          # '↦'

# Letters:

unicode_mathbbR = '\u211D'         # 'ℝ'

unicode_mathbbC = '\u2102'         # 'ℂ'


