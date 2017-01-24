from sage.misc.superseded import deprecation
deprecation(20205, '''pari/decl.pxi is deprecated, use "from sage.libs.cypari2.paridecl cimport *" instead''')

from sage.libs.cypari2.paridecl cimport *
