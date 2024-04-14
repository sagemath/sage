# sage_setup: distribution = sagemath-symbolics

cpdef _repr_Expression(x)
cpdef _latex_Expression(x)
cpdef new_Expression(parent, x)
cpdef new_Expression_from_pyobject(parent, x, bint force=?, bint recursive=?)
cpdef new_Expression_wild(parent, unsigned int n=?)
cpdef new_Expression_symbol(parent, name=?, latex_name=?, domain=?)
