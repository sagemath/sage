from cpython.ref import PyTypeObject
from sage.libs.gmp.types import mpz_ptr

def Integer_AS_MPZ(x: object) -> mpz_ptr:
    ...

def is_Integer(x: object) -> bool:
    ...
