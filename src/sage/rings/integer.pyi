from sage.libs.gmp.types import __mpz_struct, mpz_t, mpz_ptr
from sage.libs.gmp.mpz import mpz_set

from sage.structure.element import EuclideanDomainElement, RingElement
from sage.categories.morphism import Morphism

class Integer(EuclideanDomainElement):
    value: __mpz_struct[1]

    def set_from_mpz(self, value: mpz_t) -> None:
        ...

    def hash_c(self) -> int:
        ...

    def __pari__(self) -> object:
        ...

    def _shift_helper(self, y: object, sign: int) -> object:
        ...

    def _add_(self, other: object) -> object:
        ...

    def _mul_(self, other: object) -> object:
        ...

    def _pow_(self, other: object) -> object:
        ...

    def _and(self, other: 'Integer') -> 'Integer':
        ...

    def _or(self, other: 'Integer') -> 'Integer':
        ...

    def _xor(self, other: 'Integer') -> 'Integer':
        ...

    def _exact_log_log2_iter(self, m: 'Integer') -> int:
        ...

    def _exact_log_mpfi_log(self, m: object) -> int:
        ...

    def _valuation(self, p: 'Integer') -> RingElement:
        ...

    def _val_unit(self, p: 'Integer') -> object:
        ...

    def _divide_knowing_divisible_by(self, right: 'Integer') -> 'Integer':
        ...

    def _is_power_of(self, n: 'Integer') -> bool:
        ...

    def _pseudoprime_is_prime(self, proof: object) -> bool:
        ...

def mpz_set_str_python(z: mpz_ptr, s: str, base: int) -> int:
    ...

def smallInteger(value: int) -> 'Integer':
    ...

_small_primes_table: list[bool]

def _Integer_from_mpz(e: mpz_t) -> 'Integer':
    z = Integer.__new__(Integer)
    mpz_set(z.value, e)
    return z
