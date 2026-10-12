from typing import final

from sage.rings.integer import Integer

def prime_range(start: Integer | int, stop: Integer | int | None = None, step: Integer | int | None = None, algorithm: str | None = None, py_ints: bool = False) -> list[int]:
    ...

@final
class arith_int:
    def abs_int(self, x: int) -> int:
        ...

    def sign_int(self, n: int) -> int:
        ...

    def c_gcd_int(self, a: int, b: int) -> int:
        ...

    def gcd_int(self, a: int, b: int) -> int:
        ...

    def c_xgcd_int(self, a: int, b: int) -> tuple[int, int, int]:
        ...

    def xgcd_int(self, a: int, b: int) -> tuple[int, int, int]:
        ...

    def c_inverse_mod_int(self, a: int, m: int) -> int:
        ...

    def inverse_mod_int(self, a: int, m: int) -> int:
        ...

    def c_rational_recon_int(self, a: int, m: int) -> tuple[int, int]:
        ...

    def rational_recon_int(self, a: int, m: int) -> tuple[int, int]:
        ...

@final
class arith_llong:
    def abs_longlong(self, x: int) -> int:
        ...

    def sign_longlong(self, n: int) -> int:
        ...

    def c_gcd_longlong(self, a: int, b: int) -> int:
        ...

    def gcd_longlong(self, a: int, b: int) -> int:
        ...

    def c_xgcd_longlong(self, a: int, b: int) -> tuple[int, int, int]:
        ...

    def xgcd_longlong(self, a: int, b: int) -> tuple[int, int, int]:
        ...

    def c_inverse_mod_longlong(self, a: int, m: int) -> int:
        ...

    def inverse_mod_longlong(self, a: int, m: int) -> int:
        ...

    def c_rational_recon_longlong(self, a: int, m: int) -> tuple[int, int]:
        ...

    def rational_recon_longlong(self, a: int, m: int) -> tuple[int, int]:
        ...
