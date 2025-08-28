from typing import Any, List, Tuple, Union
from sage.matrix.matrix import Matrix
from sage.modules.free_module_element import FreeModuleElement
from sage.rings.integer import Integer

def listcat(l: List[List[Any]]) -> List[Any]:
    ...

class KhuriMakdisi_base:
    wL: Matrix
    w0: Matrix
    d0: int
    g: int

    def mu_image(self, wd: Matrix, we: Matrix, mu_mat: Matrix, expected_dim: int = 0) -> Matrix:
        ...

    def mu_preimage(self, we: Matrix, wde: Matrix, mu_mat: Matrix, expected_codim: int = 0) -> Matrix:
        ...

    def negate(self, wd: Matrix) -> Matrix:
        ...

    def add(self, wd1: Matrix, wd2: Matrix) -> Matrix:
        ...

    def subtract(self, wd1: Matrix, wd2: Matrix) -> Matrix:
        ...

    def multiple(self, wd: Matrix, n: int) -> Matrix:
        ...

    def zero_divisor(self) -> Matrix:
        ...

class KhuriMakdisi_large(KhuriMakdisi_base):
    mu_mat33: Matrix

    def __init__(self, V: Any, mu: Any, w0: Matrix, d0: int, g: int) -> None:
        ...

    def equal(self, wd: Matrix, we: Matrix) -> bool:
        ...

    def _add(self, wd: Matrix, we: Matrix) -> Matrix:
        ...

    def _flip(self, wd: Matrix) -> Matrix:
        ...

    def addflip(self, wd1: Matrix, wd2: Matrix) -> Matrix:
        ...

    def add_divisor(self, wd1: Matrix, wd2: Matrix, d1: int, d2: int) -> Matrix:
        ...

class KhuriMakdisi_medium(KhuriMakdisi_base):
    wV1: Matrix
    wV2: Matrix
    wV3: Matrix
    mu_mat22: Matrix
    mu_mat23: Matrix
    mu_mat31: Matrix
    mu_mat32: Matrix

    def __init__(self, V: Any, mu: Any, w0: Matrix, d0: int, g: int) -> None:
        ...

    def equal(self, wd: Matrix, we: Matrix) -> bool:
        ...

    def addflip(self, wd1: Matrix, wd2: Matrix) -> Matrix:
        ...

    def add_divisor(self, wd1: Matrix, wd2: Matrix, d1: int, d2: int) -> Matrix:
        ...

class KhuriMakdisi_small(KhuriMakdisi_base):
    wV2: Matrix
    wV3: Matrix
    wV4: Matrix
    mu_mat22: Matrix
    mu_mat23: Matrix
    mu_mat24: Matrix
    mu_mat32: Matrix
    mu_mat33: Matrix
    mu_mat34: Matrix
    mu_mat42: Matrix
    mu_mat43: Matrix

    def __init__(self, V: Any, mu: Any, w0: Matrix, d0: int, g: int) -> None:
        ...

    def equal(self, wd: Matrix, we: Matrix) -> bool:
        ...

    def addflip(self, wd1: Matrix, wd2: Matrix) -> Matrix:
        ...

    def negate(self, wd: Matrix) -> Matrix:
        ...

    def add_divisor(self, wd1: Matrix, wd2: Matrix, d1: int, d2: int) -> Matrix:
        ...
