
from sage.config import get_include_dirs


def test_cython_metaclass_header_found():
    dirs = get_include_dirs()
    assert any(
        (dir / "sage" / "cpython" / "cython_metaclass.h").is_file() for dir in dirs
    )
