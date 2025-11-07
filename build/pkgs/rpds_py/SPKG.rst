rpds_py: Python bindings to Rust's persistent data structures
==============================================================

Description
-----------

Python bindings to the Rust rpds crate for persistent data structures.

rpds-py provides efficient, immutable data structures including:

* ``HashTrieMap`` - Persistent hash map
* ``HashTrieSet`` - Persistent hash set  
* ``List`` - Persistent list with efficient operations

These data structures are backed by Rust implementations for high performance
while maintaining a Pythonic API. They are particularly useful for functional
programming patterns and situations requiring immutable, persistent collections.

The library is used by projects like the referencing library (part of the
Python JSON Schema ecosystem) as a faster alternative to pyrsistent.

License
-------

MIT License

Upstream Contact
----------------

- Author: Julian Berman <Julian+rpds@GrayVines.com>
- Home page: https://github.com/crate-py/rpds
- PyPI: https://pypi.org/project/rpds-py/
- Documentation: https://rpds.readthedocs.io/
- Upstream Rust crate: https://github.com/orium/rpds

Dependencies
------------

Python (>= 3.10)

Build dependencies: Rust toolchain (automatically handled by pip when
installing from source)

Special Notes
-------------

This package provides platform-specific binary wheels for multiple Python
versions and platforms:

* Python 3.11, 3.12, 3.13, 3.14 (including free-threaded 3.13t and 3.14t)
* Linux (x86_64, aarch64, musllinux)
* macOS (x86_64, arm64)
* Windows (win32, win_amd64, win_arm64)

The Sage build system automatically selects and downloads the appropriate
wheel for your platform and Python version using pip's auto-detection.

When building from source, a Rust toolchain is required as rpds-py contains
Rust extensions for performance.
