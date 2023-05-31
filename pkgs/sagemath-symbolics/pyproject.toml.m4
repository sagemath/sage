[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        sage_setup     \
        sage_conf      \
        pkgconfig      \
        sagemath_environment \
        sagemath_categories \
        sagemath_polyhedra \
        cython         \
        gmpy2          \
        numpy          \
        mpmath         \
        cypari         \
        cysignals      \
                    ')]
build-backend = "setuptools.build_meta"
