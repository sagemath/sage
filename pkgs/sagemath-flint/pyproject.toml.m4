[build-system]
# Minimum requirements for the build system to execute.
# sage_conf is needed for library name of the flint library.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        wheel          \
        sage_setup     \
        sage_conf      \
        sagemath_environment \
        sagemath_categories \
        cython         \
        gmpy2          \
        cypari         \
        cysignals      \
                    ')]
build-backend = "setuptools.build_meta"
