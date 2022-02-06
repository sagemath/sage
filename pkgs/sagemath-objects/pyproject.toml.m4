[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        wheel          \
        sage_setup     \
        cython         \
        gmpy2          \
        cysignals      \
                    ')]
build-backend = "setuptools.build_meta"
