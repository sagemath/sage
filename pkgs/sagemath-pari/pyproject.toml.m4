[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        sage_setup     \
        sagemath_environment \
        sagemath_categories \
        cython         \
        cysignals      \
        cypari \
        memory_allocator \
                    ')]
build-backend = "setuptools.build_meta"
