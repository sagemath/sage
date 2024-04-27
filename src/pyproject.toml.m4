[build-system]
# Minimum requirements for the build system to execute.
requires = [
     esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools \
        wheel          \
        cypari         \
        cysignals      \
        cython         \
        gmpy2          \
        jinja2         \
        jupyter_core   \
        numpy          \
        pkgconfig      \
        pplpy          \
        memory_allocator \
                    ')]
build-backend = "setuptools.build_meta"

[tool.conda-lock]
platforms = [
    'osx-64', 'linux-64', 'linux-aarch64', 'osx-arm64'
]
