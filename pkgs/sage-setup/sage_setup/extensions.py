
def create_extension(template, kwds):
    from Cython.Build.Dependencies import default_create_extension
    from sage.env import sage_include_directories

    # Add numpy and source folder to the include search path used by the compiler
    # This is a workaround for https://github.com/cython/cython/issues/1480
    include_dirs = kwds.get('include_dirs', []) + sage_include_directories(use_sources=True)
    kwds['include_dirs'] = include_dirs
    return default_create_extension(template, kwds)
