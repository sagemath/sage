#!/usr/bin/env sage-bootstrap-python


TEMPLATE = """
SAGE_LOCAL={SAGE_LOCAL}
prefix=${{SAGE_LOCAL}}
libdir=${{prefix}}/lib
includedir=${{prefix}}/include
Name: {target}
Version: {version}
Description: {target} for SageMath, provided by the OpenBLAS package.
Cflags: -I${{includedir}}
Libs: -L${{libdir}} {libflags}
""".lstrip()


import os


try:
    SAGE_LOCAL = os.environ['SAGE_LOCAL']
except KeyError:
    raise RuntimeError('must be run in a sage shell')

SAGE_DESTDIR = os.environ.get('SAGE_DESTDIR', '')


pkgconfigdir = os.path.join(SAGE_DESTDIR, SAGE_LOCAL.lstrip('/'), 'lib',
                            'pkgconfig')
if not os.path.isdir(pkgconfigdir):
    os.makedirs(pkgconfigdir)


with open('package-version.txt') as f:
    package_version = f.read()


def write_pc_file(target, libs, version):
    filename = os.path.join(pkgconfigdir, '{0}.pc'.format(target))
    libflags=' '.join('-l{0}'.format(lib) for lib in libs)
    pc_file = TEMPLATE.format(
        SAGE_LOCAL=SAGE_LOCAL,
        target=target,
        libflags=libflags,
        version=version,
    )
    with open(filename, 'w') as f:
        f.write(pc_file)
    print('Wrote {0}'.format(filename))


write_pc_file(  'blas', libs=['openblas'], version=package_version)
write_pc_file( 'cblas', libs=['openblas'], version=package_version)
write_pc_file('lapack', libs=['openblas'], version=package_version)
