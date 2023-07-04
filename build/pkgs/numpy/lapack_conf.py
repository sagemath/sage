#!/usr/bin/env python3

import pkgconfig, os

conf_file=open('site.cfg', 'w')

pc_blas   = pkgconfig.parse('cblas blas')
pc_lapack = pkgconfig.parse('lapack')

if os.environ['UNAME'] == 'Darwin':
    # on macOS, if openblas is available (the default), use it instead of the
    # macOS Accelerate framework
    if 'openblas' in pc_blas.get('libraries', []):
        conf_file.write('[openblas]\n')
        conf_file.write('libraries    = '+', '.join(pc_blas['libraries'])+'\n')
        inc_dir = pc_blas['include_dirs']
        if len(inc_dir) > 0:
            conf_file.write('include_dirs = '+ ':'.join(inc_dir)+'\n')
        lib_dir = pc_blas['library_dirs']
        if len(lib_dir) > 0:
            conf_file.write('library_dirs = '+ ':'.join(lib_dir)+'\n')
else:
    conf_file.write('[blas]\n')
    inc_dir = pc_blas['include_dirs']
    if len(inc_dir) > 0:
        conf_file.write('include_dirs = '+ ':'.join(inc_dir)+'\n')
    lib_dir = pc_blas['library_dirs']
    if len(lib_dir) > 0:
        conf_file.write('library_dirs = '+ ':'.join(lib_dir)+'\n')
    conf_file.write('blas_libs    = '+', '.join(pc_blas['libraries'])+'\n')
    conf_file.write('[lapack]\n')
    lib_dir = pc_lapack['library_dirs']
    if len(lib_dir) > 0:
        conf_file.write('library_dirs = '+ ':'.join(lib_dir)+'\n')
    conf_file.write('lapack_libs  = '+', '.join(pc_lapack['libraries'])+'\n')

conf_file.close()
