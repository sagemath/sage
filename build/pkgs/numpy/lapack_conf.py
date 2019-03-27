#!/usr/bin/env python

import pkgconfig, os

conf_file=open('site.cfg', 'w')

conf_file.write('[ALL]\n')
conf_file.write('library_dirs = '+ os.environ['SAGE_LOCAL']+ '/lib\n')
conf_file.write('include_dirs = '+ os.environ['SAGE_LOCAL']+ '/include\n')

pc_blas   = pkgconfig.parse('cblas blas')
pc_lapack = pkgconfig.parse('lapack')

if not (os.environ['UNAME'] == 'Darwin'):
    conf_file.write('[blas]\n')
    inc_dir = pc_blas['include_dirs']
    if len(inc_dir) > 0 :
        conf_file.write('include_dirs = '+ ':'.join(inc_dir)+'\n')
    lib_dir = pc_blas['library_dirs']
    if len(lib_dir) > 0 :
        conf_file.write('library_dirs = '+ ':'.join(lib_dir)+'\n')
    conf_file.write('blas_libs    = '+', '.join(pc_blas['libraries'])+'\n')
    conf_file.write('[lapack]\n')
    lib_dir = pc_lapack['library_dirs']
    if len(lib_dir) > 0 :
        conf_file.write('library_dirs = '+ ':'.join(lib_dir)+'\n')
    conf_file.write('lapack_libs  = '+', '.join(pc_lapack['libraries'])+'\n')

conf_file.close()
