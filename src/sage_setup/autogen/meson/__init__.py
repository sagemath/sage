# Small script that generates a meson.build file in the given folder.
# The generated build file contains all python files as `install_sources` and all cython files as `extension_module`

import sys
from pathlib import Path
from types import SimpleNamespace


def run(folder: Path, dry_run=False, force=False):

    if not folder.is_dir():
        print(f'Error: {folder} is not a directory')
        return
    folder_rel_to_src = folder.relative_to('src')

    python_files = sorted(list(folder.glob('*.py')) + list(folder.glob('*.pxd')) + list(folder.glob('*.h')))
    cython_files = sorted(list(folder.glob('*.pyx')))


    if not python_files and not cython_files:
        print(f'Error: {folder} does not contain any python or cython files')
        return

    def get_metadata(path: Path):
        with open(path, 'r') as file:
            metadata = SimpleNamespace()
            metadata.path = path
            metadata.libraries = []
            for line in file:
                if line.startswith('# distutils: libraries ='):
                    libraries = line.split('libraries =')[1].strip().split()
                    libraries = [
                        library
                        .replace('ARB_LIBRARY', 'arb')
                        .replace('NTL_LIBRARIES', 'ntl')
                        .replace('SINGULAR_LIBRARIES', 'singular')
                        .replace('LINBOX_LIBRARIES', 'linbox')
                        .replace('FFLASFFPACK_LIBRARIES', 'fflas')
                        .replace('GSL_LIBRARIES', 'gsl')
                        .replace('M4RI_LIBRARIES', 'm4ri')
                        .replace('GDLIB_LIBRARIES', 'gd')
                        .replace('LIBPNG_LIBRARIES', 'png')
                        .replace('CBLAS_LIBRARIES', 'cblas')
                        .replace('ZLIB_LIBRARIES', 'zlib')
                        .replace('Lfunction', 'lfunction')
                          for library in libraries]
                    try:
                        libraries.remove('CYGWIN_SQLITE3_LIBS')
                    except ValueError:
                        pass
                    metadata.libraries += libraries

            metadata.inc_dirs = []
            c_file = path.with_suffix('.c')
            cpp_file = path.with_suffix('.cpp')
            if cpp_file.exists():
                metadata.is_cpp = True
                c_file = cpp_file
            else:
                metadata.is_cpp = False
            if c_file.exists():
                metadata.not_yet_on_conda = False
                with open(c_file, 'r') as c_file:
                    contents = c_file.read()
                    known_inc_dirs = {
                        'sage/cpython/': 'inc_cpython',
                        'sage/rings': 'inc_rings',
                        'sage/rings/finite_rings': 'inc_rings_finite',
                        'sage/libs/flint': 'inc_flint',
                        'sage/libs/gsl': 'inc_gsl',
                        'sage/libs/ntl': 'inc_ntl',
                        'sage/libs/arb': 'inc_arb',
                        'sage/data_structures': 'inc_data_structures',
                        'sage/ext/': 'inc_ext',
                        'numpy/core/include/': 'inc_numpy',
                        'sage/symbolic/ginac/': 'inc_ginac',
                        'sage/symbolic/pynac_wrap.h': 'inc_pynac',
                        'sage/groups/perm_gps/partn_ref2/': 'inc_partn_ref2',
                        'sage/ext/interpreters/': 'inc_interpreters',
                    }
                    for known_inc_dir in known_inc_dirs:
                        if known_inc_dir in contents:
                            metadata.inc_dirs.append(known_inc_dirs[known_inc_dir])
                    known_libs = {
                        'cypari2/': 'cypari2',
                        'cysignals/': 'cysignals',
                        '/gmp.h': 'gmp',
                        '/iml.h': 'iml',
                        '/m4ri/': 'm4ri',
                        '/pari/': 'pari',
                        '/flint/': 'flint',
                        '/fflas-ffpack/': 'fflas',
                        '/givaro/': 'givaro',
                        '/gmp++/': 'gmpxx',
                        '/linbox/': 'linbox',
                        '/gsl/': 'gsl',
                        'mpfr.h': 'mpfr',
                        'sage/symbolic/ginac/': 'ginac',
                        'arb.h': 'arb',
                        'mpfi.h': 'mpfi',
                        'mpc.h': 'mpc',
                        'gmpy2/': 'gmpy2',
                    }
                    for known_lib in known_libs:
                        if known_lib in contents:
                            metadata.libraries.append(known_libs[known_lib])
            else:
                metadata.not_yet_on_conda = metadata.is_cpp is False

            return metadata

    cython_files = [get_metadata(file) for file in cython_files]
    cython_c_files = [file for file in cython_files if not file.is_cpp]
    cython_cpp_files = [file for file in cython_files if file.is_cpp]
    all_libraries = sorted(set(library for file in cython_files for library in file.libraries) | {'gmp'})
    all_inc_dirs = sorted(set(inc_dir for file in cython_files for inc_dir in file.inc_dirs))
    subdirs = sorted(list(folder.glob('*/')))

    meson_build_path = folder / 'meson.build'
    if not dry_run and not force and meson_build_path.exists():
        print(f'Error: {meson_build_path} already exists, use --force to overwrite')
        return

    with open(meson_build_path, 'w') if not dry_run else sys.stdout as meson_build:
        meson_build.write('py.install_sources(\n')
        for file in python_files:
            meson_build.write(f"    '{file.name}',\n")
        meson_build.write(f"    subdir: '{folder_rel_to_src}',\n")
        meson_build.write(')\n')

        if cython_c_files:
            meson_build.write('\n')
            meson_build.write('extension_data = {\n')
            for file in cython_c_files:
                if file.not_yet_on_conda:
                    meson_build.write(f"    # '{file.path.stem}': files('{file.path.name}'), # not yet on conda\n")
                else:
                    meson_build.write(f"    '{file.path.stem}': files('{file.path.name}'),\n")
            meson_build.write('}\n\n')

            meson_build.write('foreach name, pyx : extension_data\n')
            meson_build.write("    py.extension_module(name,\n")
            meson_build.write("        sources: pyx,\n")
            meson_build.write(f"        subdir: '{folder_rel_to_src}',\n")
            meson_build.write('        install: true,\n')
            meson_build.write(f"        include_directories: [{', '.join(all_inc_dirs)}],\n")
            meson_build.write(f"        dependencies: [py_dep{', ' if all_libraries else ''}{', '.join(all_libraries)}],\n")
            meson_build.write('    )\n')
            meson_build.write('endforeach\n')
        
        if cython_cpp_files:
            meson_build.write('\n')
            meson_build.write('extension_data_cpp = {\n')
            for file in cython_cpp_files:
                if file.not_yet_on_conda:
                    meson_build.write(f"    # '{file.path.stem}': files('{file.path.name}'), # not yet on conda\n")
                else:
                    meson_build.write(f"    '{file.path.stem}': files('{file.path.name}'),\n")
            meson_build.write('}\n\n')

            meson_build.write('foreach name, pyx : extension_data_cpp\n')
            meson_build.write("    py.extension_module(name,\n")
            meson_build.write("        sources: pyx,\n")
            meson_build.write(f"        subdir: '{folder_rel_to_src}',\n")
            meson_build.write('        install: true,\n')
            meson_build.write('        override_options : [\'cython_language=cpp\'],\n')
            meson_build.write(f"        include_directories: [{', '.join(all_inc_dirs)}],\n")
            meson_build.write(f"        dependencies: [py_dep{', ' if all_libraries else ''}{', '.join(all_libraries)}],\n")
            meson_build.write('    )\n')
            meson_build.write('endforeach\n')

        meson_build.write('\n')
        for subdir in subdirs:
            if subdir.name.startswith('_') or subdir.name.startswith('.'):
                continue
            if not list(subdir.glob('*.py*')):
                continue
            
            if not list(subdir.glob('*.pyx')):
                meson_build.write(f"install_subdir('{subdir.name}', install_dir: sage_install_dir / '{folder_rel_to_src.relative_to('sage')}')\n")
            else:
                meson_build.write(f"subdir('{subdir.name}')\n")
                run(subdir)


def generate_meson():
    import argparse

    parser = argparse.ArgumentParser(description='Generate meson.build file for a given folder.')
    parser.add_argument('folder', type=str, nargs='?', default='.',
                        help='folder for which the meson.build file will be generated')
    parser.add_argument('--dry-run', '-n', action='store_true',
                        help='do not write any files, just print the output')
    parser.add_argument('--force', '-f', action='store_true',
                        help='overwrite existing meson.build file')

    args = parser.parse_args()

    folder = Path(args.folder)
    run(folder, dry_run=args.dry_run, force=args.force)
    return 0
