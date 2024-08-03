# Small script that generates a meson.build file in the given folder.
# The generated build file contains all python files as `install_sources` and all cython files as `extension_module`

import os
import sys

from collections import defaultdict
from difflib import unified_diff
from pathlib import Path
from tempfile import NamedTemporaryFile
from types import SimpleNamespace


def distribution_condition(distribution: str) -> str:
    distribution_variable = 'distribution_' + distribution.replace('-', '_')
    return f"get_variable('{distribution_variable}', false)"


def distributions_shorthand(distributions: set) -> str:

    sagemath_distributions = sorted(distribution.removeprefix('sagemath-')
                                    for distribution in distributions
                                    if distribution.startswith('sagemath-'))
    distributions = sorted(distribution or "''"  # the catch-all distribution
                           for distribution in distributions
                           if not distribution.startswith('sagemath-'))
    if len(sagemath_distributions) > 1:
        distributions.insert(0, 'sagemath-{' + ','.join(sagemath_distributions) + '}')
    elif sagemath_distributions:
        distributions.insert(0, 'sagemath-' + sagemath_distributions[0])
    if distributions:
        return ', '.join(distributions)
    return 'no distribution'


def log_update(operation, meson_build_path, distributions, has_cython):
    line = f'{operation:12}{meson_build_path}'
    line = f'{line:62}  # {distributions_shorthand(distributions)}'
    if len(distributions) == 1 and has_cython:
        line += ', has Cython modules'
    print(line, file=sys.stderr)


def run(folder: Path, output_dir: Path, folder_rel_to_src=None, dry_run=False, force=False, monolithic=False):

    if not folder.is_dir():
        print(f'Error: {folder} is not a directory')
        return
    try:
        from sage.env import SAGE_SRC
    except ImportError:
        folder_rel_to_src = folder.relative_to('src')
    else:
        try:
            folder_rel_to_src = folder.relative_to(SAGE_SRC)
        except ValueError:
            folder_rel_to_src = folder.resolve().relative_to(Path(SAGE_SRC).resolve())

    if monolithic:
        def read_distribution(path):
            return ''
    else:
        from sage.misc.package_dir import read_distribution

    recurse_subdirs = {}
    by_distribution = defaultdict(lambda: SimpleNamespace(python_files=[],
                                                          cython_c_files=[],
                                                          cython_cpp_files=[],
                                                          install_subdirs=[]))
    distributions = set()
    subdirs = sorted(list(folder.glob('*/')))
    has_cython = False

    for subdir in subdirs:
        if subdir.name.startswith('_') or subdir.name.startswith('.'):
            continue
        subdir_distributions, subdir_has_cython = run(subdir, output_dir / subdir.name,
                                                      dry_run=dry_run, force=force, monolithic=monolithic)
        has_cython = has_cython or subdir_has_cython
        if not subdir_distributions:
            pass
        elif not monolithic or len(subdir_distributions) > 1 or subdir_has_cython:
            recurse_subdirs[subdir] = subdir_distributions
        else:
            by_distribution[list(subdir_distributions)[0]].install_subdirs.append(subdir)
        distributions.update(subdir_distributions)

    python_files = sorted(list(folder.glob('*.py')) + list(folder.glob('*.pxd')) + list(folder.glob('*.h')))
    cython_files = sorted(list(folder.glob('*.pyx')))

    def get_metadata(path: Path):
        with open(path, 'r') as file:
            metadata = SimpleNamespace()
            metadata.path = path
            metadata.distribution = read_distribution(path)

            if path.suffix not in [".pxd", ".pxi", ".pyx"]:
                return metadata

            metadata.libraries = []
            for line in file:
                if not line:
                    pass
                elif line.startswith('# distutils: libraries ='):
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
                elif line.startswith('#'):
                    pass
                else:
                    break

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

    python_files = [get_metadata(file) for file in python_files]
    cython_files = [get_metadata(file) for file in cython_files]
    cython_c_files = [file for file in cython_files if not file.is_cpp]
    cython_cpp_files = [file for file in cython_files if file.is_cpp]
    all_libraries = sorted(set(library for file in cython_files for library in file.libraries) | {'gmp'})
    all_inc_dirs = sorted(set(inc_dir for file in cython_files for inc_dir in file.inc_dirs))

    template_path = folder / 'meson.build.in'
    meson_build_path = output_dir / 'meson.build'
    if not dry_run and not force and meson_build_path.exists():
        print(f'Error: {meson_build_path} already exists, use --force to overwrite')
        return

    for python_file in python_files:
        by_distribution[python_file.distribution].python_files.append(python_file)
    for cython_c_file in cython_c_files:
        by_distribution[cython_c_file.distribution].cython_c_files.append(cython_c_file)
    for cython_cpp_file in cython_cpp_files:
        by_distribution[cython_cpp_file.distribution].cython_cpp_files.append(cython_cpp_file)

    distributions.update(by_distribution)
    if cython_files:
        has_cython = True

    if monolithic and not has_cython and len(distributions) <= 1:
        # No need for a meson.build file
        try:
            meson_build_path.unlink()
        except FileNotFoundError:
            log_update('Not writing', meson_build_path, distributions, has_cython)
        else:
            log_update('Removing', meson_build_path, distributions, has_cython)
        return distributions, has_cython

    os.makedirs(meson_build_path.parent, exist_ok=True)

    with NamedTemporaryFile(mode='w', prefix=meson_build_path.name,
                            dir=meson_build_path.parent, delete=False) if not dry_run else sys.stdout as meson_build:

        if not monolithic:
            meson_build.write(f"# Automatically generated by sage-generate-meson; do not edit manually\n\n")

        try:
            with open(template_path, "r") as template:
                if not monolithic:
                    meson_build.write(f"#line 1 {template_path}\n")
                for line in template:
                    meson_build.write(line)
            if not monolithic:
                meson_build.write(f"\n# Automatically generated from here\n\n")
        except FileNotFoundError:
            pass

        for distribution, files in by_distribution.items():

            if len(distributions) > 1:
                meson_build.write(f"if {distribution_condition(distribution)}\n")

            if files.python_files:
                if monolithic and len(files.python_files) <= 4:
                    meson_build.write('py.install_sources(')
                    for file in files.python_files:
                        meson_build.write(f"'{file.path.name}', ")
                    meson_build.write(f"subdir : '{folder_rel_to_src}'")
                    meson_build.write(')\n')
                else:
                    meson_build.write('py.install_sources(\n')
                    for file in files.python_files:
                        meson_build.write(f"  '{file.path.name}',\n")
                    meson_build.write(f"  subdir : '{folder_rel_to_src}'\n")
                    meson_build.write(')\n')

            if files.cython_c_files:
                meson_build.write('\n')
                meson_build.write('extension_data = {\n')
                items = []
                for file in files.cython_c_files:
                    if file.not_yet_on_conda:
                        items.append(f"  # '{file.path.stem}' : files('{file.path.name}'), # not yet on conda")
                    else:
                        items.append(f"  '{file.path.stem}' : files('{file.path.name}')")
                meson_build.write(',\n'.join(items))
                meson_build.write('\n}\n\n')

                meson_build.write('foreach name, pyx : extension_data\n')
                meson_build.write("    py.extension_module(name,\n")
                meson_build.write("        sources: pyx,\n")
                meson_build.write(f"        subdir: '{folder_rel_to_src}',\n")
                meson_build.write('        install: true,\n')
                meson_build.write(f"        include_directories: [{', '.join(all_inc_dirs)}],\n")
                meson_build.write(f"        dependencies: [py_dep{', ' if all_libraries else ''}{', '.join(all_libraries)}],\n")
                meson_build.write('    )\n')
                meson_build.write('endforeach\n')

            if files.cython_cpp_files:
                meson_build.write('\n')
                meson_build.write('extension_data_cpp = {\n')
                for file in files.cython_cpp_files:
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

            if not monolithic and files.install_subdirs:
                for subdir in files.install_subdirs:
                    meson_build.write(f"install_subdir('{subdir.name}', install_dir: sage_install_dir / '{folder_rel_to_src.relative_to('sage')}')\n")
                meson_build.write('\n')

            if len(distributions) > 1:
                meson_build.write('endif ########################################################################\n')

            meson_build.write('\n')

        if monolithic:
            for subdir in sorted(set(recurse_subdirs).union(files.install_subdirs)):
                if subdir in recurse_subdirs:
                    meson_build.write(f"subdir('{subdir.name}')\n")
                else:
                    meson_build.write(f"install_subdir('{subdir.name}', install_dir: sage_install_dir / '{folder_rel_to_src.relative_to('sage')}')\n")
        else:
            for subdir, subdir_distributions in recurse_subdirs.items():
                condition = " or ".join(distribution_condition(distribution)
                                        for distribution in sorted(subdir_distributions))
                meson_build.write(f"if {condition}\n    subdir('{subdir.name}')\nendif\n")

    meson_build_tmp_path = Path(meson_build.name)

    try:
        with open(meson_build_path, "r") as old:
            old_lines = list(old)
    except FileNotFoundError:
        log_update('Creating', meson_build_path, distributions, has_cython)
        meson_build_tmp_path.rename(meson_build_path)
    else:
        with open(meson_build_tmp_path, "r") as new:
            new_lines = list(new)
        changes = list(unified_diff(old_lines, new_lines, str(meson_build_path), str(os.path.relpath(meson_build_tmp_path))))
        if changes:
            meson_build_tmp_path.replace(meson_build_path)
            log_update('Updating', meson_build_path, distributions, has_cython)
            sys.stdout.writelines(changes)
        else:
            log_update('Unchanged', meson_build_path, distributions, has_cython)
            meson_build_tmp_path.unlink()

    return distributions, has_cython


def generate_meson():
    import argparse

    parser = argparse.ArgumentParser(description='Generate meson.build file for a given folder.')
    parser.add_argument('folder', type=str, nargs='?', default='.',
                        help='folder for which the meson.build file will be generated')
    parser.add_argument('--dry-run', '-n', action='store_true',
                        help='do not write any files, just print the output')
    parser.add_argument('--force', '-f', action='store_true',
                        help='overwrite existing meson.build files')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='folder where the generated meson.build files will be generated')
    parser.add_argument('--monolithic', action='store_true',
                        help='ignore distributions')

    args = parser.parse_args()

    folder = Path(args.folder)
    if args.output_dir is None:
        output_dir = folder
    else:
        output_dir = Path(args.output_dir)

    run(folder, output_dir, dry_run=args.dry_run, force=args.force, monolithic=args.monolithic)
    return 0
