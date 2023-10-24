# Small script that generates a meson.build file in the given folder.
# The generated build file contains all python files as `install_sources` and all cython files as `extension_module`

import argparse
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description='Generate meson.build file for a given folder.')
parser.add_argument('folder', type=str, nargs='?', default='.',
                    help='folder for which the meson.build file will be generated')
parser.add_argument('--dry-run', '-n', action='store_true',
                    help='do not write any files, just print the output')
parser.add_argument('--force', '-f', action='store_true',
                    help='overwrite existing meson.build file')

args = parser.parse_args()

folder = Path(args.folder)
dry_run = args.dry_run
force = args.force

if not folder.is_dir():
    print(f'Error: {folder} is not a directory')
    sys.exit(1)
folder_rel_to_src = folder.relative_to('src')

python_files = list(folder.glob('*.py'))
cython_files = list(folder.glob('*.pyx'))
python_files.sort()
cython_files.sort()

if not python_files and not cython_files:
    print(f'Error: {folder} does not contain any python or cython files')
    sys.exit(1)

meson_build_path = folder / 'meson.build'
if not dry_run and not force and meson_build_path.exists():
    print(f'Error: {meson_build_path} already exists, use --force to overwrite')
    sys.exit(1)

with open(meson_build_path, 'w') if not dry_run else sys.stdout as meson_build:
    meson_build.write('py.install_sources(\n')
    for file in python_files:
        meson_build.write(f"    '{file.name}',\n")
    meson_build.write(f"    subdir: '{folder_rel_to_src}',\n")
    meson_build.write(')\n\n')

    meson_build.write('extension_data = {\n')
    for file in cython_files:
        meson_build.write(f"    '{file.stem}': files('{file.name}'),\n")
    meson_build.write('}\n\n')
    meson_build.write('foreach name, pyx : extension_data\n')
    meson_build.write("    py.extension_module(name,\n")
    meson_build.write("        sources: pyx,\n")
    meson_build.write(f"        subdir: '{folder_rel_to_src}',\n")
    meson_build.write('        install: true,\n')
    meson_build.write('        dependencies: py_dep,\n')
    meson_build.write('    )\n')
    meson_build.write('endforeach\n')
