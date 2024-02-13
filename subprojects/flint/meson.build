# Thin wrapper around the CMake project
# Once https://github.com/mesonbuild/meson/issues/10298 is fixed, we can remove this file and set "method=cmake" in flint.wrap

project('flint', 'cpp', version : '3.0.1')

cmake = import('cmake')

subproj_flint = cmake.subproject('flint-cmake')
py = import('python')
# Copy all src/.h files from the subproject to the build directory
# We cannot use meson's fs module due to sandboxing
run_command(py.find_installation('python3'), '-c', 'import shutil; shutil.copytree("' + meson.current_source_dir() + '/../flint-3.0.1/src", "' + meson.current_build_dir() + '/include/flint", dirs_exist_ok=True)', check: true)
flint_dep = declare_dependency(dependencies: [subproj_flint.dependency('flint')], include_directories: include_directories('include'))
meson.override_dependency('flint', flint_dep)
