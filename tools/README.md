# Tools Directory

This folder contains various command-line tools that are used to facilitate different development tasks. Below is a brief description of each command available in this directory.

## Update Meson Build Files

This command is used to updates the Meson build files in the project. It automatically adds new source files (py, pyx) to the Meson files and removes deleted source files. This command is useful when adding or removing source files from the project.

Within an active virtual environment where Meson is installed, run the following command:

```bash
tools/update_meson.py
```

## Update Version Number

Increments the version number in the project. This command is useful when releasing a new version of the project.

Set `SAGE_ROOT` to the root directory of the Sage project and run the following command:

```bash
tools/update_version
```
