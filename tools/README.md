# Tools Directory

This folder contains various command-line tools that are used to facilitate different development tasks. Below is a brief description of each command available in this directory.

## Update Meson Build Files

This command is used to updates the Meson build files in the project. It automatically adds new source files (py, pyx) to the Meson files and removes deleted source files. This command is useful when adding or removing source files from the project.

Within an active virtual environment where `meson` is installed, run the following command:

```bash
tools/update_meson.py
```

## Update Lock Files

This command automatically updates the lock files for the project dependencies (the ones in the `./requirements/` directory). This command is useful when adding or removing dependencies from the project, or when updating the dependencies.

Within an active virtual environment where `pip-tools` and `toml` are installed, run the following command:

```bash
tools/update-lock.py
```
