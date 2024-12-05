# Tools Directory

This folder contains various command-line tools that are used to facilitate different development tasks. Below is a brief description of each command available in this directory.

## Update Meson Build Files

This command is used to updates the Meson build files in the project. It automatically adds new source files (py, pyx) to the Meson files and removes deleted source files. This command is useful when adding or removing source files from the project.

You can use [uv](https://docs.astral.sh/uv/) to run the command:

```bash
uv run tools/update-meson.py
```

Alternatively, within an active virtual environment where `meson` is installed, run the following command:

```bash
tools/update_meson.py
```
