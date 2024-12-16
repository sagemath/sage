# Tools Directory

This folder contains various command-line tools that are used to facilitate different development tasks. Below is a brief description of each command available in this directory.

## Update Conda Environment Files

This command is used to update the Conda environment files in the project. It automatically adds new dependencies to the Conda files, removes deleted dependencies, and updates the version of existing dependencies. The source of the dependencies is the `pyproject.toml` file, which specifies the following dependencies:

- `build-system.requires`: Python dependencies required for building
- `project.dependencies`: Python dependencies required for running
- `external.build-requires`: External dependencies required for building
- `external.host-requires`: External dependencies required for running


Within an active virtual environment where `grayskull`, `conda-lock` and `toml` are installed, run the following command:

```bash
tools/update-conda.py
```

## Update Meson Build Files

This command is used to updates the Meson build files in the project. It automatically adds new source files (py, pyx) to the Meson files and removes deleted source files. This command is useful when adding or removing source files from the project.

You can use [uv](https://docs.astral.sh/uv/) to run the command:

```bash
uv run tools/update-meson.py
```

Alternatively, within an active virtual environment where `meson` is installed, run the following command:

```bash
tools/update-meson.py
```
