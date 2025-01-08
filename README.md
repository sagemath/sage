<div>
<a href="https://sagemath.org">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="src/doc/common/static/logo_sagemath_white.svg">
    <img src="src/doc/common/static/logo_sagemath_black.svg" height="60" align="left">
  </picture>
</a>
   <em>"Creating a Viable Open Source Alternative to
   Magma, Maple, Mathematica, and MATLAB"</em>
</div>

#

Sage is open source mathematical software released under the GNU General Public
Licence GPLv2+, and includes packages that have [compatible software licenses](./COPYING.txt).
[People all around the globe](https://www.sagemath.org/development-map.html) have contributed to the
development of Sage. [Full documentation](https://doc.sagemath.org/html/en/index.html) is available online.

Flag algebras
-----------------

This repository is a copy of the official SageMath project with additional functionality to handle 
flag algebraic calculations. Explore the capabilities and functionalities of the package related to flag algebras by visiting the [tutorial notebook](https://github.com/bodnalev/sage/blob/theory-builder/flag_tutorial.ipynb).

Instructions to Build from Source
---------------------------------

Sage attempts to support all major Linux distributions, recent versions of
macOS, and Windows (using Windows Subsystem for Linux or
virtualization). The additional software and packages needed for flag algebraic calculations
are only tested on a few Debian based distributions. This guide might not work on other operating
systems or Linux distributions.

1.  Prepare the environment:
    - At least 10 GB of free space is required.
    - Pick a path for the folder that will contain the source files. Note it **can not contain spaces**.
    - After starting the build, you cannot move the source/build
      directory without breaking things.
    - For a minimal installation, it is enough to have the following:
        - Compilers: `gcc`, `gfortran`, `g++`
        - Build tools: GNU `make`, GNU `m4`, `perl` (including `ExtUtils::MakeMaker`), `ranlib`, `git`, `tar`, `bc`.
    - For a complete installation (recommended), see the linux system packages [in this guide](https://doc.sagemath.org/html/en/installation/source.html)
2.  Download Sage:
    - Open a terminal at the target folder and clone the Sage files there: `git clone https://github.com/bodnalev/sage.git`
    - Move inside the sage source files with the command `cd sage`
3. Optional, but highly recommended: Set some environment variables to
    customize the build.
    
    For example, the `MAKE` environment variable controls whether to
    run several jobs in parallel.  On a machine with 4 processors, say,
    typing `export MAKE="make -j4"` will configure the build script to
    perform a parallel compilation of Sage using 4 jobs. On some
    powerful machines, you might even consider `-j16`, as building with
    more jobs than CPU cores can speed things up further.
    
    To reduce the terminal output during the build, type `export V=0`.
    (`V` stands for "verbosity".)
4. Configure the build process with the command `make configure`. Then type `make build`.
5. Type `./sage` to try it out. In Sage, try for example `2 + 2`,
    `plot(x^2)`, `plot3d(lambda x, y: x*y, (-1, 1), (-1, 1))`
    to test a simple computation and plotting in 2D and 3D.
    Type <kbd>Ctrl</kbd>+<kbd>D</kbd> or `quit` to quit Sage.
6. Optional: Create a symlink to the installed `sage` script in a
    directory in your `PATH`, for example `/usr/local`. This will
    allow you to start Sage by typing `sage` from anywhere rather than
    having to either type the full path or navigate to the Sage
    directory and type `./sage`. This can be done by running:

        $ sudo ln -s $(./sage -sh -c 'ls $SAGE_ROOT/venv/bin/sage') /usr/local/bin

Changes to Included Software
----------------------------

All software included with Sage is copyrighted by the respective authors
and released under an open source license that is __GPL version 3 or
later__ compatible. See [COPYING.txt](./COPYING.txt) for more details.

Sources are in unmodified (as far as possible) tarballs in the
`upstream/` directory. The remaining description, version
information, patches, and build scripts are in the accompanying
`build/pkgs/<packagename>` directory. This directory is
part of the Sage git repository.

<p align="center">
   Copyright (C) 2005-2024 The Sage Development Team
</p>
<p align="center">
   https://www.sagemath.org
</p>

