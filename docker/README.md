[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://github.com/sagemath/sage/COPYING.txt)

# Supported tags

* `latest` — the stable `master` branch [![GitHub last commit (branch)](https://img.shields.io/github/last-commit/sagemath/sage/master.svg)](https://github.com/sagemath/sage/commits/master) [![GitLab CI](https://gitlab.com/sagemath/sage/badges/master/pipeline.svg)](https://gitlab.com/sagemath/sage/commits/master)
* `x.x` — all stable releases of Sage are tagged with their version number.
* `x.x.{beta,rc}x` - betas and release candidates of Sage as [tagged in our git repository](https://github.com/sagemath/sage/tags); since docker images are built and uploaded manually at the moment, there are often no recent betas here anymore.
* `develop` — the current development version of Sage which gets merged into the `master` branch when a new version of Sage is released [![GitHub last commit (branch)](https://img.shields.io/github/last-commit/sagemath/sage/develop.svg)](https://github.com/sagemath/sage/commits/develop) [![GitLab CI](https://gitlab.com/sagemath/sage/badges/develop/pipeline.svg)](https://gitlab.com/sagemath/sage/commits/develop); again, there is often no recent develop build here anymore.
* `-py3` - until Sage 9.1, we provided Python 2 builds (with no suffix) and Python 3 builds (with the `-py3` suffix). From Sage 9.2.beta0 on, all images we provide are based on Python 3 and the `-py3` suffix survives only for historical reasons: with or without it, you get Python 3.

# What is SageMath

SageMath is a free open-source mathematics software system licensed under the GPL. It builds on top of many existing open-source packages: NumPy, SciPy, matplotlib, Sympy, Maxima, GAP, FLINT, R and many more. Access their combined power through a common, Python-based language or directly via interfaces or wrappers. 

**Mission**: *Creating a viable free open source alternative to Magma, Maple, Mathematica and Matlab.*

# What's in this image

There are several flavours of this image.

* [`sagemath/sagemath`![image size](https://img.shields.io/microbadger/image-size/sagemath/sagemath/latest.svg)](https://hub.docker.com/r/sagemath/sagemath) contains everything necessary to run Sage on the command line. Run it with:
    ```
    docker run -it sagemath/sagemath:latest
    ```
    You can start a graphical [Jupyter Notebook](https://jupyter.org) at http://localhost:8888 instead. To use the notebook, follow the instructions printed when you run:
    ```
    docker run -p 8888:8888 sagemath/sagemath:latest sage-jupyter
    ```
* [`sagemath/sagemath-dev`![image size](https://img.shields.io/microbadger/image-size/sagemath/sagemath-dev.svg)](https://hub.docker.com/r/sagemath/sagemath-dev) contains all the build artifacts to rebuild Sage quickly (currently, this is broken, see [#34241](https://github.com/sagemath/sage/issues/34241).) This version is probably only relevant for Sage developers. Run this image with:
    ```
    docker run -it sagemath/sagemath-dev:develop
    ```
    This triggers a rebuild and drops you in a shell afterwards. Note that the git repository has been emptied to save space. If you want to use git, fetch from your git repository with `git fetch trac` and go to the commit that was used to create this image with
    ```
    git reset $(cat docker/.commit)
    ```

# How to build your own SageMath images

Run `docker build -f docker/Dockerfile --build-arg ARTIFACT_BASE=sagemath/sagemath-dev:develop --target TARGET .` in the Sage repository with `TARGET` one of `sagemath` or `sagemath-dev`.

# How these images get updated

Currently, these images are updated manually after every release. Instructions for this are at the top of `docker/Dockerfile`.

Every push to our [GitLab repository](https://gitlab.com/sagemath/sage) used to trigger a pipeline in GitLab CI which pushed the actual images to Docker Hub. These pipelines are not running anymore since nobody is maintaining them or providing the CI resources for it. (Please reach out to [sage-devel](https://groups.google.com/forum/#!forum/sage-devel) if you want to help!)

Have a look at `.gitlab-ci.yml` if you want to setup GitLab CI for your own fork of the SageMath repository.

# Report bugs and issues

Please tell us of any bugs or omissions at our [Issue Tracker](https://github.com/sagemath/sage/issues) or contact us through the [sage-support](https://groups.google.com/forum/#!forum/sage-support) or the [sage-devel](https://groups.google.com/forum/#!forum/sage-devel) mailing lists.

# License

The whole Sage software distribution is licensed under the General Public License, version 3. More details can be found in our [COPYING.txt](https://github.com/sagemath/sage/blob/master/COPYING.txt)

[//]: # (Please don't break long lines in this files as dockerhub then gets the formatting of this file wrong.)
