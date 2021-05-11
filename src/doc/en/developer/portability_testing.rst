.. nodoctest

.. highlight:: shell-session

.. _chapter-portability_testing:

=============================
Testing on multiple platforms
=============================

Sage is intended to build and run on a variety of platforms,
including all major Linux distributions, as well as MacOS, and
Windows (with Cygwin and WSL).

There is considerable variation among these platforms.
To ensure that Sage continues to build correctly on users'
machines, it is crucial to test changes to Sage, in particular
when external packages are added or upgraded, on a wide
spectrum of platforms.

Sage patchbots
==============

The `Sage patchbots <https://wiki.sagemath.org/patchbot>`_ will
automatically test your Trac ticket by attempting an incremental build
of Sage and running doctests.

Sage buildbots
==============

The `Sage Release buildbot <https://wiki.sagemath.org/buildbot>`_
builds entire tarballs (e.g., all the development releases) on a
variety of machines.

Developers' and users' tests on sage-release
============================================

Sage developers and users are encouraged to contribute to testing
releases that are announced on `Sage Release
<https://groups.google.com/forum/#!forum/sage-release>`_ on their
machines and to report test results (success and failures) by
responding to the announcements.

Testing Sage on a different platform using Docker
=================================================

`Docker <https://www.docker.com>`_ is a popular virtualization
software, running Linux operating system images ("Docker images") in
containers on a shared Linux kernel.  These containers can be run
using a Docker client on your Linux, Mac, or Windows box, as well as
on various cloud services.

To get started, you need to install a `Docker client
<https://docs.docker.com/install/>`_.  The clients are available for
Linux, Mac, and Windows.  The clients for the latter are known as
"Docker Desktop".

All examples in this section were obtained using Docker Desktop for
Mac; but the `command-line user interface
<https://docs.docker.com/engine/reference/commandline/cli/>`_ for the
other platforms is identical.

All major Linux distributions provide ready-to-use Docker images,
which are published via `Docker Hub <https://hub.docker.com>`_.  For
example, to run the current stable (LTS) version of Ubuntu
interactively, you can use the shell command::

  [mkoeppe@sage sage]$ docker run -it ubuntu:latest
  root@9f3398da43c2:/# 

Here ``ubuntu`` is referred to as the "image (name)" and ``latest`` as
the "tag".  Other releases of Ubuntu are available under different
tags, such as ``xenial`` or ``devel``.

The above command drops you in a root shell on the container::

  root@9f3398da43c2:/# uname -a
  Linux 9f3398da43c2 4.19.76-linuxkit #1 SMP Thu Oct 17 19:31:58 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux
  root@9f3398da43c2:/# df -h
  Filesystem      Size  Used Avail Use% Mounted on
  overlay         181G  116G   56G  68% /
  tmpfs            64M     0   64M   0% /dev
  tmpfs           2.7G     0  2.7G   0% /sys/fs/cgroup
  shm              64M     0   64M   0% /dev/shm
  /dev/sda1       181G  116G   56G  68% /etc/hosts
  tmpfs           2.7G     0  2.7G   0% /proc/acpi
  tmpfs           2.7G     0  2.7G   0% /sys/firmware

Exiting the shell terminates the container::

  root@9f3398da43c2:/# ^D
  [mkoeppe@sage sage]$

Let us work with a distclean Sage source tree.  If you are using git,
a good way to get one (without losing a precious installation in
``SAGE_LOCAL``) is by creating a new worktree::

  [mkoeppe@sage sage] git worktree add worktree-ubuntu-latest
  [mkoeppe@sage sage] cd worktree-ubuntu-latest
  [mkoeppe@sage worktree-ubuntu-latest] ls
  COPYING.txt ... Makefile ... configure.ac ... src tox.ini

This is not bootstrapped (``configure`` is missing), so let's bootstrap it::

  [mkoeppe@sage worktree-ubuntu-latest] make configure
  ...

We can start a container again with same image, ``ubuntu:latest``, but
this time let's mount the current directory into it::

  [mkoeppe@sage worktree-ubuntu-latest]$ docker run -it --mount type=bind,source=$(pwd),target=/sage ubuntu:latest
  root@39d693b2a75d:/# mount | grep sage
  osxfs on /sage type fuse.osxfs (rw,nosuid,nodev,relatime,user_id=0,group_id=0,allow_other,max_read=1048576)
  root@39d693b2a75d:/# cd sage
  root@39d693b2a75d:/sage# ls
  COPYING.txt ... Makefile ... config configure configure.ac ... src tox.ini
  
Typical Docker images provide minimal installations of packages only::
  
  root@39d693b2a75d:/sage# command -v python
  root@39d693b2a75d:/sage# command -v gcc
  root@39d693b2a75d:/sage#

As you can see above, the image ``ubuntu:latest`` has neither a Python nor
a GCC installed, which are among the build prerequisites of Sage.  We
need to install them using the distribution's package manager first.
 
Sage facilitates testing various distributions on Docker as follows.

Discovering the system's package system
---------------------------------------

::

  root@39d693b2a75d:/sage# build/bin/sage-guess-package-system 
  debian

Let's install gcc, hoping that the Ubuntu package providing it is
simply named ``gcc``.  If we forgot what the package manager on
Debian-derived distributions is called, we can ask Sage for a
reminder::

  root@39d693b2a75d:/sage# build/bin/sage-print-system-package-command debian install gcc
  apt-get install gcc

We remember that we need to fetch the current package lists
from the server first::

  root@39d693b2a75d:/sage# apt-get update
  root@39d693b2a75d:/sage# apt-get install gcc

Using Sage's database of distribution prerequisites
---------------------------------------------------

The source code of the Sage distribution contains a database of
package names in various distributions' package managers.  For
example, the file ``build/pkgs/_prereq/distros/debian.txt`` contains the following

.. code-block:: yaml

  # This file, build/pkgs/_prereq/distros/debian.txt, contains names
  # of Debian/Ubuntu packages needed for installation of Sage from source.
  #
  # In addition, the files build/pkgs/SPKG/distros/debian.txt contain the names
  # of packages that provide the equivalent of SPKG.
  #
  # Everything on a line after a # character is ignored.
  binutils
  make
  m4
  perl
  # python3-minimal is not enough on debian buster, ubuntu bionic - it does not have urllib
  python3    # system python for bootstrapping the build
  tar
  bc
  gcc
  # On debian buster, need C++ even to survive 'configure'. Otherwise:
  # checking how to run the C++ preprocessor... /lib/cpp
  # configure: error: in `/sage':
  # configure: error: C++ preprocessor "/lib/cpp" fails sanity check
  g++
  # Needed if we download some packages from a https upstream URL
  ca-certificates

From this information, we know that we can use the following command
on our container to install the necessary build prerequisites::

  root@39d693b2a75d:/sage# apt-get install binutils make m4 perl python3 tar bc gcc g++ ca-certificates
  Reading package lists... Done
  Building dependency tree       
  Reading state information... Done
  tar is already the newest version (1.29b-2ubuntu0.1).
  The following additional packages will be installed:
  ...
  Done.

(The Sage `Installation Guide <../installation/index.html>`_ also
provides such command lines for some distributions; these are
automatically generated from the database of package names.)

Now we can start the build::

  root@39d693b2a75d:/sage# ./configure 
  checking for a BSD-compatible install... /usr/bin/install -c
  checking for root user... yes
  configure: error: You cannot build Sage as root, switch to an unprivileged user.  (If building in a container, use --enable-build-as-root.)

Let's just follow this helpful hint::

  root@39d693b2a75d:/sage# ./configure --enable-build-as-root
  checking for a BSD-compatible install... /usr/bin/install -c
  ...


.. _section-equiv-distro-packages:

Using Sage's database of equivalent distribution packages
---------------------------------------------------------

At the end of the ``./configure`` run, Sage issued a message like the
following::

  configure: notice: the following SPKGs did not find equivalent system packages: arb boost boost_cropped bzip2 ... yasm zeromq zlib
  checking for the package system in use... debian
  configure: hint: installing the following system packages is recommended and may avoid building some of the above SPKGs from source:
  configure:   $ sudo apt-get install libflint-arb-dev ... yasm libzmq3-dev libz-dev
  configure: After installation, re-run configure using:
  configure:   $ ./config.status --recheck && ./config.status

This information comes from Sage's database of equivalent distribution
packages.  For example::

  root@39d693b2a75d:/sage# ls build/pkgs/arb/distros/
  arch.txt	conda.txt	debian.txt	gentoo.txt
  root@39d693b2a75d:/sage# cat build/pkgs/arb/distros/debian.txt 
  libflint-arb-dev

Note that these package equivalencies are based on a current stable or
testing version of the distribution; the packages are not guaranteed
to exist in every release or derivative distribution.

The Sage distribution is intended to build correctly no matter what
superset of the set of packages forming the minimal build
prerequisites is installed on the system.  If it does not, this is a
bug of the Sage distribution and should be reported and fixed on a
ticket.  Crucial part of a bug report is the configuration of the
system, in particular a list of installed packages and their versions.

Let us install a subset of these packages::

  root@39d693b2a75d:/sage# apt-get install libbz2-dev bzip2 yasm libz-dev
  Reading package lists... Done
  ...
  Setting up zlib1g-dev:amd64 (1:1.2.11.dfsg-0ubuntu2) ...
  root@39d693b2a75d:/sage#

  
Committing a container to disk
------------------------------

After terminating the container, we can create a new image
corresponding to its current state::

  root@39d693b2a75d:/sage# ^D
  [mkoeppe@sage worktree-ubuntu-latest]$ docker ps -a | head -n3
  CONTAINER ID        IMAGE                   COMMAND                   CREATED             STATUS
  39d693b2a75d        ubuntu:latest           "/bin/bash"               8 minutes ago       Exited (0) 6 seconds ago
  9f3398da43c2        ubuntu:latest           "/bin/bash"               8 minutes ago       Exited (0) 8 minutes ago
  [mkoeppe@sage worktree-ubuntu-latest]$ docker commit 39d693b2a75d ubuntu-latest-minimal-17
  sha256:4151c5ca4476660f6181cdb13923da8fe44082222b984c377fb4fd6cc05415c1

Here, ``39d693b2a75d`` was the container id (which appeared in the
shell prompts and in the output of ``docker ps``), and
``ubuntu-latest-minimal-17`` is an arbitrary symbolic name for the new
image.  The output of the command is the id of the new image.  We can
use either the symbolic name or the id to refer to the new image.

We can run the image and get a new container with the same state as
the one that we terminated.  Again we want to mount our worktree into
it; otherwise, because we did not make a copy, the new container will
have no access to the worktree::

  [mkoeppe@sage worktree-ubuntu-latest]$ docker run -it \
    --mount type=bind,source=$(pwd),target=/sage ubuntu-latest-minimal-17
  root@73987568712c:/# cd sage
  root@73987568712c:/sage# command -v gcc
  /usr/bin/gcc
  root@73987568712c:/sage# command -v yasm
  /usr/bin/yasm
  root@73987568712c:/sage# ^D
  [mkoeppe@sage worktree-ubuntu-latest]$

The image ``ubuntu-latest-minimal-17`` can be run in as many
containers as we want and can also be shared with other users or
developers so that they can run it in a container on their machine.
(See the Docker documentation on how to `share images on Docker Hub
<https://docs.docker.com/get-started/part3/>`_ or to `save images to a
tar archive
<https://docs.docker.com/engine/reference/commandline/save/>`_.)

This facilitates collaboration on fixing portability bugs of the Sage
distribution.  After reproducing a portability bug on a container,
several developers can work on fixing the bug using containers running
on their respective machines.


Generating Dockerfiles
----------------------

Sage also provides a script for generating a ``Dockerfile``, which is
a recipe for automatically building a new image::

  [mkoeppe@sage sage]$ build/bin/write-dockerfile.sh debian ":standard: :optional:" > Dockerfile

(The second argument is passed to ``sage -package list`` to find packages for the listed package types.)

.. this interface should be improved obviously. See #29146 - Refactor tox.ini and build/bin/write_dockerfile.sh

The ``Dockerfile`` instructs the command ``docker build`` to build a
new Docker image.  Let us take a quick look at the generated file;
this is slightly simplified::

  [mkoeppe@sage sage]$ cat Dockerfile 
  # Automatically generated by SAGE_ROOT/build/bin/write-dockerfile.sh
  # the :comments: separate the generated file into sections
  # to simplify writing scripts that customize this file
  ...

First, it instructs ``docker build`` to start from an existing base
image...::

  ...
  ARG BASE_IMAGE=ubuntu:latest
  FROM ${BASE_IMAGE}
  ...

Then, to install system packages...::

  ...
  RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -qqq --no-install-recommends --yes binutils make m4 perl python3 ... yasm libzmq3-dev libz-dev && apt-get clean

Then, to bootstrap and configure...::

  RUN mkdir -p /sage
  WORKDIR /sage
  ADD Makefile VERSION.txt README.md bootstrap configure.ac sage ./
  ADD src/doc/bootstrap src/doc/bootstrap
  ADD m4 ./m4
  ADD build ./build
  RUN ./bootstrap
  ADD src/bin src/bin
  ARG EXTRA_CONFIGURE_ARGS=""
  RUN ./configure --enable-build-as-root ${EXTRA_CONFIGURE_ARGS} || (cat config.log; exit 1)

Finally, to build and test...::

  ARG NUMPROC=8
  ENV MAKE="make -j${NUMPROC}"
  ARG USE_MAKEFLAGS="-k"
  RUN make ${USE_MAKEFLAGS} base-toolchain
  ARG TARGETS_PRE="sagelib-build-deps"
  RUN make ${USE_MAKEFLAGS} ${TARGETS_PRE}
  ADD src src
  ARG TARGETS="build ptest"
  RUN make ${USE_MAKEFLAGS} ${TARGETS}

You can customize the image build process by passing build arguments to the
command ``docker build``.  For example::

  [mkoeppe@sage sage]$ docker build . -f Dockerfile \
    --build-arg BASE_IMAGE=ubuntu:latest \
    --build-arg NUMPROC=4 \
    --build-arg EXTRA_CONFIGURE_ARGS="--with-python=/usr/bin/python3.42"

These arguments (and their default values) are defined using ``ARG``
commands in the ``Dockerfile``.

The above command will build Sage from scratch and will therefore take
quite long.  Let us instead just do a partial build, consisting of one
small package, by setting the arguments ``TARGETS_PRE`` and
``TARGETS``.  We use a silent build (``make V=0``)::

  [mkoeppe@sage sage]$ docker build . -f Dockerfile \
    --build-arg TARGETS_PRE=ratpoints \
    --build-arg TARGETS=ratpoints \
    --build-arg USE_MAKEFLAGS="V=0"
  Sending build context to Docker daemon    285MB
  Step 1/28 : ARG BASE_IMAGE=ubuntu:latest
  ...
  Step 2/28 : FROM ${BASE_IMAGE}
   ---> 549b9b86cb8d
  ...
  Step 25/28 : RUN make SAGE_SPKG="sage-spkg -y -o" ${USE_MAKEFLAGS} ${TARGETS_PRE}
  ...
  make[1]: Entering directory '/sage/build/make'
  sage-logger -p 'sage-spkg -y -o  ratpoints-2.1.3.p5' '/sage/logs/pkgs/ratpoints-2.1.3.p5.log'
  [ratpoints-2.1.3.p5] installing. Log file: /sage/logs/pkgs/ratpoints-2.1.3.p5.log
    [ratpoints-2.1.3.p5] successfully installed.
  make[1]: Leaving directory '/sage/build/make'

  real	0m18.886s
  user	0m1.779s
  sys	0m0.314s
  Sage build/upgrade complete!
  ...
  ---> 2d06689d39fa
  Successfully built 2d06689d39fa

We can now start a container using the image id shown in the last step::

  [mkoeppe@sage sage]$ docker run -it 2d06689d39fa bash
  root@fab59e09a641:/sage# ls -l logs/pkgs/
  total 236
  -rw-r--r-- 1 root root 231169 Mar 26 22:07 config.log
  -rw-r--r-- 1 root root   6025 Mar 26 22:27 ratpoints-2.1.3.p5.log
  root@fab59e09a641:/sage# ls -l local/lib/*rat*
  -rw-r--r-- 1 root root 177256 Mar 26 22:27 local/lib/libratpoints.a
  
You can customize the image build process further by editing the
``Dockerfile``.  For example, by default, the generated ``Dockerfile``
configures, builds, and tests Sage.  By deleting or commenting out the
commands for the latter, you can adjust the Dockerfile to stop after
the ``configure`` phase, for example.

``Dockerfile`` is the default filename for Dockerfiles.  You can
change it to any other name, but it is recommended to use
``Dockerfile`` as a prefix, such as ``Dockerfile-debian-standard``.
It should be placed within the tree rooted at the current directory
(``.``); if you want to put it elsewhere, you need to learn about
details of "Docker build contexts".

Note that in contrast to the workflow described in the above sections,
the ``Dockerfile`` **copies** a snapshot of your Sage worktree into
the build container, using ``ADD`` commands, instead of mounting the
directory into it.  This copying is subject to the exclusions in the
``.gitignore`` file (via a symbolic link from ``.dockerignore``).
Therefore, only the sources are copied, but not your configuration
(such as the file ``config.status``), nor the ``$SAGE_LOCAL`` tree,
nor any other build artefacts.

Because of this, you can build a Docker image using the generated
``Dockerfile`` from your main Sage development tree.  It does not have
to be distclean to start, and the build will not write into it at all.
Hence, you can continue editing and compiling your Sage development
tree even while Docker builds are running.


Debugging a portability bug using Docker
----------------------------------------

Let us do another partial build.  We choose a package that we suspect
might not work on all platforms, ``surf``, which was marked as
"experimental" in 2017::

  [mkoeppe@sage sage]$ docker build . -f Dockerfile \
    --build-arg BASE_IMAGE=ubuntu:latest \
    --build-arg NUMPROC=4 \
    --build-arg TARGETS_PRE=surf \
    --build-arg TARGETS=surf
  Sending build context to Docker daemon    285MB
  Step 1/28 : ARG BASE_IMAGE=ubuntu:latest
  Step 2/28 : FROM ${BASE_IMAGE}
   ---> 549b9b86cb8d
  ...
  Step 24/28 : ARG TARGETS_PRE="sagelib-build-deps"
   ---> Running in 17d0ddb5ad7b
  Removing intermediate container 17d0ddb5ad7b
   ---> 7b51411520c3
  Step 25/28 : RUN make SAGE_SPKG="sage-spkg -y -o" ${USE_MAKEFLAGS} ${TARGETS_PRE}
   ---> Running in 61833bea6a6d
  make -j4 build/make/Makefile --stop
  ...
  [surf-1.0.6-gcc6] Attempting to download package surf-1.0.6-gcc6.tar.gz from mirrors
  ...
  [surf-1.0.6-gcc6] http://mirrors.mit.edu/sage/spkg/upstream/surf/surf-1.0.6-gcc6.tar.gz
  ...
  [surf-1.0.6-gcc6] Setting up build directory for surf-1.0.6-gcc6
  ...
  [surf-1.0.6-gcc6] /usr/bin/ld: cannot find -lfl
  [surf-1.0.6-gcc6] collect2: error: ld returned 1 exit status
  [surf-1.0.6-gcc6] Makefile:504: recipe for target 'surf' failed
  [surf-1.0.6-gcc6] make[3]: *** [surf] Error 1
  ...
  [surf-1.0.6-gcc6] Error installing package surf-1.0.6-gcc6
  ...
  Makefile:2088: recipe for target '/sage/local/var/lib/sage/installed/surf-1.0.6-gcc6' failed
  make[1]: *** [/sage/local/var/lib/sage/installed/surf-1.0.6-gcc6] Error 1
  make[1]: Target 'surf' not remade because of errors.
  make[1]: Leaving directory '/sage/build/make'
  ...
  Error building Sage.

  The following package(s) may have failed to build (not necessarily
  during this run of 'make surf'):

  * package:         surf-1.0.6-gcc6
    last build time: Mar 26 22:07
    log file:        /sage/logs/pkgs/surf-1.0.6-gcc6.log
    build directory: /sage/local/var/tmp/sage/build/surf-1.0.6-gcc6

  ...
  Makefile:31: recipe for target 'surf' failed
  make: *** [surf] Error 1
  The command '/bin/sh -c make SAGE_SPKG="sage-spkg -y -o" ${USE_MAKEFLAGS} ${TARGETS_PRE}' returned a non-zero code: 2

Note that no image id is shown at the end; the build failed, and no
image is created.  However, the container in which the last step of
the build was attempted exists::

  [mkoeppe@sage sage]$ docker ps -a |head -n3
  CONTAINER ID        IMAGE                      COMMAND                   CREATED             STATUS
  61833bea6a6d        7b51411520c3               "/bin/sh -c 'make SA…"    9 minutes ago       Exited (2) 1 minute ago
  73987568712c        ubuntu-latest-minimal-17   "/bin/bash"               24 hours ago        Exited (0) 23 hours ago

We can copy the build directory from the container for inspection::

  [mkoeppe@sage sage]$ docker cp 61833bea6a6d:/sage/local/var/tmp/sage/build ubuntu-build
  [mkoeppe@sage sage]$ ls ubuntu-build/surf*/src
  AUTHORS         TODO            curve           misc
  COPYING         acinclude.m4    debug           missing
  ChangeLog       aclocal.m4      dither          mkinstalldirs
  INSTALL         background.pic  docs            mt
  Makefile        config.guess    draw            src
  Makefile.am     config.log      drawfunc        surf.1
  Makefile.global config.status   examples        surf.xpm
  Makefile.in     config.sub      gtkgui          yaccsrc
  NEWS            configure       image-formats
  README          configure.in    install-sh

Alternatively, we can use ``docker commit`` as explained earlier to
create an image from the container::

  [mkoeppe@sage sage]$ docker commit 61833bea6a6d
  sha256:003fbd511016fe305bd8494bb1747f0fbf4cb2c788b4e755e9099d9f2014a60d
  [mkoeppe@sage sage]$ docker run -it 003fbd511 bash
  root@2d9ac65f4572:/sage# (cd /sage/local/var/tmp/sage/build/surf* && /sage/sage --buildsh)

  Starting subshell with Sage environment variables set.  Don't forget
  to exit when you are done.
  ...
  Note: SAGE_ROOT=/sage
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ ls /usr/lib/libfl*
  /usr/lib/libflint-2.5.2.so  /usr/lib/libflint-2.5.2.so.13.5.2  /usr/lib/libflint.a  /usr/lib/libflint.so
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ apt-get update && apt-get install apt-file    
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ apt-file update
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ apt-file search "/usr/lib/libfl.a"
  flex-old: /usr/lib/libfl.a
  freebsd-buildutils: /usr/lib/libfl.a
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ apt-get install flex-old
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ ./spkg-install 
  checking for a BSD-compatible install... /usr/bin/install -c
  checking whether build environment is sane... yes
  ...
    /usr/bin/install -c  surf /sage/local/bin/surf
   /usr/bin/install -c -m 644 ./surf.1 /sage/local/share/man/man1/surf.1
  ...
  make[1]: Leaving directory '/sage/local/var/tmp/sage/build/surf-1.0.6-gcc6/src'
  (sage-buildsh) root@2d9ac65f4572:surf-1.0.6-gcc6$ exit
  root@2d9ac65f4572:/sage# exit
  [mkoeppe@sage sage]$

A standard case of bitrot.
  

Automatic Docker-based build testing using tox
----------------------------------------------

`tox <https://tox.readthedocs.io/en/latest/>`_ is a Python package that
is widely used for automating tests of Python projects.

Install ``tox`` for use with your system Python, for example using::

  [mkoeppe@sage sage]$ pip install --user tox

A tox "environment" is identified by a symbolic name composed of
several `Tox "factors"
<https://tox.readthedocs.io/en/latest/config.html#complex-factor-conditions>`_,
which are defined in the file ``$SAGE_ROOT/tox.ini``.

The **technology** factor describes how the environment is run:

- ``docker`` builds a Docker image as described above.

- ``local`` runs testing on the host OS instead.  We explain this
  technology in a later section.

The next two factors determine the host system configuration: The
**system factor** describes a base operating system image.

- Examples are ``ubuntu-focal``, ``debian-buster``,
  ``archlinux-latest``, ``fedora-30``, ``slackware-14.2``,
  ``centos-7-i386``, and ``ubuntu-bionic-arm64``.

- See ``$SAGE_ROOT/tox.ini`` for a complete list, and to which images
  on Docker hub they correspond.

The **packages factor** describes a list of system packages to be
installed on the system before building Sage:

- ``minimal`` installs the system packages known to Sage to provide
  minimal prerequisites for bootstrapping and building the Sage
  distribution.  This corresponds to the packages ``_bootstrap`` and
  ``_prereq``.

- ``standard`` additionally installs all known system packages that
  are equivalent to standard packages of the Sage distribution, for
  which the mechanism ``spkg-configure.m4`` is implemented.
  This corresponds to the packages listed by::

    [mkoeppe@sage sage]$ sage --package list --has-file=spkg-configure.m4 :standard:

- ``maximal`` does the same for all standard and optional packages.
  This corresponds to the packages listed by::

    [mkoeppe@sage sage]$ sage --package list :standard: :optional:

The factors are connected by a hyphen to name a system configuration,
such as ``debian-buster-standard`` and ``centos-7-i386-minimal``.

Finally, the **configuration** factor (which is allowed to be empty)
controls how the ``configure`` script is run.

The factors are connected by a hyphen to name a tox environment.  (The
order of the factors does not matter; however, for consistency and
because the ordered name is used for caching purposes, we recommend to
use the factors in the listed order.)

To run an environment::

  [mkoeppe@sage sage]$ tox -e docker-slackware-14.2-minimal
  [mkoeppe@sage sage]$ tox -e docker-ubuntu-bionic-standard
  
Arbitrary extra arguments to ``docker build`` can be supplied through
the environment variable ``EXTRA_DOCKER_BUILD_ARGS``.  For example,
for a non-silent build (``make V=1``), use::
  
  [mkoeppe@sage sage]$ EXTRA_DOCKER_BUILD_ARGS="--build-arg USE_MAKEFLAGS=\"V=1\"" \
    tox -e docker-ubuntu-bionic-standard

By default, tox uses ``TARGETS_PRE=sagelib-build-deps`` and
``TARGETS=build``, leading to a complete build of Sage without the
documentation.  If you pass positional arguments to tox (separated
from tox options by ``--``), then both ``TARGETS_PRE`` and ``TARGETS``
are set to these arguments.  In this way, you can build some specific
packages instead of all of Sage, for example::

  [mkoeppe@sage sage]$ tox -e docker-centos-8-standard -- ratpoints

If the build succeeds, this will create a new image named
``sage-docker-centos-8-standard-with-targets:9.1.beta9-431-gca4b5b2f33-dirty``,
where

- the image name is derived from the tox environment name and the
  suffix ``with-targets`` expresses that the ``make`` targets given in
  ``TARGETS`` have been built;

- the tag name describes the git revision of the source tree as per
  ``git describe --dirty``.

You can ask for tox to create named intermediate images as well.  For
example, to create the images corresponding to the state of the OS
after installing all system packages (``with-system-packages``) and
the one just after running the ``configure`` script (``configured``)::

  [mkoeppe@sage sage]$ DOCKER_TARGETS="with-system-packages configured with-targets" \
    tox -e docker-centos-8-standard -- ratpoints
  ...
  Sending build context to Docker daemon ...
  Step 1/109 : ARG BASE_IMAGE=fedora:latest
  Step 2/109 : FROM ${BASE_IMAGE} as with-system-packages
  ...
  Step 109/109 : RUN yum install -y zlib-devel || echo "(ignoring error)"
  ...
  Successfully built 4bb14c3d5646
  Successfully tagged sage-docker-centos-8-standard-with-system-packages:9.1.beta9-435-g861ba33bbc-dirty
  Sending build context to Docker daemon ...
  ...
  Successfully tagged sage-docker-centos-8-standard-configured:9.1.beta9-435-g861ba33bbc-dirty
  ...
  Sending build context to Docker daemon ...
  ...
  Successfully tagged sage-docker-centos-8-standard-with-targets:9.1.beta9-435-g861ba33bbc-dirty

Let's verify that the images are available::

  (base) egret:~/s/sage/sage-rebasing/worktree-algebraic-2018-spring (mkoeppe *$%>)$ docker images | head
  REPOSITORY                                                TAG                               IMAGE ID
  sage-docker-centos-8-standard-with-targets                9.1.beta9-435-g861ba33bbc-dirty   7ecfa86fceab
  sage-docker-centos-8-standard-configured                  9.1.beta9-435-g861ba33bbc-dirty   4314929e2b4c
  sage-docker-centos-8-standard-with-system-packages        9.1.beta9-435-g861ba33bbc-dirty   4bb14c3d5646
  ...


Automatic build testing on the host OS using tox -e local-direct
----------------------------------------------------------------

The ``local`` technology runs testing on the host OS instead.

In contrast to the ``docker`` technology, it does not make a copy of
the source tree.  It is most straightforward to run it from a
separate, distclean git worktree.

Let us try a first variant of the ``local`` technology, the tox
environment called ``local-direct``.  Because all builds with tox
begin by bootstrapping the source tree, you will need autotools and
other prerequisites installed in your system.  See
``build/pkgs/_bootstrap/distros/*.txt`` for a list of system packages that
provide these prerequisites.

We start by creating a fresh (distclean) git worktree.

  [mkoeppe@sage sage] git worktree add worktree-local
  [mkoeppe@sage sage] cd worktree-local
  [mkoeppe@sage worktree-local] ls
  COPYING.txt ... Makefile ... configure.ac ... src tox.ini

Again we build only a small package.  Build targets can be passed as
positional arguments (separated from tox options by ``--``)::

  [mkoeppe@sage worktree-local] tox -e local-direct -- ratpoints
  local-direct create: /Users/mkoeppe/.../worktree-local/.tox/local-direct
  local-direct run-test-pre: PYTHONHASHSEED='2211987514'
  ...
  src/doc/bootstrap:48: installing src/doc/en/installation/debian.txt...
  bootstrap:69: installing 'config/config.rpath'
  configure.ac:328: installing 'config/compile'
  configure.ac:113: installing 'config/config.guess'
  ...
  checking for a BSD-compatible install... /usr/bin/install -c
  checking whether build environment is sane... yes
  ...
  sage-logger -p 'sage-spkg -y -o  ratpoints-2.1.3.p5' '.../worktree-local/logs/pkgs/ratpoints-2.1.3.p5.log'
  [ratpoints-2.1.3.p5] installing. Log file: .../worktree-local/logs/pkgs/ratpoints-2.1.3.p5.log
    [ratpoints-2.1.3.p5] successfully installed.
  ...
    local-direct: commands succeeded
    congratulations :)

Let's investigate what happened here::

  [mkoeppe@sage worktree-local]$ ls -la
  total 2576
  drwxr-xr-x  35 mkoeppe  staff    1120 Mar 26 22:20 .
  drwxr-xr-x  63 mkoeppe  staff    2016 Mar 27 09:35 ..
  ...
  lrwxr-xr-x   1 mkoeppe  staff      10 Mar 26 20:34 .dockerignore -> .gitignore
  -rw-r--r--   1 mkoeppe  staff      74 Mar 26 20:34 .git
  ...
  -rw-r--r--   1 mkoeppe  staff    1212 Mar 26 20:41 .gitignore
  ...
  drwxr-xr-x   7 mkoeppe  staff     224 Mar 26 22:11 .tox
  ...
  -rw-r--r--   1 mkoeppe  staff    7542 Mar 26 20:41 Makefile
  ...
  lrwxr-xr-x   1 mkoeppe  staff     114 Mar 26 20:45 config.log -> .tox/local-direct/log/config.log
  -rwxr-xr-x   1 mkoeppe  staff   90411 Mar 26 20:46 config.status
  -rwxr-xr-x   1 mkoeppe  staff  887180 Mar 26 20:45 configure
  -rw-r--r--   1 mkoeppe  staff   17070 Mar 26 20:41 configure.ac
  ...
  lrwxr-xr-x   1 mkoeppe  staff     103 Mar 26 20:45 logs -> .tox/local-direct/log
  drwxr-xr-x  24 mkoeppe  staff     768 Mar 26 20:45 m4
  lrwxr-xr-x   1 mkoeppe  staff     105 Mar 26 20:45 prefix -> .tox/local-direct/local
  -rwxr-xr-x   1 mkoeppe  staff    4868 Mar 26 20:34 sage
  drwxr-xr-x  16 mkoeppe  staff     512 Mar 26 20:46 src
  -rw-r--r--   1 mkoeppe  staff   13478 Mar 26 20:41 tox.ini
  drwxr-xr-x   4 mkoeppe  staff     128 Mar 26 20:46 upstream

There is no ``local`` subdirectory.  This is part of a strategy to
keep the source tree clean to the extent possible. In particular:

- ``tox`` configured the build to use a separate ``$SAGE_LOCAL``
  hierarchy in a directory under the tox environment directory
  ``.tox/local-direct``.  It created a symbolic link ``prefix`` that
  points there, for convenience::

    [mkoeppe@sage worktree-local]$ ls -l prefix/lib/*rat*
    -rw-r--r--  1 mkoeppe  staff  165968 Mar 26 20:46 prefix/lib/libratpoints.a

- Likewise, it created a separate ``logs`` directory, again under the
  tox environment directory, and a symbolic link.

This makes it possible for advanced users to test several ``local``
tox environments (such as ``local-direct``) out of one worktree.  However, because a
build still writes configuration scripts and build artefacts (such as
``config.status``) into the worktree, only one ``local`` build can run
at a time in a given worktree.

The tox environment directory will be reused for the next ``tox`` run,
which will therefore do an incremental build.  To start a fresh build,
you can use the ``-r`` option.

Automatic build testing on the host OS with best-effort isolation using tox -e local
------------------------------------------------------------------------------------

``tox -e local`` (without ``-direct``) attempts a best-effort
isolation from the user's environment as follows:

- All environment variables are set to standard values; with the
  exception of ``MAKE`` and ``EXTRA_CONFIGURE_ARGS``.  In particular,
  ``PATH`` is set to just ``/usr/bin:/bin:/usr/sbin:/sbin``; it does
  not include ``/usr/local/bin``.


Note, however, that various packages have build scripts that use
``/usr/local`` or other popular file system locations such as
``/opt/sfw/``.  Therefore, the isolation is not complete.  Using
``/usr/local`` is considered standard behavior.  On the other hand, we
consider a package build script that inspects other file system
locations to be a bug of the Sage distribution, which should be
reported and fixed on a ticket.


Automatic build testing on macOS with a best-effort isolated installation of Homebrew
-------------------------------------------------------------------------------------

XCode on macOS does not provide the prerequisites for bootstrapping
the Sage distribution.  A good way to install them is using the
Homebrew package manager.

In fact, Sage provides a tox environment that automatically installs
an isolated copy of Homebrew with all prerequisites for bootstrapping::

  [mkoeppe@sage worktree-local]$ tox -e local-homebrew-macos-minimal -- lrslib
  local-homebrew-macos-minimal create: .../worktree-local/.tox/local-homebrew-macos-minimal
  local-homebrew-macos-minimal run-test-pre: PYTHONHASHSEED='4246149402'
  ...
  Initialized empty Git repository in .../worktree-local/.tox/local-homebrew-macos-minimal/homebrew/.git/
  ...
  Tapped 2 commands and 4942 formulae (5,205 files, 310.7MB).
  ==> Downloading https://ftp.gnu.org/gnu/gettext/gettext-0.20.1.tar.xz
  ...
  ==> Pouring autoconf-2.69.catalina.bottle.4.tar.gz
  ...
  ==> Pouring pkg-config-0.29.2.catalina.bottle.1.tar.gz
    .../worktree-local/.tox/local-homebrew-macos-minimal/homebrew/Cellar/pkg-config/0.29.2: 11 files, 623.4KB
  ==> Caveats
  ==> gettext
  gettext is keg-only, which means it was not symlinked into .../worktree-local/.tox/local-homebrew-macos-minimal/homebrew,
  because macOS provides the BSD gettext library & some software gets confused if both are in the library path.

  If you need to have gettext first in your PATH run:
    echo 'export PATH=".../worktree-local/.tox/local-homebrew-macos-minimal/homebrew/opt/gettext/bin:$PATH"' >> ~/.bash_profile

  For compilers to find gettext you may need to set:
    export LDFLAGS="-L.../worktree-local/.tox/local-homebrew-macos-minimal/homebrew/opt/gettext/lib"
    export CPPFLAGS="-I.../worktree-local/.tox/local-homebrew-macos-minimal/homebrew/opt/gettext/include"
  ...
  local-homebrew-macos-minimal run-test: commands[0] | bash -c 'export PATH=.../worktree-local/.tox/local-homebrew-macos-minimal/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin && . .homebrew-build-env && ./bootstrap && ./configure --prefix=.../worktree-local/.tox/local-homebrew-macos-minimal/local    && make -k V=0 ... lrslib'
  ...
  bootstrap:69: installing 'config/config.rpath'
  ...
  checking for a BSD-compatible install... /usr/bin/install -c
  checking whether build environment is sane... yes
  ...
  configure: notice: the following SPKGs did not find equivalent system packages: arb cbc cliquer ... tachyon xz yasm zeromq
  checking for the package system in use... homebrew
  configure: hint: installing the following system packages is recommended and may avoid building some of the above SPKGs from source:
  configure:   $ brew install cmake gcc gsl mpfi ninja openblas gpatch r readline xz yasm zeromq
  ...
  sage-logger -p 'sage-spkg -y -o  lrslib-062+autotools-2017-03-03.p1' '.../worktree-local/logs/pkgs/lrslib-062+autotools-2017-03-03.p1.log'
  [lrslib-062+autotools-2017-03-03.p1] installing. Log file: .../worktree-local/logs/pkgs/lrslib-062+autotools-2017-03-03.p1.log
    [lrslib-062+autotools-2017-03-03.p1] successfully installed.
  ...
    local-homebrew-macos-minimal: commands succeeded
    congratulations :)
  
The tox environment uses the subdirectory ``homebrew`` of the
environment directory ``.tox/local-homebrew-macos-minimal`` as the
Homebrew prefix.  This installation does not interact in any way with
a Homebrew installation in ``/usr/local`` that you may have.

The test script sets the ``PATH`` to the ``bin`` directory of the
Homebrew prefix, followed by ``/usr/bin:/bin:/usr/sbin:/sbin``.  It
then uses the script ``$SAGE_ROOT/.homebrew-build-env`` to set
environment variables so that Sage's build scripts will find
"keg-only" packages such as ``gettext``.

The ``local-homebrew-macos-minimal`` environment does not install
Homebrew's ``python3`` package.  It uses XCode's ``/usr/bin/python3``
as system python.  However, because various packages are missing
that Sage considers as dependencies, Sage builds its own copy of
these packages and of ``python3``.

The ``local-homebrew-macos-standard`` environment additionally
installs (in its separate isolated copy of Homebrew) all Homebrew
packages known to Sage for which the ``spkg-configure.m4`` mechanism
is implemented; this is similar to the ``docker-standard`` tox
environments described earlier.  In particular it installs and uses
Homebrew's ``python3`` package.

By using configuration factors, more variants can be tested.
The ``local-homebrew-macos-standard-python3_xcode`` environment
installs the same packages, but uses XCode's ``/usr/bin/python3``.

The ``local-homebrew-macos-standard-python3_pythonorg`` expects an
installation of Python 3.7 in
``/Library/Frameworks/Python.framework``; this is where the binary
packages provided by python.org install themselves.


Automatic build testing with a best-effort isolated installation of Conda
-------------------------------------------------------------------------

Sage provides environments ``local-conda-forge-standard`` and
``local-conda-forge-minimal`` that create isolated installations of
Miniconda in the subdirectory ``conda`` of the environment directory.
They do not interact in any way with other installations of Anaconda
or Miniconda that you may have on your system.

The environments use the conda-forge channel and use the ``python``
package and the compilers from this channel.


Options for build testing with the local technology
---------------------------------------------------

The environments using the ``local`` technology can be customized
by setting environment variables.

 - If ``SKIP_SYSTEM_PKG_INSTALL`` is set to ``1`` (or ``yes``),
   then all steps of installing system packages are skipped in this run.
   When reusing a previously created tox environment, this option can
   save time and also give developers more control for experiments
   with system packages.

 - If ``SKIP_BOOTSTRAP`` is set to ``1`` (or ``yes``), then the
   bootstrapping phase is skipped.  When reusing a previously created
   tox environment, this option can save time.

 - If ``SKIP_CONFIGURE`` is set to ``1`` (or ``yes``), then the
   ``configure`` script is not run explicitly.  When reusing a
   previously created tox environment, this option can save time.
   (The ``Makefile`` may still rerun configuration using
   ``config.status --recheck``.)

The ``local`` technology also defines a special target ``bash``:
Instead of building anything with ``make``, it just starts an
interactive shell.  For example, in combination with the above
options::

  [mkoeppe@sage worktree-local]$ SKIP_SYSTEM_PKG_INSTALL=yes SKIP_BOOTSTRAP=1 SKIP_CONFIGURE=1 tox -e local-homebrew-macos-minimal -- bash


Automatic parallel tox runs on GitHub Actions
---------------------------------------------

The Sage source tree includes a default configuration for GitHub
Actions that runs tox on a multitude of platforms on every pull
request and on every push of a tag (but not of a branch) to a
repository for which GitHub Actions are enabled.

This is defined in the file ``$SAGE_ROOT/.github/workflows/tox.yml``.

An additional GitHub Actions workflow for testing on Cygwin, not based
on tox, is defined in the file
``$SAGE_ROOT/.github/workflows/ci-cygwin.yml``.

GitHub Actions runs these build jobs on 2-core machines with 7 GB of
RAM memory and 14 GB of SSD disk space, cf.
`here <https://help.github.com/en/actions/reference/virtual-environments-for-github-hosted-runners#supported-runners-and-hardware-resources>`_,
and has a time limit of 6h per job. This is just barely enough for a
typical ``minimal`` build followed by ``make ptest`` to succeed; and
plenty of time for a typical ``standard`` build to succeed.

Build logs become available as "artifacts" when all jobs of the
workflow have finished.  Each job generates one tarball.
"Annotations" highlight certain top-level errors or warnings issued
during the build.

The following procedure triggers a run of tests with the default set of
system configurations.

- Push your changes to trac.
- Go to the `Actions page on the GitHub mirror <https://github.com/sagemath/sagetrac-mirror/actions>`_ and select the workflow you would like to run.
- Click on "Run workflow" above the list of workflow runs and select the branch where the workflow will run.

For more information, see the `GitHub documentation <https://docs.github.com/en/free-pro-team@latest/actions/managing-workflow-runs/manually-running-a-workflow>`_.

Alternatively, you can create and push a custom tag in order to trigger a run of tests as follows. 
Let's assume that ``github`` is the name of
the remote corresponding to your GitHub fork of the Sage repository::

  $ git remote -v | grep /my-github
  my-github      https://github.com/mkoeppe/sage.git (fetch)
  my-github      https://github.com/mkoeppe/sage.git (push)

- Create a ("lightweight", not "annotated") tag with an arbitrary
  name, say ``ci`` (for "Continuous Integration")::

    git tag -f ci

- Then push the tag to your GitHub repository::

    git push -f my-github ci

(In both commands, the "force" option (``-f``) allows overwriting a
previous tag of that name.)

For testing branches against a custom set of system configurations
during development, the following procedure seems to work well.  It
avoids changing the CI configuration on your development branch:

- Create a branch from a recent beta release that contains the default
  GitHub Actions configuration; name it ``TESTER``, say.

- Edit ``$SAGE_ROOT/.github/workflows/tox.yml`` to include the system
  config you wish to test.

- Commit and push the branch to your GitHub fork of sage.

- Push your development branch to your GitHub repository and create a
  pull request against the ``TESTER`` branch. This will trigger the
  GitHub Actions workflow.

You will find a workflow status page in the "Actions" tab of your
repository.

Here is how to read it.  Each of the items in the left pane represents
a full build of Sage on a particular system configuration.  A test
item in the left pane is marked with a green checkmark in the left
pane if ``make build doc-html`` finished without error.  (It also runs
package testsuites and the Sage doctests but failures in these are not
in reflected in the left pane; see below.)

The right pane ("Artifacts") offers archives of the logs for download.

Scrolling down in the right pane shows "Annotations":

* Red "check failure" annotations appear for each log file that
  contains a build error. For example, you might see::

    docker (fedora-28, standard)
    artifacts/logs-commit-8ca1c2df8f1fb4c6d54b44b34b4d8320ebecb164-tox-docker-fedora-28-standard/logs/pkgs/sagetex-3.4.log#L1
    ==== ERROR IN LOG FILE artifacts/logs-commit-8ca1c2df8f1fb4c6d54b44b34b4d8320ebecb164-tox-docker-fedora-28-standard/logs/pkgs/sagetex-3.4.log ====

* Yellow "check warning" annotations. There are 2 types of these:

  a) Package testsuite or Sage doctest failures, like the following::

       docker (fedora-30, standard)
       artifacts/logs-commit-8ca1c2df8f1fb4c6d54b44b34b4d8320ebecb164-tox-docker-fedora-30-standard/logs/ptest.log#L1
       ==== TESTSUITE FAILURE IN LOG FILE artifacts/logs-commit-8ca1c2df8f1fb4c6d54b44b34b4d8320ebecb164-tox-docker-fedora-30-standard/logs/ptest.log ====

  b) Notices from ./configure about not finding equivalent system
     packages, like the following::

       docker (fedora-31, standard)
       artifacts/logs-commit-8ca1c2df8f1fb4c6d54b44b34b4d8320ebecb164-tox-docker-fedora-31-standard/config.log#L1
       configure: notice: the following SPKGs did not find equivalent system packages: arb cbc cddlib cmake eclib ecm fflas_ffpack flint flintqs fplll givaro gp

Clicking on the annotations does not take you to a very useful
place. To view details, click on one of the items in the pane. This
changes the right pane to a log viewer.

The ``docker`` workflows automatically push images to
``docker.pkg.github.com``.  You find them in the Packages tab of your
GitHub repository.

In order to pull them for use on your computer, you need to first
generate a Personal Access Token providing the ``read:packages`` scope
as follows.  Visit https://github.com/settings/tokens/new (this may
prompt you for your GitHub password).  As "Note", type "Access
docker.pkg.github.com"; then in "Select scopes", select the checkbox
for ``read:packages``.  Finally, push the "Generate token" button at
the bottom.  This will lead to a page showing your token, such as
``de1ec7ab1ec0ffee5ca1dedbaff1ed0ddba11``.  Copy this token and paste
it to the command line::

  $ echo de1ec7ab1ec0ffee5ca1dedbaff1ed0ddba11 | docker login docker.pkg.github.com --username YOUR-GITHUB-USERNAME

where you replace the token by your token, of course, and
``YOUR-GITHUB-USERNAME`` by your GitHub username.

Now you can pull the image and run it::

  $ docker pull docker.pkg.github.com/YOUR-GITHUB-USERNAME/sage/sage-docker-fedora-31-standard-configured:f4bd671
  $ docker run -it docker.pkg.github.com/YOUR-GITHUB-USERNAME/sage/sage-docker-fedora-31-standard-configured:f4bd671 bash
