.. _sec-installation-from-binaries:

Install from Pre-built Binaries
===============================

Installation from a pre-built binary tarball is an easy and
fast way to install Sage. Note that on GNU/Linux a preferred
way is to use your package manager (e.g. apt, pacman, yum).

In all cases, we assume that you have a computer with at least
4 GB of free disk space.

Download Guide
--------------

Not sure what to download? Just follow these steps.

- Determine your operating system (Windows, Linux or macOS).

- According to your operating system, go to the appropriate Download
  section of the `SageMath website <http://www.sagemath.org/>`_.

- Choose a download server (aka mirror) that is close to your location.

- Download the binary that is appropriate to your system. Depending on your
  operating system you might need additional information such as your CPU
  type (e.g. 64 bits or 32 bits) and your operating system version. If you
  use macOS you will have the choice between a tarball (whose names ends with
  ``tar.bz2``) and two kinds of mountable disk images (whose names end with
  ``app.dmg`` and simply ``.dmg``). Except for Windows, the naming scheme of
  the files is always ``sage-VERSION-OS-CPU.EXTENSION`` where ``EXTENSION``
  can be ``tar.gz``, ``tar.bz2``, ``dmg`` or ``app.dmg``.
 
- Then choose the appropriate section below corresponding to your situation.

Linux
-----

Make sure that you have an SSL library installed
(OpenSSL recommended).

It is highly recommended that you have LaTeX installed. If you want
to view animations, you should install either ImageMagick or ffmpeg.
ImageMagick or dvipng is also used for displaying some LaTeX output
in the notebooks.

Choose an appropriate directory where to install Sage. If you have
administrator rights on your computer a good choice is ``/opt``
otherwise it can be anywhere in your home directory. Avoid spaces and
Unicode characters in the path name.

Next, download the latest binary tarball available
(see "Download Guide" above). The tarball name should end
with ``.tar.gz`` or ``.tar.bz2``. If you want to use the ``.dmg``
or ``.app.dmg`` for macOS switch to the next section.

Unpack the tarball where you intend to install Sage. This is done
from the command line using the ``tar`` program. Next, to launch
Sage, go to the ``SageMath`` directory and run the program that
is called ``sage`` (via ``./sage`` on the command line).

The first time you run Sage, you will see a message like

.. CODE-BLOCK:: text

   Rewriting paths for your new installation directory
   ===================================================

   This might take a few minutes but only has to be done once.

   patching ...  (long list of files)

At this point, you can no longer move your Sage installation and
expect Sage to function.

Once you are able to launch Sage you might want to create a shortcut
so that ``sage`` just works from the command line. To do so simply use
the ``ln`` program from the command line:

.. CODE-BLOCK:: shell-session

    $ sudo ln -s /path/to/SageMath/sage /usr/local/bin/sage

where ``/path/to/SageMath/sage`` is the actual path to your SageMath
installation.

macOS
-----

On macOS there are two possible binaries for each version. They can be
recognized by their suffixes, but their actual contents are identical.

- ``tar.bz2``: a binary tarball
- ``dmg``: a compressed image of the binary

This section explains how to install from ``dmg``. For the
installation of the binary tarball ``tar.bz2`` just follow the steps
of the Linux installation.

After downloading the file, double click on the dmg file to mount it,
which will take some time.  Then drag the folder ``SageMath`` that
just appeared to ``/Applications/``.  You might want to have shortcuts
so that ``sage`` in the console simply works out of the box.  For that
purpose, follows the steps at the end of the section "Linux".

Alternative macOS binaries are available `here
<https://github.com/3-manifolds/Sage_macOS/releases/>`_.  These
have been signed and notarized, eliminating various errors caused by
Apple's gatekeeper antimalware protections.


Microsoft Windows (Cygwin)
--------------------------

SageMath on Windows requires a 64-bit Windows (which is likely to be the case
on a modern computer). If you happen to have a 32-bit Windows, you can consider
the alternatives mentioned at the end of :ref:`installation-guide`.

To install SageMath on Windows, just download the installer (see the above
"Download Guide" section) and run it.
