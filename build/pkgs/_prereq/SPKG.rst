_prereq: Represents system packages required for installing SageMath from source
================================================================================

Description
-----------

This dummy package represents the minimal requirements (system packages)
for installing SageMath from source.

In addition to standard :wikipedia:`POSIX <POSIX>` utilities
and the :wikipedia:`bash <Bash_(Unix_shell)>` shell,
the following standard command-line development tools must be installed on your
computer:

- **make**: GNU make, version 3.80 or later. Version 3.82 or later is recommended.
- **m4**: GNU m4 1.4.2 or later (non-GNU or older versions might also work).
- **perl**: version 5.8.0 or later.
- **ar** and **ranlib**: can be obtained as part of GNU binutils.
- **tar**: GNU tar version 1.17 or later, or BSD tar (as provided on macOS).
- **python**: Python 3.4 or later, or Python 2.7.
  (This range of versions is a minimal requirement for internal purposes of the SageMath
  build system, which is referred to as ``sage-bootstrap-python``.)

Other versions of these may work, but they are untested.

On macOS, suitable versions of all of these tools are provided
by the Xcode Command Line Tools.  To install them, open a terminal
window and run ``xcode-select --install``; then click "Install" in the
pop-up window.  If the Xcode Command Line Tools are already installed,
you may want to check if they need to be updated by typing
``softwareupdate -l``.

On Linux, ``ar`` and ``ranlib`` are in the `binutils
<https://www.gnu.org/software/binutils/>`_ package.  The other
programs are usually located in packages with their respective names.

On Redhat-derived systems not all perl components are installed by
default and you might have to install the ``perl-ExtUtils-MakeMaker``
package.

To check if you have the above prerequisites installed, for example ``perl``,
type::

    $ command -v perl

or::

    $ which perl

on the command line. If it gives an error (or returns nothing), then
either ``perl`` is not installed, or it is installed but not in your
:wikipedia:`PATH <PATH_%28variable%29>`.
