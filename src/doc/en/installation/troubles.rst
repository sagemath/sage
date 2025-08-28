.. _sec-troubles:

Troubleshooting
===============

If no binary version is available for your system, you can fallback to
the :ref:`sec-installation-from-sources` or use one of the alternatives
proposed at the end of :ref:`installation-guide`.

If you have any problems building or running Sage, please take a look
at the `release tour
<https://github.com/sagemath/sage/releases>`_ corresponding to the version
that you are installing.  It may offer version-specific installation
help that has become available after the release was made and is
therefore not covered by this manual.

Also please do not hesitate to ask for help in the `SageMath forum
<https://ask.sagemath.org/questions/>`_ or the sage-support mailing
list at https://groups.google.com/forum/#!forum/sage-support.

Also note the following. Each separate component of Sage is
contained in an SPKG; these are stored in ``build/pkgs/``. As each one
is built, a build log is stored in ``logs/pkgs/``, so you can browse these
to find error messages. If an SPKG fails to build, the whole build
process will stop soon after, so check the most recent log files
first, or run::

       grep -li "^Error" logs/pkgs/*

from the top-level Sage directory to find log files with error
messages in them.  Send the file ``config.log`` as well as the
log file(s) of the packages that have failed to build
in their entirety to the sage-support mailing list
at https://groups.google.com/group/sage-support; probably someone
there will have some helpful suggestions.
