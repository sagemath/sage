.. _sec-troubles:

Troubleshooting
===============

If you have any problems installing or running Sage, please take a look at the
`Release Tour <https://github.com/sagemath/sage/releases>`_ corresponding to
the version that you are installing. It may offer version-specific installation
help that has become available after the release was made and is therefore not
covered by this guide.

You may ask for help in
`sage-support <https://groups.google.com/forum/#!forum/sage-support>`_ or
`sage-develop <https://groups.google.com/forum/#!forum/sage-develop>`_.
When you ask for help, report on your system and be prepared to send log
files generated during the build process. Note that ``./configure`` generates
``config.log`` file in the Sage root directory. Each component of Sage, called Sage
package or SPKG, leaves a build log in ``logs/pkgs/``. Hence if ``make``
fails, you can browse these to find error messages. Check the most recent log
files by runnning::

    grep -li "^Error" logs/pkgs/*

from the Sage root directory to find relevant log files. Send the file
``config.log`` as well as the log files of the failed packages in their
entirety so that we can investigate problems on your platform.
