.. _sec-launching:

Launching SageMath
==================

Now we assume that you installed SageMath properly on your system. This
section quickly explains how to start the Sage console and the Jupyter
Notebook from the command line.

If you did install the Windows version or the macOS application you
should have icons available on your desktops or launching menus. Otherwise
you are strongly advised to create shortcuts for Sage as indicated at the end
of the "Linux" Section in :ref:`sec-installation-from-binaries`. Assuming
that you have this shortcut, running

.. CODE-BLOCK:: bash

    sage

in a console starts a Sage session.  To quit the session enter ``quit`` and
then press ``<Enter>``.

To start a Jupyter Notebook instead of a Sage console, run the command

.. CODE-BLOCK:: bash

    sage -n jupyter

instead of just ``sage``. To quit the Jupyter Notebook press ``<Ctrl> + <c>``
twice in the console where you launched the command.

Using a Jupyter Notebook remotely
---------------------------------

If Sage is installed on a remote machine to which you have ``ssh`` access, you
can launch a Jupyter Notebook using a command such as

.. CODE-BLOCK:: bash

    ssh -L localhost:8888:localhost:8888 -t USER@REMOTE sage -n jupyter --no-browser --port=8888

where ``USER@REMOTE`` needs to be replaced by the login details to the remote
machine. This uses local port forwarding to connect your local machine to the
remote one. The command will print a URL to the console which you can copy and
paste in a web browser.

Note that this assumes that a firewall which might be present between server
and client allows connections on port 8888. See details on port forwarding on
the internet, e.g. https://www.ssh.com/ssh/tunneling/example.

------------------------------------------------------------------------

For further reading you can have a look at the other documents in the
SageMath documentation at http://doc.sagemath.org/.

