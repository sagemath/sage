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

.. _sec-launching-wsl-post-installation:

WSL Post-installation steps
---------------------------

If you've installed Sage Math from source on WSL, there are a couple of extra steps you can do to make your life easier:

Create a notebook launch script
"""""""""""""""""""""""""""""""

If you plan to use JupyterLab, install that first.

Now create a script called ``~/sage_nb.sh`` containing the following lines, and fill in the correct paths for your desired starting directory and ``SAGE_ROOT``


.. CODE-BLOCK:: bash

    #!/bin/bash
    # Switch to desired windows directory
    cd /mnt/c/path/to/desired/starting/directory
    # Start the Jupyter notebook
    SAGE_ROOT/sage --notebook
    # Alternatively you can run JupyterLab - delete the line above, and uncomment the line below
    #SAGE_ROOT/sage --notebook jupyterlab

Make it executable:

.. CODE-BLOCK:: bash

    chmod ug+x ~/sage_nb.sh

Run it to test:

.. CODE-BLOCK:: bash

    cd ~
    ./sage_nb.sh

The Jupyter(Lab) server should start in the terminal window, and you windows browser should open a page showing the Jupyter or JupyterLab starting page, at the directory you specified.

Create a shortcut
"""""""""""""""""

This is a final nicety that lets you start the Jupyter or JupyterLab server in one click:

* Open Windows explorer, and type ``%APPDATA%\Microsoft\Windows\Start Menu\Programs`` in the address bar and press enter. This is the folder that contains you start menu shortcuts. If you want the sage shortcut somewhere else (like your desktop), open that folder instead.
* Open a separate window and go to ``%LOCALAPPDATA%\Microsoft\WindowsApps\``
* Right-click-drag the ``ubuntu.exe`` icon from the second window into the first, then choose ``Create shortcuts here`` from the context menu when you drop it. 
* To customize this shortcut, right-click on it and choose properties.

  * On the General tab:

    * Change the name to whatever you want, e.g. "Sage 9.2 JupyterLab"
  
  * On the Shortcut tab:

    * Change Target to: ``ubuntu.exe run ~/sage_nb.sh``
    * Change Start in to: ``%USERPROFILE%``
    * Change Run to: Minimised
    * Change the icon if you want

Now hit the start button or key and type the name you gave it. it should appear in the list, and should load the server and fire up your browser when you click on it.

------------------------------------------------------------------------

For further reading you can have a look at the other documents in the
SageMath documentation at http://doc.sagemath.org/.


Setting up SageMath as a Jupyter kernel in an existing Jupyter notebook or JupyterLab installation
--------------------------------------------------------------------------------------------------

You may already have a global installation of Jupyter.  For added
convenience, it is possible to link your installation of SageMath into
your Jupyter installation, adding it to the list of available kernels
that can be selected in the notebook or JupyterLab interface.

If ``$SAGE_LOCAL`` is the installation prefix of your Sage
installation (the default is ``$SAGE_ROOT/local``) and you can start
the Jupyter notebook by typing ``jupyter notebook``, then the
following command will install SageMath as a new kernel.

.. CODE-BLOCK:: bash

    jupyter kernelspec install --user $SAGE_LOCAL/share/jupyter/kernels/sagemath

This installs the kernel under the name ``sagemath``.  If you wish to
rename it to something more specific in order to distinguish between
different installations of SageMath, you can use the additional option
``--name``, for example

.. CODE-BLOCK:: bash

    jupyter kernelspec install --user $SAGE_LOCAL/share/jupyter/kernels/sagemath --name sagemath-dev-worktree

The ``jupyter kernelspec`` approach by default does lead to about 2Gb of
sagemath documentation being copied into your personal jupyter configuration
directory. You can avoid that by instead putting a symlink in the relevant spot.
Assuming that sagemath is properly installed, you can use

.. CODE-BLOCK:: bash

    sage -sh -c 'ls -d $SAGE_LOCAL/share/jupyter/kernels/sagemath'

to find location of the sagemath kernel description and

.. CODE-BLOCK:: bash

    jupyter --paths

to find valid data directories for your jupyter installation.
A command along the lines of

.. CODE-BLOCK:: bash

    ln -s `sage -sh -c 'ls -d $SAGE_LOCAL/share/jupyter/kernels/sagemath'` $HOME/.local/share/jupyter

can then be used to create a symlink to the sagemath kernel description
in a location where your own ``jupyter`` can find it.

To get the full functionality of the SageMath kernel in your global
Jupyter installation, the following Notebook Extension packages also
need to be installed (or linked) in the environment from which the
Jupyter installation runs.

You can check the presence of some of these packages using the command
``jupyter nbextension list``.

 - For the Sage interacts, you will need the package
   ``widgetsnbextension`` installed in the Python environment of the
   Jupyter installation.  If your Jupyter installation is coming from
   the system package manager, it is best to install
   ``widgetsnbextension`` in the same way.  Otherwise, install it
   using ``pip``.

   To verify that interacts work correctly, you can evaluate the following code
   in the notebook::

     @interact
     def _(k=slider(vmin=-1.0, vmax= 3.0, step_size=0.1, default=0), auto_update=True):
     plot([lambda u:u^2-1, lambda u:u+k], (-2,2),
          ymin=-1, ymax=3, fill={1:[0]}, fillalpha=0.5).show()

 - For 3D graphics using Three.js, by default, internet connectivity
   is needed, as SageMath's custom build of the Javascript package
   Three.js is retrieved from a content delivery network.

   To verify that online 3D graphics with Three.js works correctly,
   you can evaluate the following code in the notebook::

     plot3d(lambda u,v:(u^2+v^2)/4-2,(-2,2),(-2,2)).show()

   However, it is possible to configure graphics with Three.js for
   offline use.  In this case, the Three.js installation from the Sage
   distribution needs to be made available in the environment of the
   Jupyter installation.  This can be done by copying or symlinking.
   The Three.js installation in the environment of the Jupyter
   installation must exactly match the version that comes from the
   Sage distribution.  It is not supported to use several Jupyter
   kernels corresponding to different versions of the Sage distribution.

   To verify that offline 3D graphics with Three.js works correctly,
   you can evaluate the following code in the notebook::

     plot3d(lambda u,v:(u^2+v^2)/4-2,(-2,2),(-2,2), online=False).show()

 - For 3D graphics using jsmol, you will need the package
   ``jupyter-jsmol`` installed in the Python environment of the
   Jupyter installation. You can install it using ``pip``.
   (Alternatively, you can copy or symlink it.)

   To verify that jsmol graphics work correctly, you can evaluate the
   following code in the notebook::

     plot3d(lambda u,v:(u^2+v^2)/4-2,(-2,2),(-2,2)).show(viewer="jmol")
