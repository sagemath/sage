.. highlight:: shell-session

.. _chapter-workspace-setup:

=========================
Setting up Your Workspace
=========================

.. _section-ide:

Text editors and IDEs for use with Sage
=======================================

In meta :issue:`30500` we are collecting instructions how to configure
various text editors and IDEs for use with Sage.


.. _section-gitpod:

Gitpod
======

`Gitpod <https://www.gitpod.io>`_ is a free service that will let you build and
run Sage from an online development environment based on VS Code.
Without needing to install anything on your computer, Gitpod creates a virtual
fully-functional workspace with all the dependencies and tools pre-installed.

To get started, `go to Gitpod
<https://gitpod.io/#https://github.com/sagemath/sage>`_ and log in using your
GitHub or GitLab account.  Wait while Gitpod creates a workspace.  The first
time, it may take some time to build Sage.

You can now run and edit Sage's code. Contributing your changes follows the normal
:ref:`Git workflow <chapter-git-basic>`.
For this to work, you first have to authorize Gitpod with GitHub:

1. In the running Gitpod workspace, generate a new SSH key pair by ``ssh-keygen -f tempkey``.
2. Save the private key as a secure environment variable in Gitpod using ``gp
   env PRIVATE_SSH_KEY="$(<tempkey)"``, or by using the `Gitpod UI
   <https://www.gitpod.io/docs/environment-variables#using-the-account-settings>`_.
3. Register the public key with GitHub following the instructions in :ref:`section-github-ssh-key`.
4. Close this Gitpod workspace.

After following this procedure, every new Gitpod workspace will have a
working ``origin`` remote to which you can push your changes.

You can also `use your VS Code Desktop <https://www.gitpod.io/docs/develop/vscode-desktop-support>`_ to keep
your local IDE configuration while still benefiting from Gitpod’s high-spec servers and automated prebuilds.
