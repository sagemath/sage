.. _sagenb:

=============================
Sage Notebook Developer Guide
=============================


Development of the Sage notebook currently occurs on Github using the
Git revision control system. The development model for the `Sage
Notebook`_ project is a `git <http://git-scm.com>`_ and github_
workflow.

To update to the latest development source, run the commands below,
where ``SAGE_ROOT`` is the root directory of the Sage installation,
and where ``hackdir`` is a directory you create for working on code
changes (it need not have the name or location given below).

.. warning:: This will create a new sagenb repository ignoring any
   changes you have made to the files.

::

    mkdir ~/hackdir
    cd ~/hackdir
    git clone git://github.com/sagemath/sagenb.git sagenb-git
    cd SAGE_ROOT/devel
    rm sagenb
    ln -s ~/hackdir/sagenb sagenb
    cd sagenb
    ../../sage setup.py develop

What this has done is to create a new directory, move to that
directory, and create a clone of the most up-to-date version of
the upstream notebook sources there.  Then we remove a symbolic
link ``sagenb`` in the Sage folder and replace it with a link
to your clone of upstream, finally making sure that the notebook
has the correct dependencies.

An advantage of having the separate directory for sagenb is that
you would later be able to keep it and do development work in it
even when you upgrade Sage, or even if you accidentally destroy your
Sage installation somehow.

The rest of these instructions is some very generic documentation,
slightly adapted to help develop the notebook using Git and Github.

The most important section involves how to update your new sagenb
source repository and create a "fork" of the master copy, so that you
will be able to request your changes to be merged in the Sage notebook,
called a "pull request"; see :ref:`github-development`.


.. toctree::
   :maxdepth: 2

   following_latest
   patching
   github_development


.. include:: links.inc
