Installing Visualization Software on OS X
=========================================

The first thing we need to do is rebuild Python to use OSX's
frameworks, so that it can create graphical windows. To do this
first from the terminal do

::

    cd $SAGE_ROOT/local/lib
    rm libpng*.dylib

where ``$SAGE_ROOT`` is the directory of your
Sage install. Next from within Sage,

.. skip

::

    sage: install_package('python-2.5.1-framework')

Next we will build vtk, this will take a while

.. skip

::

    sage: install_package('vtk-5.0.3.p1')

Finally

.. skip

::

    sage: install_package('MayaVi-1.5')
    sage: install_package('scitools++')
    sage: install_package('PyVTK-0.4.74')
