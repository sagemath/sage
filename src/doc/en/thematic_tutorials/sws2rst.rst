.. _sws2rst:

=====================================================
Creating a Tutorial from an old Sage Worksheet (.sws)
=====================================================

A lot of pedagogical material has been written in the Sage Worksheet format, which is no longer supported by Sage after the transition to Python 3 and the removal of the SageNB package.

However, it is possible to convert Sage Worksheet files.
Once you have created a worksheet and are satisfied with the text and
computations, download it to a directory.

We will assume here that the worksheet is called ``Tutorial.sws``
and the directory is called ``make_tutorial``.  We also assume that
``sage`` is your Sage command; if it is not in your ``PATH`` then replace
this with the path to your Sage installation, such as
``/Applications/Sage-9.2.app/Contents/Resources/sage/sage`` if you are
using the Mac app and have placed it in your Applications directory.

* Next, you will need an optional package to convert your worksheet.  Use the
  command:

  .. CODE-BLOCK:: shell-session

      $ sage -i sage_sws2rst

  to install it (or, in the Mac app, use the ``Terminal Session`` advanced
  menu).

* Then we will use the ``sws2rst`` script to turn the worksheet into
  a document in the `ReStructuredText <http://sphinx-doc.org/rest.html>`_
  format.  Be sure you are in the same directory as the worksheet:

  .. CODE-BLOCK:: shell-session

      $ sage --sws2rst Tutorial.sws

  This will create an ``.rst`` file along with a subdirectory of image
  files (which may be empty if there are no images).

  You can find help for ``sws2rst`` with the command
  ``sage --sws2rst -h`` once you have installed beautifulsoup4.

* In principle, such a file could be added directly to Sage's documentation (see
  the `developer's manual <../developer/index.html>`_). However, you probably
  want to check whether it looks right first. So next we will compile this file
  to html documentation.

  * Follow the instructions of ``sage --sws2rst --sphinxify``.  First,
    we will open a Sage shell session, where all appropriate Sage
    references already work properly:

    .. CODE-BLOCK:: shell-session

        $ sage --sh

    From here, you should be able to just type:

    .. CODE-BLOCK:: shell-session

        $ sphinx-quickstart

    and then respond to prompts for turning your ``.rst`` file into
    documentation.  For most of them you can just hit enter/return to
    accept the defaults.  However, you will probably want to

    * Enter a name for the project
    * Enter a name for you
    * Type ``y`` for the question about using MathJax

    Keep note of the instructions; the main other thing to do is add
    your file's name to ``index.rst``, and then just do:

    .. CODE-BLOCK:: shell-session

        $ make html

    and wait while magic happens.  To see the results, open the file
    ``make_tutorial/_build/html/Tutorial.html`` with a browser, or
    use your graphical file system to navigate to the same place.

* Now you can modify the ``.rst`` file more and repeat the steps
  of compiling it until it is ready for inclusion, or just for distribution
  among other Sage users as an HTML file.  (Do ``make pdf`` for a PDF
  version.)
