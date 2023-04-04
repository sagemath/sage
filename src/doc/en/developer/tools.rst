.. nodoctest

.. highlight:: shell-session

.. _chapter-tools:

========================================
Additional development and testing tools
========================================

.. _section-tools-tox:

Tox
===

`Tox <https://tox.readthedocs.io/en/latest/>`_ is a popular package that is
used by a large number of Python projects as the standard entry point
for testing and linting.

Sage includes tox as a standard package and uses it for three purposes:

- For portability testing of the Sage distribution, as we explain in
  :ref:`chapter-portability_testing`.  This is configured in the file
  ``SAGE_ROOT/tox.ini``.

- For testing modularized distributions of the Sage library. This is configured
  in ``tox.ini`` files in subdirectories of ``SAGE_ROOT/pkgs/``, such as
  ``SAGE_ROOT/pkgs/sagemath-standard/tox.ini``. Each distribution's configuration
  defines tox environments for testing the distribution with different Python
  versions and different ways how the dependencies are provided.
  We explain this in :ref:`chapter-modularization`.

- As an entry point for testing and linting of the Sage library, as we describe below.
  This is configured in the file ``SAGE_ROOT/src/tox.ini``.

The tox configuration ``SAGE_ROOT/src/tox.ini`` can be invoked by using the command
``./sage --tox``.  (If ``tox`` is available in your system installation,
you can just type ``tox`` instead.)

This configuration provides an entry point for various testing/linting methods,
known as "tox environments".  We can type ``./sage --advanced`` to see what is
available::

  $ ./sage --advanced
  SageMath version 9.2
  ...
  Testing files:
  ...
  --tox [options] <files|dirs> -- general entry point for testing
                                  and linting of the Sage library
     -e <envlist>     -- run specific test environments
                         (default: run all except full pycodestyle)
        doctest                -- run the Sage doctester
                                  (same as "sage -t")
        coverage               -- give information about doctest coverage of files
                                  (same as "sage --coverage[all]")
        startuptime            -- display how long each component of Sage takes to start up
                                  (same as "sage --startuptime")
        pycodestyle-minimal    -- check against Sage's minimal style conventions
        relint                 -- check whether some forbidden patterns appear
                                  (includes all patchbot pattern-exclusion plugins)
        rst                    -- validate Python docstrings markup as reStructuredText
        codespell              -- check for misspelled words in source code
        cython-lint            -- Check Cython files for code style
        pyright                -- run the static typing checker pyright
        pycodestyle            -- check against the Python style conventions of PEP8
     -p auto          -- run test environments in parallel
     --help           -- show tox help


Doctest
=======

The command ``./sage -tox -e doctest`` requires that Sage has been
built already.  ``doctest`` is a special tox environment that runs the
Sage doctester in the normal Sage environment.  This is equivalent to
using the command ``./sage -t``; see :ref:`chapter-doctesting`.


.. _section-tools-coverage:

Coverage
========

The command ``./sage -tox -e coverage`` requires that Sage has been
built already.  ``coverage`` is a special tox environment that is
equivalent to using the command ``./sage --coverageall`` (if no
arguments are provided) or ``./sage --coverage`` (if arguments are
provided).


.. _section-tools-startuptime:

Startuptime
===========

The command ``./sage -tox -e startuptime`` requires that Sage has been
built already.  ``startuptime`` is a special tox environment that is
equivalent to using the command ``./sage --startuptime``.


.. _section-tools-pycodestyle:

Pycodestyle
===========
`Pycodestyle <https://pycodestyle.pycqa.org/en/latest/>`_ (formerly known as pep8)
checks Python code against the style conventions of `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_.
Tox automatically installs pycodestyle in a separate virtual environment
on the first use.

Sage defines two configurations for pycodestyle.  The command ``./sage -tox -e pycodestyle-minimal`` uses
pycodestyle in a minimal configuration.
As of Sage 9.5, the entire Sage library conforms to this configuration::

  $ ./sage -tox -e pycodestyle-minimal -- src/sage/
  pycodestyle-minimal installed: pycodestyle==2.8.0
  pycodestyle-minimal run-test-pre: PYTHONHASHSEED='28778046'
  pycodestyle-minimal run-test: commands[0] | pycodestyle --select E401,E70,W605,E711,E712,E721 sage
  ___________ summary ____________
    pycodestyle-minimal: commands succeeded
    congratulations :)

When preparing a branch for a Sage ticket, developers should verify that ``./sage -tox -e
pycodestyle-minimal`` passes.  When the Sage patchbot runs on the ticket, it will perform similar
coding style checks; but running the check locally reduces the turnaround time from hours
to seconds.

The second configuration is used with the command ``./sage -tox -e pycodestyle`` and runs a
more thorough check::

  $ ./sage -tox -e pycodestyle -- src/sage/quadratic_forms/quadratic_form.py
  pycodestyle installed: pycodestyle==2.8.0
  pycodestyle run-test-pre: PYTHONHASHSEED='2520226550'
  pycodestyle run-test: commands[0] | pycodestyle sage/quadratic_forms/quadratic_form.py
  sage/quadratic_forms/quadratic_form.py:135:9: E225 missing whitespace around operator
  sage/quadratic_forms/quadratic_form.py:163:64: E225 missing whitespace around operator
  sage/quadratic_forms/quadratic_form.py:165:52: E225 missing whitespace around operator
  sage/quadratic_forms/quadratic_form.py:173:42: E228 missing whitespace around modulo operator
  ...
  sage/quadratic_forms/quadratic_form.py:1620:9: E266 too many leading '#' for block comment
  sage/quadratic_forms/quadratic_form.py:1621:9: E266 too many leading '#' for block comment
  25      E111 indentation is not a multiple of 4
  2       E117 over-indented
  129     E127 continuation line over-indented for visual indent
  1       E128 continuation line under-indented for visual indent
  4       E201 whitespace after '['
  4       E202 whitespace before ']'
  2       E222 multiple spaces after operator
  7       E225 missing whitespace around operator
  1       E228 missing whitespace around modulo operator
  25      E231 missing whitespace after ','
  1       E262 inline comment should start with '# '
  3       E265 block comment should start with '# '
  62      E266 too many leading '#' for block comment
  2       E272 multiple spaces before keyword
  2       E301 expected 1 blank line, found 0
  17      E303 too many blank lines (2)
  ERROR: InvocationError for command .../pycodestyle sage/quadratic_forms/quadratic_form.py (exited with code 1)
  ___________ summary ____________
  ERROR:   pycodestyle: commands failed

When preparing a branch for a Sage ticket that adds new code,
developers should verify that ``./sage -tox -e pycodestyle`` does not
issue warnings for the added code.  This will avoid later cleanup
tickets as the Sage codebase is moving toward full PEP 8 compliance.

On the other hand, it is usually not advisable to mix coding-style
fixes with productive changes on the same ticket because this would
makes it harder for reviewers to evaluate the changes.

By passing the options ``--count -qq`` we can reduce the output to
only show the number of style violation warnings.  This can be helpful
for planning work on coding-style clean-up tickets that focus on one
or a few related issues::

  $ ./sage -tox -e pycodestyle -- --count -qq src/sage
  pycodestyle installed: pycodestyle==2.8.0
  pycodestyle run-test-pre: PYTHONHASHSEED='3166223974'
  pycodestyle run-test: commands[0] | pycodestyle --count -qq sage
  557     E111 indentation is not a multiple of 4
  1       E112 expected an indented block
  194     E114 indentation is not a multiple of 4 (comment)
  ...
  7       E743 ambiguous function definition 'l'
  335     W291 trailing whitespace
  4       W292 no newline at end of file
  229     W293 blank line contains whitespace
  459     W391 blank line at end of file
  97797
  ERROR: InvocationError for command .../pycodestyle --count -qq sage (exited with code 1)
  ___________ summary ____________
  ERROR:   pycodestyle: commands failed

*Installation:* (for manual use:) ``pip install -U pycodestyle --user``

*Usage:*

- With tox: See above.

- Manual: Run ``pycodestyle path/to/the/file.py``.

- VS Code: The minimal version of pycodestyle is activated by default in
  ``SAGE_ROOT/.vscode/settings.json`` (the corresponding setting is
  ``"python.linting.pycodestyleEnabled": true``). Note that the
  ``settings.json`` file is not ignored by git so be aware to keep it in sync
  with the trac repo. For further details, see the
  `official VS Code documentation <https://code.visualstudio.com/docs/python/linting>`__.

*Configuration:* ``[pycodestyle]`` block in ``SAGE_ROOT/src/tox.ini``

*Documentation:* https://pycodestyle.pycqa.org/en/latest/index.html


.. _section-tools-cython-lint:

Cython-lint
===========

`Cython-lint <https://pypi.org/project/cython-lint/>`_ checks Cython source files
for coding style.

.. _section-tools-relint:

Relint
======

`Relint <https://pypi.org/project/relint/>`_ checks all source files for forbidden
text patterns specified by regular expressions.

Our configuration of relint flags some outdated Python constructions, plain TeX
commands when equivalent LaTeX commands are available, common mistakes in
documentation markup, and modularization anti-patterns.

*Configuration:* ``SAGE_ROOT/src/.relint.yml``

*Documentation:* https://pypi.org/project/relint/


.. _section-tools-codespell:

Codespell
=========
`Codespell <https://pypi.org/project/codespell/>`_ uses a dictionary to check for
misspelled words in source code.

Sage defines a configuration for codespell::

  $ ./sage -tox -e codespell -- src/sage/homology/
  codespell installed: codespell==2.1.0
  codespell run-test-pre: PYTHONHASHSEED='1285169064'
  codespell run-test: commands[0] | codespell '--skip=*.png,*.jpg,*.JPG,*.inv,*.dia,*.pdf,*.ico,*#*,*~*,*.bak,*.orig,*.log,*.sobj,*.tar,*.gz,*.pyc,*.o,*.sws,*.so,*.a,.DS_Store' --skip=doc/ca,doc/de,doc/es,doc/hu,doc/ja,doc/ru,doc/fr,doc/it,doc/pt,doc/tr --skip=src/doc/ca,src/doc/de,src/doc/es,src/doc/hu,src/doc/ja,src/doc/ru,src/doc/fr,src/doc/it,src/doc/pt,src/doc/tr '--skip=.git,.tox,worktree*,dist,upstream,logs,local,cythonized,scripts-3,autom4te.cache,tmp,lib.*,*.egg-info' --dictionary=- --dictionary=/Users/mkoeppe/s/sage/sage-rebasing/src/.codespell-dictionary.txt --ignore-words=/Users/mkoeppe/s/sage/sage-rebasing/src/.codespell-ignore.txt sage/homology
  sage/homology/hochschild_complex.py:271: mone ==> mono, money, none
  sage/homology/hochschild_complex.py:277: mone ==> mono, money, none
  sage/homology/hochschild_complex.py:280: mone ==> mono, money, none
  sage/homology/chain_complex.py:2185: mor ==> more
  sage/homology/chain_complex.py:2204: mor ==> more
  sage/homology/chain_complex.py:2210: mor ==> more
  sage/homology/chain_complex.py:2211: mor ==> more
  sage/homology/chain_complex.py:2214: mor ==> more
  sage/homology/chain_complex.py:2215: mor ==> more
  ERROR: InvocationError for command .../codespell '--skip=*.png,...' --dictionary=- --dictionary=/Users/mkoeppe/s/sage/sage-rebasing/src/.codespell-dictionary.txt --ignore-words=/Users/mkoeppe/s/sage/sage-rebasing/src/.codespell-ignore.txt sage/homology (exited with code 65)
  ___________ summary ____________
  ERROR:   codespell: commands failed

*Configuration:*

- ``[testenv:codespell]`` block in ``SAGE_ROOT/src/tox.ini``

- ``SAGE_ROOT/src/.codespell-dictionary.txt`` and ``SAGE_ROOT/src/.codespell-ignore.txt``


.. _section-tools-pytest:

Pytest
======
`Pytest <https://docs.pytest.org/en/stable/>`_ is a testing framework.
It is included in the Sage distribution as an optional package.

Currently, Sage only makes very limited use of pytest, for testing the
package :mod:`sage.numerical.backends` and some modules in
:mod:`sage.manifolds`.

*Installation:*

- ``./sage -i pytest pytest_xdist``.

*Usage:*

- Tox, Sage doctester: At the end of ``./sage -t`` (or ``./sage --tox -e doctest``), Pytest is automatically invoked.

- Manual: Run ``./sage -pytest path/to/the/test_file.py`` or ``./sage -pytest``
  to run all tests. The additional argument ``-n`` can be used to
  distribute tests across multiple CPUs to speed up test execution.
  For example, ``./sage -pytest -n 4`` will run 4 tests in parallel, while
  ``./sage -pytest -n auto`` will spawn a number of workers processes equal
  to the number of available CPUs.

- VS Code: Install the `Python extension <https://marketplace.visualstudio.com/items?itemName=ms-python.python>`_ and follow the `offical VS Code documentation <https://code.visualstudio.com/docs/python/testing>`__.

*Configuration:* ``SAGE_ROOT/src/conftest.py``

*Documentation:* https://docs.pytest.org/en/stable/index.html


.. _section-tools-pyright:

Pyright
=======
`Pyright <https://github.com/microsoft/pyright>`_ is static type checker.

*Installation:*

- (for manual use:) ``npm install -g pyright``, see `documentation <https://github.com/microsoft/pyright#installation>`__ for details.

*Usage:*

- Tox: Run ``./sage -tox -e pyright path/to/the/file.py``

- Manual: Run ``pyright path/to/the/file.py``

- VS Code: Install the `Pylance <https://marketplace.visualstudio.com/items?itemName=ms-python.vscode-pylance>`__ extension.

*Configuration:* ``SAGE_ROOT/pyrightconfig.json``

*Documentation:* https://github.com/microsoft/pyright#documentation


.. _section-tools-pyflakes:

Pyflakes
========
`Pyflakes <https://github.com/PyCQA/pyflakes>`_ checks for common coding errors.
