.. nodoctest

.. highlight:: shell-session

.. _chapter-tools:

========================================
Additional development and testing tools
========================================

Pytest
===============================
`Pytest <https://docs.pytest.org/en/stable/>`_ is a testing framework.
At the moment, Sage is not yet using any tests based on pytest.

*Installation:* ``pip install -U pytest``, see `documentation <https://docs.pytest.org/en/stable/getting-started.html#installation-and-getting-started>`__ for details.
*Usage:*
- Manual: Run ``pytest path/to/the/test_file.py`` or ``pytest`` to run all tests (from a virtual environment with Sage installed)
- VS Code: Install the `Python <https://marketplace.visualstudio.com/items?itemName=ms-python.python>`_ extension and follow the `offical VS Code documentation <https://code.visualstudio.com/docs/python/testing>`__.
*Configuration:* ``conftest.py`` in the source folder
*Documentation:* https://docs.pytest.org/en/stable/index.html

Pyright 
===============================
`Pyright <https://github.com/microsoft/pyright>`_ is static type checker.

*Installation:* ``npm install -g pyright``, see `documentation <https://github.com/microsoft/pyright#installation>`__ for details.
*Usage:*
- Manual: Run ``pyright path/to/the/file.py``
- VS Code: Install the `Pylance <https://marketplace.visualstudio.com/items?itemName=ms-python.vscode-pylance>`__ extension.
*Configuration:* ``pyrightconfig.json`` in the root
*Note*: Currently, only the manifolds package is checked. Further packages can be added in the ``pyrightconfig.json`` file.
*Documentation:* https://github.com/microsoft/pyright#documentation

Pycodestyle
===============================
`Pycodestyle <https://pycodestyle.pycqa.org/en/latest/>`_ checks against the style conventions of PEP8 Python.

*Installation:* ``pip install -U pycodestyle --user``
*Usage:*
- Manual: Run ``pycodestyle path/to/the/file.py``
- VS Code: Activate by adding the setting ``"python.linting.pycodestyleEnabled": true``, see `official VS Code documentation <https://code.visualstudio.com/docs/python/linting>`__ for details.
*Configuration:* ``[pycodestyle]`` block in ``src/tox.ini``
*Documentation:* https://pycodestyle.pycqa.org/en/latest/index.html

Pyflakes
===============================
`Pyflakes <https://github.com/PyCQA/pyflakes>`_ checks for common coding errors.
