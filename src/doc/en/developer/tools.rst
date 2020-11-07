.. nodoctest

.. highlight:: shell-session

.. _chapter-tools:

========================================
Additional development and testing tools
========================================

Pyright 
===============================
`Pyright <https://github.com/microsoft/pyright>` is static type checker.

*Installation:* ``npm install -g pyright``, see `documentation <https://github.com/microsoft/pyright#installation>`_ for details.
*Usage:*
- Manual: Run ``pyright path/to/the/file.py``
- VS Code: Install the `Pylance <https://marketplace.visualstudio.com/items?itemName=ms-python.vscode-pylance>`_ extension.
*Configuration:* ``pyrightconfig.json`` in the root
*Note*: Currently, only the manifolds package is checked. Further packages can be added in the ``pyrightconfig.json`` file.
*Documentation:* https://github.com/microsoft/pyright#documentation

Pycodestyle
===============================
`Pycodestyle <https://pycodestyle.pycqa.org/en/latest/>`_ checks against the style conventions of PEP8 Python.

*Installation:* ``pip install -U pycodestyle --user``
*Usage:*
- Manual: Run ``pycodestyle path/to/the/file.py``
- VS Code: Activate by adding the setting ``"python.linting.pycodestyleEnabled": true``, see `official VS Code documentation <https://code.visualstudio.com/docs/python/linting>`_ for details.
*Configuration:* ``[pycodestyle]`` block in ``src/tox.ini``
*Documentation:* https://pycodestyle.pycqa.org/en/latest/index.html

Pyflakes
===============================
`Pyflakes <https://github.com/PyCQA/pyflakes>`_ checks for common coding errors.
