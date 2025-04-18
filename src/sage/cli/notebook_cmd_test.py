import argparse

import pytest

from sage.cli.notebook_cmd import JupyterNotebookCmd


def test_jupyter_as_default():
    parser = argparse.ArgumentParser()
    JupyterNotebookCmd.extend_parser(parser)
    args = parser.parse_args(["--notebook"])
    assert args.notebook == "jupyter"


def test_jupyter_explicitly():
    parser = argparse.ArgumentParser()
    JupyterNotebookCmd.extend_parser(parser)
    args = parser.parse_args(["--notebook", "jupyter"])
    assert args.notebook == "jupyter"


def test_jupyterlab_explicitly():
    parser = argparse.ArgumentParser()
    JupyterNotebookCmd.extend_parser(parser)
    args = parser.parse_args(["--notebook", "jupyterlab"])
    assert args.notebook == "jupyterlab"


def test_invalid_notebook_choice():
    parser = argparse.ArgumentParser()
    JupyterNotebookCmd.extend_parser(parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["--notebook", "invalid"])


def test_help():
    parser = argparse.ArgumentParser()
    JupyterNotebookCmd.extend_parser(parser)
    assert parser.format_usage() == "usage: pytest [-h] [-n [{jupyter,jupyterlab}]]\n"
