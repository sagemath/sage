import argparse

from sage.cli.jupyter_kernel_cmd import JupyterKernelCmd


def test_no_jupyter_kernel_by_default():
    parser = argparse.ArgumentParser()
    JupyterKernelCmd.extend_parser(parser)
    args = parser.parse_args([])
    assert args.jupyter_kernel is None


def test_install_subcommand_collected():
    parser = argparse.ArgumentParser()
    JupyterKernelCmd.extend_parser(parser)
    args = parser.parse_args(["--jupyter-kernel", "install"])
    assert args.jupyter_kernel == ["install"]


def test_remaining_options_passed_through():
    parser = argparse.ArgumentParser()
    JupyterKernelCmd.extend_parser(parser)
    args = parser.parse_args(["--jupyter-kernel", "install", "--user", "--portable"])
    assert args.jupyter_kernel == ["install", "--user", "--portable"]
