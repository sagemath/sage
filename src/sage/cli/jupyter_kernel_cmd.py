import argparse

from sage.cli.options import CliOptions


class JupyterKernelCmd:
    @staticmethod
    def extend_parser(parser: argparse.ArgumentParser):
        r"""
        Extend the parser with the Jupyter kernel installation command.

        INPUT:

        - ``parser`` -- the parser to extend

        OUTPUT: the extended parser
        """
        parser.add_argument(
            "--jupyter-kernel",
            nargs=argparse.REMAINDER,
            metavar="...",
            help="install the SageMath Jupyter kernel spec; run "
            "'sage --jupyter-kernel install --help' for options",
        )

    def __init__(self, options: CliOptions):
        r"""
        Initialize the command.
        """
        self.options = options

    def run(self) -> int:
        r"""
        Install the SageMath Jupyter kernel spec.
        """
        from sage.repl.ipython_kernel.install import main

        main(self.options.jupyter_kernel)
        return 0
