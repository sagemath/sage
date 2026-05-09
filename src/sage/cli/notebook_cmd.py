import argparse

from sage.cli.options import CliOptions


class JupyterNotebookCmd:
    @staticmethod
    def extend_parser(parser: argparse.ArgumentParser):
        r"""
        Extend the parser with the Jupyter notebook command.

        INPUT:

        - ``parsers`` -- the parsers to extend.

        OUTPUT:

        - the extended parser.
        """
        parser.add_argument(
            "-n",
            "--notebook",
            nargs="?",
            const="jupyter",
            choices=["jupyter", "jupyterlab"],
            help="start the Jupyter notebook server (default: jupyter)",
        )

    def __init__(self, options: CliOptions):
        r"""
        Initialize the command.
        """
        self.options = options

    def run(self) -> int:
        r"""
        Start the Jupyter notebook server.
        """
        if self.options.notebook == "jupyter":
            try:
                # notebook 6
                from notebook.notebookapp import main
            except ImportError:
                # notebook 7
                from notebook.app import main
        elif self.options.notebook == "jupyterlab":
            from jupyterlab.labapp import main
        else:
            raise ValueError(f"Unknown notebook type: {self.options.notebook}")

        return main([])
