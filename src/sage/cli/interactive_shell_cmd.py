from sage.cli.options import CliOptions


class InteractiveShellCmd:
    def __init__(self, options: CliOptions):
        r"""
        Initialize the command.
        """
        self.options = options

    def run(self) -> int:
        r"""
        Start the interactive shell.
        """
        # Display startup banner. Do this before anything else to give the user
        # early feedback that Sage is starting.
        from sage.misc.banner import banner

        banner()

        from sage.repl.interpreter import SageTerminalApp

        app = SageTerminalApp.instance()
        app.initialize([])
        return app.start()  # type: ignore
