from dataclasses import dataclass


@dataclass
class CliOptions:
    """
    A TypedDict for command-line interface options.
    """

    """Indicates whether verbose output is enabled."""
    verbose: bool = False

    """The notebook type to start."""
    notebook: str = "jupyter"

    """The command to execute."""
    command: str | None = None
