from dataclasses import dataclass


@dataclass
class CliOptions:
    """
    A TypedDict for command-line interface options.
    """

    """Indicates whether verbose output is enabled."""
    verbose: bool = False

    """Indicates whether the banner should be displayed."""
    quiet: bool = False

    """Indicates whether the IPython simple prompt should be used."""
    simple_prompt: bool = False

    """The notebook type to start."""
    notebook: str = "jupyter"

    """The command to execute."""
    command: str | None = None

    """The file to execute."""
    file: str | None = None
