import argparse
import logging
import os
import sys


def _source_sagerc() -> None:
    """
    Read the shell startup file before importing Sage.

    Restart through Bash so that both exported variables and shell settings
    (such as the working directory and umask) apply to the Python process.
    The process ID identifies this restart; a child Sage process must still
    read its own startup file.
    """
    if os.environ.pop('_SAGE_RC_SOURCED', None) == str(os.getpid()):
        return
    # Match the environment marker used by src/bin/sage-env. An inherited
    # marker for a different environment does not mean this one is ready.
    sage_env_version = '6:' + ':'.join(
        os.environ.get(name, '') for name in ('SAGE_LOCAL', 'SAGE_VENV', 'SAGE_SRC')
    )
    if os.environ.get('SAGE_ENV_SOURCED') == sage_env_version:
        # The legacy Bash launcher has already read the startup file.
        return

    dot_sage = os.environ.get('DOT_SAGE') or os.path.expanduser('~/.sage')
    rc_file = os.environ.get('SAGE_RC_FILE') or os.path.join(dot_sage, 'sagerc')
    if not os.access(rc_file, os.R_OK):
        return

    # Keep interpreter options as well as CLI arguments. A relative console
    # script must remain reachable even if sagerc changes the working directory.
    argv = sys.orig_argv.copy()
    script_index = len(argv) - len(sys.argv)
    if argv[script_index] == sys.argv[0] and sys.argv[0] not in ('-c', '-'):
        argv[script_index] = os.path.abspath(sys.argv[0])

    environment = os.environ.copy()
    environment['DOT_SAGE'] = dot_sage
    os.execvpe('bash', [
        'bash', '--noprofile', '--norc', '-c', '''
_sage_source_rc() {
    SAGE_RC_FILE=${SAGE_RC_FILE:-$DOT_SAGE/sagerc}
    source "$1"
}
# A function keeps any "set --" in sagerc from replacing the CLI arguments.
_sage_source_rc "$1" >&2
if [ $? -ne 0 ]; then
    printf 'Error sourcing %s\\n' "$1" >&2
    exit 1
fi
shift
export _SAGE_RC_SOURCED=$$ || exit 1
exec "$@"
''', 'sage-sagerc', os.path.abspath(rc_file), sys.executable, *argv[1:]
    ], environment)


def main() -> int:
    _source_sagerc()

    # These imports eventually load sage.env and sage.all, which must see the
    # environment set by sagerc.
    from sage.cli.eval_cmd import EvalCmd
    from sage.cli.interactive_shell_cmd import InteractiveShellCmd
    from sage.cli.notebook_cmd import JupyterNotebookCmd
    from sage.cli.options import CliOptions
    from sage.cli.run_file_cmd import RunFileCmd
    from sage.cli.version_cmd import VersionCmd

    input_args = sys.argv[1:]
    parser = argparse.ArgumentParser(
        prog="sage",
        description="If no command is given, starts the interactive interpreter where you can enter statements and expressions, immediately execute them and see their results.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="print additional information",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        default=False,
        help="do not display the banner",
    )
    parser.add_argument(
        "--simple-prompt",
        action="store_true",
        default=False,
        help="use simple prompt IPython mode",
    )

    VersionCmd.extend_parser(parser)
    JupyterNotebookCmd.extend_parser(parser)
    EvalCmd.extend_parser(parser)
    RunFileCmd.extend_parser(parser)

    if not input_args:
        return InteractiveShellCmd(CliOptions()).run()

    args = parser.parse_args(input_args)
    options = CliOptions(**vars(args))

    logging.basicConfig(level=logging.DEBUG if options.verbose else logging.INFO)

    if args.file:
        return RunFileCmd(options).run()
    if args.command:
        return EvalCmd(options).run()
    if args.notebook:
        return JupyterNotebookCmd(options).run()
    return InteractiveShellCmd(options).run()
