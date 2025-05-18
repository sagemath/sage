from sage.cli.eval_cmd import EvalCmd
from sage.cli.options import CliOptions


def test_eval_cmd_print(capsys):
    options = CliOptions(command="print(3^33)")
    eval_cmd = EvalCmd(options)

    result = eval_cmd.run()
    captured = capsys.readouterr()
    assert captured.out == "5559060566555523\n"
    assert result == 0


def test_eval_cmd_invalid_command(capsys):
    options = CliOptions(command="invalid_command")
    eval_cmd = EvalCmd(options)

    result = eval_cmd.run()
    captured = capsys.readouterr()
    assert (
        "An error occurred while executing the command: name 'invalid_command' is not defined"
        in captured.out
    )
    assert result == 1
