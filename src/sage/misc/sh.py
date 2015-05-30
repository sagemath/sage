"Evaluating shell scripts"

import os


class Sh:
    r"""
    Evaluates a shell script and returns the output.

    To use this from the notebook type ``sh`` at the beginning of
    the input cell.  The working directory is then the (usually
    temporary) directory where the Sage worksheet process is
    executing.
    """
    def eval(self, code, globals=None, locals=None):
        r"""
        This is difficult to test because the output goes to the
        screen rather than being captured by the doctest program, so
        the following really only tests that the command doesn't bomb,
        not that it gives the right output::

            sage: sh.eval('''echo "Hello there"\nif [ $? -eq 0 ]; then\necho "good"\nfi''') # random output
        """
        # Print out the current absolute path, which is where the code
        # will be evaluated.    Evidently, users find this comforting,
        # though I personally find it to be a bit much (William Stein).
        print(os.path.abspath('.'))
        # Evaluate the input code block.  Fortunately, os.system works
        # fine with multiline input (in contrast to subprocess.Popen).
        os.system(str(code))
        # Return '' so nothing extra (for example an unsightly None)
        # gets printed when doing %sh in the notebook.
        return ''

# Create the sh object, so that %sh mode works in the notebook.
sh = Sh()
