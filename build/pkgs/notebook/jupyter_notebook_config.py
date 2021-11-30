# Configuration file for Sage's builtin Jupyter notebook server

# Disable XSRF checking to fix Thebe. See
# * https://trac.sagemath.org/ticket/22458
# * https://github.com/oreillymedia/thebe/issues/93

c.NotebookApp.disable_check_xsrf = True

# send2trash sometimes doesn't work, so disable that for now
# See https://github.com/jupyter/notebook/issues/3249

c.FileContentsManager.delete_to_trash = False
