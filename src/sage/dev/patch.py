r"""
Handling of deprecated hg patch files.

This file contains the methods of :class:`sagedev.SageDev` which are related to
mercurial patches.

AUTHORS:

- David Roe, Frej Drejhammar, Julian Rueth, Martin Raum, Nicolas M. Thiery, R.
  Andrew Ohana, Robert Bradshaw, Timo Kluck: initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 David Roe <roed.math@gmail.com>
#                          Frej Drejhammar <frej.drejhammar@gmail.com>
#                          Julian Rueth <julian.rueth@fsfe.org>
#                          Martin Raum <martin@raum-brothers.eu>
#                          Nicolas M. Thiery <Nicolas.Thiery@u-psud.fr>
#                          R. Andrew Ohana <andrew.ohana@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          Timo Kluck <tkluck@infty.nl>
#                          Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import re

# regular expressions to parse mercurial patches
HG_HEADER_REGEX = re.compile(r"^# HG changeset patch$")
HG_USER_REGEX = re.compile(r"^# User (.*)$")
HG_DATE_REGEX = re.compile(r"^# Date (\d+) (-?\d+)$")
HG_NODE_REGEX = re.compile(r"^# Node ID ([0-9a-f]+)$")
HG_PARENT_REGEX = re.compile(r"^# Parent +([0-9a-f]+)$")
HG_DIFF_REGEX = re.compile(r"^diff (?:-r [0-9a-f]+ ){1,2}(.*)$")
PM_DIFF_REGEX = re.compile(r"^(?:(?:\+\+\+)|(?:---)) [ab]/([^ ]*)(?: .*)?$")
MV_DIFF_REGEX = re.compile(r"^rename (?:(?:to)|(?:from)) (.*)$")
CP_DIFF_REGEX = re.compile(r"^copy (?:(?:to)|(?:from)) (.*)$")

# regular expressions to parse git patches -- at least those created by us
GIT_FROM_REGEX = re.compile(r"^From: (.*)$")
GIT_SUBJECT_REGEX = re.compile(r"^Subject: (.*)$")
GIT_DATE_REGEX = re.compile(r"^Date: (.*)$")
GIT_DIFF_REGEX = re.compile(r"^diff --git a/(.*) b/(.*)$") # this regex should work for our patches
                                                           # since we do not have spaces in file names

# regular expressions to determine whether a path was written for the new git
# repository or for the old hg repositoryw
HG_PATH_REGEX = re.compile(r"^(?=sage/)|(?=doc/)|(?=module_list\.py)|(?=setup\.py)|(?=c_lib/)")
GIT_PATH_REGEX = re.compile(r"^(?=src/)")

from user_interface_error import OperationCancelledError
from git_error import GitError


class MercurialPatchMixin(object):

    def import_patch(self, patchname=None, url=None, local_file=None,
                     diff_format=None, header_format=None, path_format=None):
        r"""
        Legacy support: Import a patch into the current branch.

        If no arguments are given, then all patches from the ticket are
        downloaded and applied using :meth:`download_patch`.

        If ``local_file`` is specified, apply the file it points to.

        Otherwise, download the patch using :meth:`download_patch` and apply
        it.

        INPUT:

        - ``patchname`` -- a string or ``None`` (default: ``None``), passed on
          to :meth:`download_patch`

        - ``url`` -- a string or ``None`` (default: ``None``), passed on to
          :meth:`download_patch`

        - ``local_file`` -- a string or ``None`` (default: ``None``), if
          specified, ``url`` and ``patchname`` must be ``None``; instead of
          downloading the patch, apply this patch file.

        - ``diff_format`` -- a string or ``None`` (default: ``None``), per
          default the format of the patch file is autodetected; it can be
          specified explicitly with this parameter

        - ``header_format`` -- a string or ``None`` (default: ``None``), per
          default the format of the patch header is autodetected; it can be
          specified explicitly with this parameter

        - ``path_format`` -- a string or ``None`` (default: ``None``), per
          default the format of the paths is autodetected; it can be specified
          explicitly with this parameter

        .. NOTE::

            This method calls :meth:`_rewrite_patch` if necessary to rewrite
            patches which were created for sage before the move to git
            happened. In other words, this is not just a simple wrapper for
            ``git am``.

        .. SEEALSO::

            - :meth:`download_patch` -- download a patch to a local file.

            - :meth:`download` -- merges in changes from a git branch
              rather than a patch.

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create a patch::

            sage: open("tracked", "w").close()
            sage: open("tracked2", "w").close()
            sage: patchfile = os.path.join(dev._sagedev.tmp_dir,"tracked.patch")
            sage: dev.git.silent.add("tracked", "tracked2")
            sage: with open(patchfile, "w") as f: f.write(dev.git.diff(cached=True))
            sage: dev.git.silent.reset()

        Applying this patch fails since we are not in a sage repository::

            sage: dev.import_patch(local_file=patchfile, path_format="new")
            There are untracked files in your working directory:
            tracked
            tracked2
            The patch cannot be imported unless these files are removed.

        After moving away ``tracked`` and ``tracked2``, this works::

            sage: os.unlink("tracked")
            sage: os.unlink("tracked2")
            sage: dev.import_patch(local_file=patchfile, path_format="new")
            Applying: No Subject. Modified: tracked, tracked2

         We create a patch which does not apply::

            sage: with open("tracked", "w") as f: f.write("foo")
            sage: dev.git.silent.add("tracked")
            sage: with open("tracked", "w") as f: f.write("boo")
            sage: with open("tracked2", "w") as f: f.write("boo")
            sage: with open(patchfile, "w") as f: f.write(dev.git.diff())
            sage: dev.git.clean_wrapper()
            sage: open("tracked").read()
            ''

         The import fails::

            sage: UI.append("abort")
            sage: UI.append("y")
            sage: dev.import_patch(local_file=patchfile, path_format="new")
            Applying: No Subject. Modified: tracked, tracked2
            Patch failed at 0001 No Subject. Modified: tracked, tracked2
            The copy of the patch that failed is found in:
               .../rebase-apply/patch
            <BLANKLINE>
            The patch does not apply cleanly. Reject files will be created for the parts
            that do not apply if you proceed.
            Apply it anyway? [yes/No] y
            The patch did not apply cleanly. Please integrate the `.rej` files that were
            created and resolve conflicts. After you do, type `resolved`. If you want to
            abort this process, type `abort`. [resolved/abort] abort
            Removing tracked.rej
            sage: open("tracked").read()
            ''

            sage: UI.append("resolved")
            sage: UI.append("y")
            sage: dev.import_patch(local_file=patchfile, path_format="new")
            Applying: No Subject. Modified: tracked, tracked2
            Patch failed at 0001 No Subject. Modified: tracked, tracked2
            The copy of the patch that failed is found in:
               .../rebase-apply/patch
            <BLANKLINE>
            The patch does not apply cleanly. Reject files will be created for the parts
            that do not apply if you proceed.
            Apply it anyway? [yes/No] y
            The patch did not apply cleanly. Please integrate the `.rej` files that were
            created and resolve conflicts. After you do, type `resolved`. If you want to
            abort this process, type `abort`. [resolved/abort] resolved
            Removing tracked.rej
            sage: open("tracked").read() # we did not actually incorporate the .rej files in this doctest, so nothing has changed
            ''
            sage: open("tracked2").read()
            'boo'
        """
        try:
            self.reset_to_clean_state()
            self.clean()
        except OperationCancelledError:
            self._UI.error("Cannot import patch. Your working directory is not in a clean state.")
            raise

        untracked = self.git.untracked_files()
        # do not exclude .patch files here: they would be deleted by clean() later
        if untracked:
            self._UI.error("There are untracked files in your working directory:\n{0}\nThe patch cannot be imported unless these files are removed.".format("\n".join(untracked)))
            raise OperationCancelledError("untracked files make import impossible")

        if not local_file:
            local_files = self.download_patch(patchname=patchname, url=url)
            try:
                for local_file in local_files:
                    lines = open(local_file).read().splitlines()
                    if not self._UI.confirm("I will apply a patch which starts with the following lines:\n{0}\n\nIt modifies the following files:\n{1}\nIs this what you want?".format("\n".join(lines[:8]),", ".join(self._detect_patch_modified_files(lines))), default=True):
                        self._UI.show("Skipping this patch.")
                    else:
                        self.import_patch(
                            local_file=local_file,
                            diff_format=diff_format, header_format=header_format, path_format=path_format)
            finally:
                for local_file in local_files:
                    self._UI.debug("Deleting {0}.".format(local_file))
                    os.unlink(local_file)
        elif patchname or url:
            from sagedev import SageDevValueError
            raise SageDevValueError("if local_file is specified, patchname and url must not be specified")
        else:
            lines = open(local_file).read().splitlines()
            lines = self._rewrite_patch(lines, to_header_format="git",
                    to_path_format="new", from_diff_format=diff_format,
                    from_header_format=header_format,
                    from_path_format=path_format)

            from sage.dev.misc import tmp_filename
            outfile = tmp_filename()
            with open(outfile, 'w') as f:
                f.write("\n".join(lines)+"\n")

            self._UI.debug("Trying to apply reformatted patch `%s`"%outfile)
            # am, apply and add need to be in the root directory
            curdir = os.getcwd()
            os.chdir(self.git._src)
            try:
                try:
                    self.git.echo.am(outfile, "--resolvemsg= ", ignore_whitespace=True)
                except GitError as err:
                    self._UI.error([err.stdout, ''])
                    self._UI.warning("The patch does not apply cleanly. Reject files will be"
                                     " created for the parts that do not apply if you proceed.")
                    if not self._UI.confirm("Apply it anyway?", default=False):
                        self._UI.debug("Not applying patch.")
                        self.git.reset_to_clean_state()
                        self.git.clean_wrapper(remove_untracked_files=True)
                        raise OperationCancelledError("User requested to cancel the apply.")

                    try:
                        try:
                            self.git.silent.apply(outfile, ignore_whitespace=True, reject=True)
                        except GitError:
                            if self._UI.select("The patch did not apply cleanly. Please integrate the `.rej` files that were created and resolve conflicts. After you do, type `resolved`. If you want to abort this process, type `abort`.", ("resolved","abort")) == "abort":
                                raise OperationCancelledError("User requested to cancel the apply.")
                        else:
                            self._UI.show("It seemed that the patch would not apply, but in fact it did.")

                        self.git.super_silent.add(update=True)
                        untracked = [fname for fname in self.git.untracked_files() if not fname.endswith(".rej")]
                        if untracked and self._UI.confirm("The patch will introduce the following new files to the repository:\n{0}\nIs this correct?".format("\n".join(untracked)), default=True):
                            self.git.super_silent.add(*untracked)
                        self.git.am('--resolvemsg= ', resolved=True)
                        self._UI.debug("A commit on the current branch has been created from the patch.")
                    finally:
                        self.git.reset_to_clean_state()
                        self.git.clean_wrapper(remove_untracked_files=True, remove_untracked_directories=True)
            finally:
                os.chdir(curdir)

    def download_patch(self, ticket=None, patchname=None, url=None):
        r"""
        Legacy support: Download a patch to a temporary directory.

        If only ``ticket`` is specified, then try to make sense of the
        ``apply`` statements in the comments on the ticket to download the
        tickets in the right order just like the patchbot would do.

        If no ``ticket`` is specified, use the current ticket.

        If ``ticket`` and ``patchname`` are specified, download the
        patch ``patchname`` attached to ``ticket``.

        If ``url`` is specified, download ``url``.

        Raise an error on any other combination of parameters.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``)

        - ``patchname`` -- a string or ``None`` (default: ``None``)

        - ``url`` -- a string or ``None`` (default: ``None``)

        OUTPUT:

        Otherwise, returns a tuple of the names of the local patch files that
        have been downloaded.

        .. SEEALSO::

            - :meth:`import_patch` -- also creates a commit on the current
              branch from the patch.

        EXAMPLES::

            sage: from sage.env import SAGE_ROOT
            sage: os.chdir(SAGE_ROOT) # silence possible warnings about not being in SAGE_ROOT
            sage: dev.download_patch(ticket=14882,        # optional: internet
            ....:        patchname='trac_14882-backtrack_longtime-dg.patch')
            Downloading "https://trac.sagemath.org/raw-attachment/ticket/14882/trac_14882
            -backtrack_longtime-dg.patch"...
            Downloaded "https://trac.sagemath.org/raw-attachment/ticket/14882/trac_14882
            -backtrack_longtime-dg.patch" to "/...patch".
            ('/...patch',)

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.sagedev import single_user_setup
            sage: dev, config, UI, server = single_user_setup()

        Create a new ticket::

            sage: UI.append("Summary: summary1\ndescription")
            sage: dev.create_ticket()
            Created ticket #1 at https://trac.sagemath.org/1.
            <BLANKLINE>
            #  (use "sage --dev checkout --ticket=1" to create a new local branch)
            1

        There are no attachments to download yet::

            sage: dev.download_patch(ticket=1)     # optional - internet
            Ticket #1 has no attachments.

        After adding one attachment, this works::

            sage: server.tickets[1].attachments['first.patch'] = ''
            sage: dev.download_patch(ticket=1)    # not tested, just an example

        After adding another attachment, this does not work anymore, one needs
        to specify which attachment should be downloaded::

            sage: server.tickets[1].attachments['second.patch'] = ''
            sage: dev.download_patch(ticket=1)   # optional - internet
            I could not understand the comments on ticket #1. To apply use one of the
            patches on the ticket, set the parameter `patchname` to one of: first.patch,
            second.patch
            sage: dev.download_patch(ticket=1,     # not tested, just an example
            ....:     patchname='second.patch')

        It is an error not to specify any parameters if not on a ticket::

            sage: dev.vanilla()
            sage: dev.download_patch()
            ticket or url must be specified if not currently on a ticket

        Check that the parser for the rss stream works::

            sage: UI.append("n")
            sage: dev._sagedev.trac = sage.all.dev.trac
            sage: dev.download_patch(ticket=12415)      # optional: internet
            It seems that the following patches have to be applied in this order:
            12415_spkg_bin_sage.patch
            12415_script.patch
            12415_framework.patch
            12415_doctest_review.patch
            12415_script_review.patch
            12415_review_review.patch
            12415_doc.patch
            12415_test.patch
            12415_review.patch
            12415_review3.patch
            12415_doctest_fixes.patch
            12415_manifest.patch
            12415_rebase_58.patch
            Should I download these patches? [Yes/no] n
            Ticket #12415 has more than one attachment but you chose not to download
            them in the proposed order. To use only one of these patches set the
            parameter `patchname` to one of: 12415_doc.patch, 12415_doctest_fixes.patch,
            12415_doctest_review.patch, 12415_framework.patch, 12415_manifest.patch,
            12415_rebase_58.patch, 12415_review.patch, 12415_review3.patch,
            12415_review_review.patch, 12415_script.patch, 12415_script_review.patch,
            12415_spkg_bin_sage.patch, 12415_test.patch
        """
        if url is not None:
            if ticket or patchname:
                raise ValueError('If "url" is specifed then neither "ticket" nor "patchname"'
                                 ' may be specified.')
            import urllib
            self._UI.show('Downloading "{0}"...', url)
            ret = urllib.urlretrieve(url)[0]
            self._UI.show('Downloaded "{0}" to "{1}".', url, ret)
            return (ret,)

        if ticket is None:
            ticket = self._current_ticket()

        if ticket is None:
            from sagedev import SageDevValueError
            raise SageDevValueError("ticket or url must be specified if not currently on a ticket")

        ticket = self._ticket_from_ticket_name(ticket)

        if patchname:
            from sage.env import TRAC_SERVER_URI
            url = TRAC_SERVER_URI+"/raw-attachment/ticket/%s/%s"%(ticket,patchname)
            if url.startswith("https://"):
                try:
                    import ssl
                except ImportError:
                    # python is not build with ssl support by default. to make
                    # downloading patches work even if ssl is not present, we try
                    # to access trac through http
                    url = url.replace("https","http",1)
            return self.download_patch(url = url)
        else:
            attachments = self.trac.attachment_names(ticket)
            if len(attachments) == 0:
                from sagedev import SageDevValueError
                raise SageDevValueError("Ticket #%s has no attachments."%ticket)
            if len(attachments) == 1:
                ret = self.download_patch(ticket = ticket, patchname = attachments[0])
                self._UI.show('Attachment "{0}" for ticket #{1} has been downloaded to "{2}".'
                              .format(attachments[0], ticket, ret[0]))
                return ret
            else:
                from sage.env import TRAC_SERVER_URI
                rss = TRAC_SERVER_URI+"/ticket/%s?format=rss"%ticket
                self._UI.debug('There is more than one attachment on ticket #{0}. '
                               'Reading "{1}" to try to find out in which order they must be applied.',
                               ticket, rss)
                import urllib2
                rss = urllib2.urlopen(rss).read()

                # the following block has been copied from the patchbot
                all_patches = []
                patches = []
                import re
                folded_regex = re.compile('all.*(folded|combined|merged)')
                subsequent_regex = re.compile('second|third|fourth|next|on top|after')
                attachment_regex = re.compile(r"<strong>attachment</strong>\s*set to <em>(.*)</em>", re.M)
                rebased_regex = re.compile(r"([-.]?rebased?)|(-v\d)")
                def extract_tag(sgml, tag):
                    """
                    Find the first occurance of the tag start (including
                    attributes) and return the contents of that tag (really, up
                    until the next end tag of that type).

                    Crude but fast.
                    """
                    tag_name = tag[1:-1]
                    if ' ' in tag_name:
                        tag_name = tag_name[:tag_name.index(' ')]
                    end = "</%s>" % tag_name
                    start_ix = sgml.find(tag)
                    if start_ix == -1:
                        return None
                    end_ix = sgml.find(end, start_ix)
                    if end_ix == -1:
                        return None
                    return sgml[start_ix + len(tag) : end_ix].strip()
                for item in rss.split('<item>'):
                    description = extract_tag(item, '<description>').replace('&lt;', '<').replace('&gt;', '>')
                    m = attachment_regex.search(description)
                    comments = description[description.find('</ul>') + 1:]
                    # Look for apply... followed by patch names
                    for line in comments.split('\n'):
                        if 'apply' in line.lower():
                            new_patches = []
                            for p in line[line.lower().index('apply') + 5:].split(','):
                                for pp in p.strip().split():
                                    if pp in all_patches:
                                        new_patches.append(pp)
                            if new_patches or (m and not subsequent_regex.search(line)):
                                patches = new_patches
                        elif m and folded_regex.search(line):
                            patches = [] # will add this patch below
                    if m is not None:
                        attachment = m.group(1)
                        import os.path
                        base, ext = os.path.splitext(attachment)
                        if '.' in base:
                            try:
                                base2, ext2 = os.path.splitext(base)
                                count = int(ext2[1:])
                                for i in range(count):
                                    if i:
                                        older = "%s.%s%s" % (base2, i, ext)
                                    else:
                                        older = "%s%s" % (base2, ext)
                                    if older in patches:
                                        patches.remove(older)
                            except Exception:
                                pass
                        if rebased_regex.search(attachment):
                            older = rebased_regex.sub('', attachment)
                            if older in patches:
                                patches.remove(older)
                        if ext in ('.patch', '.diff'):
                            all_patches.append(attachment)
                            patches.append(attachment)

                # now patches contains the list of patches to apply
                if patches:
                    if self._UI.confirm("It seems that the following patches have to be applied in this order: \n{0}\nShould I download these patches?".format("\n".join(patches)),default=True):
                        ret = []
                        for patch in patches:
                            ret.extend(self.download_patch(ticket=ticket, patchname=patch))
                        return ret
                    else:
                        self._UI.error("Ticket #{0} has more than one attachment but you chose not to download them in the proposed order. To use only one of these patches set the parameter `patchname` to one of: {1}".format(ticket, ", ".join(sorted(attachments))))
                        raise OperationCancelledError("user requested")
                else:
                    self._UI.error("I could not understand the comments on ticket #{0}. To apply use one of the patches on the ticket, set the parameter `patchname` to one of: {1}".format(ticket, ", ".join(sorted(attachments))))
                    raise OperationCancelledError("user requested")

    def _detect_patch_diff_format(self, lines):
        r"""
        Determine the format of the ``diff`` lines in ``lines``.

        INPUT:

        - ``lines`` -- a list of strings

        OUTPUT:

        Either ``git`` (for ``diff --git`` lines) or ``hg`` (for ``diff -r`` lines).

        .. NOTE::

            Most Sage developpers have configured mercurial to export
            patches in git format.

        TESTS::

            sage: dev = dev._sagedev
            sage: dev._detect_patch_diff_format(
            ....:     ["diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py"])
            'hg'
            sage: dev._detect_patch_diff_format(
            ....:     ["diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi"])
            'git'

            sage: from sage.env import SAGE_SRC
            sage: dev._detect_patch_diff_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines())
            'git'
            sage: dev._detect_patch_diff_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","diff.patch"
            ....:         )).read().splitlines())
            'hg'

            sage: dev._detect_patch_diff_format(["# HG changeset patch"])
            Traceback (most recent call last):
            ...
            NotImplementedError: Failed to detect diff format.
            sage: dev._detect_patch_diff_format(
            ... ["diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py",
            ...  "diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi"])
            Traceback (most recent call last):
            ...
            SageDevValueError: File appears to have mixed diff formats.
        """
        format = None
        regexs = { "hg" : HG_DIFF_REGEX, "git" : GIT_DIFF_REGEX }

        for line in lines:
            for name,regex in regexs.items():
                if regex.match(line):
                    if format is None:
                        format = name
                    if format != name:
                        from sagedev import SageDevValueError
                        raise SageDevValueError("File appears to have mixed diff formats.")

        if format is None:
            raise NotImplementedError("Failed to detect diff format.")
        else:
            return format

    def _detect_patch_path_format(self, lines, diff_format=None):
        r"""
        Determine the format of the paths in the patch given in ``lines``.

        INPUT:

        - ``lines`` -- a list (or iterable) of strings

        - ``diff_format`` -- ``'hg'``,``'git'``, or ``None`` (default:
          ``None``), the format of the ``diff`` lines in the patch. If
          ``None``, the format will be determined by
          :meth:`_detect_patch_diff_format`.

        OUTPUT:

        A string, ``'new'`` (new repository layout) or ``'old'`` (old
        repository layout).

        EXAMPLES::

            sage: dev._wrap("_detect_patch_path_format")
            sage: dev._detect_patch_path_format(
            ....:     ["diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py"])
            'old'
            sage: dev._detect_patch_path_format(
            ....:     ["diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py"],
            ....:     diff_format="git")
            Traceback (most recent call last):
            ...
            NotImplementedError: Failed to detect path format.
            sage: dev._detect_patch_path_format(
            ....:     ["diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi"])
            'old'
            sage: dev._detect_patch_path_format(
            ....:     ["diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi"])
            'new'
            sage: dev._detect_patch_path_format(
            ....:     ["rename to sage/rings/number_field/totallyreal.pyx"], diff_format='hg')
            'old'
            sage: dev._detect_patch_path_format(
            ....:     ["rename from src/sage/rings/number_field/totalyreal.pyx"], diff_format='git')
            'new'

            sage: from sage.env import SAGE_SRC
            sage: dev._detect_patch_path_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines())
            'old'
        """
        lines = list(lines)
        if diff_format is None:
            diff_format = self._detect_patch_diff_format(lines)

        path_format = None

        if diff_format == "git":
            diff_regexs = (GIT_DIFF_REGEX, PM_DIFF_REGEX, MV_DIFF_REGEX, CP_DIFF_REGEX)
        elif diff_format == "hg":
            diff_regexs = (HG_DIFF_REGEX, PM_DIFF_REGEX, MV_DIFF_REGEX, CP_DIFF_REGEX)
        else:
            raise NotImplementedError(diff_format)

        regexs = { "old" : HG_PATH_REGEX, "new" : GIT_PATH_REGEX }

        for line in lines:
            for regex in diff_regexs:
                match = regex.match(line)
                if match:
                    for group in match.groups():
                        for name, regex in regexs.items():
                            if regex.match(group):
                                if path_format is None:
                                    path_format = name
                                if path_format != name:
                                    from sagedev import SageDevValueError
                                    raise SageDevValueError("File appears to have mixed path formats.")

        if path_format is None:
            raise NotImplementedError("Failed to detect path format.")
        else:
           return path_format

    def _rewrite_patch_diff_paths(self, lines, to_format, from_format=None, diff_format=None):
        r"""
        Rewrite the ``diff`` lines in ``lines`` to use ``to_format``.

        INPUT:

        - ``lines`` -- a list or iterable of strings

        - ``to_format`` -- ``'old'`` or ``'new'``

        - ``from_format`` -- ``'old'``, ``'new'``, or ``None`` (default:
          ``None``), the current formatting of the paths; detected
          automatically if ``None``

        - ``diff_format`` -- ``'git'``, ``'hg'``, or ``None`` (default:
          ``None``), the format of the ``diff`` lines; detected automatically
          if ``None``

        OUTPUT:

        A list of string, ``lines`` rewritten to conform to ``lines``.

        EXAMPLES:

        Paths in the old format::

            sage: dev._wrap("_rewrite_patch_diff_paths")
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py'],
            ....:     to_format="old")
            ['diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi'],
            ....:     to_format="old")
            ['diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['--- a/sage/rings/padics/pow_computer_ext.pxd',
            ....:      '+++ b/sage/rings/padics/pow_computer_ext.pxd'],
            ....:     to_format="old", diff_format="git")
            ['--- a/sage/rings/padics/pow_computer_ext.pxd',
             '+++ b/sage/rings/padics/pow_computer_ext.pxd']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py'],
            ....:     to_format="new")
            ['diff -r 1492e39aff50 -r 5803166c5b11 src/sage/schemes/elliptic_curves/ell_rational_field.py']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi'],
            ....:     to_format="new")
            ['diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['--- a/sage/rings/padics/pow_computer_ext.pxd',
            ....:      '+++ b/sage/rings/padics/pow_computer_ext.pxd'],
            ....:     to_format="new", diff_format="git")
            ['--- a/src/sage/rings/padics/pow_computer_ext.pxd',
             '+++ b/src/sage/rings/padics/pow_computer_ext.pxd']

        Paths in the new format::

            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff -r 1492e39aff50 -r 5803166c5b11 src/sage/schemes/elliptic_curves/ell_rational_field.py'],
            ....:     to_format="old")
            ['diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi'],
            ....:     to_format="old")
            ['diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['--- a/src/sage/rings/padics/pow_computer_ext.pxd',
            ....:      '+++ b/src/sage/rings/padics/pow_computer_ext.pxd'],
            ....:     to_format="old", diff_format="git")
            ['--- a/sage/rings/padics/pow_computer_ext.pxd',
             '+++ b/sage/rings/padics/pow_computer_ext.pxd']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff -r 1492e39aff50 -r 5803166c5b11 src/sage/schemes/elliptic_curves/ell_rational_field.py'],
            ....:     to_format="new")
            ['diff -r 1492e39aff50 -r 5803166c5b11 src/sage/schemes/elliptic_curves/ell_rational_field.py']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi'],
            ....:     to_format="new")
            ['diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['--- a/src/sage/rings/padics/pow_computer_ext.pxd',
            ....:      '+++ b/src/sage/rings/padics/pow_computer_ext.pxd'],
            ....:     to_format="new", diff_format="git")
            ['--- a/src/sage/rings/padics/pow_computer_ext.pxd',
             '+++ b/src/sage/rings/padics/pow_computer_ext.pxd']

            sage: dev._rewrite_patch_diff_paths(
            ....:     ['rename from sage/combinat/crystals/letters.py',
            ....:      'rename to sage/combinat/crystals/letters.pyx'],
            ....:     to_format="new", diff_format="hg")
            ['rename from src/sage/combinat/crystals/letters.py',
             'rename to src/sage/combinat/crystals/letters.pyx']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['rename from src/sage/combinat/crystals/letters.py',
            ....:      'rename to src/sage/combinat/crystals/letters.pyx'],
            ....:     to_format="old", diff_format="git")
            ['rename from sage/combinat/crystals/letters.py',
             'rename to sage/combinat/crystals/letters.pyx']

            sage: from sage.env import SAGE_SRC
            sage: result = dev._rewrite_patch_diff_paths(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines(),
            ....:     to_format="new", diff_format="git")
            sage: len(result)
            2980
            sage: result[0]
            '#8703: Enumerated sets and data structure for ordered and binary trees'
            sage: result[12]
            'diff --git a/src/doc/en/reference/combinat/index.rst b/src/doc/en/reference/combinat/index.rst'
        """
        lines = list(lines)
        if diff_format is None:
            diff_format = self._detect_patch_diff_format(lines)

        if from_format is None:
            from_format = self._detect_patch_path_format(lines, diff_format=diff_format)

        if to_format == from_format:
            return lines

        def hg_path_to_git_path(path):
            if any([path.startswith(p) for p in
                    ("module_list.py", "setup.py", "c_lib/", "sage/", "doc/")]):
                return "src/%s"%path
            else:
                raise NotImplementedError('mapping hg path "%s"'%path)

        def git_path_to_hg_path(path):
            if any([path.startswith(p) for p in
                    ("src/module_list.py", "src/setup.py", "src/c_lib/", "src/sage/", "src/doc/")]):
                return path[4:]
            else:
                raise NotImplementedError('mapping git path "%s"'%path)

        def apply_replacements(lines, diff_regexs, replacement):
            ret = []
            for line in lines:
                for diff_regex in diff_regexs:
                    m = diff_regex.match(line)
                    if m:
                        line = line[:m.start(1)] + \
                               ("".join([ line[m.end(i-1):m.start(i)]+replacement(m.group(i))
                                          for i in range(1,m.lastindex+1) ])) + \
                               line[m.end(m.lastindex):]
                ret.append(line)
            return ret

        diff_regex = None
        if diff_format == "hg":
            diff_regex = (HG_DIFF_REGEX, PM_DIFF_REGEX, MV_DIFF_REGEX, CP_DIFF_REGEX)
        elif diff_format == "git":
            diff_regex = (GIT_DIFF_REGEX, PM_DIFF_REGEX, MV_DIFF_REGEX, CP_DIFF_REGEX)
        else:
            raise NotImplementedError(diff_format)

        if from_format == "old":
            return self._rewrite_patch_diff_paths(
                apply_replacements(lines, diff_regex, hg_path_to_git_path),
                from_format="new", to_format=to_format, diff_format=diff_format)
        elif from_format == "new":
            if to_format == "old":
                return apply_replacements(lines, diff_regex, git_path_to_hg_path)
            else:
                raise NotImplementedError(to_format)
        else:
            raise NotImplementedError(from_format)

    def _detect_patch_header_format(self, lines):
        r"""
        Detect the format of the patch header in ``lines``.

        INPUT:

        - ``lines`` -- a list (or iterable) of strings

        OUTPUT:

        A string, ``'hg-export'`` (mercurial export header), ``'hg'``
        (mercurial header), ``'git'`` (git mailbox header), ``'diff'`` (no
        header)

        EXAMPLES::

            sage: dev._wrap("_detect_patch_header_format")
            sage: dev._detect_patch_header_format(
            ... ['# HG changeset patch','# Parent 05fca316b08fe56c8eec85151d9a6dde6f435d46'])
            'hg'
            sage: dev._detect_patch_header_format(
            ... ['# HG changeset patch','# User foo@bar.com'])
            'hg-export'
            sage: dev._detect_patch_header_format(
            ... ['From: foo@bar'])
            'git'

            sage: from sage.env import SAGE_SRC
            sage: dev._detect_patch_header_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines())
            'diff'
            sage: dev._detect_patch_header_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","diff.patch"
            ....:         )).read().splitlines())
            'diff'
        """
        lines = list(lines)
        if not lines:
            from sagedev import SageDevValueError
            raise SageDevValueError("patch is empty")

        if HG_HEADER_REGEX.match(lines[0]):
            if HG_USER_REGEX.match(lines[1]):
                return "hg-export"
            elif HG_PARENT_REGEX.match(lines[1]):
                return "hg"
        elif GIT_FROM_REGEX.match(lines[0]):
            return "git"
        return "diff"

    def _detect_patch_modified_files(self, lines, diff_format = None):
        r"""
        Return a list of files which are modified by the patch in ``lines``.

        TESTS::

            sage: dev._wrap("_detect_patch_modified_files")
            sage: from sage.env import SAGE_SRC
            sage: dev._detect_patch_modified_files(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines())
            ['ordered_tree.py', 'binary_tree.pyx', 'list_clone.pyx', 'permutation.py',
             'index.rst', 'abstract_tree.py', 'all.py', 'binary_tree.py']
        """
        if diff_format is None:
            diff_format = self._detect_patch_diff_format(lines)

        if diff_format == "hg":
            regex = HG_DIFF_REGEX
        elif diff_format == "git":
            regex = GIT_DIFF_REGEX
        else:
            raise NotImplementedError(diff_format)

        ret = set()
        for line in lines:
            m = regex.match(line)
            if m:
                for group in m.groups():
                    split = group.split('/')
                    if split:
                        ret.add(split[-1])
        return list(ret)

    def _rewrite_patch_header(self, lines, to_format, from_format = None, diff_format = None):
        r"""
        Rewrite ``lines`` to match ``to_format``.

        INPUT:

        - ``lines`` -- a list of strings, the lines of the patch file

        - ``to_format`` -- one of ``'hg'``, ``'hg-export'``, ``'diff'``,
          ``'git'``, the format of the resulting patch file.

        - ``from_format`` -- one of ``None``, ``'hg'``, ``'hg-export'``, ``'diff'``, ``'git'``
          (default: ``None``), the format of the patch file.  The format is
          determined automatically if ``format`` is ``None``.

        OUTPUT:

        A list of lines, in the format specified by ``to_format``.

        Some sample patch files are in data/, in hg and git
        format. Since the translation is not perfect, the resulting
        file is also put there for comparison.

        EXAMPLES::

            sage: from sage.env import SAGE_SRC
            sage: hg_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "test", "data", "hg.patch")
            ....:     ).read().splitlines()
            sage: hg_output_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "test", "data", "hg-output.patch")
            ....:     ).read().splitlines()
            sage: git_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "test", "data", "git.patch")
            ....:     ).read().splitlines()
            sage: git_output_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "test", "data", "git-output.patch")
            ....:     ).read().splitlines()

            sage: dev._wrap("_rewrite_patch_header")
            sage: dev._rewrite_patch_header(git_lines, 'git') == git_lines
            True
            sage: dev._rewrite_patch_header(hg_lines, 'hg-export') == hg_lines
            True

            sage: dev._rewrite_patch_header(git_lines, 'hg-export') == hg_output_lines
            True
            sage: dev._rewrite_patch_header(hg_lines, 'git') == git_output_lines
            True

            sage: dev._rewrite_patch_header(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines(), 'git')[:5]
            ['From: "Unknown User" <unknown@sagemath.org>',
            'Subject: #8703: Enumerated sets and data structure for ordered and binary trees',
            'Date: ...',
            '',
            '- The Class Abstract[Labelled]Tree allows for inheritance from different']
        """
        import email.utils, time

        lines = list(lines)
        if not lines:
            from sagedev import SageDevValueError
            raise SageDevValueError("empty patch file")

        if from_format is None:
            from_format = self._detect_patch_header_format(lines)

        if from_format == to_format:
            return lines

        def parse_header(lines, regexs, mandatory=False):
            header = {}
            i = 0
            for (key, regex) in regexs:
                if i > len(lines):
                    if mandatory:
                        from sagedev import SageDevValueError
                        raise SageDevValueError('Malformed patch. Missing line for regular expression "%s".'
                                                %(regex.pattern))
                    else:
                        return
                match = regex.match(lines[i])
                if match is not None:
                    if len(match.groups()) > 0:
                        header[key] = match.groups()[0]
                    i += 1
                elif mandatory:
                    from sagedev import SageDevValueError
                    raise SageDevValueError('Malformed patch. Line "%s" does not match regular expression "%s".'
                                            %(lines[i],regex.pattern))

            message = []
            for i in range(i,len(lines)):
                if lines[i].startswith("diff -"):
                    break
                else:
                    message.append(lines[i])

            header["message"] = message
            return header, lines[i:]

        if from_format == "git":
            header, diff = parse_header(lines,
                    (("user", GIT_FROM_REGEX), ("subject", GIT_SUBJECT_REGEX), ("date", GIT_DATE_REGEX)),
                    mandatory=True)

            if to_format == "hg-export":
                ret = []
                ret.append('# HG changeset patch')
                ret.append('# User %s'%(header["user"]))
                old_TZ = os.environ.get('TZ')
                try:
                    os.environ['TZ'] = 'UTC'
                    time.tzset()
                    t = time.mktime(email.utils.parsedate(header["date"]))
                    ret.append('# Date %s 00000'%int(t)) # this is not portable
                finally:
                    if old_TZ:
                        os.environ['TZ'] = old_TZ
                    else:
                        del os.environ['TZ']
                    time.tzset()
                ret.append('# Node ID 0000000000000000000000000000000000000000')
                ret.append('# Parent  0000000000000000000000000000000000000000')
                ret.append(header["subject"])
                ret.extend(header["message"])
                ret.extend(diff)
                return ret
            else:
                raise NotImplementedError(to_format)
        elif from_format in ["hg", "diff", "hg-export"]:
            header, diff = parse_header(lines,
                                        (("hg_header", HG_HEADER_REGEX),
                                         ("user", HG_USER_REGEX),
                                         ("date", HG_DATE_REGEX),
                                         ("node", HG_NODE_REGEX),
                                         ("parent", HG_PARENT_REGEX)))
            user    = header.get("user", '"Unknown User" <unknown@sagemath.org>')
            date    = email.utils.formatdate(int(header.get("date", time.time())))
            message = header.get("message", [])
            if message:
                subject = message[0]
                message = message[1:]
            else:
                subject = 'No Subject. Modified: %s'%(", ".join(
                    sorted(self._detect_patch_modified_files(lines))))
            ret = []
            ret.append('From: %s'%user)
            ret.append('Subject: %s'%subject)
            ret.append('Date: %s'%date)
            ret.append('')
            if message and message != ['']: # avoid a double empty line
                ret.extend(message)
            ret.extend(diff)
            return self._rewrite_patch_header(ret, to_format=to_format, from_format="git",
                                              diff_format=diff_format)
        else:
            raise NotImplementedError(from_format)

    def _rewrite_patch(self, lines, to_path_format, to_header_format, from_diff_format=None,
                       from_path_format=None, from_header_format=None):
        r"""
        Rewrite the patch in ``lines`` to the path format given in
        ``to_path_format`` and the header format given in ``to_header_format``.

        TESTS::

            sage: dev._wrap("_rewrite_patch")
            sage: from sage.env import SAGE_SRC
            sage: git_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "test", "data", "git.patch")
            ....:     ).read().splitlines()
            sage: dev._rewrite_patch(git_lines, "old", "git") == git_lines
            True
        """
        return self._rewrite_patch_diff_paths(
            self._rewrite_patch_header(lines,
                                       to_format=to_header_format,
                                       from_format=from_header_format,
                                       diff_format=from_diff_format),
            to_format=to_path_format,
            diff_format=from_diff_format,
            from_format=from_path_format)
