"""
The documentation builder

It is the starting point for building documentation, and is
responsible to figure out what to build and with which options. The
actual documentation build for each individual document is then done
in a subprocess call to sphinx, see :func:`builder_helper`.

* The builder can be configured in build_options.py
* The sphinx subprocesses are configured in conf.py
"""

from __future__ import absolute_import
from __future__ import print_function
from six.moves import range

import logging, optparse, os, shutil, subprocess, sys, re

import sphinx.cmdline
import sphinx.util.console
import sphinx.ext.intersphinx

import sage.all
from sage.misc.cachefunc import cached_method
from sage.misc.misc import sage_makedirs
from sage.env import SAGE_DOC_SRC, SAGE_DOC, SAGE_SRC

from .build_options import (LANGUAGES, SPHINXOPTS, PAPER, OMIT,
     PAPEROPTS, ALLSPHINXOPTS, NUM_THREADS, WEBSITESPHINXOPTS,
     INCREMENTAL_BUILD, ABORT_ON_ERROR)


def excepthook(*exc_info):
    """
    Print docbuild error and hint how to solve it
    """
    logger.error('Error building the documentation.', exc_info=exc_info)
    if INCREMENTAL_BUILD:
        logger.error('''
Note: incremental documentation builds sometimes cause spurious
error messages. To be certain that these are real errors, run
"make doc-clean" first and try again.''')


def delete_empty_directories(root_dir):
    """
    Delete all empty directories found under ``root_dir``

    INPUT:

    - ``root_dir`` -- string. A valid directory name.
    """
    for dirpath, dirnames, filenames in os.walk(root_dir, topdown=False):
        if not dirnames + filenames:
            logger.warning('Deleting empty directory {0}'.format(dirpath))
            os.rmdir(dirpath)


##########################################
#      Parallel Building Ref Manual      #
##########################################
def build_ref_doc(args):
    doc = args[0]
    lang = args[1]
    format = args[2]
    kwds = args[3]
    args = args[4:]
    if format == 'inventory':  # you must not use the inventory to build the inventory
        kwds['use_multidoc_inventory'] = False
    getattr(ReferenceSubBuilder(doc, lang), format)(*args, **kwds)

##########################################
#             Builders                   #
##########################################

def builder_helper(type):
    """
    Returns a function which builds the documentation for
    output type type.
    """
    def f(self, *args, **kwds):
        output_dir = self._output_dir(type)

        options = ALLSPHINXOPTS

        if self.name == 'website':
            # WEBSITESPHINXOPTS is either empty or " -A hide_pdf_links=1 "
            options += WEBSITESPHINXOPTS

        if kwds.get('use_multidoc_inventory', True):
            options += ' -D multidoc_first_pass=0'
        else:
            options += ' -D multidoc_first_pass=1'

        build_command = '-b %s -d %s %s %s %s'%(type, self._doctrees_dir(),
                                                  options, self.dir,
                                                  output_dir)
        logger.debug(build_command)

        # Run Sphinx with Sage's special logger
        sys.argv = ["sphinx-build"] + build_command.split()
        from .sphinxbuild import runsphinx
        try:
            runsphinx()
        except Exception:
            if ABORT_ON_ERROR:
                raise

        if "/latex" in output_dir:
            logger.warning("LaTeX file written to {}".format(output_dir))
        else:
            logger.warning(
                "Build finished.  The built documents can be found in {}".
                format(output_dir))

    f.is_output_format = True
    return f


class DocBuilder(object):
    def __init__(self, name, lang='en'):
        """
        INPUT:

        - ``name`` - the name of a subdirectory in SAGE_DOC_SRC, such as
          'tutorial' or 'bordeaux_2008'

        - ``lang`` - (default "en") the language of the document.
        """
        doc = name.split(os.path.sep)

        if doc[0] in LANGUAGES:
            lang = doc[0]
            doc.pop(0)

        self.name = os.path.join(*doc)
        self.lang = lang
        self.dir = os.path.join(SAGE_DOC_SRC, self.lang, self.name)

    def _output_dir(self, type):
        """
        Returns the directory where the output of type ``type`` is stored.
        If the directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: from sage_setup.docbuild import DocBuilder
            sage: b = DocBuilder('tutorial')
            sage: b._output_dir('html')
            '.../html/en/tutorial'
        """
        d = os.path.join(SAGE_DOC, type, self.lang, self.name)
        sage_makedirs(d)
        return d

    def _doctrees_dir(self):
        """
        Returns the directory where the doctrees are stored.  If the
        directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: from sage_setup.docbuild import DocBuilder
            sage: b = DocBuilder('tutorial')
            sage: b._doctrees_dir()
            '.../doctrees/en/tutorial'
        """
        d = os.path.join(SAGE_DOC, 'doctrees', self.lang, self.name)
        sage_makedirs(d)
        return d

    def _output_formats(self):
        """
        Returns a list of the possible output formats.

        EXAMPLES::

            sage: from sage_setup.docbuild import DocBuilder
            sage: b = DocBuilder('tutorial')
            sage: b._output_formats()
            ['changes', 'html', 'htmlhelp', 'inventory', 'json', 'latex', 'linkcheck', 'pickle', 'web']

        """
        #Go through all the attributes of self and check to
        #see which ones have an 'is_output_format' attribute.  These
        #are the ones created with builder_helper.
        output_formats = []
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'is_output_format'):
                output_formats.append(attr)
        output_formats.sort()
        return output_formats

    def pdf(self):
        """
        Builds the PDF files for this document.  This is done by first
        (re)-building the LaTeX output, going into that LaTeX
        directory, and running 'make all-pdf' (or for the special case of
        the ja docs, 'all-pdf-ja(ex,to run platex)' there.

        EXAMPLES::

            sage: from sage_setup.docbuild import DocBuilder
            sage: b = DocBuilder('tutorial')
            sage: b.pdf() #not tested
        """
        self.latex()
        tex_dir = self._output_dir('latex')
        pdf_dir = self._output_dir('pdf')
        make_target = "cd '%s' && $MAKE %s && mv -f *.pdf '%s'"
        error_message = "failed to run $MAKE %s in %s"
        MB_LANG = {'ja': 'all-pdf-ja'} # language name : the modified target

        # Replace the command for languages that require special processing
        if self.lang in MB_LANG:
            command = MB_LANG[self.lang]
        else:
            command = 'all-pdf'

        if subprocess.call(make_target%(tex_dir, command, pdf_dir), shell=True):
            raise RuntimeError(error_message%(command, tex_dir))
        logger.warning("Build finished.  The built documents can be found in %s", pdf_dir)

    def clean(self, *args):
        shutil.rmtree(self._doctrees_dir())
        output_formats = list(args) if args else self._output_formats()
        for format in output_formats:
            shutil.rmtree(self._output_dir(format), ignore_errors=True)

    html = builder_helper('html')
    pickle = builder_helper('pickle')
    web = pickle
    json = builder_helper('json')
    htmlhelp = builder_helper('htmlhelp')
    latex = builder_helper('latex')
    changes = builder_helper('changes')
    linkcheck = builder_helper('linkcheck')
    # import the customized builder for object.inv files
    inventory = builder_helper('inventory')


if NUM_THREADS > 1:
    def build_many(target, args):
        from multiprocessing import Pool
        pool = Pool(NUM_THREADS, maxtasksperchild=1)
        # map_async handles KeyboardInterrupt correctly. Plain map and
        # apply_async does not, so don't use it.
        x = pool.map_async(target, args, 1)
        try:
            ret = x.get(99999)
            pool.close()
            pool.join()
        except Exception:
            pool.terminate()
            if ABORT_ON_ERROR:
                raise
        return ret
else:
    def build_many(target, args):
        results = []

        for arg in args:
            try:
                results.append(target(arg))
            except Exception:
                if ABORT_ON_ERROR:
                    raise

        return results


##########################################
#      Parallel Building Ref Manual      #
##########################################
def build_other_doc(args):
    document = args[0]
    name = args[1]
    kwds = args[2]
    args = args[3:]
    logger.warning("\nBuilding %s.\n" % document)
    getattr(get_builder(document), name)(*args, **kwds)

class AllBuilder(object):
    """
    A class used to build all of the documentation.
    """
    def __getattr__(self, attr):
        """
        For any attributes not explicitly defined, we just go through
        all of the documents and call their attr.  For example,
        'AllBuilder().json()' will go through all of the documents
        and call the json() method on their builders.
        """
        from functools import partial
        return partial(self._wrapper, attr)

    def _wrapper(self, name, *args, **kwds):
        """
        This is the function which goes through all of the documents
        and does the actual building.
        """
        import time
        start = time.time()
        docs = self.get_all_documents()
        refs = [x for x in docs if x.endswith('reference')]
        others = [x for x in docs if not x.endswith('reference')]

        # Build the reference manual twice to resolve references.  That is,
        # build once with the inventory builder to construct the intersphinx
        # inventory files, and then build the second time for real.  So the
        # first build should be as fast as possible;
        logger.warning("\nBuilding reference manual, first pass.\n")
        for document in refs:
            getattr(get_builder(document), 'inventory')(*args, **kwds)

        logger.warning("Building reference manual, second pass.\n")
        for document in refs:
            getattr(get_builder(document), name)(*args, **kwds)

        # build the other documents in parallel
        L = [(doc, name, kwds) + args for doc in others]
        build_many(build_other_doc, L)
        logger.warning("Elapsed time: %.1f seconds."%(time.time()-start))
        logger.warning("Done building the documentation!")

    def get_all_documents(self):
        """
        Returns a list of all of the documents. A document is a directory within one of
        the language subdirectories of SAGE_DOC_SRC specified by the global LANGUAGES
        variable.

        EXAMPLES::

            sage: from sage_setup.docbuild import AllBuilder
            sage: documents = AllBuilder().get_all_documents()
            sage: 'en/tutorial' in documents
            True
            sage: documents[0] == 'en/reference'
            True
        """
        documents = []
        for lang in LANGUAGES:
            for document in os.listdir(os.path.join(SAGE_DOC_SRC, lang)):
                if (document not in OMIT
                    and os.path.isdir(os.path.join(SAGE_DOC_SRC, lang, document))):
                    documents.append(os.path.join(lang, document))

        # Ensure that the reference guide is compiled first so that links from
        # the other documents to it are correctly resolved.
        if 'en/reference' in documents:
            documents.remove('en/reference')
        documents.insert(0, 'en/reference')

        return documents


class WebsiteBuilder(DocBuilder):
    def html(self):
        """
        After we've finished building the website index page, we copy
        everything one directory up.  Then we call
        :meth:`create_html_redirects`.
        """
        DocBuilder.html(self)
        html_output_dir = self._output_dir('html')
        for f in os.listdir(html_output_dir):
            src = os.path.join(html_output_dir, f)
            dst = os.path.join(html_output_dir, '..', f)
            if os.path.isdir(src):
                shutil.rmtree(dst, ignore_errors=True)
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)
        self.create_html_redirects()

    def create_html_redirects(self):
        """
        Writes a number of small HTML files; these are files which
        used to contain the main content of the reference manual
        before before splitting the manual into multiple
        documents. After the split, those files have moved, so in each
        old location, write a file which redirects to the new version.
        (This is so old URLs to pieces of the reference manual still
        open the correct files.)
        """
        # The simple html template which will cause a redirect to the
        # correct file
        html_template = """<html><head>
            <meta HTTP-EQUIV="REFRESH" content="0; url=%s">
            </head><body></body></html>"""

        reference_dir = os.path.abspath(os.path.join(self._output_dir('html'),
                                                     '..', 'reference'))
        reference_builder = ReferenceBuilder('reference')
        refdir = os.path.join(os.environ['SAGE_DOC_SRC'], 'en', 'reference')
        for document in reference_builder.get_all_documents(refdir):
            #path is the directory above reference dir
            path = os.path.abspath(os.path.join(reference_dir, '..'))

            # the name of the subdocument
            document_name = document.split('/')[1]

            # the sage directory within a subdocument, for example
            # local/share/doc/sage/html/en/reference/algebras/sage
            sage_directory = os.path.join(path, document, 'sage')

            # Walk through all of the files in the sage_directory
            for dirpath, dirnames, filenames in os.walk(sage_directory):
                # a string like reference/algebras/sage/algebras
                short_path = dirpath[len(path)+1:]

                # a string like sage/algebras
                shorter_path = os.path.join(*short_path.split(os.sep)[2:])

                #Make the shorter path directory
                try:
                    os.makedirs(os.path.join(reference_dir, shorter_path))
                except OSError:
                    pass

                for filename in filenames:
                    if not filename.endswith('html'):
                        continue

                    # the name of the html file we are going to create
                    redirect_filename = os.path.join(reference_dir, shorter_path, filename)

                    # the number of levels up we need to use in the relative url
                    levels_up = len(shorter_path.split(os.sep))

                    # the relative url that we will redirect to
                    redirect_url = "/".join(['..']*levels_up + [document_name, shorter_path, filename])

                    # write the html file which performs the redirect
                    with open(redirect_filename, 'w') as f:
                        f.write(html_template % redirect_url)


    def clean(self):
        """
        When we clean the output for the website index, we need to
        remove all of the HTML that were placed in the parent
        directory.
        """
        html_output_dir = self._output_dir('html')
        parent_dir = os.path.realpath(os.path.join(html_output_dir, '..'))
        for filename in os.listdir(html_output_dir):
            parent_filename = os.path.join(parent_dir, filename)
            if not os.path.exists(parent_filename):
                continue
            if os.path.isdir(parent_filename):
                shutil.rmtree(parent_filename, ignore_errors=True)
            else:
                os.unlink(parent_filename)

        DocBuilder.clean(self)


class ReferenceBuilder(AllBuilder):
    """
    This class builds the reference manual.  It uses DocBuilder to
    build the top-level page and ReferenceSubBuilder for each
    sub-component.
    """
    def __init__(self, name, lang='en'):
        """
        Records the reference manual's name, in case it's not
        identical to 'reference'.
        """
        AllBuilder.__init__(self)
        doc = name.split(os.path.sep)

        if doc[0] in LANGUAGES:
            lang = doc[0]
            doc.pop(0)

        self.name = doc[0]
        self.lang = lang

    def _output_dir(self, type, lang='en'):
        """
        Returns the directory where the output of type type is stored.
        If the directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: from sage_setup.docbuild import ReferenceBuilder
            sage: b = ReferenceBuilder('reference')
            sage: b._output_dir('html')
            '.../html/en/reference'
        """
        d = os.path.join(SAGE_DOC, type, lang, self.name)
        sage_makedirs(d)
        return d

    def _wrapper(self, format, *args, **kwds):
        """
        Builds reference manuals.  For each language, it builds the
        top-level document and its components.
        """
        for lang in LANGUAGES:
            refdir = os.path.join(SAGE_DOC_SRC, lang, self.name)
            if not os.path.exists(refdir):
                continue
            output_dir = self._output_dir(format, lang)
            L = [(doc, lang, format, kwds) + args for doc in self.get_all_documents(refdir)]
            build_many(build_ref_doc, L)
            # The html refman must be build at the end to ensure correct
            # merging of indexes and inventories.
            # Sphinx is run here in the current process (not in a
            # subprocess) and the IntersphinxCache gets populated to be
            # used for the second pass of the reference manual and for
            # the other documents.
            getattr(DocBuilder(self.name, lang), format)(*args, **kwds)

            # PDF: we need to build master index file which lists all
            # of the PDF file.  So we create an html file, based on
            # the file index.html from the "website" target.
            if format == 'pdf':
                # First build the website page.  (This only takes a
                # few seconds.)
                getattr(get_builder('website'), 'html')()
                # Copy the relevant pieces of
                # SAGE_DOC/html/en/website/_static to output_dir.
                # (Don't copy all of _static to save some space: we
                # don't need all of the MathJax stuff, and in
                # particular we don't need the fonts.)
                website_dir = os.path.join(SAGE_DOC, 'html',
                                           'en', 'website')
                static_files = ['COPYING.txt', 'basic.css', 'blank.gif',
                         'default.css', 'doctools.js', 'favicon.ico',
                         'file.png', 'jquery.js', 'minus.png',
                         'pdf.png', 'plus.png', 'pygments.css',
                         'sage.css', 'sageicon.png', 'sagelogo.png',
                         'searchtools.js', 'sidebar.js', 'underscore.js']
                sage_makedirs(os.path.join(output_dir, '_static'))
                for f in static_files:
                    try:
                        shutil.copyfile(os.path.join(website_dir, '_static', f),
                                        os.path.join(output_dir, '_static', f))
                    except IOError: # original file does not exist
                        pass
                # Now modify website's index.html page and write it
                # to output_dir.
                f = open(os.path.join(website_dir, 'index.html'))
                html = f.read().replace('Documentation', 'Reference')
                f.close()
                html_output_dir = os.path.dirname(website_dir)
                html = html.replace('http://www.sagemath.org',
                                    os.path.join(html_output_dir, 'index.html'))
                # From index.html, we want the preamble and the tail.
                html_end_preamble = html.find('<h1>Sage Reference')
                html_bottom = html.rfind('</table>') + len('</table>')
                # For the content, we modify doc/en/reference/index.rst,
                # which has two parts: the body and the table of contents.
                f = open(os.path.join(SAGE_DOC_SRC, lang, 'reference', 'index.rst'))
                rst = f.read()
                f.close()
                # Replace rst links with html links.  There are two forms:
                #
                #   `blah`__    followed by __ LINK
                #
                #   :doc:`blah <module/index>`
                #
                # Change the first form to
                #
                #   <a href="LINK">blah</a>
                #
                # Change the second form to
                #
                #   <a href="module/module.pdf">blah <img src="_static/pdf.png" /></a>
                #
                rst = re.sub('`([^`]*)`__\.\n\n__ (.*)',
                                  r'<a href="\2">\1</a>.', rst)
                rst = re.sub(r':doc:`([^<]*?)\s+<(.*)/index>`',
                             r'<a href="\2/\2.pdf">\1 <img src="_static/pdf.png" /></a>',
                             rst)
                # Get rid of todolist and miscellaneous rst markup.
                rst = rst.replace('.. toctree::', '')
                rst = rst.replace(':maxdepth: 2', '')
                rst = rst.replace('todolist', '')
                start = rst.find('=\n') + 1
                end = rst.find('Table of Contents')
                # Body: add paragraph <p> markup.
                rst_body = rst[start:end]
                rst_body = rst_body.replace('\n\n', '</p>\n<p>')
                start = rst.find('Table of Contents') + 2*len('Table of Contents') + 1
                # Don't include the indices.
                end = rst.find('Indices and Tables')
                # TOC: change * to <li>, change rst headers to html headers.
                rst_toc = rst[start:end]
                rst_toc = rst_toc.replace('*', '<li>')
                rst_toc = re.sub('\n([A-Z][a-zA-Z, ]*)\n-*\n',
                             '</ul>\n\n\n<h2>\\1</h2>\n\n<ul>\n', rst_toc)
                # Now write the file.
                new_index = open(os.path.join(output_dir, 'index.html'), 'w')
                new_index.write(html[:html_end_preamble])
                new_index.write('<h1>' + rst[:rst.find('\n')] +
                                ' (PDF version)'+ '</h1>')
                new_index.write(rst_body)
                new_index.write('<h2>Table of Contents</h2>\n\n<ul>')
                new_index.write(rst_toc)
                new_index.write('</ul>\n\n')
                new_index.write(html[html_bottom:])
                new_index.close()
                logger.warning('''
PDF documents have been created in subdirectories of

  %s

Alternatively, you can open

  %s

for a webpage listing all of the documents.''' % (output_dir,
                                                 os.path.join(output_dir,
                                                              'index.html')))

    def get_all_documents(self, refdir):
        """
        Returns a list of all reference manual components to build.
        We add a component name if it's a subdirectory of the manual's
        directory and contains a file named 'index.rst'.

        We return the largest component (most subdirectory entries)
        first since they will take the longest to build.

        EXAMPLES::

            sage: from sage_setup.docbuild import ReferenceBuilder
            sage: b = ReferenceBuilder('reference')
            sage: refdir = os.path.join(os.environ['SAGE_DOC_SRC'], 'en', b.name)
            sage: sorted(b.get_all_documents(refdir))
            ['reference/algebras',
             'reference/arithgroup',
             ...,
             'reference/tensor_free_modules']
        """
        documents = []

        for doc in os.listdir(refdir):
            directory = os.path.join(refdir, doc)
            if os.path.exists(os.path.join(directory, 'index.rst')):
                n = len(os.listdir(directory))
                documents.append((-n, os.path.join(self.name, doc)))

        return [ doc[1] for doc in sorted(documents) ]


class ReferenceSubBuilder(DocBuilder):
    """
    This class builds sub-components of the reference manual.  It is
    resposible for making sure the auto generated ReST files for the
    Sage library are up to date.

    When building any output, we must first go through and check
    to see if we need to update any of the autogenerated ReST
    files.  There are two cases where this would happen:

    1. A new module gets added to one of the toctrees.

    2. The actual module gets updated and possibly contains a new
       title.
    """
    def __init__(self, *args, **kwds):
        DocBuilder.__init__(self, *args, **kwds)
        self._wrap_builder_helpers()

    def _wrap_builder_helpers(self):
        from functools import partial, update_wrapper
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'is_output_format'):
                f = partial(self._wrapper, attr)
                f.is_output_format = True
                update_wrapper(f, getattr(self, attr))
                setattr(self, attr, f)

    def _wrapper(self, build_type, *args, **kwds):
        """
        This is the wrapper around the builder_helper methods that
        goes through and makes sure things are up to date.
        """
        # Delete the auto-generated .rst files, if the inherited
        # and/or underscored members options have changed.
        global options
        inherit_prev = self.get_cache().get('option_inherited')
        underscore_prev = self.get_cache().get('option_underscore')
        if (inherit_prev is None or inherit_prev != options.inherited or
            underscore_prev is None or underscore_prev != options.underscore):
            logger.info("Detected change(s) in inherited and/or underscored members option(s).")
            self.clean_auto()
            self.get_cache.clear_cache()

        # After "sage -clone", refresh the .rst file mtimes in
        # environment.pickle.
        if options.update_mtimes:
            logger.info("Checking for .rst file mtimes to update...")
            self.update_mtimes()

        #Update the .rst files for modified Python modules
        logger.info("Updating .rst files with modified modules...")
        for module_name in self.get_modified_modules():
            self.write_auto_rest_file(module_name.replace(os.path.sep, '.'))

        #Write the .rst files for newly included modules
        logger.info("Writing .rst files for newly-included modules...")
        for module_name in self.get_newly_included_modules(save=True):
            self.write_auto_rest_file(module_name)

        #Copy over the custom .rst files from _sage
        _sage = os.path.join(self.dir, '_sage')
        if os.path.exists(_sage):
            logger.info("Copying over custom .rst files from %s ...", _sage)
            shutil.copytree(_sage, os.path.join(self.dir, 'sage'))

        getattr(DocBuilder, build_type)(self, *args, **kwds)

    def cache_filename(self):
        """
        Returns the filename where the pickle of the dictionary of
        already generated ReST files is stored.
        """
        return os.path.join(self._doctrees_dir(), 'reference.pickle')

    @cached_method
    def get_cache(self):
        """
        Retrieve the cache of already generated ReST files.  If it
        doesn't exist, then we just return an empty dictionary.  If it
        is corrupted, return an empty dictionary.
        """
        filename = self.cache_filename()
        if not os.path.exists(filename):
            return {}
        from six.moves import cPickle
        file = open(self.cache_filename(), 'rb')
        try:
            cache = cPickle.load(file)
        except Exception:
            logger.debug("Cache file '%s' is corrupted; ignoring it..."% filename)
            cache = {}
        else:
            logger.debug("Loaded .rst file cache: %s", filename)
        finally:
            file.close()
        return cache

    def save_cache(self):
        """
        Save the cache of already generated ReST files.
        """
        cache = self.get_cache()

        global options
        cache['option_inherited'] = options.inherited
        cache['option_underscore'] = options.underscore

        from six.moves import cPickle
        file = open(self.cache_filename(), 'wb')
        cPickle.dump(cache, file)
        file.close()
        logger.debug("Saved .rst file cache: %s", self.cache_filename())

    def get_sphinx_environment(self):
        """
        Returns the Sphinx environment for this project.
        """
        from sphinx.environment import BuildEnvironment
        class Foo(object):
            pass
        config = Foo()
        config.values = []

        env_pickle = os.path.join(self._doctrees_dir(), 'environment.pickle')
        try:
            env = BuildEnvironment.frompickle(self.dir, config, env_pickle)
            logger.debug("Opened Sphinx environment: %s", env_pickle)
            return env
        except IOError as err:
            logger.debug("Failed to open Sphinx environment: %s", err)
            pass

    def update_mtimes(self):
        """
        Updates the modification times for ReST files in the Sphinx
        environment for this project.
        """
        env = self.get_sphinx_environment()
        if env is not None:
            import time
            for doc in env.all_docs:
                env.all_docs[doc] = time.time()
            logger.info("Updated %d .rst file mtimes", len(env.all_docs))
            # This is the only place we need to save (as opposed to
            # load) Sphinx's pickle, so we do it right here.
            env_pickle = os.path.join(self._doctrees_dir(),
                                      'environment.pickle')

            # When cloning a new branch (see
            # SAGE_LOCAL/bin/sage-clone), we hard link the doc output.
            # To avoid making unlinked, potentially inconsistent
            # copies of the environment, we *don't* use
            # env.topickle(env_pickle), which first writes a temporary
            # file.  We adapt sphinx.environment's
            # BuildEnvironment.topickle:
            from six.moves import cPickle
            import types

            # remove unpicklable attributes
            env.set_warnfunc(None)
            del env.config.values
            picklefile = open(env_pickle, 'wb')
            # remove potentially pickling-problematic values from config
            for key, val in vars(env.config).items():
                if key.startswith('_') or isinstance(val, (types.ModuleType,
                                                           types.FunctionType,
                                                           type)):
                    del env.config[key]
            try:
                cPickle.dump(env, picklefile, cPickle.HIGHEST_PROTOCOL)
            finally:
                picklefile.close()

            logger.debug("Saved Sphinx environment: %s", env_pickle)

    def get_modified_modules(self):
        """
        Returns an iterator for all the modules that have been modified
        since the documentation was last built.
        """
        env = self.get_sphinx_environment()
        if env is None:
            logger.debug("Stopped check for modified modules.")
            return
        try:
            added, changed, removed = env.get_outdated_files(False)
            logger.info("Sphinx found %d modified modules", len(changed))
        except OSError as err:
            logger.debug("Sphinx failed to determine modified modules: %s", err)
            self.clean_auto()
            return
        for name in changed:
            # Only pay attention to files in a directory sage/... In
            # particular, don't treat a file like 'sagetex.rst' in
            # doc/en/reference/misc as an autogenerated file: see
            # #14199.
            if name.startswith('sage' + os.sep):
                yield name

    def print_modified_modules(self):
        """
        Prints a list of all the modules that have been modified since
        the documentation was last built.
        """
        for module_name in self.get_modified_modules():
            print(module_name)

    def get_all_rst_files(self, exclude_sage=True):
        """
        Returns an iterator for all rst files which are not
        autogenerated.
        """
        for directory, subdirs, files in os.walk(self.dir):
            if exclude_sage and directory.startswith(os.path.join(self.dir, 'sage')):
                continue
            for filename in files:
                if not filename.endswith('.rst'):
                    continue
                yield os.path.join(directory, filename)

    def get_all_included_modules(self):
        """
        Returns an iterator for all modules which are included in the
        reference manual.
        """
        for filename in self.get_all_rst_files():
            for module in self.get_modules(filename):
                yield module

    def get_newly_included_modules(self, save=False):
        """
        Returns an iterator for all modules that appear in the
        toctrees that don't appear in the cache.
        """
        cache = self.get_cache()
        new_modules = 0
        for module in self.get_all_included_modules():
            if module not in cache:
                cache[module] = True
                new_modules += 1
                yield module
        logger.info("Found %d newly included modules", new_modules)
        if save:
            self.save_cache()

    def print_newly_included_modules(self):
        """
        Prints all of the modules that appear in the toctrees that
        don't appear in the cache.
        """
        for module_name in self.get_newly_included_modules():
            print(module_name)

    def get_modules(self, filename):
        """
        Given a filename for a ReST file, return an iterator for
        all of the autogenerated ReST files that it includes.
        """
        #Create the regular expression used to detect an autogenerated file
        auto_re = re.compile('^\s*(..\/)*(sage(nb)?\/[\w\/]*)\s*$')

        #Read the lines
        f = open(filename)
        lines = f.readlines()
        f.close()

        for line in lines:
            match = auto_re.match(line)
            if match:
                yield match.group(2).replace(os.path.sep, '.')

    def get_module_docstring_title(self, module_name):
        """
        Returns the title of the module from its docstring.
        """
        #Try to import the module
        try:
            __import__(module_name)
        except ImportError as err:
            logger.error("Warning: Could not import %s %s", module_name, err)
            return "UNABLE TO IMPORT MODULE"
        module = sys.modules[module_name]

        #Get the docstring
        doc = module.__doc__
        if doc is None:
            doc = module.doc if hasattr(module, 'doc') else ""

        #Extract the title
        i = doc.find('\n')
        if i != -1:
            return doc[i+1:].lstrip().splitlines()[0]
        else:
            return doc

    def auto_rest_filename(self, module_name):
        """
        Returns the name of the file associated to a given module

        EXAMPLES::

            sage: from sage_setup.docbuild import ReferenceSubBuilder
            sage: ReferenceSubBuilder("reference").auto_rest_filename("sage.combinat.partition")
            '.../doc/en/reference/sage/combinat/partition.rst'
        """
        return self.dir + os.path.sep + module_name.replace('.',os.path.sep) + '.rst'

    def write_auto_rest_file(self, module_name):
        """
        Writes the autogenerated ReST file for module_name.
        """
        if not module_name.startswith('sage'):
            return
        filename = self.auto_rest_filename(module_name)
        sage_makedirs(os.path.dirname(filename))

        outfile = open(filename, 'w')

        title = self.get_module_docstring_title(module_name)

        if title == '':
            logger.error("Warning: Missing title for %s", module_name)
            title = "MISSING TITLE"

        # Don't doctest the autogenerated file.
        outfile.write(".. nodoctest\n\n")
        # Now write the actual content.
        outfile.write(".. _%s:\n\n"%(module_name.replace(".__init__","")))
        outfile.write(title + '\n')
        outfile.write('='*len(title) + "\n\n")
        outfile.write('.. This file has been autogenerated.\n\n')

        global options
        inherited = ':inherited-members:' if options.inherited else ''

        automodule = '''
.. automodule:: %s
   :members:
   :undoc-members:
   :show-inheritance:
   %s

'''
        outfile.write(automodule % (module_name, inherited))

        outfile.close()

    def clean_auto(self):
        """
        Remove the cache file for the autogenerated files as well as
        the files themselves.
        """
        if os.path.exists(self.cache_filename()):
            os.unlink(self.cache_filename())
            logger.debug("Deleted .rst cache file: %s", self.cache_filename())

        import shutil
        try:
            shutil.rmtree(os.path.join(self.dir, 'sage'))
            logger.debug("Deleted auto-generated .rst files in: %s",
                         os.path.join(self.dir, 'sage'))
        except OSError:
            pass

    def get_unincluded_modules(self):
        """
        Returns an iterator for all the modules in the Sage library
        which are not included in the reference manual.
        """
        #Make a dictionary of the included modules
        included_modules = {}
        for module_name in self.get_all_included_modules():
            included_modules[module_name] = True

        base_path = os.path.join(SAGE_SRC, 'sage')
        for directory, subdirs, files in os.walk(base_path):
            for filename in files:
                if not (filename.endswith('.py') or
                        filename.endswith('.pyx')):
                    continue

                path = os.path.join(directory, filename)

                #Create the module name
                module_name = path[len(base_path):].replace(os.path.sep, '.')
                module_name = 'sage' + module_name
                module_name = module_name[:-4] if module_name.endswith('pyx') else module_name[:-3]

                #Exclude some ones  -- we don't want init the manual
                if module_name.endswith('__init__') or module_name.endswith('all'):
                    continue

                if module_name not in included_modules:
                    yield module_name

    def print_unincluded_modules(self):
        """
        Prints all of the modules which are not included in the Sage
        reference manual.
        """
        for module_name in self.get_unincluded_modules():
            print(module_name)

    def print_included_modules(self):
        """
        Prints all of the modules that are included in the Sage reference
        manual.
        """
        for module_name in self.get_all_included_modules():
            print(module_name)


class SingleFileBuilder(DocBuilder):
    """
    This is the class used to build the documentation for a single
    user-specified file. If the file is called 'foo.py', then the
    documentation is built in ``DIR/foo/`` if the user passes the
    command line option "-o DIR", or in ``DOT_SAGE/docbuild/foo/``
    otherwise.
    """
    def __init__(self, path):
        """
        INPUT:

        - ``path`` - the path to the file for which documentation
          should be built
        """
        global options

        self.lang = 'en'
        self.name = 'single_file'
        path = os.path.abspath(path)

        # Create docbuild and relevant subdirectories, e.g.,
        # the static and templates directories in the output directory.
        # By default, this is DOT_SAGE/docbuild/MODULE_NAME, but can
        # also be specified at the command line.
        module_name = os.path.splitext(os.path.basename(path))[0]
        latex_name = module_name.replace('_', r'\_')

        if options.output_dir:
            base_dir = os.path.join(options.output_dir, module_name)
            if os.path.exists(base_dir):
                logger.warning('Warning: Directory %s exists. It is safer to build in a new directory.' % base_dir)
        else:
            DOT_SAGE = os.environ['DOT_SAGE']
            base_dir = os.path.join(DOT_SAGE, 'docbuild', module_name)
            try:
                shutil.rmtree(base_dir)
            except OSError:
                pass
        self.dir = os.path.join(base_dir, 'source')

        sage_makedirs(os.path.join(self.dir, "static"))
        sage_makedirs(os.path.join(self.dir, "templates"))
        # Write self.dir/conf.py
        conf = """# -*- coding: utf-8 -*-
# This file is automatically generated by {}, do not edit!

import sys, os
sys.path.append(os.environ['SAGE_DOC_SRC'])
sys.path.append({!r})
from common.conf import *
project = u'Documentation for {}'
release = 'unknown'
name = {!r}
html_title = project
html_short_title = project
htmlhelp_basename = name

latex_domain_indices = False
latex_documents = [
  ('index', name + '.tex', u'Documentation for {}',
   u'unknown', 'manual'),
]
""".format(__file__, self.dir, module_name, module_name, latex_name)

        if 'SAGE_DOC_UNDERSCORE' in os.environ:
            conf +="""
def setup(app):
    app.connect('autodoc-skip-member', skip_member)
"""

        with open(os.path.join(self.dir, 'conf.py'), 'w') as conffile:
            conffile.write(conf)

        # Write self.dir/index.rst
        title = 'Docs for file %s' % path
        heading = title + "\n" + ("=" * len(title))
        index = """{}

.. This file is automatically generated by {}, do not edit!

.. automodule:: {}
   :members:
   :undoc-members:
   :show-inheritance:

""".format(heading, __file__, module_name)
        with open(os.path.join(self.dir, 'index.rst'), 'w') as indexfile:
            indexfile.write(index)

        # Create link from original file to self.dir. Note that we
        # append self.dir to sys.path in conf.py. This is reasonably
        # safe (but not perfect), since we just created self.dir.
        try:
            os.symlink(path, os.path.join(self.dir, os.path.basename(path)))
        except OSError:
            pass

    def _output_dir(self, type):
        """
        Returns the directory where the output of type type is stored.
        If the directory does not exist, then it will automatically be
        created.
        """
        base_dir = os.path.split(self.dir)[0]
        d = os.path.join(base_dir, "output", type)
        sage_makedirs(d)
        return d

    def _doctrees_dir(self):
        """
        Returns the directory where the doctrees are stored.  If the
        directory does not exist, then it will automatically be
        created.
        """
        return self._output_dir('doctrees')


def get_builder(name):
    """
    Returns an appropriate *Builder object for the document ``name``.
    DocBuilder and its subclasses do all the real work in building the
    documentation.
    """
    if name == 'all':
        return AllBuilder()
    elif name.endswith('reference'):
        return ReferenceBuilder(name)
    elif 'reference' in name and os.path.exists(os.path.join(SAGE_DOC_SRC, 'en', name)):
        return ReferenceSubBuilder(name)
    elif name.endswith('website'):
        return WebsiteBuilder(name)
    elif name.startswith('file='):
        path = name[5:]
        if path.endswith('.sage') or path.endswith('.pyx'):
            raise NotImplementedError('Building documentation for a single file only works for Python files.')
        return SingleFileBuilder(path)
    elif name in get_documents() or name in AllBuilder().get_all_documents():
        return DocBuilder(name)
    else:
        print("'%s' is not a recognized document. Type 'sage --docbuild -D' for a list"%name)
        print("of documents, or 'sage --docbuild --help' for more help.")
        sys.exit(1)

def format_columns(lst, align='<', cols=None, indent=4, pad=3, width=80):
    """
    Utility function that formats a list as a simple table and returns
    a Unicode string representation.  The number of columns is
    computed from the other options, unless it's passed as a keyword
    argument.  For help on Python's string formatter, see

    http://docs.python.org/library/string.html#format-string-syntax
    """
    # Can we generalize this (efficiently) to other / multiple inputs
    # and generators?
    size = max(map(len, lst)) + pad
    if cols is None:
        import math
        cols = math.trunc((width - indent) / size)
    s = " " * indent
    for i in range(len(lst)):
        if i != 0 and i % cols == 0:
            s += "\n" + " " * indent
        s += "{0:{1}{2}}".format(lst[i], align, size)
    s += "\n"
    return unicode(s)

def help_usage(s=u"", compact=False):
    """
    Appends and returns a brief usage message for the Sage
    documentation builder.  If 'compact' is False, the function adds a
    final newline character.
    """
    s += "sage --docbuild [OPTIONS] DOCUMENT (FORMAT | COMMAND)"
    if not compact:
        s += "\n"
    return s

def help_description(s=u"", compact=False):
    """
    Appends and returns a brief description of the Sage documentation
    builder.  If 'compact' is False, the function adds a final newline
    character.
    """
    s += "Build or return information about Sage documentation.\n\n"
    s += "    DOCUMENT    name of the document to build\n"
    s += "    FORMAT      document output format\n"
    s += "    COMMAND     document-specific command\n\n"
    s += "Note that DOCUMENT may have the form 'file=/path/to/FILE',\n"
    s += "which builds the documentation for the specified file.\n\n"
    s += "A DOCUMENT and either a FORMAT or a COMMAND are required,\n"
    s += "unless a list of one or more of these is requested."
    if not compact:
        s += "\n"
    return s

def help_examples(s=u""):
    """
    Appends and returns some usage examples for the Sage documentation
    builder.
    """
    s += "Examples:\n"
    s += "    sage --docbuild -FDC all\n"
    s += "    sage --docbuild constructions pdf\n"
    s += "    sage --docbuild reference html -jv3\n"
    s += "    sage --docbuild --mathjax tutorial html\n"
    s += "    sage --docbuild reference print_unincluded_modules\n"
    s += "    sage --docbuild developer -j html --sphinx-opts -q,-aE --verbose 2"
    return s

def get_documents():
    """
    Returns a list of document names the Sage documentation builder
    will accept as command-line arguments.
    """
    all_b = AllBuilder()
    docs = all_b.get_all_documents()
    docs = [(d[3:] if d[0:3] == 'en/' else d) for d in docs]
    return docs

def help_documents(s=u""):
    """
    Appends and returns a tabular list of documents, including a
    shortcut 'all' for all documents, available to the Sage
    documentation builder.
    """
    docs = get_documents()
    s += "DOCUMENTs:\n"
    s += format_columns(docs + ['all  (!)'])
    s += "(!) Builds everything.\n\n"
    if 'reference' in docs:
        s+= "Other valid document names take the form 'reference/DIR', where\n"
        s+= "DIR is a subdirectory of SAGE_DOC_SRC/en/reference/.\n"
        s+= "This builds just the specified part of the reference manual.\n"
    s += "DOCUMENT may also have the form 'file=/path/to/FILE', which builds\n"
    s += "the documentation for the specified file.\n"
    return s

def get_formats():
    """
    Returns a list of output formats the Sage documentation builder
    will accept on the command-line.
    """
    tut_b = DocBuilder('en/tutorial')
    formats = tut_b._output_formats()
    formats.remove('html')
    return ['html', 'pdf'] + formats

def help_formats(s=u""):
    """
    Appends and returns a tabular list of output formats available to
    the Sage documentation builder.
    """
    s += "FORMATs:\n"
    s += format_columns(get_formats())
    return s

def help_commands(name='all', s=u""):
    """
    Appends and returns a tabular list of commands, if any, the Sage
    documentation builder can run on the indicated document.  The
    default is to list all commands for all documents.
    """
    # To do: Generate the lists dynamically, using class attributes,
    # as with the Builders above.
    command_dict = { 'reference' : [
        'print_included_modules',   'print_modified_modules       (*)',
        'print_unincluded_modules', 'print_newly_included_modules (*)',
        ] }
    for doc in command_dict:
        if name == 'all' or doc == name:
            s += "COMMANDs for the DOCUMENT '" + doc + "':\n"
            s += format_columns(command_dict[doc])
            s += "(*) Since the last build.\n"
    return s

def help_message_long(option, opt_str, value, parser):
    """
    Prints an extended help message for the Sage documentation builder
    and exits.
    """
    help_funcs = [ help_usage, help_description, help_documents,
                   help_formats, help_commands, parser.format_option_help,
                   help_examples ]
    for f in help_funcs:
        print(f())
    sys.exit(0)

def help_message_short(option=None, opt_str=None, value=None, parser=None,
                       error=False):
    """
    Prints a help message for the Sage documentation builder.  The
    message includes command-line usage and a list of options.  The
    message is printed only on the first call.  If error is True
    during this call, the message is printed only if the user hasn't
    requested a list (e.g., documents, formats, commands).
    """
    if not hasattr(parser.values, 'printed_help'):
        if error == True:
            if not hasattr(parser.values, 'printed_list'):
                parser.print_help()
        else:
            parser.print_help()
        setattr(parser.values, 'printed_help', 1)

def help_wrapper(option, opt_str, value, parser):
    """
    A helper wrapper for command-line options to the Sage
    documentation builder that print lists, such as document names,
    formats, and document-specific commands.
    """
    if option.dest == 'commands':
        print(help_commands(value), end="")
    if option.dest == 'documents':
        print(help_documents(), end="")
    if option.dest == 'formats':
        print(help_formats(), end="")
    setattr(parser.values, 'printed_list', 1)


class IndentedHelpFormatter2(optparse.IndentedHelpFormatter, object):
    """
    Custom help formatter class for optparse's OptionParser.
    """
    def format_description(self, description):
        """
        Returns a formatted description, preserving any original
        explicit new line characters.
        """
        if description:
            lines_in = description.split('\n')
            lines_out = [self._format_text(line) for line in lines_in]
            return "\n".join(lines_out) + "\n"
        else:
            return ""

    def format_heading(self, heading):
        """
        Returns a formatted heading using the superclass' formatter.
        If the heading is 'options', up to case, the function converts
        it to ALL CAPS.  This allows us to match the heading 'OPTIONS' with
        the same token in the builder's usage message.
        """
        if heading.lower() == 'options':
            heading = "OPTIONS"
        return super(IndentedHelpFormatter2, self).format_heading(heading)

def setup_parser():
    """
    Sets up and returns a command-line OptionParser instance for the
    Sage documentation builder.
    """
    # Documentation: http://docs.python.org/library/optparse.html
    parser = optparse.OptionParser(add_help_option=False,
                                   usage=help_usage(compact=True),
                                   formatter=IndentedHelpFormatter2(),
                                   description=help_description(compact=True))

    # Standard options. Note: We use explicit option.dest names
    # to avoid ambiguity.
    standard = optparse.OptionGroup(parser, "Standard")
    standard.add_option("-h", "--help",
                        action="callback", callback=help_message_short,
                        help="show a help message and exit")
    standard.add_option("-H", "--help-all",
                        action="callback", callback=help_message_long,
                        help="show an extended help message and exit")
    standard.add_option("-D", "--documents", dest="documents",
                        action="callback", callback=help_wrapper,
                        help="list all available DOCUMENTs")
    standard.add_option("-F", "--formats", dest="formats",
                        action="callback", callback=help_wrapper,
                        help="list all output FORMATs")
    standard.add_option("-C", "--commands", dest="commands",
                        type="string", metavar="DOC",
                        action="callback", callback=help_wrapper,
                        help="list all COMMANDs for DOCUMENT DOC; use 'all' to list all")

    standard.add_option("-i", "--inherited", dest="inherited",
                        default=False, action="store_true",
                        help="include inherited members in reference manual; may be slow, may fail for PDF output")
    standard.add_option("-u", "--underscore", dest="underscore",
                        default=False, action="store_true",
                        help="include variables prefixed with '_' in reference manual; may be slow, may fail for PDF output")

    standard.add_option("-j", "--mathjax", "--jsmath", dest="mathjax",
                        action="store_true",
                        help="render math using MathJax; FORMATs: html, json, pickle, web")
    standard.add_option("--no-plot", dest="no_plot",
                        action="store_true",
                        help="do not include graphics auto-generated using the '.. plot' markup")
    standard.add_option("--no-tests", dest="skip_tests", default=False,
                        action="store_true",
                        help="do not include TESTS blocks in the reference manual")
    standard.add_option("--no-pdf-links", dest="no_pdf_links",
                        action="store_true",
                        help="do not include PDF links in DOCUMENT 'website'; FORMATs: html, json, pickle, web")
    standard.add_option("--warn-links", dest="warn_links",
                        default=False, action="store_true",
                        help="issue a warning whenever a link is not properly resolved; equivalent to '--sphinx-opts -n' (sphinx option: nitpicky)")
    standard.add_option("--check-nested", dest="check_nested",
                        action="store_true",
                        help="check picklability of nested classes in DOCUMENT 'reference'")
    standard.add_option("-N", "--no-colors", dest="color", default=True,
                        action="store_false",
                        help="do not color output; does not affect children")
    standard.add_option("-q", "--quiet", dest="verbose",
                        action="store_const", const=0,
                        help="work quietly; same as --verbose=0")
    standard.add_option("-v", "--verbose", dest="verbose",
                        type="int", default=1, metavar="LEVEL",
                        action="store",
                        help="report progress at LEVEL=0 (quiet), 1 (normal), 2 (info), or 3 (debug); does not affect children")
    standard.add_option("-o", "--output", dest="output_dir", default=None,
                        metavar="DIR", action="store",
                        help="if DOCUMENT is a single file ('file=...'), write output to this directory")
    parser.add_option_group(standard)

    # Advanced options.
    advanced = optparse.OptionGroup(parser, "Advanced",
                                    "Use these options with care.")
    advanced.add_option("-S", "--sphinx-opts", dest="sphinx_opts",
                        type="string", metavar="OPTS",
                        action="store",
                        help="pass comma-separated OPTS to sphinx-build")
    advanced.add_option("-U", "--update-mtimes", dest="update_mtimes",
                        default=False, action="store_true",
                        help="before building reference manual, update modification times for auto-generated ReST files")
    advanced.add_option("-k", "--keep-going", dest="keep_going",
                        default=False, action="store_true",
                        help="Do not abort on errors but continue as much as possible after an error")
    parser.add_option_group(advanced)

    return parser

def setup_logger(verbose=1, color=True):
    """
    Sets up and returns a Python Logger instance for the Sage
    documentation builder.  The optional argument sets logger's level
    and message format.
    """
    # Set up colors. Adapted from sphinx.cmdline.
    import sphinx.util.console as c
    if not color or not sys.stdout.isatty() or not c.color_terminal():
        c.nocolor()

    # Available colors: black, darkgray, (dark)red, dark(green),
    # brown, yellow, (dark)blue, purple, fuchsia, turquoise, teal,
    # lightgray, white.  Available styles: reset, bold, faint,
    # standout, underline, blink.

    # Set up log record formats.
    format_std = "%(message)s"
    formatter = logging.Formatter(format_std)

    # format_debug = "%(module)s #%(lineno)s %(funcName)s() %(message)s"
    fields = ['%(module)s', '#%(lineno)s', '%(funcName)s()', '%(message)s']
    colors = ['darkblue', 'darkred', 'brown', 'reset']
    styles = ['reset', 'reset', 'reset', 'reset']
    format_debug = ""
    for i in range(len(fields)):
        format_debug += c.colorize(styles[i], c.colorize(colors[i], fields[i]))
        if i != len(fields):
            format_debug += " "

    # Documentation:  http://docs.python.org/library/logging.html
    logger = logging.getLogger('docbuild')

    # Note: There's also Handler.setLevel().  The argument is the
    # lowest severity message that the respective logger or handler
    # will pass on.  The default levels are DEBUG, INFO, WARNING,
    # ERROR, and CRITICAL.  We use "WARNING" for normal verbosity and
    # "ERROR" for quiet operation.  It's possible to define custom
    # levels.  See the documentation for details.
    if verbose == 0:
        logger.setLevel(logging.ERROR)
    if verbose == 1:
        logger.setLevel(logging.WARNING)
    if verbose == 2:
        logger.setLevel(logging.INFO)
    if verbose == 3:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(format_debug)

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


class IntersphinxCache:
    """
    Replace sphinx.ext.intersphinx.fetch_inventory by an in-memory
    cached version.
    """
    def __init__(self):
        self.inventories = {}
        self.real_fetch_inventory = sphinx.ext.intersphinx.fetch_inventory
        sphinx.ext.intersphinx.fetch_inventory = self.fetch_inventory

    def fetch_inventory(self, app, uri, inv):
        """
        Return the result of ``sphinx.ext.intersphinx.fetch_inventory()``
        from a cache if possible. Otherwise, call
        ``sphinx.ext.intersphinx.fetch_inventory()`` and cache the result.
        """
        t = (uri, inv)
        try:
            return self.inventories[t]
        except KeyError:
            i = self.real_fetch_inventory(app, uri, inv)
            self.inventories[t] = i
            return i


def patch_domain_init():
    """
    Applies a monkey-patch to the __init__ method of the Domain class in
    Sphinx, in order to work around a bug.

    See https://trac.sagemath.org/ticket/21044 as well as
    https://github.com/sphinx-doc/sphinx/pull/2816 for details about that
    bug.
    """

    from sphinx.domains import Domain
    import copy

    orig_init = Domain.__init__

    def __init__(self, *args, **kwargs):
        orig_init(self, *args, **kwargs)

        # Replace the original initial_data class attribute with a new
        # deep-copy of itself, since the bug will cause the original
        # initial_data to be modified in-place
        self.__class__.initial_data = copy.deepcopy(self.initial_data)

    Domain.__init__ = __init__


def main():
    # Parse the command-line.
    parser = setup_parser()
    global options
    options, args = parser.parse_args()

    # Get the name and type (target format) of the document we are
    # trying to build.
    try:
        name, type = args
    except ValueError:
        help_message_short(parser=parser, error=True)
        sys.exit(1)

    # Set up module-wide logging.
    global logger
    logger = setup_logger(options.verbose, options.color)

    sys.excepthook = excepthook

    # Process selected options.
    #
    # MathJax: this check usually has no practical effect, since
    # SAGE_DOC_MATHJAX is set to "True" by the script sage-env.
    # To disable MathJax, set SAGE_DOC_MATHJAX to "no" or "False".
    if options.mathjax or (os.environ.get('SAGE_DOC_MATHJAX', 'no') != 'no'
                           and os.environ.get('SAGE_DOC_MATHJAX', 'no') != 'False'):
        os.environ['SAGE_DOC_MATHJAX'] = 'True'

    if options.check_nested:
        os.environ['SAGE_CHECK_NESTED'] = 'True'

    if options.underscore:
        os.environ['SAGE_DOC_UNDERSCORE'] = "True"

    global ALLSPHINXOPTS, WEBSITESPHINXOPTS, ABORT_ON_ERROR
    if options.sphinx_opts:
        ALLSPHINXOPTS += options.sphinx_opts.replace(',', ' ') + " "
    if options.no_pdf_links:
        WEBSITESPHINXOPTS = " -A hide_pdf_links=1 "
    if options.warn_links:
        ALLSPHINXOPTS += "-n "
    if options.no_plot:
        os.environ['SAGE_SKIP_PLOT_DIRECTIVE'] = 'yes'
    if options.skip_tests:
        os.environ['SAGE_SKIP_TESTS_BLOCKS'] = 'True'

    ABORT_ON_ERROR = not options.keep_going

    patch_domain_init()

    # Delete empty directories. This is needed in particular for empty
    # directories due to "git checkout" which never deletes empty
    # directories it leaves behind. See Trac #20010.
    delete_empty_directories(SAGE_DOC_SRC)

    # Set up Intersphinx cache
    C = IntersphinxCache()

    builder = getattr(get_builder(name), type)
    builder()
