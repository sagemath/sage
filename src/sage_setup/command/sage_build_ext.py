import os
import errno

# Import setuptools before importing distutils, so that setuptools
# can replace distutils by its own vendored copy.
import setuptools

from distutils import log
from setuptools.command.build_ext import build_ext
from distutils.dep_util import newer_group
from distutils.errors import DistutilsSetupError
from sage_setup.run_parallel import execute_list_of_commands

class sage_build_ext(build_ext):
    def finalize_options(self):
        build_ext.finalize_options(self)
        self.check_flags()

    def run(self):
        # Always run the Cythonize command before building extensions
        self.run_command('build_cython')
        build_ext.run(self)

    def check_flags(self):
        """
        Sanity check the compiler flags used to build the extensions
        """
        forbidden = None
        if os.environ.get("SAGE_FAT_BINARY") == "yes":
            # When building with SAGE_FAT_BINARY=yes, we should not
            # enable CPU features which do not exist on every CPU.
            # Such flags usually come from other libraries adding the
            # flags to the pkgconfig configuration.  So if you hit these
            # errors, the problem is most likely with some external
            # library and not with Sage.
            import re
            forbidden = re.compile(r"-march=|-mpcu=|-msse3|-msse4|-mpopcnt|-mavx")

        if forbidden is not None:
            errors = 0
            for ext in self.extensions:
                flags = ext.extra_compile_args
                for flag in flags:
                    if forbidden.match(flag):
                        log.error("%s uses forbidden flag '%s'", ext.name, flag)
                        errors += 1
            if errors:
                raise RuntimeError("forbidden flags used")

    def build_extensions(self):

        from distutils.debug import DEBUG

        if DEBUG:
            print("self.compiler.compiler:")
            print(self.compiler.compiler)
            print("self.compiler.compiler_cxx:")
            print(self.compiler.compiler_cxx) # currently not used
            print("self.compiler.compiler_so:")
            print(self.compiler.compiler_so)
            print("self.compiler.linker_so:")
            print(self.compiler.linker_so)
            # There are further interesting variables...

        if DEBUG:
            print("self.compiler.linker_so (after fixing library dirs):")
            print(self.compiler.linker_so)


        # First, sanity-check the 'extensions' list
        self.check_extensions_list(self.extensions)

        import time
        t = time.time()

        compile_commands = []
        for ext in self.extensions:
            need_to_compile, p = self.prepare_extension(ext)
            if need_to_compile:
                compile_commands.append((self.build_extension, p))

        execute_list_of_commands(compile_commands)

        print("Total time spent compiling C/C++ extensions: %.2f seconds." % (time.time() - t))

    def prepare_extension(self, ext):
        sources = ext.sources
        if sources is None or not isinstance(sources, (list, tuple)):
            raise DistutilsSetupError(("in 'ext_modules' option (extension '%s'), " +
                   "'sources' must be present and must be " +
                   "a list of source filenames") % ext.name)
        sources = list(sources)

        fullname = self.get_ext_fullname(ext.name)
        if self.inplace:
            # ignore build-lib -- put the compiled extension into
            # the source tree along with pure Python modules

            modpath = fullname.split('.')
            package = '.'.join(modpath[0:-1])
            base = modpath[-1]

            build_py = self.get_finalized_command('build_py')
            package_dir = build_py.get_package_dir(package)
            ext_filename = os.path.join(package_dir,
                                        self.get_ext_filename(base))
            relative_ext_filename = self.get_ext_filename(base)
        else:
            ext_filename = os.path.join(self.build_lib,
                                        self.get_ext_filename(fullname))
            relative_ext_filename = self.get_ext_filename(fullname)

        # while dispatching the calls to gcc in parallel, we sometimes
        # hit a race condition where two separate build_ext objects
        # try to create a given directory at the same time; whoever
        # loses the race then seems to throw an error, saying that
        # the directory already exists. so, instead of fighting to
        # fix the race condition, we simply make sure the entire
        # directory tree exists now, while we're processing the
        # extensions in serial.
        relative_ext_dir = os.path.split(relative_ext_filename)[0]
        prefixes = ['', self.build_lib, self.build_temp]
        for prefix in prefixes:
            path = os.path.join(prefix, relative_ext_dir)
            try:
                os.makedirs(path)
            except OSError as e:
                assert e.errno==errno.EEXIST, 'Cannot create %s.' % path
        depends = sources + ext.depends
        if not (self.force or newer_group(depends, ext_filename, 'newer')):
            log.debug("skipping '%s' extension (up-to-date)", ext.name)
            need_to_compile = False
        else:
            log.info("building '%s' extension", ext.name)
            need_to_compile = True

        return need_to_compile, (sources, ext, ext_filename)

    def build_extension(self, p):

        sources, ext, ext_filename = p

        # First, scan the sources for SWIG definition files (.i), run
        # SWIG on 'em to create .c files, and modify the sources list
        # accordingly.
        sources = self.swig_sources(sources, ext)

        # Next, compile the source code to object files.

        # XXX not honouring 'define_macros' or 'undef_macros' -- the
        # CCompiler API needs to change to accommodate this, and I
        # want to do one thing at a time!

        # Two possible sources for extra compiler arguments:
        #   - 'extra_compile_args' in Extension object
        #   - CFLAGS environment variable (not particularly
        #     elegant, but people seem to expect it and I
        #     guess it's useful)
        # The environment variable should take precedence, and
        # any sensible compiler will give precedence to later
        # command line args.  Hence we combine them in order:
        extra_args = ext.extra_compile_args or []

        macros = ext.define_macros[:]
        for undef in ext.undef_macros:
            macros.append((undef,))

        objects = self.compiler.compile(sources,
                                        output_dir=self.build_temp,
                                        macros=macros,
                                        include_dirs=ext.include_dirs,
                                        debug=self.debug,
                                        extra_postargs=extra_args,
                                        depends=ext.depends)

        # XXX -- this is a Vile HACK!
        #
        # The setup.py script for Python on Unix needs to be able to
        # get this list so it can perform all the clean up needed to
        # avoid keeping object files around when cleaning out a failed
        # build of an extension module.  Since Distutils does not
        # track dependencies, we have to get rid of intermediates to
        # ensure all the intermediates will be properly re-built.
        #
        self._built_objects = objects[:]

        # Now link the object files together into a "shared object" --
        # of course, first we have to figure out all the other things
        # that go into the mix.
        if ext.extra_objects:
            objects.extend(ext.extra_objects)
        extra_args = ext.extra_link_args or []

        # Detect target language, if not provided
        language = ext.language or self.compiler.detect_language(sources)

        self.compiler.link_shared_object(
            objects, ext_filename,
            libraries=self.get_libraries(ext),
            library_dirs=ext.library_dirs,
            runtime_library_dirs=ext.runtime_library_dirs,
            extra_postargs=extra_args,
            export_symbols=self.get_export_symbols(ext),
            debug=self.debug,
            build_temp=self.build_temp,
            target_lang=language)
