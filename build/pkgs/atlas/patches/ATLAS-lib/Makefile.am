ACLOCAL_AMFLAGS = -I m4

# Upstream doesn't version libraries, so this is kind of arbitrary
SO_VERSION=3:0

# You shouldn't have to customize anything from here on

ATLAS_SERIAL_LIBS=atlas cblas f77blas lapack

ATLAS_PARALLEL_LIBS=ptcblas ptf77blas ptlapack

all: all_parallel

all_serial: $(ATLAS_SERIAL_LIBS)

all_parallel: all_serial $(ATLAS_PARALLEL_LIBS)

clean:
	$(RM) -r -f .libs *.la *.dll *.dll.a *-obj

install: install_parallel

install_serial: $(ATLAS_SERIAL_LIBS)
	$(MKDIR_P) @libdir@
	$(MKDIR_P) @bindir@
	for module in $(ATLAS_SERIAL_LIBS); do \
		@LTINSTALL@; \
	done;
	@LIBTOOL@ --finish @libdir@
	@LIBTOOL@ --finish @bindir@

install_parallel: install_serial $(ATLAS_PARALLEL_LIBS)
	$(MKDIR_P) @libdir@
	$(MKDIR_P) @bindir@
	for module in $(ATLAS_PARALLEL_LIBS); do \
		@LTINSTALL@; \
	done;
	@LIBTOOL@ --finish @libdir@
	@LIBTOOL@ --finish @bindir@

# In presence of multiple definitions of thread related functions, use the ones from *_mut.o
# Works with all nm output formats (BSD/POSIX/System V)
atlas: libatlas.a
	-$(RM) -r -f libatlas-obj
	$(MKDIR) libatlas-obj
	cd libatlas-obj && $(AR) x ../$< && cd ..
	if [ `$(NM) -g $< | $(EGREP) -w 'ATL_(Set|Reset|Free|Dec)AtomicCount' | $(GREP) -w T | wc -l` -gt 4 ]; then \
	  $(RM) `ls -1 libatlas-obj/ATL_{Set,Reset,Free,Dec}AtomicCount_*.o | $(GREP) -v '_mut.o$$'`; \
	fi
	@LTCLINK@ \
	    -o @LTPREFIX@$@@LTSUFFIX@ libatlas-obj/*.o \
	    @PTHREAD_LIB@ -lm \
	    @LTLINKFLAGS@

cblas: libcblas.a @LTPREFIX@atlas@LTSUFFIX@
	-$(RM) -r -f libcblas-obj
	$(MKDIR) libcblas-obj
	cd libcblas-obj && $(AR) x ../$< && cd ..
	@LTCLINK@ \
	    -o @LTPREFIX@$@@LTSUFFIX@ libcblas-obj/*.o \
	    -latlas \
	    @LTLINKFLAGS@

ptcblas: libptcblas.a @LTPREFIX@atlas@LTSUFFIX@
	-$(RM) -r -f libptcblas-obj
	$(MKDIR) libptcblas-obj
	cd libptcblas-obj && $(AR) x ../$< && cd ..
	@LTCLINK@ \
	    -o @LTPREFIX@$@@LTSUFFIX@ libptcblas-obj/*.o \
	    @PTHREAD_LIB@ -latlas \
	    @LTLINKFLAGS@

f77blas: libf77blas.a @LTPREFIX@atlas@LTSUFFIX@
	-$(RM) -r -f libf77blas-obj
	$(MKDIR) libf77blas-obj
	cd libf77blas-obj && $(AR) x ../$< && cd ..
	@LTF77LINK@ \
	    -o @LTPREFIX@$@@LTSUFFIX@ libf77blas-obj/*.o \
	    -latlas \
	    @LTLINKFLAGS@

ptf77blas: libptf77blas.a @LTPREFIX@atlas@LTSUFFIX@
	-$(RM) -r -f libptf77blas-obj
	$(MKDIR) libptf77blas-obj
	cd libptf77blas-obj && $(AR) x ../$< && cd ..
	@LTF77LINK@ \
	    -o @LTPREFIX@$@@LTSUFFIX@ libptf77blas-obj/*.o \
	    @PTHREAD_LIB@ -latlas \
	    @LTLINKFLAGS@

f77refblas: libf77refblas.a
	-$(RM) -r -f libf77refblas-obj
	$(MKDIR) libf77refblas-obj
	cd libf77refblas-obj && $(AR) x ../$< && cd ..
	@LTF77LINK@ \
	    -o @LTPREFIX@$@@LTSUFFIX@ libf77refblas-obj/*.o \
	    -latlas \
	    @LTLINKFLAGS@

lapack: liblapack.a @LTPREFIX@atlas@LTSUFFIX@ @LTPREFIX@cblas@LTSUFFIX@ @LTPREFIX@f77blas@LTSUFFIX@
	-$(RM) -r -f liblapack-obj
	$(MKDIR) liblapack-obj
	cd liblapack-obj && $(AR) x ../$< && cd ..
	@LTF77LINK@ \
	    -o @LTPREFIX@$@@LTSUFFIX@ liblapack-obj/*.o \
	    -lcblas -lf77blas -latlas -lm \
	    @LTLINKFLAGS@

ptlapack: libptlapack.a @LTPREFIX@atlas@LTSUFFIX@ @LTPREFIX@ptcblas@LTSUFFIX@ @LTPREFIX@ptf77blas@LTSUFFIX@
	-$(RM) -r -f libptlapack-obj
	$(MKDIR) libptlapack-obj
	cd libptlapack-obj && $(AR) x ../$< && cd ..
	@LTF77LINK@ \
	    -o @LTPREFIX@$@@LTSUFFIX@ libptlapack-obj/*.o \
	    @PTHREAD_LIB@ -lptcblas -lptf77blas -latlas -lm \
	    @LTLINKFLAGS@

.PHONY: all all_serial all_parallel build install install_serial install_parallel clean $(ATLAS_SERIAL_LIBS) $(ATLAS_PARALLEL_LIBS)
