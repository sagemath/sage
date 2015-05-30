lib_LTLIBRARIES = libbz2.la

libbz2_la_SOURCES = \
	blocksort.c \
	huffman.c \
	crctable.c \
	randtable.c \
	compress.c \
	decompress.c \
	bzlib.c

libbz2_la_LDFLAGS = \
	-version-info $(BZIP2_LT_CURRENT):$(BZIP2_LT_REVISION):$(BZIP2_LT_AGE) \
	-no-undefined

include_HEADERS = bzlib.h

noinst_HEADERS = bzlib_private.h

bin_PROGRAMS = bzip2 bzip2recover

bzip2_SOURCES = bzip2.c
bzip2_LDADD = libbz2.la

bzip2recover_SOURCES = bzip2recover.c
bzip2recover_LDADD = libbz2.la

bin_SCRIPTS = bzgrep bzmore bzdiff

man_MANS = bzip2.1 bzgrep.1 bzmore.1 bzdiff.1

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = bzip2.pc

$(pkgconfig_DATA): $(srcdir)/bzip2.pc.in config.status

install-exec-hook:
	rm -f $(DESTDIR)$(bindir)/`echo "bunzip2" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzcat" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzegrep" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzfgrep" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzless" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzcmp" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	$(LN_S) `echo "bzip2" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'` $(DESTDIR)$(bindir)/`echo "bunzip2" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	$(LN_S) `echo "bzip2" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'` $(DESTDIR)$(bindir)/`echo "bzcat" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	$(LN_S) `echo "bzgrep" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'` $(DESTDIR)$(bindir)/`echo "bzegrep" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	$(LN_S) `echo "bzgrep" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'` $(DESTDIR)$(bindir)/`echo "bzfgrep" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	$(LN_S) `echo "bzmore" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'` $(DESTDIR)$(bindir)/`echo "bzless" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	$(LN_S) `echo "bzdiff" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'` $(DESTDIR)$(bindir)/`echo "bzcmp" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`

install-data-hook:
	echo ".so man1/`echo "bzip2" | sed 's,^.*/,,;$(transform)'`.1" >$(DESTDIR)$(mandir)/man1/`echo "bunzip2" | sed 's,^.*/,,;$(transform)'`.1
	echo ".so man1/`echo "bzip2" | sed 's,^.*/,,;$(transform)'`.1" >$(DESTDIR)$(mandir)/man1/`echo "bzcat" | sed 's,^.*/,,;$(transform)'`.1
	echo ".so man1/`echo "bzgrep" | sed 's,^.*/,,;$(transform)'`.1" >$(DESTDIR)$(mandir)/man1/`echo "bzegrep" | sed 's,^.*/,,;$(transform)'`.1
	echo ".so man1/`echo "bzgrep" | sed 's,^.*/,,;$(transform)'`.1" >$(DESTDIR)$(mandir)/man1/`echo "bzfgrep" | sed 's,^.*/,,;$(transform)'`.1
	echo ".so man1/`echo "bzmore" | sed 's,^.*/,,;$(transform)'`.1" >$(DESTDIR)$(mandir)/man1/`echo "bzless" | sed 's,^.*/,,;$(transform)'`.1
	echo ".so man1/`echo "bzdiff" | sed 's,^.*/,,;$(transform)'`.1" >$(DESTDIR)$(mandir)/man1/`echo "bzcmp" | sed 's,^.*/,,;$(transform)'`.1

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/`echo "bunzip2" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzcat" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzegrep" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzfgrep" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzless" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(bindir)/`echo "bzcmp" | sed 's,^.*/,,;$(transform);s/$$/$(EXEEXT)/'`
	rm -f $(DESTDIR)$(mandir)/man1/`echo "bunzip2" | sed 's,^.*/,,;$(transform)'`.1
	rm -f $(DESTDIR)$(mandir)/man1/`echo "bzcat" | sed 's,^.*/,,;$(transform)'`.1
	rm -f $(DESTDIR)$(mandir)/man1/`echo "bzegrep" | sed 's,^.*/,,;$(transform)'`.1
	rm -f $(DESTDIR)$(mandir)/man1/`echo "bzfgrep" | sed 's,^.*/,,;$(transform)'`.1
	rm -f $(DESTDIR)$(mandir)/man1/`echo "bzless" | sed 's,^.*/,,;$(transform)'`.1
	rm -f $(DESTDIR)$(mandir)/man1/`echo "bzcmp" | sed 's,^.*/,,;$(transform)'`.1

test: bzip2
	@cat $(srcdir)/words1
	./bzip2 -1  <$(srcdir)/sample1.ref >sample1.rb2
	./bzip2 -2  <$(srcdir)/sample2.ref >sample2.rb2
	./bzip2 -3  <$(srcdir)/sample3.ref >sample3.rb2
	./bzip2 -d  <$(srcdir)/sample1.bz2 >sample1.tst
	./bzip2 -d  <$(srcdir)/sample2.bz2 >sample2.tst
	./bzip2 -ds <$(srcdir)/sample3.bz2 >sample3.tst
	cmp $(srcdir)/sample1.bz2 sample1.rb2
	cmp $(srcdir)/sample2.bz2 sample2.rb2
	cmp $(srcdir)/sample3.bz2 sample3.rb2
	cmp sample1.tst $(srcdir)/sample1.ref
	cmp sample2.tst $(srcdir)/sample2.ref
	cmp sample3.tst $(srcdir)/sample3.ref
	@cat $(srcdir)/words3

manual: $(srcdir)/manual.html $(srcdir)/manual.ps $(srcdir)/manual.pdf

manual.ps: $(MANUAL_SRCS)
	cd $(srcdir); ./xmlproc.sh -ps manual.xml

manual.pdf: $(MANUAL_SRCS)
	cd $(srcdir); ./xmlproc.sh -pdf manual.xml

manual.html: $(MANUAL_SRCS)
	cd $(srcdir); ./xmlproc.sh -html manual.xml

EXTRA_DIST = \
	$(bin_SCRIPTS) \
	$(man_MANS) \
	README.autotools \
	README.XML.STUFF \
	bz-common.xsl \
	bz-fo.xsl \
	bz-html.xsl \
	bzip.css \
	bzip2.1.preformatted \
	bzip2.pc.in \
	bzip2.txt \
	dlltest.c \
	dlltest.dsp \
	entities.xml \
	format.pl \
	libbz2.def \
	libbz2.dsp \
	makefile.msc \
	manual.html \
	manual.pdf \
	manual.ps \
	manual.xml \
	mk251.c \
	sample1.bz2 \
	sample1.ref \
	sample2.bz2 \
	sample2.ref \
	sample3.bz2 \
	sample3.ref \
	spewG.c \
	unzcrash.c \
	words0 \
	words1 \
	words2 \
	words3 \
	xmlproc.sh

ACLOCAL_AMFLAGS = -I m4
