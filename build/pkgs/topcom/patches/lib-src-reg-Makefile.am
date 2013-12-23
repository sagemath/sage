lib_LTLIBRARIES = libCHECKREG.la

libCHECKREG_la_SOURCES = \
SPXinterface.cc \
LPinterface.cc \
RegularityCheck.cc

include_HEADERS = \
SPXinterface.hh \
RegularityCheck.hh \
LPinterface.hh

libCHECKREG_la_CPPFLAGS = \
	-DVERBOSE -DGMPRATIONAL -ftemplate-depth-30 \
	$(GMP_CFLAGS) $(CDD_CFLAGS) \
	-I../lib-src -I../wrap-gmp-gmpxx

libCHECKREG_la_LIBADD = $(CDD_LIBS)
