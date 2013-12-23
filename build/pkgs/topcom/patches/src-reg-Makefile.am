bin_PROGRAMS = checkregularity

checkregularity_SOURCES = checkregularity.cc

LDADD           = ../lib-src-reg/libCHECKREG.la \
                  ../lib-src/libTOPCOM.la  

AM_CPPFLAGS     = $(GMP_CFLAGS) $(CDD_CFLAGS)
AM_CPPFLAGS    += -DVERBOSE -DGMPRATIONAL -ftemplate-depth-30
AM_CPPFLAGS    += -I../wrap-gmp-gmpxx
AM_CPPFLAGS    += -I../lib-src
AM_CPPFLAGS    += -I../lib-src-reg
AM_LDFLAGS      = $(CDD_LIBS) $(GMP_LIBS)
