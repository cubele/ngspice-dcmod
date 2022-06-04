## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libsparse.la libgmres.la libgraphops.la

libsparse_la_SOURCES = \
	spalloc.c	\
	spbuild.c	\
	spconfig.h	\
	spdefs.h	\
	spextra.c	\
	spfactor.c	\
	spoutput.c	\
	spsmp.c		\
	spsolve.c	\
	sputils.c

libgmres_la_SOURCES = \
	gmres.c	\
	gmresutils.h	\
	gmresutils.c	\
	newiter.c

libgraphops_la_SOURCES = \
	graphops.hpp	\
	graphops.cpp	\
	feGRASS.cpp	\
	feGRASS.hpp	\
	graphwrapper.cpp	\
	trialmodel.hpp

AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src/include -O2
AM_CXXFLAGS = -I$(top_srcdir)/src/include -std=c++11 -O2
AM_CFLAGS = $(STATIC)

MAINTAINERCLEANFILES = Makefile.in
