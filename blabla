# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

#include Makefile.arch

#pour une architecture standard avec command use root
include $(ROOTSYS)/etc/Makefile.arch 

#------------------------------------------------------------------------------

#ANALYSEFO       = TAnalyseFASTER.$(ObjSuf)
#ANALYSEFS       = TAnalyseFASTER.$(SrcSuf)

#DATAMANIPO		= TDataManip.$(ObjSuf)
#DATAMANIPS		= TDataManip.$(SrcSuf)

#DETECTORO		= TDetecteur.$(ObjSuf)
#DETECTORS   	= TDetecteur.$(SrcSuf)

EXPERIMENT0     = CExperiment.$(ObjSuf)
EXPERIMENTS     = CExperiment.$(SrcSuf)

ifeq ($(PLATFORM),win32)
CLASSLIB      = libClasslib.lib
else
CLASSLIB      = $(CLASSSO)
endif

MAINO    = Analyse_Main.$(ObjSuf) 
MAINS    = Analyse_Main.$(SrcSuf)
MAIN     = myAnalyse$(ExeSuf)

OBJS          = $(EXPERIMENT0) $(MAINO)

PROGRAMS      = $(MAIN)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    Aclock Hello Tetris

all:            $(PROGRAMS)

$(CLASSSO):     $(CLASSO) $(EXPERIMENT0)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) -g $(SOFLAGS) $^ $(OutPutOpt) $@
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -g -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@) -lSpectrum
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) -g $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
else
		$(LD) -g $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
endif
		@echo "$@ done"


$(MAIN):        $(MAINO) $(EXPERIMENT0) $(CLASSSO)
		$(LD) $(LDFLAGS) -g $(MAINO) $(EXPERIMENT0) $(CLASSLIB) $(LIBS) \
		  -lSpectrum $(OutPutOpt)$(MAIN) 
		@echo "$@ done"

clean:
		@rm -f $(OBJS) $(CLASSSO) *Dict* core

.SUFFIXES: .$(SrcSuf)

###

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

###



