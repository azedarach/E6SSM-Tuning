# Switch between local machine and cloud
IS_LOCAL        := yes

ABSBASEDIR      := /data/dylan/essmScans
ifeq ($(IS_LOCAL),yes)
ABSBASEDIR      := /home/dylan/Documents/Postgraduate/E6SSM-Tuning
endif

# Switches
ENABLE_FFLITE    := no
ENABLE_LOOPTOOLS := no
ENABLE_THREADS   := no
ifeq ($(IS_LOCAL),yes)
ENABLE_THREADS   := yes
endif

LIBEXT          := .a

MODEL           := genericE6SSM
MODULES         := src test

# locations of spectrum generator libraries
MODELDIR        := $(ABSBASEDIR)/spectrum/models/$(MODEL)
CONFIGDIR       := $(ABSBASEDIR)/spectrum/config
FLEXIDIR        := $(ABSBASEDIR)/spectrum/src
LEGACYDIR       := $(ABSBASEDIR)/spectrum/legacy
FFLITEDIR       := $(ABSBASEDIR)/spectrum/fflite

INCMODEL        := -I$(MODELDIR)
INCCONFIG       := -I$(CONFIGDIR)
INCFLEXI        := -I$(FLEXIDIR)
INCLEGACY       := -I$(LEGACYDIR)
INCFFLITE       := -I$(FFLITEDIR)
LIBMODEL        := $(MODELDIR)/lib$(MODEL)$(LIBEXT)
LIBFLEXI        := $(FLEXIDIR)/libflexisusy$(LIBEXT)
LIBLEGACY       := $(LEGACYDIR)/liblegacy$(LIBEXT)
LIBFFLITE       := $(FFLITEDIR)/libfflite$(LIBEXT)

# Variables for compilation
ifeq ($(IS_LOCAL),yes)
CXX                := g++
CPPFLAGS           := -I. $(INCCONFIG) $(INCFLEXI) $(INCLEGACY) \
                      $(INCMODEL)  
CXXFLAGS           := -std=c++11 -O2
CXX_DEP_GEN        := g++
CXXFLAGS_DEP_GEN   := -std=c++11
FC                 := gfortran
FFLAGS             := -O2 -frecursive
FLIBS              := -L/usr/lib/gcc/x86_64-linux-gnu/4.8/ -lgfortran -lm
FOR_DEP_GEN        := gfortran
MAKELIB            := ar cru
LIBEXT             := .a
BOOSTTESTLIBS      := -lboost_unit_test_framework
BOOSTTHREADLIBS    := 
BOOSTFLAGS         := 
EIGENFLAGS         := -I/usr/include/eigen3
GSLLIBS            := -L/usr/lib -lgsl -lgslcblas -lm
GSLFLAGS           := -I/usr/include
LAPACKLIBS         := -llapack
LOOPTOOLSFLAGS     := 
LOOPTOOLSLIBS      := 
THREADLIBS         := -lpthread
else
CXX             := g++
CPPFLAGS        :=  -I. $(INCCONFIG) $(INCFLEXI) $(INCLEGACY) \
                   $(INCMODEL)
CXXFLAGS        := -std=c++0x -O2
CXX_DEP_GEN     := g++
CXXFLAGS_DEP_GEN:= -std=c++0x
FC              := gfortran
FFLAGS          := -O2 -frecursive
FLIBS           := -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/ -lgfortran -lm
FOR_DEP_GEN     := gfortran
BOOSTTESTLIBS   := -L/data/dylan/lib -lboost_unit_test_framework
BOOSTTHREADLIBS := 
BOOSTFLAGS      := -I/data/dylan/include
EIGENFLAGS      := -I/data/dylan/include/eigen3
GSLLIBS         := -lgsl -lgslcblas -lm
GSLFLAGS        := -I/usr/include
LAPACKLIBS      := -llapack
LOOPTOOLSFLAGS  := 
LOOPTOOLSLIBS   := 
THREADLIBS      := -lpthread
endif

ifeq ($(ENABLE_LOOPTOOLS),yes)
LOOPFUNCFLAGS	   := $(LOOPTOOLSFLAGS)
LOOPFUNCLIBS	   := $(LOOPTOOLSLIBS)
endif
ifeq ($(ENABLE_FFLITE),yes)
LOOPFUNCFLAGS	   :=
LOOPFUNCLIBS	    = $(LIBFFLITE)
endif

# the modules add their dependency files to this variable
ALLDEP   :=
# the modules add source files to be created to this variable
ALLSRC   :=
# the modules add executables to this variable
ALLEXE   :=
# the modules add test executables to this variable
ALLTST   :=

# returns file name with absolute path, taking whitespace in directory
# names into account
abspathx        = $(foreach name,$(1),\
		$(shell echo '$(abspath $(name))' | sed s/\[\[:space:\]\]/\\\\\&/g))

.PHONY:         all allsrc allexec alltest clean clean-executables clean-dep \
		clean-obj showbuild

all: allexec

include $(patsubst %, %/module.mk, $(MODULES))

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),install-src)
ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),depend)
ifneq ($(MAKECMDGOALS),doc)
ifeq ($(ENABLE_COMPILE),yes)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
-include $(ALLDEP)
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

allsrc:   $(ALLSRC)
allexec:  $(ALLEXE)
alltest:  $(ALLTST)

clean-dep:
	-rm -f $(ALLDEP)

depend:  clean-dep
depend:  $(ALLDEP)

%.d: %.cpp
# -MT '$*.o' ensures that the target contains the full path
	$(CXX_DEP_GEN) $(CPPFLAGS) $(CXXFLAGS_DEP_GEN) -MM -MP -MG -o $@ -MT '$*.o' $^

%.d: %.f
# the sed script ensures that the target contains the full path
	$(FOR_DEP_GEN) $(CPPFLAGS) -cpp -MM -MP -MG $^ -MT '$*.o' | \
	sed 's|.*\.o:|$*.o:|' > $@

%.d: %.F
# the sed script ensures that the target contains the full path
	$(FC) $(CPPFLAGS) -cpp -MM -MP -MG $^ -MT '$*.o' | \
	sed 's|.*\.o:|$*.o:|' > $@

clean-executables:
	-rm -f $(ALLEXE)

showbuild:
	@echo "PKGNAME            = $(PKGNAME)"
	@echo "VERSION            = $(FLEXIBLESUSY_VERSION)"
	@echo "ABSBASEDIR         = $(ABSBASEDIR)"
	@echo "INSTALL_DIR        = $(INSTALL_DIR)"
	@echo ""
	@echo "MATH               = $(MATH)"
	@echo "MODELS             = $(MODELS)"
	@echo "ALGORITHMS         = $(ALGORITHMS)"
	@echo ""
	@echo "CXX                = $(CXX)"
	@echo "CPPFLAGS           = $(CPPFLAGS)"
	@echo "CXXFLAGS           = $(CXXFLAGS)"
	@echo "CXX_DEP_GEN        = $(CXX_DEP_GEN)"
	@echo "CXXFLAGS_DEP_GEN   = $(CXXFLAGS_DEP_GEN)"
	@echo "FC                 = $(FC)"
	@echo "FFLAGS             = $(FFLAGS)"
	@echo "FLIBS              = $(FLIBS)"
	@echo "FOR_DEP_GEN        = $(FOR_DEP_GEN)"
	@echo "MAKELIB            = $(MAKELIB)"
	@echo "LIBEXT             = $(LIBEXT)"
	@echo "BOOSTTESTLIBS      = $(BOOSTTESTLIBS)"
	@echo "BOOSTTHREADLIBS    = $(BOOSTTHREADLIBS)"
	@echo "BOOSTFLAGS         = $(BOOSTFLAGS)"
	@echo "EIGENFLAGS         = $(EIGENFLAGS)"
	@echo "GSLLIBS            = $(GSLLIBS)"
	@echo "GSLFLAGS           = $(GSLFLAGS)"
	@echo "LAPACKLIBS         = $(LAPACKLIBS)"
	@echo "LOOPFUNCFLAGS      = $(LOOPFUNCFLAGS)"
	@echo "LOOPFUNCLIBS       = $(LOOPFUNCLIBS)"
	@echo "THREADLIBS         = $(THREADLIBS)"
	@echo ""
	@echo "ENABLE_COMPILE     = $(ENABLE_COMPILE)"
	@echo "ENABLE_FFLITE      = $(ENABLE_FFLITE)"
	@echo "ENABLE_LOOPTOOLS   = $(ENABLE_LOOPTOOLS)"
	@echo "ENABLE_META        = $(ENABLE_META)"
	@echo "ENABLE_STATIC_LIBS = $(ENABLE_STATIC_LIBS)"
	@echo "ENABLE_THREADS     = $(ENABLE_THREADS)"
	@echo ""
	@echo "The list of modules to be built:"
	@echo "--------------------------------"
	@echo "$(MODULES)"
