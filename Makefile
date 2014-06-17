# Switch between local machine and cloud
IS_LOCAL        := yes

# Package information
PKGNAME         := FlexibleSUSY
VERSION         := 1.0.1
# Needs to be changed to correspond to the fine tuning
# executable, not FlexibleSUSY
ABSBASEDIR      := /data/dylan/FlexibleSUSY
ifeq ($(IS_LOCAL),yes)
ABSBASEDIR      := /home/dylan/Documents/Postgraduate/E6SSM-Tuning
SRCDIR          := $(ABSBASEDIR)/src
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

MODELDIR        := $(SRCDIR)/E6SSM_Spectrum_Generators/models/$(MODEL)
CONFIGDIR       := $(SRCDIR)/E6SSM_Spectrum_Generators/config
FLEXIDIR        := $(SRCDIR)/E6SSM_Spectrum_Generators/src
LEGACYDIR       := $(SRCDIR)/E6SSM_Spectrum_Generators/legacy
FFLITEDIR       := $(SRCDIR)/E6SSM_Spectrum_Generators/fflite

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

ifeq ($(IS_LOCAL),yes)
CXX                := g++
CPPFLAGS           :=  -I. $(INCCONFIG) $(INCFLEXI) $(INCLEGACY) \
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
endif

ifeq ($(ENABLE_LOOPTOOLS),yes)
LOOPFUNCFLAGS	   := $(LOOPTOOLSFLAGS)
LOOPFUNCLIBS	   := $(LOOPTOOLSLIBS)
endif
ifeq ($(ENABLE_FFLITE),yes)
LOOPFUNCFLAGS	   :=
LOOPFUNCLIBS	    = $(LIBFFLITE)
endif

TUNING_HDR := \
                $(SRCDIR)/essmtuningutils.h \
                $(SRCDIR)/tuningutils.h \
                $(SRCDIR)/flags.h \
                $(SRCDIR)/tuningnumerics.h

TUNING_SRC := \
		$(SRCDIR)/essmScanner.cpp \
                $(SRCDIR)/essmtuningutils.cpp \
	        $(SRCDIR)/flags.cpp \
                $(SRCDIR)/tuningnumerics.cpp

TUNING_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TUNING_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(TUNING_SRC)))

TUNING_DEP := \
		$(TUNING_OBJ:.o=.d)

TUNING_EXE := \
		essmScanner.x

# returns file name with absolute path, taking whitespace in directory
# names into account
abspathx        = $(foreach name,$(1),\
		$(shell echo '$(abspath $(name))' | sed s/\[\[:space:\]\]/\\\\\&/g))

.PHONY:         all clean clean-dep clean-obj showbuild

all: $(TUNING_EXE)

clean-dep:
		-rm -f $(TUNING_DEP)

clean-obj:
		-rm -f $(TUNING_OBJ)

clean: clean-dep clean-obj
		-rm -f $(TUNING_EXE)

$(TUNING_DEP) $(TUNING_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(TUNING_DEP) $(TUNING_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(TUNING_EXE): $(TUNING_OBJ) $(LIBMODEL) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(CPPFLAGS)

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),clean-dep)
ifneq ($(MAKECMDGOALS),clean-obj)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),showbuild)
-include $(TUNING_DEP)
endif
endif
endif
endif
endif

%.d: %.cpp
# -MT '$*.o' ensures that the target contains the full path
	$(CXX_DEP_GEN) $(CPPFLAGS) $(CXXFLAGS_DEP_GEN) -MM -MP -MG -o $@ -MT '$*.o' $^

%.d: %.f
# the sed script ensures that the target contains the full path
	$(FOR_DEP_GEN) $(CPPFLAGS) -cpp -MM -MP -MG $^ -MT '$*.o' | \
	sed 's|.*\.o:|$*.o:|' > $@

showbuild:
	@echo "# package information"
	@echo "PKGNAME            = $(PKGNAME)"
	@echo "VERSION            = $(VERSION)"
	@echo "ABSBASEDIR         = $(ABSBASEDIR)"
	@echo ""
	@echo "# linked FlexibleSUSY libraries"
	@echo "MODEL              = $(MODEL)"
	@echo "MODELDIR           = $(MODELDIR)"
	@echo "FLEXIDIR           = $(FLEXIDIR)"
	@echo "LEGACYDIR          = $(LEGACYDIR)"
	@echo "LIBMODEL           = $(LIBMODEL)"
	@echo "LIBFLEXI           = $(LIBFLEXI)"
	@echo "LIBLEGACY          = $(LIBLEGACY)"
	@echo ""
	@echo "# compilation information"
	@echo "CXX                = $(CXX)"
	@echo "CPPFLAGS           = $(CPPFLAGS)"
	@echo "CXXFLAGS           = $(CXXFLAGS)"
	@echo "CXX_DEP_GEN        = $(CXX_DEP_GEN)"
	@echo "CXXFLAGS_DEP_GEN   = $(CXXFLAGS_DEP_GEN)"
	@echo "FC                 = $(FC)"
	@echo "FFLAGS             = $(FFLAGS)"
	@echo "FLIBS              = $(FLIBS)"
	@echo "FOR_DEP_GEN        = $(FOR_DEP_GEN)"
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
	@echo "ENABLE_FFLITE      = $(ENABLE_FFLITE)"
	@echo "ENABLE_LOOPTOOLS   = $(ENABLE_LOOPTOOLS)"
	@echo "ENABLE_THREADS     = $(ENABLE_THREADS)"
