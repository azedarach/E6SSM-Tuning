DIR        := src
MODNAME    := essmtuning

ESSMTUNING_MK	:= \
		$(DIR)/module.mk

ESSMTUNING_HDR := \
                $(DIR)/essmtuningutils.h \
                $(DIR)/tuningutils.h \
                $(DIR)/flags.h \
                $(DIR)/tuningnumerics.h

ESSMTUNING_SRC := \
                $(DIR)/essmtuningutils.cpp \
	        $(DIR)/flags.cpp \
                $(DIR)/tuningnumerics.cpp

EXEESSMTUNING_SRC := \
		$(DIR)/essmScanner.cpp \
		$(DIR)/essmScanInputsGenerator.cpp \
		$(DIR)/essmHiggsMasses.cpp \
		$(DIR)/essmRGERunner.cpp

ESSMTUNING_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(ESSMTUNING_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(ESSMTUNING_SRC)))

ESSMTUNING_DEP := \
		$(ESSMTUNING_OBJ:.o=.d)

EXEESSMTUNING_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEESSMTUNING_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEESSMTUNING_SRC)))

EXEESSMTUNING_DEP := \
		$(EXEESSMTUNING_OBJ:.o=.d)

ESSMTUNINGSCAN_OBJ := \
			$(DIR)/essmScanner.o

ESSMTUNINGSCAN_EXE := \
			$(DIR)/essmScanner.x

ESSMGENERATOR_OBJ := \
			$(DIR)/essmScanInputsGenerator.o

ESSMGENERATOR_EXE := \
			$(DIR)/essmScanInputsGenerator.x

ESSMHIGGSSCAN_OBJ := \
			$(DIR)/essmHiggsMasses.o

ESSMHIGGSSCAN_EXE := \
			$(DIR)/essmHiggsMasses.x

ESSMRGERUNNER_OBJ := \
			$(DIR)/essmRGERunner.o

ESSMRGERUNNER_EXE := \
			$(DIR)/essmRGERunner.x

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-obj

all-$(MODNAME): $(ESSMTUNING_OBJ)

clean-$(MODNAME)-dep:
		-rm -f $(ESSMTUNING_DEP)
		-rm -f $(EXEESSMTUNING_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(ESSMTUNING_OBJ)
		-rm -f $(EXEESSMTUNING_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(ESSMTUNINGSCAN_EXE)
		-rm -f $(ESSMGENERATOR_EXE)
		-rm -f $(ESSMHIGGSSCAN_EXE)
		-rm -f $(ESSMRGERUNNER_EXE)

clean::		clean-$(MODNAME)

$(EXEESSMTUNING_DEP) $(EXEESSMTUNING_OBJ) $(ESSMTUNING_DEP) $(ESSMTUNING_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(TUNING_DEP) $(TUNING_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(ESSMTUNINGSCAN_EXE): $(ESSMTUNINGSCAN_OBJ) $(ESSMTUNING_OBJ) $(LIBMODEL) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(CPPFLAGS)

$(ESSMGENERATOR_EXE): $(ESSMGENERATOR_OBJ) $(LIBFLEXI)
		$(CXX) -o $@ $(call abspathx,$^) $(FLIBS) $(CPPFLAGS)

$(ESSMHIGGSSCAN_EXE): $(ESSMHIGGSSCAN_OBJ) $(ESSMTUNING_OBJ) $(LIBMODEL) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(CPPFLAGS)

$(ESSMRGERUNNER_EXE): $(ESSMRGERUNNER_OBJ) $(ESSMTUNING_OBJ) $(LIBMODEL) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(CPPFLAGS)

ALLDEP += $(ESSMTUNING_DEP) $(EXEESSMTUNING_DEP)
ALLSRC += $(ESSMTUNING_SRC) $(EXEESSMTUNING_SRC)
ALLEXE += $(ESSMTUNINGSCAN_EXE) $(ESSMGENERATOR_EXE) $(ESSMHIGGSSCAN_EXE) $(ESSMRGERUNNER_EXE)
