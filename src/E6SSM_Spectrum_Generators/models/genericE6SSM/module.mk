DIR          := models/genericE6SSM
MODNAME      := genericE6SSM
SARAH_MODEL  := genericE6SSM

genericE6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

genericE6SSM_MK     := \
		$(DIR)/module.mk

genericE6SSM_TWO_SCALE_MK := \
		$(DIR)/two_scale_susy.mk \
		$(DIR)/two_scale_soft.mk

genericE6SSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.genericE6SSM

genericE6SSM_GNUPLOT := \
		$(DIR)/genericE6SSM_plot_rgflow.gnuplot \
		$(DIR)/genericE6SSM_plot_spectrum.gnuplot

LIBgenericE6SSM_SRC :=
EXEgenericE6SSM_SRC :=

LIBgenericE6SSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBgenericE6SSM_SRC += \
		$(DIR)/genericE6SSM_info.cpp \
		$(DIR)/genericE6SSM_slha_io.cpp \
		$(DIR)/genericE6SSM_physical.cpp \
		$(DIR)/genericE6SSM_utilities.cpp \
		$(DIR)/genericE6SSM_two_scale_convergence_tester.cpp \
		$(DIR)/genericE6SSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/genericE6SSM_two_scale_initial_guesser.cpp \
		$(DIR)/genericE6SSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/genericE6SSM_two_scale_model.cpp \
		$(DIR)/genericE6SSM_two_scale_susy_parameters.cpp \
		$(DIR)/genericE6SSM_two_scale_soft_parameters.cpp \
		$(DIR)/genericE6SSM_two_scale_susy_scale_constraint.cpp
EXEgenericE6SSM_SRC += \
		$(DIR)/run_genericE6SSM.cpp \
		$(DIR)/scan_genericE6SSM.cpp
LIBgenericE6SSM_HDR += \
		$(DIR)/genericE6SSM_convergence_tester.hpp \
		$(DIR)/genericE6SSM_high_scale_constraint.hpp \
		$(DIR)/genericE6SSM_info.hpp \
		$(DIR)/genericE6SSM_initial_guesser.hpp \
		$(DIR)/genericE6SSM_input_parameters.hpp \
		$(DIR)/genericE6SSM_low_scale_constraint.hpp \
		$(DIR)/genericE6SSM_model.hpp \
		$(DIR)/genericE6SSM_physical.hpp \
		$(DIR)/genericE6SSM_slha_io.hpp \
		$(DIR)/genericE6SSM_spectrum_generator.hpp \
		$(DIR)/genericE6SSM_susy_scale_constraint.hpp \
		$(DIR)/genericE6SSM_utilities.hpp \
		$(DIR)/genericE6SSM_two_scale_convergence_tester.hpp \
		$(DIR)/genericE6SSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/genericE6SSM_two_scale_initial_guesser.hpp \
		$(DIR)/genericE6SSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/genericE6SSM_two_scale_model.hpp \
		$(DIR)/genericE6SSM_two_scale_soft_parameters.hpp \
		$(DIR)/genericE6SSM_two_scale_susy_parameters.hpp \
		$(DIR)/genericE6SSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(DIR)/two_scale_susy.mk
-include $(DIR)/two_scale_soft.mk
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(DIR)/two_scale_susy.mk: run-metacode-$(MODNAME)
		@true
$(DIR)/two_scale_soft.mk: run-metacode-$(MODNAME)
		@true
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

# remove duplicates in case all algorithms are used
LIBgenericE6SSM_SRC := $(sort $(LIBgenericE6SSM_SRC))
EXEgenericE6SSM_SRC := $(sort $(EXEgenericE6SSM_SRC))

LIBgenericE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBgenericE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBgenericE6SSM_SRC)))

EXEgenericE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEgenericE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEgenericE6SSM_SRC)))

LIBgenericE6SSM_DEP := \
		$(LIBgenericE6SSM_OBJ:.o=.d)

EXEgenericE6SSM_DEP := \
		$(EXEgenericE6SSM_OBJ:.o=.d)

LIBgenericE6SSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_genericE6SSM_OBJ := $(DIR)/run_genericE6SSM.o
RUN_genericE6SSM_EXE := $(DIR)/run_genericE6SSM.x

SCAN_genericE6SSM_OBJ := $(DIR)/scan_genericE6SSM.o
SCAN_genericE6SSM_EXE := $(DIR)/scan_genericE6SSM.x

METACODE_STAMP_genericE6SSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

SARAH_MODEL_FILES_genericE6SSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		distclean-$(MODNAME) run-metacode-$(MODNAME)

all-$(MODNAME): $(LIBgenericE6SSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(genericE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBgenericE6SSM_SRC) $(genericE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBgenericE6SSM_HDR) $(genericE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEgenericE6SSM_SRC) $(genericE6SSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(genericE6SSM_MK) $(genericE6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(genericE6SSM_TWO_SCALE_MK) $(genericE6SSM_INSTALL_DIR)
ifneq ($(genericE6SSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(genericE6SSM_SLHA_INPUT) $(genericE6SSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(genericE6SSM_GNUPLOT) $(genericE6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBgenericE6SSM_DEP)
		-rm -f $(EXEgenericE6SSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBgenericE6SSM_OBJ)
		-rm -f $(EXEgenericE6SSM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBgenericE6SSM)
		-rm -f $(RUN_genericE6SSM_EXE)
		-rm -f $(SCAN_genericE6SSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBgenericE6SSM_SRC) $(LIBgenericE6SSM_HDR) $(EXEgenericE6SSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_genericE6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_genericE6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_genericE6SSM)
		$(MATH) -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_genericE6SSM)"
		@echo "Note: to regenerate genericE6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_genericE6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_genericE6SSM):
		@true
endif

$(LIBgenericE6SSM_DEP) $(EXEgenericE6SSM_DEP) $(LIBgenericE6SSM_OBJ) $(EXEgenericE6SSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBgenericE6SSM_DEP) $(EXEgenericE6SSM_DEP) $(LIBgenericE6SSM_OBJ) $(EXEgenericE6SSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBgenericE6SSM): $(LIBgenericE6SSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_genericE6SSM_EXE): $(RUN_genericE6SSM_OBJ) $(LIBgenericE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS)

$(SCAN_genericE6SSM_EXE): $(SCAN_genericE6SSM_OBJ) $(LIBgenericE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS)

ALLDEP += $(LIBgenericE6SSM_DEP) $(EXEgenericE6SSM_DEP)
ALLSRC += $(LIBgenericE6SSM_SRC) $(EXEgenericE6SSM_SRC)
ALLLIB += $(LIBgenericE6SSM)
ALLEXE += $(RUN_genericE6SSM_EXE) $(SCAN_genericE6SSM_EXE)
