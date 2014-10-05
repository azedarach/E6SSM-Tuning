DIR          := models/lowE6SSM
MODNAME      := lowE6SSM
SARAH_MODEL  := genericE6SSM

lowE6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

lowE6SSM_MK     := \
		$(DIR)/module.mk

lowE6SSM_TWO_SCALE_MK := \
		$(DIR)/two_scale_susy.mk \
		$(DIR)/two_scale_soft.mk

lowE6SSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.lowE6SSM

lowE6SSM_GNUPLOT := \
		$(DIR)/lowE6SSM_plot_rgflow.gnuplot \
		$(DIR)/lowE6SSM_plot_spectrum.gnuplot

lowE6SSM_TARBALL := \
		$(MODNAME).tar.gz

LIBlowE6SSM_SRC :=
EXElowE6SSM_SRC :=

LIBlowE6SSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBlowE6SSM_SRC += \
		$(DIR)/lowE6SSM_info.cpp \
		$(DIR)/lowE6SSM_slha_io.cpp \
		$(DIR)/lowE6SSM_physical.cpp \
		$(DIR)/lowE6SSM_utilities.cpp \
		$(DIR)/lowE6SSM_two_scale_convergence_tester.cpp \
		$(DIR)/lowE6SSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/lowE6SSM_two_scale_initial_guesser.cpp \
		$(DIR)/lowE6SSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/lowE6SSM_two_scale_model.cpp \
		$(DIR)/lowE6SSM_two_scale_susy_parameters.cpp \
		$(DIR)/lowE6SSM_two_scale_soft_parameters.cpp \
		$(DIR)/lowE6SSM_two_scale_susy_scale_constraint.cpp \
		$(DIR)/lowE6SSM_two_scale_tuning_calculator.cpp
EXElowE6SSM_SRC += \
		$(DIR)/run_lowE6SSM.cpp \
		$(DIR)/scan_lowE6SSM.cpp
LIBlowE6SSM_HDR += \
		$(DIR)/lowE6SSM_convergence_tester.hpp \
		$(DIR)/lowE6SSM_high_scale_constraint.hpp \
		$(DIR)/lowE6SSM_info.hpp \
		$(DIR)/lowE6SSM_initial_guesser.hpp \
		$(DIR)/lowE6SSM_input_parameters.hpp \
		$(DIR)/lowE6SSM_low_scale_constraint.hpp \
		$(DIR)/lowE6SSM_model.hpp \
		$(DIR)/lowE6SSM_physical.hpp \
		$(DIR)/lowE6SSM_slha_io.hpp \
		$(DIR)/lowE6SSM_spectrum_generator.hpp \
		$(DIR)/lowE6SSM_susy_scale_constraint.hpp \
		$(DIR)/lowE6SSM_utilities.hpp \
		$(DIR)/lowE6SSM_two_scale_convergence_tester.hpp \
		$(DIR)/lowE6SSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/lowE6SSM_two_scale_initial_guesser.hpp \
		$(DIR)/lowE6SSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/lowE6SSM_two_scale_model.hpp \
		$(DIR)/lowE6SSM_two_scale_soft_parameters.hpp \
		$(DIR)/lowE6SSM_two_scale_susy_parameters.hpp \
		$(DIR)/lowE6SSM_two_scale_susy_scale_constraint.hpp \
		$(DIR)/lowE6SSM_two_scale_tuning_calculator.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(DIR)/two_scale_susy.mk
-include $(DIR)/two_scale_soft.mk
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(DIR)/two_scale_susy.mk: run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(DIR)/two_scale_soft.mk: run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
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

# remove duplicates in case all algorithms are used
LIBlowE6SSM_SRC := $(sort $(LIBlowE6SSM_SRC))
EXElowE6SSM_SRC := $(sort $(EXElowE6SSM_SRC))

LIBlowE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBlowE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBlowE6SSM_SRC)))

EXElowE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXElowE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXElowE6SSM_SRC)))

LIBlowE6SSM_DEP := \
		$(LIBlowE6SSM_OBJ:.o=.d)

EXElowE6SSM_DEP := \
		$(EXElowE6SSM_OBJ:.o=.d)

LIBlowE6SSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_lowE6SSM_OBJ := $(DIR)/run_lowE6SSM.o
RUN_lowE6SSM_EXE := $(DIR)/run_lowE6SSM.x

SCAN_lowE6SSM_OBJ := $(DIR)/scan_lowE6SSM.o
SCAN_lowE6SSM_EXE := $(DIR)/scan_lowE6SSM.x

METACODE_STAMP_lowE6SSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_lowE6SSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBlowE6SSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(lowE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBlowE6SSM_SRC) $(lowE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBlowE6SSM_HDR) $(lowE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXElowE6SSM_SRC) $(lowE6SSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(lowE6SSM_MK) $(lowE6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(lowE6SSM_TWO_SCALE_MK) $(lowE6SSM_INSTALL_DIR)
ifneq ($(lowE6SSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(lowE6SSM_SLHA_INPUT) $(lowE6SSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(lowE6SSM_GNUPLOT) $(lowE6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBlowE6SSM_DEP)
		-rm -f $(EXElowE6SSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBlowE6SSM_OBJ)
		-rm -f $(EXElowE6SSM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBlowE6SSM)
		-rm -f $(RUN_lowE6SSM_EXE)
		-rm -f $(SCAN_lowE6SSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(lowE6SSM_TARBALL) \
		$(LIBlowE6SSM_SRC) $(LIBlowE6SSM_HDR) \
		$(EXElowE6SSM_SRC) \
		$(lowE6SSM_MK) $(lowE6SSM_TWO_SCALE_MK) \
		$(lowE6SSM_SLHA_INPUT) $(lowE6SSM_GNUPLOT)

$(LIBlowE6SSM_SRC) $(LIBlowE6SSM_HDR) $(EXElowE6SSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_lowE6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_lowE6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_lowE6SSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_lowE6SSM)"
		@echo "Note: to regenerate lowE6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_lowE6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_lowE6SSM):
		@true
endif

$(LIBlowE6SSM_DEP) $(EXElowE6SSM_DEP) $(LIBlowE6SSM_OBJ) $(EXElowE6SSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBlowE6SSM_DEP) $(EXElowE6SSM_DEP) $(LIBlowE6SSM_OBJ) $(EXElowE6SSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBlowE6SSM): $(LIBlowE6SSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_lowE6SSM_EXE): $(RUN_lowE6SSM_OBJ) $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(SCAN_lowE6SSM_EXE): $(SCAN_lowE6SSM_OBJ) $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBlowE6SSM_DEP) $(EXElowE6SSM_DEP)
ALLSRC += $(LIBlowE6SSM_SRC) $(EXElowE6SSM_SRC)
ALLLIB += $(LIBlowE6SSM)
ALLEXE += $(RUN_lowE6SSM_EXE) $(SCAN_lowE6SSM_EXE)
