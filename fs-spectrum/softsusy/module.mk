DIR     := softsusy
MODNAME := softsusy

LIBSOFTSUSY_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

LIBSOFTSUSY_MK  := \
		   $(DIR)/module.mk

LIBSOFTSUSY_TARBALL := \
		   $(MODNAME).tar.gz

LIBSOFTSUSY_SRC := \
		   $(DIR)/flags.cpp \
		   $(DIR)/mssm_tuning_utils.cpp \
		   $(DIR)/softsusy_def.cpp \
		   $(DIR)/softsusy_essmsoftpars.cpp \
		   $(DIR)/softsusy_essmsusy.cpp \
		   $(DIR)/softsusy_flavoursoft.cpp \
		   $(DIR)/softsusy_linalg.cpp \
		   $(DIR)/softsusy_lowe.cpp \
		   $(DIR)/softsusy_numerics.cpp \
		   $(DIR)/softsusy_physpars.cpp \
		   $(DIR)/softsusy_rge.cpp \
		   $(DIR)/softsusy_rpvneut.cpp \
		   $(DIR)/softsusy_rpvsoft.cpp \
		   $(DIR)/softsusy_rpvsusypars.cpp \
		   $(DIR)/softsusy_softpars.cpp \
		   $(DIR)/softsusy_softsusy.cpp \
		   $(DIR)/softsusy_susy.cpp \
		   $(DIR)/softsusy_tensor.cpp \
		   $(DIR)/softsusy_twoloophiggs.f \
		   $(DIR)/softsusy_utils.cpp \
		   $(DIR)/tuning_numerics.cpp

EXESOFTSUSY_SRC := \
		   $(DIR)/cutoff_scan_mssm.cpp \
		   $(DIR)/softsusy_main.cpp \
		   $(DIR)/softsusy_rpvmain.cpp \
		   $(DIR)/softsusy_rpvneutmain.cpp \
		   $(DIR)/softsusy_softpoint.cpp

LIBSOFTSUSY_HDR := \
		   $(DIR)/flags.h \
		   $(DIR)/mssm_tuning_utils.h \
		   $(DIR)/softsusy_def.h \
		   $(DIR)/softsusy_essmsoftpars.h \
		   $(DIR)/softsusy_essmsusy.h \
		   $(DIR)/softsusy_flavoursoft.h \
		   $(DIR)/softsusy_linalg.h \
		   $(DIR)/softsusy_lowe.h \
		   $(DIR)/softsusy_mycomplex.h \
		   $(DIR)/softsusy_numerics.h \
		   $(DIR)/softsusy_physpars.h \
		   $(DIR)/softsusy_rge.h \
		   $(DIR)/softsusy_rpvneut.h \
		   $(DIR)/softsusy_rpvsoft.h \
		   $(DIR)/softsusy_rpvsusypars.h \
		   $(DIR)/softsusy_softpars.h \
		   $(DIR)/softsusy_softsusy.h \
		   $(DIR)/softsusy_susy.h \
		   $(DIR)/softsusy_tensor.h \
		   $(DIR)/softsusy_twoloophiggs.h \
		   $(DIR)/softsusy_utils.h \
		   $(DIR)/softsusy_xpr_base.h \
		   $(DIR)/softsusy_xpr_matrix.h \
		   $(DIR)/softsusy_xpr_vector.h \
		   $(DIR)/tuning_numerics.h \
		   $(DIR)/tuning_utils.h

EXESOFTSUSY_HDR := \
		   $(DIR)/softsusy_main.h \
		   $(DIR)/softsusy_rpvmain.h \
		   $(DIR)/softsusy_rpvneutmain.h \
		   $(DIR)/softsusy_softpoint.h

LIBSOFTSUSY_OBJ := \
		   $(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSOFTSUSY_SRC))) \
		   $(patsubst %.f, %.o, $(filter %.f, $(LIBSOFTSUSY_SRC)))

EXESOFTSUSY_OBJ := \
		   $(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESOFTSUSY_SRC))) \
		   $(patsubst %.f, %.o, $(filter %.f, $(EXESOFTSUSY_SRC)))

LIBSOFTSUSY_DEP := \
		   $(LIBSOFTSUSY_OBJ:.o=.d)

EXESOFTSUSY_DEP := \
		   $(EXESOFTSUSY_OBJ:.o=.d)

LIBSOFTSUSY	:= $(DIR)/lib$(MODNAME)$(LIBEXT)

CUTOFF_SCAN_OBJ := $(DIR)/cutoff_scan_mssm.o
CUTOFF_SCAN_EXE := $(DIR)/cutoff_scan_mssm.x

SOFTSUSY_MAIN_OBJ := $(DIR)/softsusy_main.o
SOFTSUSY_MAIN_EXE := $(DIR)/softsusy_main.x

SOFTSUSY_RPVMAIN_OBJ := $(DIR)/softsusy_rpvmain.o
SOFTSUSY_RPVMAIN_EXE := $(DIR)/softsusy_rpvmain.x

SOFTSUSY_RPVNEUTMAIN_OBJ := $(DIR)/softsusy_rpvneutmain.o
SOFTSUSY_RPVNEUTMAIN_EXE := $(DIR)/softsusy_rpvneutmain.x

SOFTSUSY_SOFTPOINT_OBJ := $(DIR)/softsusy_softpoint.o
SOFTSUSY_SOFTPOINT_EXE := $(DIR)/softsusy_softpoint.x

.PHONY: 	all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSOFTSUSY)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIBSOFTSUSY_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSOFTSUSY_SRC) $(LIBSOFTSUSY_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSOFTSUSY_HDR) $(LIBSOFTSUSY_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXESOFTSUSY_SRC) $(LIBSOFTSUSY_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXESOFTSUSY_HDR) $(LIBSOFTSUSY_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(LIBSOFTSUSY_MK) $(LIBSOFTSUSY_INSTALL_DIR) -m u=rw,g=r,o=r 
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBSOFTSUSY_DEP)
		-rm -f $(EXESOFTSUSY_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSOFTSUSY_OBJ)
		-rm -f $(EXESOFTSUSY_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSOFTSUSY)
		-rm -f $(CUTOFF_SCAN_EXE)
		-rm -f $(SOFTSUSY_RPVNEUTMAIN_EXE)
		-rm -f $(SOFTSUSY_RPVMAIN_EXE)
		-rm -f $(SOFTSUSY_MAIN_EXE)
		-rm -f $(SOFTSUSY_SOFTPOINT_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::		clean-$(MODNAME)

distclean::	distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(SOFTSUSY_TARBALL) \
		$(LIBSOFTSUSY_SRC) $(LIBSOFTSUSY_HDR) \
		$(EXESOFTSUSY_SRC) $(EXESOFTSUSY_HDR) \
		$(LIBSOFTSUSY_MK)

$(LIBSOFTSUSY): $(LIBSOFTSUSY_OBJ)
		$(MAKELIB) $@ $^

$(EXESOFTSUSY_DEP) $(EXESOFTSUSY_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

$(CUTOFF_SCAN_EXE): $(CUTOFF_SCAN_OBJ) $(LIBSOFTSUSY) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(SOFTSUSY_MAIN_EXE): $(SOFTSUSY_MAIN_OBJ) $(LIBSOFTSUSY)
		$(CXX) -o $@ $(call abspathx,$^) $(FLIBS)

$(SOFTSUSY_RPVMAIN_EXE): $(SOFTSUSY_RPVMAIN_OBJ) $(LIBSOFTSUSY)
		$(CXX) -o $@ $(call abspathx,$^) $(FLIBS)

$(SOFTSUSY_RPVNEUTMAIN_EXE): $(SOFTSUSY_RPVNEUTMAIN_OBJ) $(LIBSOFTSUSY)
		$(CXX) -o $@ $(call abspathx,$^) $(FLIBS)

$(SOFTSUSY_SOFTPOINT_EXE): $(SOFTSUSY_SOFTPOINT_OBJ) $(LIBSOFTSUSY)
		$(CXX) -o $@ $(call abspathx,$^) $(FLIBS)


ALLDEP += $(LIBSOFTSUSY_DEP) $(EXESOFTSUSY_DEP)
ALLSRC += $(LIBSOFTSUSY_SRC) $(EXESOFTSUSY_SRC)
ALLLIB += $(LIBSOFTSUSY)
ALLEXE += $(CUTOFF_SCAN_EXE) $(SOFTSUSY_MAIN_EXE) $(SOFTSUSY_RPVMAIN_EXE) $(SOFTSUSY_RPVNEUTMAIN_EXE) $(SOFTSUSY_SOFTPOINT_EXE)
