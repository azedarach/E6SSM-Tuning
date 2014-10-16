DIR     := test
MODNAME := test

LIBTEST_SRC := \
		$(DIR)/stopwatch.cpp

LIBTEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBTEST_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBTEST_SRC)))

LIBTEST_DEP := \
		$(LIBTEST_OBJ:.o=.d)

LIBTEST     := $(DIR)/lib$(MODNAME)$(LIBEXT)

TEST_SRC := \
		$(DIR)/test_lowE6SSM_stop_first_derivs.cpp \
		$(DIR)/test_lowE6SSM_stop_gN_second_derivs.cpp \
		$(DIR)/test_lowE6SSM_stop_g1_second_derivs.cpp \
		$(DIR)/test_lowE6SSM_stop_g2_second_derivs.cpp \
		$(DIR)/test_lowE6SSM_stop_mq222_second_derivs.cpp \
		$(DIR)/test_lowE6SSM_stop_vd_second_derivs.cpp \
		$(DIR)/test_lowE6SSM_stop_vu_second_derivs.cpp \
		$(DIR)/test_lowE6SSM_stop_vs_second_derivs.cpp \
		$(DIR)/test_lowE6SSM_stop_Yu22_second_derivs.cpp \
		$(DIR)/test_lowE6SSM_top_first_derivs.cpp \
		$(DIR)/test_lowE6SSM_top_second_derivs.cpp

ifneq ($(findstring lattice,$(ALGORITHMS)),)
TEST_SRC +=
endif

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
TEST_SRC +=
ifneq (,$(findstring test,$(MAKECMDGOALS)))
ALLDEP += $(LIBFFLITE_DEP)
endif
endif

TEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_SRC)))

TEST_DEP := \
		$(TEST_OBJ:.o=.d)

TEST_EXE := \
		$(TEST_OBJ:.o=.x)

TEST_EXE_LOG := $(TEST_EXE:.x=.x.log)

TEST_LOG := $(TEST_EXE_LOG)

.PHONY:		all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		clean-$(MODNAME)-log \
		execute-tests execute-compiled-tests

all-$(MODNAME): $(LIBTEST) $(TEST_EXE)

clean-$(MODNAME)-dep:
		-rm -f $(TEST_DEP)
		-rm -f $(LIBTEST_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(TEST_OBJ)
		-rm -f $(LIBTEST_OBJ)

clean-$(MODNAME)-log:
		-rm -f $(TEST_LOG)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		  clean-$(MODNAME)-log
		-rm -f $(LIBTEST)
		-rm -f $(TEST_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

$(DIR)/%.x.log: $(DIR)/%.x
		@rm -f $@
		@echo "**************************************************" >> $@;
		@echo "* executing test: $< " >> $@;
		@echo "**************************************************" >> $@;
		@$< --log_level=test_suite >> $@ 2>&1; \
		if [ $$? = 0 ]; then echo "$<: OK"; else echo "$<: FAILED"; fi

execute-tests: $(TEST_LOG)

execute-compiled-tests: $(TEST_EXE_LOG)

clean:: clean-$(MODNAME)

distclean:: distclean-$(MODNAME)

$(DIR)/test_lowE6SSM_stop_first_derivs.x: $(DIR)/test_lowE6SSM_stop_first_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS) 

$(DIR)/test_lowE6SSM_stop_gN_second_derivs.x: $(DIR)/test_lowE6SSM_stop_gN_second_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_lowE6SSM_stop_g1_second_derivs.x: $(DIR)/test_lowE6SSM_stop_g1_second_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_lowE6SSM_stop_g2_second_derivs.x: $(DIR)/test_lowE6SSM_stop_g2_second_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_lowE6SSM_stop_mq222_second_derivs.x: $(DIR)/test_lowE6SSM_stop_mq222_second_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_lowE6SSM_stop_vd_second_derivs.x: $(DIR)/test_lowE6SSM_stop_vd_second_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_lowE6SSM_stop_vu_second_derivs.x: $(DIR)/test_lowE6SSM_stop_vu_second_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_lowE6SSM_stop_vs_second_derivs.x: $(DIR)/test_lowE6SSM_stop_vs_second_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_lowE6SSM_stop_Yu22_second_derivs.x: $(DIR)/test_lowE6SSM_stop_Yu22_second_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS)

$(DIR)/test_lowE6SSM_top_first_derivs.x: $(DIR)/test_lowE6SSM_top_first_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS) 

$(DIR)/test_lowE6SSM_top_second_derivs.x: $(DIR)/test_lowE6SSM_top_second_derivs.o $(LIBlowE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS) 

$(DIR)/test_%.x: $(DIR)/test_%.o
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(BOOSTTESTLIBS) $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS) $(THREADLIBS) 

$(TEST_OBJ) $(TEST_DEP): CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)

ifeq ($(ENABLE_STATIC_LIBS),yes)
$(LIBTEST): $(LIBTEST_OBJ)
		$(MAKELIB) $@ $^
else
$(LIBTEST): $(LIBTEST_OBJ)
		$(MAKELIB) $@ $^ $(BOOSTTHREADLIBS) $(GSLLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS) 
endif

ALLDEP += $(LIBTEST_DEP) $(TEST_DEP)
ALLLIB += $(LIBTEST)
ALLTST += $(TEST_EXE)
