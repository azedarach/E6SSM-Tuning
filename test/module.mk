DIR      := test
MODNAME  := test

TEST_SRC := \
		$(DIR)/test_essmBetaDerivs.cpp

TEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_SRC)))

TEST_DEP := \
		$(TEST_OBJ:.o=.d)

TEST_EXE := \
		$(TEST_OBJ:.o=.x)

TEST_EXE_LOG  := $(TEST_EXE:.x=.x.log)

.PHONY:		all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-log \
		execute-compiled-tests

all-$(MODNAME): $(TEST_EXE)

clean-$(MODNAME)-dep:
		-rm -f $(TEST_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(TEST_OBJ)

clean-$(MODNAME)-log:
		-rm -f $(TEST_EXE_LOG)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
                  clean-$(MODNAME)-log
		-rm -f $(TEST_EXE)

$(DIR)/%.x.log: $(DIR)/%.x
		@rm -f $@
		@echo "**************************************************" >> $@;
		@echo "* executing test: $< " >> $@;
		@echo "**************************************************" >> $@;
		@$< >> $@ 2>&1; \
		if [ $$? = 0 ]; then echo "$<: OK"; else echo "$<: FAILED"; fi

execute-compiled-tests: $(TEST_EXE_LOG)

clean::         clean-$(MODNAME)

$(DIR)/test_essmBetaDerivs.x: $(DIR)/test_essmBetaDerivs.o $(ESSMTUNING_OBJ) $(LIBMODEL) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(FLIBS) $(CPPFLAGS)

# general test rule which links all libraries needed for a generated model
$(DIR)/test_%.x: $(DIR)/test_%.o
		$(CXX) -o $@ $(call abspathx,$^) $(BOOSTTHREADLIBS) $(THREADLIBS) $(GSLLIBS) $(FLIBS)

# add boost and eigen flags for the test object files and dependencies
$(TEST_OBJ) $(TEST_DEP): CPPFLAGS += $(BOOSTFLAGS) $(EIGENFLAGS)

ALLDEP += $(TEST_DEP)
ALLTST += $(TEST_EXE)
