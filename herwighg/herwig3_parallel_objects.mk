# Compile the libtool objects of all Herwig directories in parallel before
# the ordinary recursive "make all", which then performs the links in the
# upstream order. Automake's recursive build otherwise visits the SUBDIRS
# one at a time, using at most a few cores.
#
# Run from the top of a configured Herwig tree, then run "make all":
#   make -f herwig3_parallel_objects.mk parallel-objects
# The generated sources of every directory are created first, and a failure
# there stops the build: some __all.cc files are written by a subdirectory
# and compiled by its parent, so no compilation may start before they are
# complete. Object compilation errors are ignored here; "make all" then
# retries those objects and reports the error. This relies on a clean tree:
# with dependency tracking disabled, objects do not depend on headers.
# The "+" prefixes pass the jobserver to sub-makes invoked via a variable.

LCG_HERWIG_HELPER := $(abspath $(lastword $(MAKEFILE_LIST)))
LCG_HERWIG_SUB = $(MAKE) --no-print-directory -f Makefile -f $(LCG_HERWIG_HELPER) LCG_HERWIG_SUBMAKE=1

ifeq ($(LCG_HERWIG_SUBMAKE),)

parallel-objects:
	@dirs=`$(LCG_HERWIG_SUB) -s lcg-subdirs` || exit 1; \
	$(MAKE) --no-print-directory -f $(LCG_HERWIG_HELPER) LCG_HERWIG_DIRS=". $$dirs" lcg-run-sources && \
	$(MAKE) --no-print-directory -f $(LCG_HERWIG_HELPER) LCG_HERWIG_DIRS=". $$dirs" lcg-run-objects

lcg-run-sources: $(addprefix lcg-sources/,$(LCG_HERWIG_DIRS))
lcg-run-objects: $(addprefix lcg-objects/,$(LCG_HERWIG_DIRS))

lcg-sources/%:
	+@$(LCG_HERWIG_SUB) -C $* lcg-built-sources
lcg-objects/%:
	+-@$(LCG_HERWIG_SUB) -C $* -k lcg-objects

.PHONY: parallel-objects lcg-run-sources lcg-run-objects

else

lcg-subdirs:
	+@for dir in $(SUBDIRS); do \
	  test "$$dir" = . && continue; \
	  echo "$(LCG_HERWIG_PREFIX)$$dir"; \
	  $(LCG_HERWIG_SUB) -s -C "$$dir" LCG_HERWIG_PREFIX="$(LCG_HERWIG_PREFIX)$$dir/" lcg-subdirs || exit 1; \
	done
lcg-built-sources: $(BUILT_SOURCES)
lcg-objects: $(foreach lib,$(LTLIBRARIES),$($(subst /,_,$(subst -,_,$(subst .,_,$(lib))))_OBJECTS))

.PHONY: lcg-subdirs lcg-built-sources lcg-objects

endif
