# Compile the libtool objects of all Herwig directories in parallel before
# the ordinary recursive "make all", which then performs the links in the
# upstream order. Automake's recursive build otherwise visits the SUBDIRS
# one at a time, using at most a few cores.
#
# Run from the top of a configured Herwig build tree, then run "make all":
#   make -j8 parallel-objects
#   make -j8
# The parallel-objects target of the top-level Makefile runs
#   make -f $(srcdir)/herwig3_parallel_objects.mk parallel-objects
# The generated sources of every directory are created first, and a failure
# there stops the build: some __all.cc files are written by a subdirectory
# and compiled by its parent, so no compilation may start before they are
# complete. Object compilation errors are ignored here, so this step
# succeeds even if an object fails; "make all" then retries those objects
# and reports the error. With dependency tracking disabled, objects do not
# depend on headers, so then run this only on a clean tree.
# All -j jobs can be compilations at once, so peak memory grows with -j.
# The "+" prefixes pass the jobserver to sub-makes invoked via a variable.

PAROBJ_HELPER := $(abspath $(lastword $(MAKEFILE_LIST)))
PAROBJ_SUB = $(MAKE) --no-print-directory -f Makefile -f $(PAROBJ_HELPER) PAROBJ_SUBMAKE=1

ifeq ($(PAROBJ_SUBMAKE),)

# Only the marked lines name directories: a sub-make that regenerates an
# out-of-date Makefile also prints the messages of config.status. The names
# are joined by spaces, because make 3.81 does not split words at newlines.
parallel-objects:
	@out=`$(PAROBJ_SUB) -s parobj-subdirs` || exit 1; \
	dirs=`printf '%s\n' "$$out" | sed -n 's/^parobj-dir //p' | tr '\n' ' '`; \
	$(MAKE) --no-print-directory -f $(PAROBJ_HELPER) PAROBJ_DIRS=". $$dirs" parobj-run-sources && \
	$(MAKE) --no-print-directory -f $(PAROBJ_HELPER) PAROBJ_DIRS=". $$dirs" parobj-run-objects

parobj-run-sources: $(addprefix parobj-sources/,$(PAROBJ_DIRS))
parobj-run-objects: $(addprefix parobj-objects/,$(PAROBJ_DIRS))

parobj-sources/%:
	+@$(PAROBJ_SUB) -C $* parobj-built-sources
parobj-objects/%:
	+-@$(PAROBJ_SUB) -C $* -k parobj-objects

.PHONY: parallel-objects parobj-run-sources parobj-run-objects

else

# The empty ":" recipes stop make from printing "Nothing to be done" for
# each directory that has nothing to generate or compile.
parobj-subdirs:
	+@for dir in $(SUBDIRS); do \
	  test "$$dir" = . && continue; \
	  echo "parobj-dir $(PAROBJ_PREFIX)$$dir"; \
	  $(PAROBJ_SUB) -s -C "$$dir" PAROBJ_PREFIX="$(PAROBJ_PREFIX)$$dir/" parobj-subdirs || exit 1; \
	done
parobj-built-sources: $(BUILT_SOURCES) ; @:
parobj-objects: $(foreach lib,$(LTLIBRARIES),$($(subst /,_,$(subst -,_,$(subst .,_,$(lib))))_OBJECTS)) ; @:

.PHONY: parobj-subdirs parobj-built-sources parobj-objects

endif
