BUILT_SOURCES = done-all-links

DIRLINKS = $(top_srcdir)/Analysis \
	   $(top_srcdir)/Decay \
           $(top_srcdir)/Hadronization \
           $(top_srcdir)/MatrixElement \
           $(top_srcdir)/Models \
           $(top_srcdir)/PDF \
           $(top_srcdir)/PDT \
           $(top_srcdir)/Shower \
           $(top_srcdir)/Sampling \
           $(top_srcdir)/UnderlyingEvent \
           $(top_srcdir)/Utilities \
           $(top_srcdir)/API

LOOPTOOLHEADERS = $(top_srcdir)/Looptools/include/clooptools.h

CLEANFILES = done-all-links 

done-all-links: $(DIRLINKS) $(LOOPTOOLHEADERS)
	mkdir -p Herwig/Config Herwig/Looptools
	$(LN_S) -f $(addprefix ../, $(DIRLINKS)) Herwig
	$(LN_S) -f $(addprefix ../../, $(LOOPTOOLHEADERS)) Herwig/Looptools
	touch done-all-links

install-data-local:
	for dir in `find Herwig -follow \( -name '*.h' -or -name '*.icc' \
	-or -name '*.tcc' -or -name '*.fh' -or -name '*.xh' \) -print | \
	sed -e 's|/[^/]*$$||' | sort -u`; do \
	  files=`find $$dir -maxdepth 1 -follow \( -name '*.h' -or -name '*.icc' \
	  -or -name '*.tcc' -or -name '*.fh' -or -name '*.xh' \) -print`; \
	  $(MKDIR_P) "$(DESTDIR)$(includedir)/$$dir" && \
	  $(INSTALL_DATA) $$files "$(DESTDIR)$(includedir)/$$dir" || exit 1; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/Herwig

clean-local:
	rm -rf Herwig
