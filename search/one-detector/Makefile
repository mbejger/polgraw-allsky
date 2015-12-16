SUBDIRS = src src-gpu

subdirs: $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) -C $@

clean:
	@for D in $(SUBDIRS) ; do \
		$(MAKE) -C $$D clean ; \
	done

.PHONY: subdirs $(SUBDIRS)
