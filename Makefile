wdls = $(wildcard wdl/*.wdl)

.PHONY : validate_wdl pkg check clean docs
gelpers_version := $(shell grep -F Version src/gelpers/DESCRIPTION | cut -d ' ' -f 2)

validate_wdl:
ifneq ($(findstring validate_wdl,$(MAKECMDGOALS)),)
  ifeq ($(strip $(WOMTOOL)),)
    $(error WOMTOOL must be set to path of womtool JAR)
  endif
endif
ifneq ($(strip $(wdls)),)
	@for wdl in $(wdls); do \
	  printf "\n\n\e[38;5;11m*** Validating %s ***\e[0m\n" "$${wdl}"; \
	  java -jar $(WOMTOOL) validate "$${wdl}" || exit 1; \
	done
endif

pkg:
	R CMD build gelpers

check:
	R CMD check gelpers_$(gelpers_version).tar.gz

docs:
	R -e 'roxygen2::roxygenize("src/gelpers")'

clean:
	rm -rf *.Rcheck
