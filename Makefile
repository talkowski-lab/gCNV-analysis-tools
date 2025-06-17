docker_name := gcnv-tools
git_tag := $(shell git tag -l --contains HEAD)
ifneq ($(.SHELLSTATUS), 0)
$(error something went wrong while getting tag)
endif

ifeq ($(strip $(git_tag)),)
  git_branch:= $(shell git rev-parse --abbrev-ref HEAD)
  ifneq ($(.SHELLSTATUS), 0)
  $(error failed to get git branch)
  endif
endif

docker_tag := $(if $(strip $(git_tag)),$(git_tag),$(git_branch))
docker_path := $(if $(DOCKER_REPO),$(DOCKER_REPO)/)$(docker_name):$(docker_tag)

wdls = $(wildcard wdl/*.wdl)

.PHONY : validate_wdl build check clean install
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

build:
	R -e 'roxygen2::roxygenize("src/gelpers")'
	R CMD build src/gelpers

check: build
	R CMD check gelpers_$(gelpers_version).tar.gz

install: build
	R CMD INSTALL gelpers_$(gelpers_version).tar.gz
clean:
	rm -rf *.Rcheck

.PHONY: docker-build
docker-build:
	docker build -t $(docker_path) .

.PHONY: docker-push
docker-push:
	docker push $(docker_path)
