# Custom Makefile for direct compilation and installation of the rDock library and binaries
# 2023-03-31, Zaragoza, Spain
# 
# How to use this Makefile:
#
# set any environment variables you need to set before (or while) invoking make with the desired target, e.g.:
#   export CXX=g++-7; make -j 4
# or
#   CXX=g++-9 PREFIX=~/.local make -j 4
#
# the following targets are available for the user:
#   build: build the library and binaries (default target)
#   build_lib: build the library
#   build_bin: build the binaries
#   test: run the tests suite
#   clean: removes the object files
#   clean_lib: removes the lib folder
#   clean_bin: removes the compiled binaries
#   veryclean: removes all the above
#   install: install the library, binaries, and headers in PREFIX folder
#   full: cleans everything and rebuilds the library and binaries from scratch,
#         then runs the tests and install everything into PREFIX folder
#
#
# the following variables can be configured by the user:
#   CONFIG: the configuration to build. Can be DEBUG or RELEASE. Default: RELEASE
#           if DEBUG is set, the library and binaries are built with debug symbols and
#           the -D_DEBUG macro is defined, making the code more verbose (but slower).
#
#   PREFIX: the folder where to install the library, binaries, and headers. Default: /usr
#
#   CXX: the c++ compiler to use. Default: g++
#
#   CXX_EXTRA_FLAGS: extra flags to pass to the compiler. Default: empty
#

PREFIX                      ?= /usr
CONFIG                      ?= RELEASE
CXX                         ?= g++
CXX_EXTRA_FLAGS				?=
CXX_BASE_FLAGS              += -pipe -std=c++14 -fPIC -fpermissive
CXX_DEBUG_CONFIG_FLAGS      += -O0 -g
CXX_RELEASE_CONFIG_FLAGS    += -O3 -ffast-math
CXX_WARNING_FLAGS           += -Wno-deprecated

DEBUG_DEFINES               := -D_DEBUG
RELEASE_DEFINES             := -D_NDEBUG

ifeq ($(CONFIG),DEBUG)
	CXX_CONFIG_FLAGS += $(CXX_DEBUG_CONFIG_FLAGS)
	DEFINES += $(DEBUG_DEFINES)
else
	CXX_CONFIG_FLAGS += $(CXX_RELEASE_CONFIG_FLAGS)
	DEFINES += $(RELEASE_DEFINES)
endif

CXX_FLAGS                   := $(CXX_BASE_FLAGS) $(CXX_CONFIG_FLAGS) $(CXX_WARNING_FLAGS) $(CXX_EXTRA_FLAGS) $(DEFINES)
LINK_FLAGS                  := -shared
LIBS                        += -lm -lpopt -lRbt
INCLUDE                     := $(addprefix -I./, $(shell find include/ -type d )) $(addprefix -I./, $(shell find import/ -type d ))
LIBRARY                     := ./lib


simplex_sources = $(shell find import/simplex/src/ -type f -name '*.cxx')
simplex_objects = $(subst import/simplex/src, obj/simplex, $(simplex_sources:.cxx=.o))

GP_sources      = $(shell find src/GP/ -type f -name 'Rbt*.cxx')
GP_objects		= $(subst src/GP/, obj/GP/, $(GP_sources:.cxx=.o))

RBT_sources     = $(shell find src/lib/ -type f -name '*.cxx')
RBT_objects     = $(subst src/lib/, obj/, $(RBT_sources:.cxx=.o))

objects = $(RBT_objects) $(simplex_objects) $(GP_objects)

objdirs = obj obj/simplex obj/GP
$(shell mkdir -p $(objdirs) ./lib ./bin)

bin_names   = rbdock rbcavity rbmoegrid rblist rbcalcgrid
bins        = $(addprefix bin/, $(bin_names))
bin_sources = $(addprefix src/exe/, $(addsuffix .cxx $(bins)))

.PHONY:	\
	install \
	target_folders build_directories \
	lib bin scripts \
	build build_lib build_bin \
	test \
	clean clean_bin clean_lib veryclean \

## User directed targets

build: build_lib build_bin scripts

install: build target_folders
	@cp -r bin/* $(PREFIX)/bin
	@cp -r lib/* $(PREFIX)/lib
	@cp -r include/* $(PREFIX)/include

full:
	$(MAKE) veryclean
	$(MAKE) build
	$(MAKE) test
	$(MAKE) install

build_bin: build_directories
	$(MAKE) $(bins)

build_lib: build_directories
	$(MAKE) lib

test:
	@echo "[test] target to be implemented"

clean:
	@rm -rf obj

clean_bin:
	@rm -f $(bins)

clean_lib:
	@rm -f lib/libRbt.so

veryclean: clean clean_bin clean_lib

## Internal targets

target_folders:
	@mkdir $(PREFIX)/{bin,lib,include} -p

build_directories:
	@mkdir -p $(objdirs) ./lib ./bin ./log ./tests/tmp

obj/%.o: src/lib/%.cxx
	@echo $(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<
	@$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

obj/simplex/%.o: import/simplex/src/%.cxx
	@echo $(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

obj/GP/%.o: src/GP/%.cxx
	@echo $(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

scripts: build_directories
	@cp -r scripts/* bin/

lib: $(objects)
	@echo $(objects)
	$(CXX) $(CXX_FLAGS) -shared -L$(LIBRARY) $^ -o lib/libRbt.so $(LIBS)

bin/%: src/exe/%.cxx lib
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -L$(LIBRARY) -o $@ $< $(LIBS)
