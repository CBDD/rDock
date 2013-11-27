# Builds the rDock library and exes for a specified
# combination of platform and configuration
#
# Platform =
#   linux-g++ (default)
#   linux-g++-64
#   linux-pathCC-64
#
# Configuration =
#   release (default)
#   debug
#
# Called by master makefile or can be called standalone
# e.g.
#   make "PLATFORM=linux-g++-64" "CONFIG=debug" "VERSION=enspiral_dev"
#
# The VERSION string affects the naming of the libRbt shared library,
# serves no other purpose.
#
# This Makefile uses the tmake, exemake, progen, buildnum.csh scripts
# that must be installed in the build directory.
#
# Generated files:
# ./Makefile.$(PLATFORM)-$(CONFIG)-lib	Custom Makefile for libRbt.so
# ./Makefile.$(PLATFORM)-$(CONFIG)-exe	Custom Makefile for rDock exes
# ./$(PLATFORM)/$(CONFIG)/		Intermediate build files
#
PROJECT		= libRbt.pro	# tmake project file to create
BINDIR		= ../bin	# Install location for exes
LIBDIR		= ../lib	# Install location for libs
CP		= cp -df
PARALLEL	= 4		# Uses make -j to build across 2 processors

# Configurations
RELEASE		= release
DEBUG		= debug
RELEASE_DEFINES	= _NDEBUG
DEBUG_DEFINES	= _DEBUG

# Default platform, configuration, and version
# PLATFORM must define a valid tmake configuration
# i.e. ./tmakelib/$(PLATFORM)/tmake.conf must exist.
PLATFORM	= linux-g++-64
CONFIG		= $(RELEASE)
VERSION		= enspiral_dev

# Platform/config dependent macros
ifeq "$(CONFIG)" "$(RELEASE)"
DEFINES = $(RELEASE_DEFINES)
else
DEFINES = $(DEBUG_DEFINES)
endif
TARGET		= $(PLATFORM)-$(CONFIG)
BUILDDIR	= $(PLATFORM)/$(CONFIG)
MAKEFILE_LIB	= Makefile.$(TARGET)-lib
MAKEFILE_EXE	= Makefile.$(TARGET)-exe

# All of the source files for libRbt.so
SRC_FILES	=	../include/*.h \
			../include/GP/*.h \
			../import/simplex/include/*.h \
			../import/tnt/include/*.h \
			../src/lib/*.cxx \
			../src/GP/Rbt*.cxx \
			../import/simplex/src/*.cxx \
##################################################
rdock: rdock-install-lib rdock-install-exe

rdock-dt: rdock-install-lib rdock-install-dt-exe

rdock-install-lib: rdock-lib
	$(CP) $(BUILDDIR)/lib/* $(LIBDIR)

rdock-install-exe: rdock-exe
	$(CP) $(BUILDDIR)/exe/* $(BINDIR)

rdock-install-dt-exe: rdock-dt-exe
	$(CP) $(BUILDDIR)/exe/* $(BINDIR)

rdock-exe: rdock-lib
	./exemake -t"$(PLATFORM)" -c"$(CONFIG)" -o"$(MAKEFILE_EXE)"
	$(MAKE) -j $(PARALLEL) -f $(MAKEFILE_EXE) exe
	$(MAKE) -j $(PARALLEL) -f $(MAKEFILE_EXE) unit_test

# Include Daylight-dependent progs, not built by default
rdock-dt-exe: rdock-lib
	./exemake -t"$(PLATFORM)" -c"$(CONFIG)" -o"$(MAKEFILE_EXE)"
	$(MAKE) -j $(PARALLEL) -f $(MAKEFILE_EXE) exe
	$(MAKE) -j $(PARALLEL) -f $(MAKEFILE_EXE) unit_test
	$(MAKE) -j $(PARALLEL) -f $(MAKEFILE_EXE) dt_exe

rdock-lib:
	./progen -t lib -o $(PROJECT)	'VERSION = $(VERSION)' \
					'CONFIG = dll warn_off $(CONFIG)' \
					'INCLUDEPATH = ../include;../include/GP;../import/simplex/include;../import/tnt/include' \
					'DEPENDPATH = $$INCLUDEPATH' \
					'DEFINES += $(DEFINES)' \
					'TARGET = Rbt' \
					'OBJECTS_DIR = $(BUILDDIR)/obj' \
					'DESTDIR = $(BUILDDIR)/lib' \
					$(SRC_FILES)
	./tmake 'TMAKEPATH = ./tmakelib/$(PLATFORM)' $(PROJECT) -o $(MAKEFILE_LIB)
	./buildnum.csh
	$(MAKE) -j $(PARALLEL) -f $(MAKEFILE_LIB)
