
BOUT_TOP	= ../..

TARGET = sd1d

DIRS = atomicpp

SOURCEC		= sd1d.cxx div_ops.cxx loadmetric.cxx radiation.cxx

# Capture the git version, to be printed in the outputs
GIT_VERSION := $(shell git describe --abbrev=40 --dirty --always --tags)
CXXFLAGS += -DSD1D_REVISION=\"$(GIT_VERSION)\"

include $(BOUT_TOP)/make.config

# Note: Need to create revision.hxx
# Modifying all, $(TARGET) or sd1d targets doesn't work,
# but target depends on makefile
makefile: revision.hxx
revision.hxx: revision.hxx.in
	cp $< $@
