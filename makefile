#-------------------------------------------------------------------------------
#
# makefile
#
# Purpose:
#
#   Linux make file for software libraries
#
# Last modified:
#
#   2022/12/09  AHA  Created
#   2024/11/12  AHA  Add build target for documentation
#   2024/11/12  AHA  Remove and re-create xsd schema file
#   2025/04/10  AHA  Avoid using rm -rf *
#
#-------------------------------------------------------------------------------

# Paths

GROOPS     = .
GROOPS_bin = $(GROOPS)/bin
GROOPS_doc = $(GROOPS)/docs
GROOPS_gui = $(GROOPS)/gui
GROOPS_bld = $(GROOPS)/source/build

# Get number of parallel build jobs

ifneq ($(shell which nproc 2> /dev/null),)
  NJOBS = $(shell nproc)
else
  NJOBS = 1
endif

# Parallel compilation on Linux and Cygwin

PMAKE = make
OS = $(shell uname -o)
ifeq ($(OS),GNU/Linux)
  PMAKE = make -j $(NJOBS)
endif
ifeq ($(OS),Cygwin)
  PMAKE = make -j $(NJOBS)
endif

# Operating system dependent settings for qmake
#
# Default: use standard qmake
QMAKE = qmake

# Linux: use qmake or qmake-qt5
ifeq ($(OS),GNU/Linux)
  ifeq ($(shell which qmake-qt5  2> /dev/null),)
    QMAKE = qmake
  else
    QMAKE = qmake-qt5  # if qmake is not available (e.g. openSuse)
  endif
endif

# Targets

all: init \
	 groops_ groopsgui_

# Create directory tree

init:
	if [ ! -d "$(GROOPS_bld)" ]; then mkdir -p $(GROOPS_bld); fi
	
groops_:
	cd $(GROOPS_bld); cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../..; $(PMAKE); $(PMAKE) install
	rm groops.xsd; groops --xsd groops.xsd
	

groopsgui_:
	cd $(GROOPS_gui); $(QMAKE); $(PMAKE)

# Documentation

doc:
	cd $(GROOPS_doc); ./makeDocumentation.sh

# Clean up

clean:
	if [ -d "$(GROOPS_bin)" ]; then cd $(GROOPS_bin); rm -f ./*; fi
	if [ -d "$(GROOPS_bld)" ]; then cd $(GROOPS_bld); $(PMAKE) clean; rm -rf ./*; fi
	cd $(GROOPS_gui); $(PMAKE) clean
