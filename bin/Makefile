#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (run-directory)
#
# !DESCRIPTION: Makefile for the GEOS-IT data processing code.  This is a
#  "router" makefile, which calls to the Makefile in the Code directory.
#\\
#\\
# !REMARKS:
# To build the program, call "make" with the following syntax:
#
#   make TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# !REVISION HISTORY: 
# See git history
#
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
SHELL = /bin/sh
DIR   = ../Code
EXE0   = GeosIt_Driver0.x
EXE1   = GeosIt_Driver1.x
EXE2   = GeosIt_Driver2.x

#==============================================================================
# Makefile targets (for now, build w/ the sample data)
#==============================================================================

.PHONY: all norun debug driver lib check clean realclean doc docclean help

all: 
	$(MAKE) -C $(DIR) all

lib:
	$(MAKE) -C $(DIR) lib

clean:
	$(MAKE) -C $(DIR) clean

realclean:
	$(MAKE) -C $(DIR) realclean

doc:
	$(MAKE) -C $(DIR) doc

docclean:
	$(MAKE) -C $(DIR) docclean

help:
	$(MAKE) -C $(DIR) help

driver: all

template:
	$(MAKE) -C $(DIR) -f Makefile_Template all

#EOC
