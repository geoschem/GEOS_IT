#EOC
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in Code subdirectory)
#
# !DESCRIPTION: Makefile for the GEOS-IT data processing code.
#\\
#\\
# !REMARKS:
# To build the program, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# BIN        Specifies the directory where the executable will be sent
# DOC        Specifies the directory where documentation will be built
# EXE        Specifies the name of the executable file
# F90        Contains the default Fortran compilation commands
# LD         Contains the command to link to libraries & make executable
# LIB        Specifies the directory where the libNcUtils.a files will be sent
# LINK_NC    Contains the comand to link to the netCDF library
# MOD        Specifies the directory where module files (*.mod) will be sent
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
#                                                                             .
# F90, LD, LINK_NC are included from "Makefile_header.mk"
#
# !REVISION HISTORY: 
#  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
SHELL   :=/bin/bash
BIN     :=../bin
DOC     :=../doc
EXE     :=GeosItDriver.x
EXE0    :=GeosItDriver0.x
EXE1    :=GeosItDriver1.x
EXE2    :=GeosItDriver2.x
LIB     :=../lib
MOD     :=../mod

# Include header file.  This returns variables F90, LD, LINK_NC
# as well as the default Makefile compilation rules for source code files.
include ./Makefile_header.mk

#==============================================================================
# Filename section
#
# NOTES:
# (1) The function $(wildcard m_*.F90) returns all files matching m_*.F90, etc.
# (2) The syntax $(SRC:.F90=.o) is a shorthand for the "patsubst" function. 
#==============================================================================

#----------------------
# Objects
#----------------------

# Source code & object files
SRC  := $(wildcard *.F) $(wildcard *.F90)            # Get these files
SRC  := $(filter-out GeosItDriver.F90  \
                     GeosItDriver0.F90 \
                     GeosItDriver1.F90 \
                     GeosItDriver2.F90 \
                     CreateTemplateFile.F90,$(SRC))  # but skip these files
TMP  := $(SRC:.F=.o)
OBJ  := $(TMP:.F90=.o)

# Object file for executable
OBJx  := $(EXE:.x=.o)
OBJx0 := $(EXE0:.x=.o)
OBJx1 := $(EXE1:.x=.o)
OBJx2 := $(EXE2:.x=.o)

#-----------------------
# Dependency listing
#-----------------------
m_do_err_out.o           : m_do_err_out.F90

m_netcdf_io_close.o      : m_netcdf_io_close.F90        \
                           m_do_err_out.o

m_netcdf_io_create.o     : m_netcdf_io_create.F90       \
                           m_do_err_out.o

m_netcdf_io_define.o     : m_netcdf_io_define.F90       \
                           m_do_err_out.o

m_netcdf_io_get_dimlen.o : m_netcdf_io_get_dimlen.F90   \
                           m_do_err_out.o

m_netcdf_io_handle_err.o : m_netcdf_io_handle_err.F90   \
                           m_do_err_out.o

m_netcdf_io_open.o       : m_netcdf_io_open.F90         \
                           m_do_err_out.o

m_netcdf_io_read.o       : m_netcdf_io_read.F90         \
                           m_do_err_out.o

m_netcdf_io_write.o      : m_netcdf_io_write.F90        \
                           m_do_err_out.o

CharpakModule.o          : CharpakModule.F

GeosItInputsModule.o    : GeosItInputsModule.F90       \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

GeosItRegridModule.o    : GeosItRegridModule.F90       \
                           GeosItInputsModule.o

GeosItUtilityModule.o   : GeosItUtilityModule.F90      \
                           CharpakModule.o

GeosItA1Module.o        : GeosItA1Module.F90           \
                           CharpakModule.o              \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItUtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

GeosItA3CldModule.o     : GeosItA3CldModule.F90        \
                           CharpakModule.o              \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItUtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

GeosItA3DynModule.o     : GeosItA3DynModule.F90        \
                           CharpakModule.o              \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItUtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

GeosItA3MstCModule.o    : GeosItA3MstCModule.F90       \
                           CharpakModule.o              \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItUtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

GeosItA3MstEModule.o    : GeosItA3MstEModule.F90       \
                           CharpakModule.o              \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItUtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

GeosItCnModule.o        : GeosItCnModule.F90           \
                           CharpakModule.o              \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItUtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

GeosItI3Module.o        : GeosItI3Module.F90           \
                           CharpakModule.o              \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItUtilityModule.o        \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o

GeosItDriver.o          : GeosItDriver.F90             \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItA1Module.o             \
                           GeosItA3CldModule.o          \
                           GeosItA3DynModule.o          \
                           GeosItA3MstCModule.o         \
                           GeosItA3MstEModule.o         \
                           GeosItCnModule.o             \
                           GeosItI3Module.o

GeosItDriver0.o         : GeosItDriver0.F90            \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItA1Module.o             \
                           GeosItCnModule.o             \
                           GeosItI3Module.o

GeosItDriver1.o         : GeosItDriver1.F90            \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItA3CldModule.o          \
                           GeosItA3DynModule.o   

GeosItDriver2.o         : GeosItDriver2.F90            \
                           GeosItInputsModule.o         \
                           GeosItRegridModule.o         \
                           GeosItA3MstCModule.o         \
                           GeosItA3MstEModule.o

#==============================================================================
# Makefile targets
#==============================================================================

.PHONY: all doc docclean clean realclean help test driver

all: lib exe

lib: $(OBJ)
	$(AR) crs libMerra2.a $(OBJ)
	mv libMerra2.a $(LIB)

exe: $(OBJx) $(OBJx) $(OBJx0) $(OBJx1) $(OBJx2)
	$(LD) $(OBJx) -L$(LIB) -lMerra2 $(LINK_NC) -o $(EXE)
	$(LD) $(OBJx0) -L$(LIB) -lMerra2 $(LINK_NC) -o $(EXE0)
	$(LD) $(OBJx1) -L$(LIB) -lMerra2 $(LINK_NC) -o $(EXE1)
	$(LD) $(OBJx2) -L$(LIB) -lMerra2 $(LINK_NC) -o $(EXE2)
	cp -f $(EXE) $(BIN)
	cp -f $(EXE0) $(BIN)
	cp -f $(EXE1) $(BIN)
	cp -f $(EXE2) $(BIN)

doc:
	$(MAKE) -C $(DOC) doc

docclean:
	$(MAKE) -C $(DOC) clean

clean:
	rm -f *.o *.nc *.x *.a *.mod *.nc

realclean:
	$(MAKE) clean
	$(MAKE) docclean
	rm -f $(BIN)/*.x
	rm -f $(LIB)/*.a
	rm -f $(MOD)/*.mod

driver: all

#==============================================================================
# Help Screen
#==============================================================================

help:
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'lib        Compiles code and creates the *.mod, *.o, and *.a files'
	@echo 'exe        Builds the executable file'
	@echo 'all        Synonym for "lib exe" (default target)'

	@echo 'check      Compiles code and runs the executable using GMAO sample data (for debugging)'
	@echo 'clean      Removes *.o *.mod and executable files in Code subdirectory'
	@echo 'realclean  Makes clean, then removes contents of bin, lib, doc subdirectories'
	@echo 'doc        Builds documentation files (*.ps, *.pdf) in doc subdirectory'
	@echo 'docclean   Removes documentation files from doc subdirectory'
	@echo 'help       Displays this help screen'
	@echo ''
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'COMPILER=___   Options: ifort pgi sun xlf (default is ifort)'
	@echo 'DEBUG=yes      Compiles GEOS-Chem for use w/ a debugger'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debugging purposes)'
	@echo 'TRACEBACK=yes  Compiles with -traceback flag (ifort only)'
#
#EOC
