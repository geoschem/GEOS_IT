#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in doc subdirectory)
#
# !DESCRIPTION: Makefile for building the documentation (in PDF and 
#  PostScript formats) for the GEOS-5.7.x data processing code, Makefiles,
#  and Perl scripts.
#\\
#\\
# !REMARKS:
# To build the documentation, call "make" with the following syntax:
#
#   make TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# You must have the LaTeX utilities (latex, dvips, dvipdf) installed
# on your system in order to build the documentation.
#
# !REVISION HISTORY: 
#  12 Jan 2012 - R. Yantosca - Initial version
#  12 Jan 2012 - R. Yantosca - Split source code doc from netCDF code doc
#  24 Sep 2013 - R. Yantosca - Now add NCL scripts in adjust to the doc list
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
SHELL   :=/bin/bash
ADJDIR  :=../adjust
BINDIR  :=../bin
CODEDIR :=../Code
MAINDIR :=../
PERLDIR :=../perl

#==============================================================================
# Default targets
#==============================================================================

.PHONY: doc clean

doc: srcdoc ncdoc makedoc perldoc adjdoc

clean:
	rm -f *.tex *.ps *.pdf

#==============================================================================
# Build documentation for source code
#==============================================================================

# Source files
SRC1 :=                                \
./src.intro                            \
$(CODEDIR)/CreateTemplateFile.F90      \
$(CODEDIR)/CharpakModule.F	       \
$(CODEDIR)/GeosItA1Module.F90	       \
$(CODEDIR)/GeosItA3CldModule.F90       \
$(CODEDIR)/GeosItA3DynModule.F90       \
$(CODEDIR)/GeosItA3MstCModule.F90      \
$(CODEDIR)/GeosItA3MstEModule.F90      \
$(CODEDIR)/GeosItCnModule.F90	       \
$(CODEDIR)/GeosItI3Module.F90	       \
$(CODEDIR)/GeosItInputsModule.F90      \
$(CODEDIR)/GeosItRegridModule.F90      \
$(CODEDIR)/GeosItUtilityModule.F90     \
$(CODEDIR)/GeosItDriver.F90            \
$(CODEDIR)/GeosItDriver0.F90	       \
$(CODEDIR)/GeosItDriver1.F90	       \
$(CODEDIR)/GeosItDriver2.F90


# LaTeX etc. files
TEX1 :=GeosItDriver.tex
DVI1 :=GeosItDriver.dvi
PDF1 :=GeosItDriver.pdf
PS1  :=GeosItDriver.ps

srcdoc:
	rm -f $(TEX1)
	protex -sf $(SRC1) > $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	dvipdf $(DVI1) $(PDF1)
	dvips $(DVI1) -o $(PS1)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for netCDF library code
#==============================================================================

# Source files
SRC4 =                                 \
./ncdf.intro                           \
$(CODEDIR)/m_do_err_out.F90	       \
$(CODEDIR)/m_netcdf_io_checks.F90      \
$(CODEDIR)/m_netcdf_io_close.F90       \
$(CODEDIR)/m_netcdf_io_create.F90      \
$(CODEDIR)/m_netcdf_io_define.F90      \
$(CODEDIR)/m_netcdf_io_get_dimlen.F90  \
$(CODEDIR)/m_netcdf_io_handle_err.F90  \
$(CODEDIR)/m_netcdf_io_open.F90	       \
$(CODEDIR)/m_netcdf_io_read.F90	       \
$(CODEDIR)/m_netcdf_io_write.F90

# LaTeX etc. files
TEX4 :=NcdfUtilities.tex
DVI4 :=NcdfUtilities.dvi
PDF4 :=NcdfUtilities.pdf
PS4  :=NcdfUtilities.ps

ncdoc:
	rm -f $(TEX4)
	protex -sf $(SRC4) > $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	dvipdf $(DVI4) $(PDF4)
	dvips $(DVI4) -o $(PS4)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for Makefiles
#==============================================================================

# Source files 
SRC2 =                        \
./make.intro                  \
$(MAINDIR)/Makefile           \
$(CODEDIR)/Makefile           \
$(CODEDIR)/Makefile_header.mk \
./Makefile

# LaTeX etc. files
TEX2 = GeosItDriver_Makefiles.tex
DVI2 = GeosItDriver_Makefiles.dvi
PDF2 = GeosItDriver_Makefiles.pdf
PS2  = GeosItDriver_Makefiles.ps

makedoc:
	rm -f $(TEX2)
	protex -fS $(SRC2) > $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	dvipdf $(DVI2) $(PDF2)
	dvips $(DVI2) -o $(PS2)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for Perl scripts
#==============================================================================

# Source files 
SRC3 =                    \
./perl.intro              \
$(PERLDIR)/doGeosIt       \
$(PERLDIR)/doGeosIt.input \
$(PERLDIR)/getGeosIt      \
$(PERLDIR)/checkGeosIt    \
$(PERLDIR)/sleepGeosIt    \
$(PERLDIR)/schedGeosIt    \
$(PERLDIR)/cleanJobs      \
$(PERLDIR)/cleanLogs      \
$(PERLDIR)/delGeosIt      \
$(PERLDIR)/Dates.pm

# LaTeX etc. files
TEX3 = GeosItDriver_Perl.tex
DVI3 = GeosItDriver_Perl.dvi
PDF3 = GeosItDriver_Perl.pdf
PS3  = GeosItDriver_Perl.ps

perldoc:
	rm -f $(TEX3)
	protex -sfS $(SRC3) > $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	dvipdf $(DVI3) $(PDF3)
	dvips $(DVI3) -o $(PS3)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for NCL scripts
#==============================================================================

# Source files 
SRC4=                                 \
./adjust.intro                        \
$(ADJDIR)/adjust_all.ncl              \
$(ADJDIR)/add_coards_global_atts.ncl  \
$(ADJDIR)/add_coards_var_atts.ncl     \
$(ADJDIR)/adjust_a1.ncl               \
$(ADJDIR)/adjust_a3cld.ncl            \
$(ADJDIR)/adjust_a3dyn.ncl            \
$(ADJDIR)/adjust_a3mstc.ncl           \
$(ADJDIR)/adjust_a3mste.ncl           \
$(ADJDIR)/adjust_cn.ncl               \
$(ADJDIR)/adjust_i3.ncl               \
$(ADJDIR)/notdir.ncl                  \
$(ADJDIR)/adjust_apbp.ncl

# LaTeX etc. files
TEX4 = AdjustForCoards.tex
DVI4 = AdjustForCoards.dvi
PDF4 = AdjustForCoards.pdf
PS4  = AdjustForCoards.ps

adjdoc:
	rm -f $(TEX4)
	protex -sfN $(SRC4) > $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	dvipdf $(DVI4) $(PDF4)
	dvips $(DVI4) -o $(PS4)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
#  Help Screen
#==============================================================================

help:
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'doc      Builds documentation (*.ps, *.pdf) for both source code and makefiles (default)'
	@echo 'srcdoc   Builds documentation (*.ps, *.pdf) only for source code files'
	@echo 'makedoc  Builds documentation (*.ps, *.pdf) only for makefiles'
	@echo 'perldoc  Builds documentation (*.ps, *.pdf) only for Perl scripts'
	@echo 'adjdoc   Builds documentation (*.ps, *.pdf) only for NCL scripts'
	@echo 'clean    Removes all files (*.tex, *.ps, *.pdf)'
	@echo 'help     Displays this help screen'

#EOC
