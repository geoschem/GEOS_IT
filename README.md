[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/geoschem/GEOS_IT/blob/master/LICENSE.txt)

This package was originally developed by Bob Yantosca for GEOS-FP data processing. It is
adapted for GEOS-IT by Lizzie Lundgren. README text is taken from GEOS-FP data processing
tutorial written by Bob Yantosca in 2014.

GEOS_IT
=======

Development of the code and scripts used to process the "raw" GEOS-IT met data for input into GEOS-Chem

Root-level subdirectories
=========================

Code/   : Contains Fortran code
adjust/ : Contains NCL code (not used much anymore)
bin/    : Contains executable and input files that specify data directories and
          options for data reprocessing
doc/    : Manual pages will be built here when you type 'make doc'
jobs/   : Contains job scripts created by the doGeosItMulti driver script
lib/    : Fortran library files (*.a) will be built here during compilation
logs/   : Log files fromt he GOEs-IT data processing jobs are sent here
mod/    : Fortran module files (*.mod) will be built here during compilation
perl/   : Contains scripts (written in Perl) to download and processing GEOS-IT data

The most important directories here to be aware of are Code, bin, logs, and perl.


Scripts in perl/ subdirectory
=============================

There are several scripts in the perl/ subdirectory that control the entire data
download and regridding process.

Scripts for downloading data:
   getGeosIt   : Downloads 1 day of GEOS-IT met data from NASA
   checkGeosIT : Checks to see if the GEOS-IT data files were downloaded

Scripts for file management:
   cleanJobs   : Removes job scripts in the jobs/ subdirectory
   cleanLogs   : Removes log files in the logs/ subdirectory
   moveGeosIt  : Moves files from temp dir to data dir (called by doGeosIt)
   delGeosIt   : Removes GeosIt "raw" met data files (for manual use)
   purgeGeosIt : Removes old GEOS-IT "raw" data (called by doGeosIt)

Scripts for regridding data:
   Dates.pm       : Perl modules containing common subroutines
   doGeosItMulti  : Main driver script; extracts/regrids 1 day of GEOS-IT data
   doGeosIt.input : Input file with settings for doGeosItMulti
   runMet         : Called by doGeosIt


Scripts for automating the data download and regridding process:
   schedGeosIt           : Calls purgeGeosIt, getGeosIt, and doGeosIt
   schedGeosItInAdvance  : Schedules regridding jobs in days advance

