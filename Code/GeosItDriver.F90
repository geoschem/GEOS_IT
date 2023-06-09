!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GeosItDriver
!
! !DESCRIPTION: Program GeosItDriver is the top-level driver for the 
!  GEOS-5.7.x regridding programs.  GeosItDriver will call routines to 
!  extract, regrid, and save the GEOS-5.7.x met data to files for 
!  input to GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
PROGRAM GeosItDriver
!
! !USES:
!
  USE GeosItA1Module
  USE GeosItA3CldModule
  USE GeosItA3DynModule
  USE GeosItA3MstCModule
  USE GeosItA3MstEModule
  USE GeosItCnModule
  USE GeosItI3Module
  USE GeosItInputsModule
  USE GeosItRegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  20 Sep 2013 - R. Yantosca - Renamed Geos57 to GeosIt in routine names
!------------------------------------------------------------------------------

! Read filenames and fields to process from an input file
CALL GeosItInitialize

! Initialize GEOS-5 regridding code
CALL GeosItRegridInit

! Create the constant data file
IF ( doMakeCn ) CALL GeosItMakeCn

! Create the 1-hour average data file
CALL GeosItMakeA1

! Create the 3-hour average data files
CALL GeosItMakeA3Cld
CALL GeosItMakeA3Dyn
CALL GeosItMakeA3MstC
CALL GeosItMakeA3MstE

! Create the 6-hour instantaneous data file
CALL GeosItMakeI3

! Cleanup and quit 
CALL GeosItCleanup

END PROGRAM GeosItDriver
!EOP
