!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GeosItDriver0
!
! !DESCRIPTION: Program GeosItDriver0 is a top-level driver for the 
!  GEOS-5.7.x regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM GeosItDriver0
!
! !USES:
!
  USE GeosItA1Module
  USE GeosItCnModule
  USE GeosItI3Module
  USE GeosItInputsModule
  USE GeosItRegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  GeosItDriver1 creates the CN (constant), A1 (1hr time average), and
!  I3 (3hr instantaneous) data files for input into GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  20 Sep 2013 - R. Yantosca - Now renamed Geos57 to GeosIT in routine names
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL GeosItInitialize

  ! Initialize GEOS-5 regridding code
  CALL GeosItRegridInit

  ! Create the constant data file
  IF ( doMakeCn ) CALL GeosItMakeCn

  ! Create the 1-hour average data file
  CALL GeosItMakeA1

  ! Create the 3-hour instantaneous files
  CALL GeosItMakeI3

  ! Cleanup and quit 
  CALL GeosItCleanup

END PROGRAM GeosItDriver0
!EOP
