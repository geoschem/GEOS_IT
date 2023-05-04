!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GeosItDriver1
!
! !DESCRIPTION: Program GeosItDriver1 is a top-level driver for the 
!  GEOS-5.7.x regridding programs. 
!\\
!\\
! !INTERFACE:
!
PROGRAM GeosItDriver1
!
! !USES:
!
  USE GeosItA3CldModule
  USE GeosItA3DynModule
  USE GeosItInputsModule
  USE GeosItRegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  GeosItDriver1 creates the A3cld (3hr time-averaged cloud parameters) and
!  A3dyn  (3hr time-averaged dynamics parameters) data files for input into 
!  GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  20 Sep 2013 - R. Yantosca - Now renamed Geos57 to GeosIt in routine names
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL GeosItInitialize

  ! Initialize GEOS-5 regridding code
  CALL GeosItRegridInit

  ! Create the 3-hour average data files
  CALL GeosItMakeA3Cld
  CALL GeosItMakeA3Dyn

  ! Cleanup and quit 
  CALL GeosItCleanup

END PROGRAM GeosItDriver1
!EOP
