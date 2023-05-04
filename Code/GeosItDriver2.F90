!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GeosItDriver2
!
! !DESCRIPTION: Program GeosItDriver2 is the top-level driver for the 
!  GEOS-5.7.x regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM GeosItDriver2
!
! !USES:
!
  USE GeosItA3MstCModule
  USE GeosItA3MstEModule
  USE GeosItInputsModule
  USE GeosItRegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  GeosItDriver1 creates the A3mstC (3hr time-averaged moist parameters, on
!  level centers) and A3MstE (3hr time-averaged moist parameters on level
!  edges) data files for input into GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  20 Sep 2013 - R. Yantosca - Rename Geos57 to GeosIt in routine names
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL GeosItInitialize

  ! Initialize GEOS-5 regridding code
  CALL GeosItRegridInit

  ! Create the 3-hour average data file
  CALL GeosItMakeA3MstC
  CALL GeosItMakeA3MstE

  ! Cleanup and quit 
  CALL GeosItCleanup

END PROGRAM GeosItDriver2
!EOP
