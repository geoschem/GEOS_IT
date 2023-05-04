!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GeosItFixA3Cld
!
! !DESCRIPTION: Program GeosItFixA3Cld is the driver routine that we must
!  use to reprocess the A3Cld files.  We need to use the CLOUD field from
!  tavg3_3d_rad_Nv instead of taking it from CFAN + CFLS (which causes
!  higher O3 and OH in GEOS-Chem).
!\\
!\\
! !INTERFACE:
!
PROGRAM GeosItFixA3Cld
!
! !USES:
!
  USE GeosItA3CldModule
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

  ! Cleanup and quit 
  CALL GeosItCleanup

END PROGRAM GeosItFixA3Cld
!EOP
