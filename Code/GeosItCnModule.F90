!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GeosItCnModule
!
! !DESCRIPTION: Module GeosItCnModule contains routines to create the
!  GEOS-Chem constant data files from the GEOS-IT raw data.
!\\
!\\
! !INTERFACE:

MODULE GeosItCnModule
!
! !USES:
!
  ! GEOS-IT data modules
  USE CharpakModule
  USE GeosItInputsModule
  USE GeosItRegridModule
  USE GeosItUtilityModule

  ! Modules for writing netCDF
  USE m_netcdf_io_create
  USE m_netcdf_io_define
  USE m_netcdf_io_write
  USE m_netcdf_io_close

  ! Modules for reading netCDF
  USE m_netcdf_io_open
  USE m_netcdf_io_close
  USE m_netcdf_io_get_dimlen
  USE m_netcdf_io_read

  IMPLICIT NONE
  PRIVATE

  ! Include files
# include "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: GeosItMakeCn
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: GetNFields
  PRIVATE :: ProcessCn2dAsmNx
!
! !REMARKS:
!  NOTE: Hardwire the constant data file to 00:00 GMT on 1998/01/01.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  25 Oct 2011 - R. Yantosca - Initial version for GEOS-FP
!  07 Jun 2023 - E. Lundgren - Adapted for GEOS-IT
!  See git history for additional revision history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcOutFileDef
!
! !DESCRIPTION: Subroutine NcOutFileDef pre-defines variable names and
!  attributes that will be added to the netCDF output files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcOutFileDef( X,        Y,           T,    &
                           xMid,     YMid,        time, &
                           gridName, outFileName, fOut )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)    :: X             ! Longitude dimension
    INTEGER,          INTENT(IN)    :: Y             ! Latitude dimension
    INTEGER,          INTENT(IN)    :: T             ! Time dimension
    REAL*4,           INTENT(IN)    :: xMid(X)       ! Array of lon centers
    REAL*4,           INTENT(IN)    :: yMid(Y)       ! Array of lat centers
    INTEGER,          INTENT(IN)    :: time(T)       ! Array of times
    CHARACTER(LEN=*), INTENT(IN)    :: gridName      ! Name of the grid
    CHARACTER(LEN=*), INTENT(IN)    :: outFileName   ! Output file name
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: fOut          ! Output netCDF file ID
!
! !REMARKS:
!  NOTE: Hardwire the constant data file to date 1998/01/01; 00:00 GMT.

! !REVISION HISTORY:
!  25 Oct 2011 - R. Yantosca - Initial version for GEOS-FP
!  07 Jun 2023 - E. Lundgren - Adapted for GEOS-IT
!  See git history for additional revision history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255) :: sysTime
    CHARACTER(LEN=255) :: lName,   units,   gamap,   DI,   DJ
    CHARACTER(LEN=255) :: delta_t, begin_d, begin_t, incr, msg,  cal
    INTEGER            :: idLon,   idLat,   idTime,  vId,  oMode
    LOGICAL            :: is_nc4
    ! Arrays
    INTEGER            :: var1(1), var3(3)

    !=========================================================================
    ! %%% BEGINNING OF NETCDF DEFINITION SECTION %%%
    !=========================================================================

    ! Echo info
    msg = '%%%%%% ENTERING ROUTINE NcOutFileDef %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Echo info
    WRITE( 6, 100 ) TRIM( gridName )
100 FORMAT ( '%%% Defining netCDF file vars & attrs for ', a' grid' )

    is_nc4 = .TRUE.

    ! Open netCDF file for writing
    CALL NcCr_Wr( fOut, TRIM( outFileName ), WRITE_NC4=is_nc4 )

    ! Turn filling off
    CALL NcSetFill( fOut, NF_NOFILL, oMode )

    !-------------------------------------------------------------------------
    ! Define global attributes and filling mode
    !-------------------------------------------------------------------------

    ! Title string
    lName = 'GEOS-IT constant parameters (CN), processed for GEOS-Chem input'
    CALL NcDef_Glob_Attributes( fOut, 'Title',                TRIM( lName ) )

    ! Contact
    lName = "GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)"
    CALL NcDef_Glob_Attributes( fOut, 'Contact',              TRIM( lName ) )

    ! References
    lName = "www.geos-chem.org; wiki.geos-chem.org"
    CALL NcDef_Glob_Attributes( fOut, 'References',           TRIM( lName ) )

    ! Filename
    lName = NotDir( outFileName )
    CALL NcDef_Glob_Attributes( fOut, 'Filename',             TRIM( lName ) )

    ! History
    sysTime = SystemTimeStamp()
    lName = 'File generated on: ' // TRIM( sysTime )
    CALL NcDef_Glob_Attributes( fOut, 'History' ,             TRIM( lName ) )
    CALL NcDef_Glob_Attributes( fOut, 'ProductionDateTime',   TRIM( lName ) )
    CALL NcDef_Glob_Attributes( fOut, 'ModificationDateTime', TRIM( lName ) )

    ! Format
    lName = "NetCDF-4" ;
    CALL NcDef_Glob_Attributes( fOut, 'Format' ,              TRIM( lName ) )

    ! Format
    lName = "global" ;
    CALL NcDef_Glob_Attributes( fOut, 'SpatialCoverage',      TRIM( lName ) )

    ! Conventions
    lName = 'COARDS'
    CALL NcDef_Glob_Attributes( fOut, 'Conventions',          TRIM( lName ) )

    ! Version
    lName = 'GEOS-IT'
    CALL NcDef_Glob_Attributes( fOut, 'Version',              TRIM( lName ) )

    ! Model
    lName = 'GEOS-5'
    CALL NcDef_Glob_Attributes( fOut, 'Model',                TRIM( lName ) )

    ! NLayers
    lName = '72'
    CALL NcDef_Glob_Attributes( fOut, 'Nlayers',              TRIM( lName ) )

    ! Start Date (hardwire to 1998/01/01)
    lName = '19980101'
    CALL NcDef_Glob_Attributes( fOut, 'Start_Date',           TRIM( lName ) )

    ! Start Time
    lName = '00:00:00.0'
    CALL NcDef_Glob_Attributes( fOut, 'Start_Time',           TRIM( lName ) )

    ! End Date (hardwire to 1998/01/01)
    lName = '19980101'
    CALL NcDef_Glob_Attributes( fOut, 'End_Date',             TRIM( lName ) )

    ! End Time
    lName = '00:00:00.0'
    CALL NcDef_Glob_Attributes( fOut, 'End_Time',             TRIM( lName ) )

    ! Delta-time
    lName = '000000'
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Time',           TRIM( lName ) )

    ! Pick DI and DJ attributes based on the grid
    SELECT CASE ( TRIM( gridName ) )
       CASE( 'nested EU 05', 'nested NA 05', 'nested AS 05', '0.5 x 0.625 global')
          DI = '0.625'
          DJ = '0.5'
       CASE( '2 x 2.5 global' )
          DI = '2.5'
          DJ = '2'
       CASE( '4 x 5 global' )
          DI = '5'
          DJ = '4'
    END SELECT

    ! Delta-lon
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lon',            TRIM( DI    ) )

    ! Delta-lat
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lat',            TRIM( DJ    ) )

    !-------------------------------------------------------------------------
    ! Define dimensions and index arrays.  NOTE: COARDS specifies that index
    ! arrays will have the same names as the dimensions that define them.
    !-------------------------------------------------------------------------

    ! netCDF dimension variables
    CALL NcDef_Dimension( fOut, 'time', 1, idTime )
    CALL NcDef_Dimension( fOut, 'lat',  Y, idLat  )
    CALL NcDef_Dimension( fOut, 'lon',  X, idLon  )

    ! Time index array (hardwire date to 1998/01/01)
    var1    = (/ idTime /)
    cal     = 'gregorian'
    lName   = 'time'
    units   = UnitsForTime( 19980101 )
    delta_t = '0000-00-00 00:00:00'
    begin_d = '19980101'
    begin_t = '000000'
    incr    = '000000'
    CALL NcDef_Variable      ( fOut, 'time', NF_INT,  1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'calendar',       TRIM( cal     )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'delta_t',        TRIM( delta_t )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_date',     TRIM( begin_d )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_time',     TRIM( begin_t )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'time_increment', TRIM( incr    )  )

    ! Latitude index array
    var1    = (/ idLat /)
    lName   = 'latitude'
    units   = 'degrees_north'
    CALL NcDef_Variable      ( fOut, 'lat', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName   )  )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    )

    ! Longitude index array
    var1    = (/ idLon /)
    lName   = 'longitude'
    units   = 'degrees_east'
    CALL NcDef_Variable      ( fOut, 'lon', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units )    )

    !-------------------------------------------------------------------------
    ! Define data arrays
    !-------------------------------------------------------------------------

    ! FRLAKE
    IF ( StrPos( 'FRLAKE', asm_const_0hr_slv_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)
       
       lName = 'Fraction of lake type in grid box'
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRLAKE', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FRLAND
    IF ( StrPos( 'FRLAND', asm_const_0hr_slv_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)
       
       lName = 'Fraction of land in grid box'
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRLAND', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FRLANDICE
    IF ( StrPos( 'FRLANDIC', asm_const_0hr_slv_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)
       
       lName = 'Fraction of land ice in grid box'
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRLANDIC', NF_FLOAT, 3, var3, vId  )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FROCEAN
    IF ( StrPos( 'FROCEAN', asm_const_0hr_slv_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)
       
       lName = 'Fraction of ocean in grid box'
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FROCEAN', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PHIS
    IF ( StrPos( 'PHIS', asm_const_0hr_slv_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)
       
       lName = 'Surface geopotential'
       units = 'm2 s-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PHIS', NF_FLOAT, 3, var3, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    !=========================================================================
    ! %%% END OF NETCDF DEFINITION SECTION %%%
    !=========================================================================

    ! End the definition section
    CALL NcEnd_def( fOut )

    ! Write index arrays
    CALL NcWr( xMid, fOut, 'lon',  (/ 1 /), (/ X /) )
    CALL NcWr( yMid, fOut, 'lat',  (/ 1 /), (/ Y /) )
    CALL NcWr( time, fOut, 'time', (/ 1 /), (/ T /) )

    ! Echo info
    msg = '%%%%%% LEAVING ROUTINE NcOutFileDef %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE NcOutFileDef
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GeosItMakeCn
!
! !DESCRIPTION: Routine GeosItMakeCn is the the driver routine for
! \begin{enumerate}
! \item Extracting constant data fields (surface values) from
!       the GEOS-IT raw data files (netCDF-4 format)
! \item Regridding the fields to GEOS-Chem data resolution, and
! \item Saving the regridded data to netCDF format.
! \end{enumerate}
! This routine is called directly from the main program GeosItDriver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosItMakeCn
!
! !REVISION HISTORY:
!  27 Jul 2010 - R. Yantosca - Initial version for GEOS-FP
!  07 Jun 2023 - E. Lundgren - Adapted for GEOS-IT
!  See git history for additional revision history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: F
    INTEGER                 :: nFields
    INTEGER                 :: nAllFields
    CHARACTER(LEN=MAX_CHAR) :: msg
    CHARACTER(LEN=MAX_CHAR) :: fName
    CHARACTER(LEN=MAX_CHAR) :: gName

    ! Arrays
    INTEGER                 :: time(1)
    CHARACTER(LEN=MAX_CHAR) :: fields(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info
    msg = '%%%%%%%%%% ENTERING ROUTINE GeosItMakeCn %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Return the list of fields and number of fields to process
    ! from each of the MERRA raw met data files
    CALL GetNFields( asm_const_0hr_slv_data, nFields, fields )

    ! Total number of fields that we will process
    nAllFields = nFields

    ! Echo info
    WRITE( IU_LOG, 100 ) TRIM( asm_const_0hr_slv_file ), nFields
    WRITE( IU_LOG, 110 ) nAllFields

    ! Formats
100 FORMAT( '%%% # of fields from ', a, ' : ',          i5 )
110 FORMAT( '%%% TOTAL # OF FIELDS TO BE REGRIDDED : ', i5 )

    !=======================================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Hours in the day
    time = (/ 0 /)

    ! Open 0.5x0.625 output file
    IF ( do05x0625 ) THEN
      fName = TRIM( tempDirTmpl05x0625 ) // TRIM( dataTmpl05x0625 )
      gName = '0.5 x 0.625 global'
      CALL ExpandDate  ( fName,     19980101,     000000      )
      CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
      CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
      CALL NcOutFileDef( I05x0625,     J05x0625,        1,    &
                         xMid_05x0625, nc_yMid_05x0625, time,      &
                         gName,     fName,        fOut05x0625    )
   ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     19980101,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.           )
       CALL NcOutFileDef( I2x25,     J2x25,        1,           &
                          xMid_2x25, nc_yMid_2x25, time,        &
                          gName,     fName,        fOut2x25    )
    ENDIF

    ! Open 4 x 5 output file
    IF ( do4x5 ) THEN
       fName = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     19980101,  000000         )
       CALL StrRepl     ( fName,     '%%%%%%',  'CN    '       )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I4x5,      J4x5,      1,              &
                          xMid_4x5,  nc_yMid_4x5,  time,        &
                          gName,     fName,     fOut4x5        )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu05 ) THEN
       fName = TRIM( tempDirTmplNestEu05 ) // TRIM( dataTmplNestEu05 )
       gName = 'nested EU 05'
       CALL ExpandDate  ( fName,     19980101,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestEu05,  J_NestEu05,     1,       &
                          xMid_05x0625(I0_eu05:I1_eu05),          &
                          yMid_05x0625(J0_eu05:J1_eu05),          &
                          time,    gName,        fName,       &
                          fOut05NestEu                           )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa05 ) THEN
       fName = TRIM( tempDirTmplNestNa05 ) // TRIM( dataTmplNestNa05 )
       gName = 'nested NA 05'
       CALL ExpandDate  ( fName,     19980101,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I_NestNa05,  J_NestNa05,  1,              &
                          xMid_05x0625(I0_na05:I1_na05),          &
                          yMid_05x0625(J0_na05:J1_na05),          &
                          time,      gName,        fName,       &
                          fOut05NestNa                            )
    ENDIF

    ! Open nested AS output file
    IF ( doNestAs05 ) THEN
       fName = TRIM( tempDirTmplNestAs05 ) // TRIM( dataTmplNestAs05 )
       gName = 'nested AS 05'
       CALL ExpandDate  ( fName,     19980101,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I_NestAs05,  J_NestAs05,  1,              &
                          xMid_05x0625(I0_as05:I1_as05),          &
                          yMid_05x0625(J0_as05:J1_as05),          &
                          time,      gName,        fName,       &
                          fOut05NestAs                            )
    ENDIF

    !=======================================================================
    ! Process data
    !=======================================================================

    ! Regrid fields from the various raw data files
    CALL ProcessCn2dAsmNx( nFields, fields )

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing CN output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( do05x0625  ) CALL NcCl( fOut05x0625  )
    IF ( do2x25     ) CALL NcCl( fOut2x25     )
    IF ( do4x5      ) CALL NcCl( fOut4x5      )
    IF ( doNestEu05 ) CALL NcCl( fOut05NestEu )
    IF ( doNestNa05 ) CALL NcCl( fOut05NestNa )
    IF ( doNestAs05 ) CALL NcCl( fOut05NestAs )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE GeosItMakeCn %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE GeosItMakeCn
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ProcessCn2dAsmNx
!
! !DESCRIPTION: Subroutine ProcessCn2dAsmNx regrids the GEOS-IT met fields
!  from the "asm\_const\_0hr\_slv" file and saves output to netCDF format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ProcessCn2dAsmNx( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY:
!  04 Jan 2012 - R. Yantosca - Initial version for GEOS-FP
!  07 Jun 2023 - E. Lundgren - Adapted for GEOS-IT
!  See git history for additional revision history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop variables
    INTEGER                 :: F

    ! Variables for netCDF I/O
    INTEGER                 :: X,        Y,        T

    INTEGER                 :: X05x0625, Y05x0625, T05x0625
    INTEGER                 :: X2x25,    Y2x25,    T2x25
    INTEGER                 :: X4x5,     Y4x5,     T4x5

    INTEGER                 :: XNestEu05,  YNestEu05,  TNestEu05
    INTEGER                 :: XNestNa05,  YNestNa05,  TNestNa05
    INTEGER                 :: XNestAs05,  YNestAs05,  TNestAs05

    INTEGER                 :: st2d(2),  st3d(3)
    INTEGER                 :: ct2d(2),  ct3d(3)

    ! Data arrays
    REAL*4, TARGET          :: Q    ( I05x0625, J05x0625, 1     )
    REAL*4                  :: Q2x25( I2x25,    J2x25           )
    REAL*4                  :: Q4x5 ( I4x5,     J4x5            )

    ! Pointers
    REAL*4,  POINTER        :: QNest(:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name8
    CHARACTER(LEN=9       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

   ! 0.5x0.625 global grid
    IF ( do05x0625 ) THEN
      CALL NcGet_DimLen( fOut05x0625,   'lon',  X05x0625   )
      CALL NcGet_DimLen( fOut05x0625,   'lat',  Y05x0625   )
      CALL NcGet_DimLen( fOut05x0625,   'time', T05x0625   )
    ENDIF

    ! 2x2.5 global grid
    IF ( do2x25 ) THEN
       CALL NcGet_DimLen( fOut2x25,   'lon',  X2x25   )
       CALL NcGet_DimLen( fOut2x25,   'lat',  Y2x25   )
       CALL NcGet_DimLen( fOut2x25,   'time', T2x25   )
    ENDIF

    ! 4x5 global grid
    IF ( do4x5 ) THEN
       CALL NcGet_DimLen( fOut4x5,    'lon',  X4x5    )
       CALL NcGet_DimLen( fOut4x5,    'lat',  Y4x5    )
       CALL NcGet_DimLen( fOut4x5,    'time', T4x5    )
    ENDIF

    ! Nested EU grid 0625
    IF ( doNestEu05 ) THEN
       CALL NcGet_DimLen( fOut05NestEu, 'lon',  XNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'lat',  YNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'time', TNestEu05 )
    ENDIF

    ! Nested NA grid 0625
    IF ( doNestNa05 ) THEN
       CALL NcGet_DimLen( fOut05NestNa, 'lon',  XNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'lat',  YNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'time', TNestNa05 )
    ENDIF

    ! Nested AS grid 0625
    IF ( doNestAs05 ) THEN
       CALL NcGet_DimLen( fOut05NestAs, 'lon',  XNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'lat',  YNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'time', TNestAs05 )
    ENDIF

    !=======================================================================
    ! Open input file
    ! NOTE: For constant file, hardwire date to 1998/01/01
    !=======================================================================

    ! Echo info
    msg = '%%%%%% ENTERING ROUTINE ProcessCn2dAsmNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Create input filename from the template
    fNameInput = TRIM( inputDataDir ) // TRIM( asm_const_0hr_slv_file )
    CALL expandDate( fNameInput, yyyymmdd, 000000 )

    ! Echo info
    msg = '%%% Opening ' // TRIM( fNameInput )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Open the netCDF4 file for input
    CALL NcOp_Rd( fIn, TRIM( fNameInput ) )

    ! Get the dimensions from the netCDF file
    CALL NcGet_DimLen( fIn, 'lon',  X )
    CALL NcGet_DimLen( fIn, 'lat',  Y )
    CALL NcGet_DimLen( fIn, 'time', T )

    !=======================================================================
    ! Process data
    !=======================================================================

    ! Loop over data fields
    DO F = 1, nFields

       ! Save field name into an 9-char variable.
       ! This will truncate field names longer than 8 chars.
       name  = TRIM( fields(F) )

       ! Skip if fieldname is empty
       IF ( name == '' ) CYCLE

       ! Zero data arrays
       Q     = 0e0
       Q2x25 = 0e0
       Q4x5  = 0e0

       !-----------------------------
       ! Read data
       !-----------------------------
       msg = '%%% Reading     ' // name
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Start and count index arrays for netCDF
       ! (There is only one data block per file)
       st3d = (/ 1, 1, 1 /)
       ct3d = (/ X, Y, T /)

       ! Read data from file
       CALL NcRd( Q, fIn, TRIM( name ), st3d, ct3d )

       ! Replace missing values with zeroes
       WHERE( Q == FILL_VALUE ) Q = 0e0

       ! Do not need to flip in vertical since 2d fields

       !-----------------------------
       ! Regrid data
       !-----------------------------
       msg = '%%% Regridding  ' // name
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Regrid to 2 x 2.5
       IF ( do2x25 ) THEN
          CALL RegridGeosItto2x25( 0, Q(:,:,1), Q2x25 )
       ENDIF

       ! Regrid to 4x5
       IF ( do4x5 ) THEN
          CALL RegridGeosItTo4x5( 0, Q(:,:,1), Q4x5 )
       ENDIF

       !-----------------------------
       ! Write netCDF output
       !-----------------------------

       msg = '%%% Archiving   ' // name
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Special handing
       IF ( TRIM( name ) == 'FRLANDICE' ) name ='FRLANDIC'

       ! Write 0.5x0.625 data
       IF ( do05x0625 ) THEN
          st3d = (/ 1,     1,     1     /)
          ct3d = (/ X05x0625, Y05x0625, T05x0625 /)
          CALL NcWr( Q, fOut05x0625, TRIM( name ), st3d, ct3d )
       ENDIF

       ! Write 2x2.5 data
       IF ( do2x25 ) THEN
          st3d = (/ 1,     1,     1     /)
          ct3d = (/ X2x25, Y2x25, T2x25 /)
          CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st3d, ct3d )
       ENDIF

       ! Write 4x5 data
       IF ( do4x5 ) THEN
          st3d = (/ 1,    1,    1    /)
          ct3d = (/ X4x5, Y4x5, T4x5 /)
          CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st3d, ct3d )
       ENDIF

       ! Nested EU
       IF ( doNestEu05 ) THEN
          QNest => Q( I0_eu05:I1_eu05, J0_eu05:J1_eu05, 1 )
          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestEu05, YNestEu05, TNestEu05 /)
          CALL NcWr( QNest, fOut05NestEu, TRIM( name ), st3d, ct3d )
          NULLIFY( QNest )
       ENDIF

       ! Nested NA
       IF ( doNestNa05 ) THEN
          QNest => Q( I0_na05:I1_na05, J0_na05:J1_na05, 1 )
          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestNa05, YNestNa05, TNestNa05 /)
          CALL NcWr( QNest, fOut05NestNa, TRIM( name ), st3d, ct3d )
          NULLIFY( QNest )
       ENDIF

       ! Nested AS
       IF ( doNestAs05 ) THEN
          QNest => Q( I0_as05:I1_as05, J0_as05:J1_as05, 1 )
          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestAs05, YNestAs05, TNestAs05 /)
          CALL NcWr( QNest, fOut05NestAs, TRIM( name ), st3d, ct3d )
          NULLIFY( QNest )
       ENDIF

    ENDDO

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close input file
    msg = '%%% Closing ' // TRIM( fNameInput )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )
    CALL NcCl( fIn )

    ! Echo info
    msg = '%%%%%% EXITING ROUTINE ProcessCn2dAsmNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE ProcessCn2dAsmNx
!EOC
END MODULE GeosItCnModule
