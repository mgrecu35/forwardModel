!===============================================================================
!
!---Function: Performs a prel;iminary check of all files used as inputs to the 
!             L2 of Combined Algorithm code for GPM; reduces possibilities
!             of errant file entries crashing program-promotes smooth exits.
!
!   status values
!      0 = opened ok
!     -1 = could not find file
!     -2 = nil file found, unavailable
!     -3 = empty granule
!     -4 = not an HDF5 file
!
!===============================================================================
SUBROUTINE file_check (igmi1, igmi2, igmi3, i2aku, i2adpr, i2akuenv, isnow,    &
                       iseaice, i2cmb, file_gmi1, file_gmi2, file_gmi3,        &
		       file_2adpr, file_2akuenv, file_2aku, file_snow,         &
		       file_seaice, file_2cmb, stat_gmi1, stat_gmi2,           &
		       stat_gmi3, stat_2akuenv, stat_2adpr, stat_2aku,         &
		       stat_snow, stat_seaice, stat_cmb, main_date, main_orbit) 
    
    IMPLICIT NONE
    
!....Calling sequence parameters  
    INTEGER*4,INTENT(in)   :: igmi1          ! length of file name for cgmi1     
    INTEGER*4,INTENT(in)   :: igmi2          ! length of file name for cgmi2
    INTEGER*4,INTENT(in)   :: igmi3          ! length of file name for cgmi3
    INTEGER*4,INTENT(in)   :: i2aku          ! length of file name for 2aku
    INTEGER*4,INTENT(in)   :: i2adpr         ! length of file name for 2adpr
    INTEGER*4,INTENT(in)   :: i2akuenv       ! length of file name for 2akuenv
    INTEGER*4,INTENT(in)   :: isnow          ! length of file name for snow
    INTEGER*4,INTENT(in)   :: iseaice        ! length of file name for seaice
    INTEGER*4,INTENT(in)   :: i2cmb          ! length of file name for cmb
    CHARACTER*1,INTENT(in) :: file_gmi1(1000)    ! full name of cgmi1 file
    CHARACTER*1,INTENT(in) :: file_gmi2(1000)    ! full name of cgmi3 file
    CHARACTER*1,INTENT(in) :: file_gmi3(1000)    ! full name of cgmi3 file
    CHARACTER*1,INTENT(in) :: file_2adpr(1000)   ! full name of 2adpr file
    CHARACTER*1,INTENT(in) :: file_2akuenv(1000) ! full name of 2akuenv file
    CHARACTER*1,INTENT(in) :: file_2aku(1000)    ! full name of 2aku file
    CHARACTER*1,INTENT(in) :: file_snow(1000)    ! directory name of 2aku file
    CHARACTER*1,INTENT(in) :: file_seaice(1000)  ! full name of 2aku file
    CHARACTER*1,INTENT(in) :: file_2cmb(1000)    ! full name of 2cmb file
    INTEGER*4,INTENT(out)  :: stat_gmi1      ! cgmi file1 status
    INTEGER*4,INTENT(out)  :: stat_gmi2      ! cgmi file2 status
    INTEGER*4,INTENT(out)  :: stat_gmi3      ! cgmi file3 status
    INTEGER*4,INTENT(out)  :: stat_2akuenv   ! 2aku environmental file status
    INTEGER*4,INTENT(out)  :: stat_2adpr     ! 2adpr file status
    INTEGER*4,INTENT(out)  :: stat_2aku      ! 2aku file status
    INTEGER*4,INTENT(out)  :: stat_snow      ! snow cover file status
    INTEGER*4,INTENT(out)  :: stat_seaice    ! reynolds seaice file status
    INTEGER*4,INTENT(out)  :: stat_cmb       ! combined file status
    INTEGER*4,INTENT(out)  :: main_date      ! principal date
    INTEGER*4,INTENT(out)  :: main_orbit     ! principal orbit
    
!...Local variables
    INTEGER*4    orbit_gmi1, orbit_gmi2, orbit_gmi3,    & ! derived orbit numbers
                 orbit_2adpr, orbit_2akuenv, orbit_2aku
    INTEGER*4    date_gmi1, date_gmi2, date_gmi3,       & ! derived date numbers
                 date_2adpr, date_2akuenv, date_2aku
    CHARACTER*10 alg_gmi1, alg_gmi2, alg_gmi3,          & ! derived alg IDs
                 alg_2adpr, alg_2akuenv, alg_2aku
    CHARACTER*3 nil_checker        ! local variable to check for "nil" files
    CHARACTER(len=1000) xfile      ! temporary file name
    INTEGER*4 i                    ! loop index
    LOGICAL*4 exist_flag           ! tests for file existence
    INTEGER*4 ix                   ! array index inside string
    CHARACTER(len=1000) yfile      ! temporary file name

!...Initialize status flags for all input files
    stat_snow    = 0
    stat_seaice  = 0

!...Check f2AKu
    CALL single_verify(i2aku,file_2aku,orbit_2aku,date_2aku,alg_2aku,          &
                       stat_2aku)
    main_date  = date_2aku
    main_orbit = orbit_2aku  
    PRINT *,'Retrieved 2AKu Orbit Number : ',main_orbit
    PRINT *,'Retrieved 2AKu Date Number  : ',main_date
    
!...Check CGMI 1 file
    CALL single_verify(igmi1,file_gmi1,orbit_gmi1,date_gmi1,alg_gmi1,          &
                       stat_gmi1)
    
!...Check CGMI 2 file
    CALL single_verify(igmi2,file_gmi2,orbit_gmi2,date_gmi2,alg_gmi2,          &
                       stat_gmi2)
    
!...Check CGMI 3 file
    CALL single_verify(igmi3,file_gmi3,orbit_gmi3,date_gmi3,alg_gmi3,          &
                       stat_gmi3)
    
!...Check f2ADPR
    CALL single_verify(i2adpr,file_2adpr,orbit_2adpr,date_2adpr,alg_2adpr,     &
                       stat_2adpr)
    
!...Check f2AKuENV
    CALL single_verify(i2akuenv,file_2akuenv,orbit_2akuenv,date_2akuenv,       &
                       alg_2akuenv,stat_2akuenv)
    
!...Check SeaIce file
    stat_seaice = 0
    nil_checker = TRANSFER(file_seaice,nil_checker)
    IF (nil_checker .NE. 'nil') THEN
        xfile = ' '
        DO i=1,iseaice
            xfile(i:i)= file_seaice(i)
        ENDDO
        INQUIRE(FILE=TRIM(xfile),EXIST=exist_flag)
	IF (.NOT. exist_flag) THEN
	    stat_seaice = -1
        ENDIF				  
    ELSE
        stat_seaice = -2
    ENDIF

!...Directory Check for CMB File
    stat_cmb = 0 
    !...extract just the directory from full name
    xfile = ' '
    DO i=1,i2cmb
        xfile(i:i)= file_2cmb(i)
    ENDDO
    ix = SCAN (xfile, '/', BACK = .TRUE.)
    yfile = ' '
    yfile(1:ix) = xfile(1:ix)
    !...verify directory exists
    IF (ix .GE. 1) THEN
        INQUIRE (DIRECTORY=yfile(1:ix), EXIST=exist_flag)
        IF (exist_flag) THEN
            stat_cmb = 0
        ELSE
            stat_cmb = -1
        ENDIF        
    ELSE
        stat_cmb = 0
    ENDIF        
    !...check for nil 
    nil_checker = TRANSFER(file_2cmb,nil_checker)
    IF (nil_checker .EQ. 'nil') stat_cmb = -2
    
!...Directory Check for Snow File
    INQUIRE (DIRECTORY=file_snow(1:isnow), EXIST=exist_flag)
    IF (exist_flag) THEN
        stat_snow = 0
    ELSE
        stat_snow = -1
    ENDIF        

!...Verify Orbit Numbers
    IF (orbit_gmi2 .NE. orbit_2aku .OR. orbit_2adpr .NE. orbit_2aku .OR.       &
        orbit_2akuenv .NE. orbit_2aku .OR. (orbit_gmi1+1) .NE. orbit_2aku .OR. &
	(orbit_gmi3-1) .NE. orbit_2aku)                                        &
        PRINT *,'WARNING: Orbit numbers not correlated ',                      &
	         orbit_2aku,orbit_gmi1,orbit_gmi2,orbit_gmi3,orbit_2adpr,      &
		 orbit_2akuenv
		 
!...Verify Date Numbers
    IF (date_gmi2 .NE. date_2aku .OR.                                          &
        date_2adpr .NE. date_2aku .OR.                                         &
        (date_gmi1 .NE. date_2aku .AND. date_gmi1 .NE. 0) .OR.                 &
        (date_gmi3 .NE. date_2aku .AND. date_gmi3 .NE. 0) .OR.                 &
        date_2akuenv .NE. date_2aku)                                           &
        PRINT *,'WARNING: date numbers not correlated ',date_2aku,             &
	         date_gmi1,date_gmi2,date_gmi3,date_2adpr,date_2akuenv
		 
!...Verify File Typing
    IF (alg_2aku .NE. '2AKu' .OR.  alg_2adpr .NE. '2ADPR' .OR.                 &
	alg_gmi1 .NE. '1CGMI' .OR. alg_gmi2 .NE. '1CGMI' .OR.                  &
	alg_gmi3 .NE. '1CGMI' .OR. alg_2akuenv .NE. '2AKuENV')                 &
        PRINT *,'WARNING: some algorithm IDs incorrect ',alg_2aku,             &
	         alg_gmi1,alg_gmi2,alg_gmi3,alg_2adpr,alg_2akuenv
	
    RETURN
    END

!===============================================================================
!
!---function: Performs file information verification and ability to process
!             for HDF5 format files
!
!===============================================================================
SUBROUTINE single_verify (ilen, file_name, orbit, date, algorithm, status)

    IMPLICIT NONE
    
!...Calling sequence parameters  
    INTEGER*4,INTENT(in)     :: ilen             ! length of file name      
    CHARACTER*1,INTENT(in)   :: file_name(1000)  ! full name of file
    INTEGER*4,INTENT(out)    :: orbit            ! internal orbit number
    INTEGER*4,INTENT(out)    :: date             ! internal date yyyymmdd
    CHARACTER*10,INTENT(out) :: algorithm        ! standard TRMM format data ID
    INTEGER*4,INTENT(out)    :: status           ! access status

!...Local variables
    CHARACTER*3 nil_checker        ! local variable to check for "nil" files
    CHARACTER(len=1000) xfile      ! temporary file name
    INTEGER*4 i                    ! loop index
    LOGICAL*4 exist_flag           ! tests for file existence

!...Set default status
    status = 0
    date = 0
    orbit = 0
    algorithm = ' '

!...Check file availability
    nil_checker = TRANSFER(file_name,nil_checker)
    IF (nil_checker .NE. 'nil') THEN
        xfile = ' '
        DO i=1,ilen
            xfile(i:i)= file_name(i)
        ENDDO
        INQUIRE(FILE=TRIM(xfile),EXIST=exist_flag)
	IF (.NOT. exist_flag) THEN
	    status = -1
        ELSE
            CALL header_getter_L2(ilen, file_name(1:ilen), orbit,              &
	                          date, algorithm, status)
            IF (status .NE. 0) THEN
	        status = -4
            ENDIF
        ENDIF				  
    ELSE
        status = -2
    ENDIF

    RETURN
    END

