!=============================================================================
!
!  NAME : header_getter    FORTRAN
!
!  PURPOSE: This is a quickie utility to retrieve the orbit number and date
!           number from an HDF5 formatted PPS file.  This information can
!           also be retrieved using Toolkit utilities; however, this is not
!           always convenient.  Note that if you don't give it an HDF5 file, 
!           it will give "H5F_OPEN" message and stop.
!
!  INPUT PARAMS:
!        input_file - HDF formatted input file
!
!  OUTPUT VALUE:
!        orbit_number  - orbit number in integer format
!        date_number   - date number in integer format, yyyymmdd
!        algorithm_ID  - type file (e.g. 1CGMI, 2AKu, etc)
!        error_code    - error output, integer; anything less than 0 ==> bad
!                        and you get a message
!
!  AUTHOR:  S.McLaughlin   03/18/2013
!
!  HISTORY: 
!    Date       Author        Change activity
!  03/18/2013 S.McLaughlin -  initial version
!  04/03/2013 S.McLaughlin -  added algorithm_ID
!  12/06/2013 S.McLaughlin -  adapted for use in L2 code (file name passing)
!
!=============================================================================

SUBROUTINE header_getter_L2(len_file, filename, orbit_number, date_number,     &
                            algorithm_ID, error_code)

USE HDF5  ! module contains necessary modules

IMPLICIT NONE

!...Calling sequence variables
    INTEGER*4, intent(in)   :: len_file
    CHARACTER(len=len_file) :: filename
    INTEGER*4, intent(out)  :: orbit_number
    INTEGER*4, intent(out)  :: date_number
    CHARACTER(len=10), intent(out) :: algorithm_ID
    INTEGER*4, intent(out)  :: error_code

!...Programm variables needed for HDF5 access
    INTEGER*4       :: error       ! Error flag from HDF5 utilities
    INTEGER(HID_T)  :: file_id     ! File identifier 
    INTEGER(HID_T)  :: attr_id     ! Attribute identifier
    INTEGER(SIZE_T), PARAMETER :: buf_size = 255   ! Buffer size, attr name
    CHARACTER(LEN=250) :: buffer   ! Buffer to hold attribute name
    INTEGER(HID_T)  :: type_id     ! Attribute datatype identifier
    INTEGER*4 :: num_attr          ! Number of attributes in file      
    INTEGER(HSIZE_T), DIMENSION(*) :: dims(1) ! dimension size(s) of data
                                              !  buffers, in this case just 1
    LOGICAL :: f_corder_valid      ! Indicates whether the creation order 
                                   !  data is valid for this attribute 
    INTEGER :: corder              ! Is a positive integer containing the 
                                   !  creation order of the attribute
    INTEGER :: cset                ! Indicates the character set used for 
                                   !  the attribute?s name
    INTEGER(HSIZE_T) :: data_size  ! Indicates the size, in the number
                                   ! of characters, of the attribute
    CHARACTER(LEN=1000) :: big_buffer   ! Buffer, holds total header string
    CHARACTER(LEN=1000) :: new_string   ! Buffer, revised headerstring

!...Local computational variables
    INTEGER*4 :: place1, place3, place2, place_length
    INTEGER*4 :: year, month, day

!...Set initial error code
    error_code = 0
    
!...Initialize FORTRAN interface.
    CALL h5open_f(error) 

!...Open file
    CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error_code)
    IF (error_code .NE. 0) THEN
        PRINT *, 'WARNING: Error on open of HDF5 file ',filename
        RETURN  
    ENDIF

!...Get number of attributes
    CALL h5aget_num_attrs_f(file_id, num_attr, error)
    IF (num_attr .LE. 0) THEN 
        PRINT *, 'No attributes this file ',filename
        error = -1
	RETURN
    ENDIF

!...Open the attribute where header is expected (at "0")
    CALL h5aopen_idx_f(file_id, 0, attr_id, error)

!...Get attribute's name, verify it is the expected header
    CALL h5aget_name_f(attr_id, buf_size, buffer, error)
    IF (buffer(1:error) .NE. 'FileHeader') THEN 
        PRINT *, 'Attribute is not FileHeader ',buffer
        error_code = -2
	RETURN
    ENDIF

!...Retrieve attribute type
    CALL h5aget_type_f(attr_id, type_id, error)

!...Retrieve size of the header character string
    CALL h5aget_info_f(attr_id, f_corder_valid, corder, cset, data_size, error)
    dims(1) = data_size

!...Acquire header string
    CALL h5aread_f(attr_id, type_id, big_buffer, dims, error)

!...Locate the orbit number in the string
    place1 = INDEX(big_buffer,'GranuleNumber=')
    IF (place1 .LE. 0) THEN
        PRINT *, 'No orbit number ',place1
        error_code = -3
	RETURN
    ENDIF

!...Check length of orbit number string; need start/end selected string
    new_string = big_buffer(place1:500)
    place2 = INDEX(new_string,'=')
    place3 = INDEX(new_string,';')
    place_length = place3 - place2
    
!...Isolate orbit number string
    IF (place_length .GT. 1) THEN
        READ(UNIT=new_string(place2+1:place3-1),fmt=*) orbit_number
    ELSE  
        orbit_number = 0
    ENDIF

!...Locate the date number in the string
    place1 = INDEX(big_buffer,'StartGranuleDateTime=')
    IF (place1 .LE. 0) THEN
        PRINT *, 'No date number ',place1
        error_code = -4
	RETURN
    ENDIF

!...Check length of date number string; need start/end selected string
    new_string = big_buffer(place1:500)
    place2 = INDEX(new_string,'=')
    place3 = INDEX(new_string,';')
    place_length = place3 - place2
    
!...Isolate date number strings
    IF (place_length .GT. 1) THEN
        READ(UNIT=new_string(place2+1:place3-1),FMT=53) year, month, day
        date_number = 10000*year + 100*month + day
    ELSE  
        date_number = 0
    ENDIF
53  FORMAT(i4,1x,i2,1x,i2)

!...Locate the algorithm ID in the string
    place1 = INDEX(big_buffer,'AlgorithmID=')
    IF (place1 .LE. 0) THEN
        PRINT *, 'No algorithm ID ',place1
        error_code = -4
	RETURN
    ENDIF

!...Check length of algorithm ID string; need start/end selected string
    new_string = big_buffer(place1:500)
    place2 = INDEX(new_string,'=')
    place3 = INDEX(new_string,';')
    place_length = place3 - place2

!...Isolate algorithm ID string
    IF (place_length .GT. 1) THEN
        READ(UNIT=new_string(place2+1:place3-1),fmt=*) algorithm_ID
    ELSE  
        algorithm_ID = ' '
    ENDIF
    algorithm_ID = trim(algorithm_ID)

!...Close attribute
    CALL h5aclose_f(attr_id, error)

!...Close file
    CALL h5fclose_f(file_id, error)
     
!...Close FORTRAN interface.
    CALL h5close_f(error)
 
RETURN
END
