!#######################################################################
subroutine remove_extra_spaces( str )
  ! remove extra embedded spaces in a character string
  use program_constants
  implicit none
  character(len=*), intent(inout) :: str
  integer(our_int) :: last, i
  do
     last = len_trim( str )
     if( last == 0 ) exit
     i = index( trim(str), "  " )
     if( i == 0 ) exit
     str( (i+1):last ) = str( (i+2):last ) // " "
  end do
end subroutine remove_extra_spaces
!#######################################################################
subroutine remove_spaces_around_equals( str )
  ! remove spaces around any equals sign in a character string
  use program_constants
  implicit none
  character(len=*), intent(inout) :: str
  integer(our_int) :: last, i
  do
     last = len_trim( str )
     if( last == 0 ) exit
     i = index( trim(str), " =" )
     if( i == 0 ) exit
     str( i:last ) = str( (i+1):last ) // " "
  end do
  do
     last = len_trim( str )
     if( last == 0 ) exit
     i = index( trim(str), "= " )
     if( i == 0 ) exit
     str( (i+1):last ) = str( (i+2):last ) // " "
  end do
end subroutine remove_spaces_around_equals
!#######################################################################
subroutine remove_trailing_comments( str )
  ! blanks everything in character string after exclamation point
  use program_constants
  implicit none
  character(len=*), intent(inout) :: str
  integer(our_int) :: last, i
  do
     last = len_trim( str )
     if( last == 0 ) exit
     i = index( trim(str), "!" )
     if( i == 0 ) exit
     str( i:last ) = ""
  end do
end subroutine remove_trailing_comments
!#######################################################################
subroutine count_program_statements( filename, iounit, line, no_lines )
  ! count every line that is not blank, that does not begin with 
  ! an exclamation point, and that does not end with an ampersand
  use program_constants
  implicit none
  character(len=*), intent(in) :: filename
  integer(our_int), intent(in) :: iounit
  character(len=*), intent(inout) :: line
  integer(our_int), intent(out) :: no_lines
  integer(our_int) :: linewidth, status
  ! begin
  linewidth = len(line)
  open( unit=iounit, file=filename, action="READ", status="OLD" )
  no_lines = 0
  do
     read( unit=iounit, fmt="(A)", iostat=status ) line
     if( status < 0 ) exit ! end of file reached
     line = adjustl( line )
     if( line(1:1) == "!" ) cycle
     line = adjustr( line )
     if( line(linewidth:linewidth) == "&" ) cycle
     no_lines = no_lines + 1
  end do
  close(iounit)
end subroutine count_program_statements
!#######################################################################
subroutine count_all_lines( filename, iounit, line, no_lines )
  ! count every line in the file
  use program_constants
  implicit none
  character(len=*), intent(in) :: filename
  integer(our_int), intent(in) :: iounit
  character(len=*), intent(inout) :: line
  integer(our_int), intent(out) :: no_lines
  integer(our_int) :: linewidth, status
  ! begin
  linewidth = len(line)
  open( unit=iounit, file=filename, action="READ", status="OLD" )
  no_lines = 0
  do
     read( unit=iounit, fmt="(A)", iostat=status ) line
     if( status < 0 ) exit ! end of file reached
     no_lines = no_lines + 1
  end do
  close(iounit)
end subroutine count_all_lines
!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################
program errcodes
  use program_constants
  use error_handler
  use dynalloc
  use tabulate
  implicit none
  integer(our_int) :: fileno, status, nlines, i, j, no_subnames, &
       no_modnames, no_comments, posn, st, fin, nlines_all_files
  integer(our_int), allocatable :: file_st(:), file_fin(:), &
       file_nlines(:)
  character(len=80) :: mod_filename
  integer(our_int), parameter :: linewidth=132, err_msg_width=72
  character(len=linewidth) :: line
  character(len=40*linewidth) :: long_line
  character(len=40*linewidth), allocatable :: program_statements(:), &
       program_statements_all_files(:)
  character(len=err_msg_width), allocatable :: subnames(:), modnames(:), &
       comments(:)
  character(len=err_msg_width), pointer :: &
       nonredundant_subnames(:) => null(), &
       nonredundant_modnames(:) => null(), &
       nonredundant_comments(:) => null()
  character(len=12) :: sInt
  type(table_type) :: table
  type(error_type) :: err
  integer(our_int) :: ijunk
  !##################
  call err_reset(err)
  ! count total no. of program statements in all *.f90 files
  allocate( file_nlines( size( modname_list ) ) )
  nlines_all_files = 0
  do fileno = 1, size( modname_list )
     mod_filename = adjustl( trim( modname_list( fileno ) ) // ".f90" )
     call count_program_statements( mod_filename, 99, line, nlines )
     nlines_all_files = nlines_all_files + nlines
     file_nlines( fileno ) = nlines
  end do
  allocate( file_st( size( modname_list ) ) )
  allocate( file_fin( size( modname_list ) ) )
  file_st(1) = 1
  file_fin(1) = file_nlines(1)
  do fileno = 2, size( modname_list )
     file_st(fileno) = file_fin( fileno-1 ) + 1
     file_fin(fileno) = file_st(fileno) + file_nlines(fileno) - 1
  end do
  ! read in all statements from all files
  allocate( program_statements_all_files( nlines_all_files ) )
  do fileno = 1, size( modname_list )
     mod_filename = adjustl( trim( modname_list( fileno ) ) // ".f90" )
     nlines = file_nlines( fileno )
     allocate( program_statements( nlines ) )
     !### read in program statements, concatenating continuation lines
     program_statements(:) = ""
     open( unit=99, file=mod_filename, action="READ", status="OLD" )
     i = 0
     do
        ! go to next line that doesn't start with exclamation point
        do
           read( unit=99, fmt="(A)", iostat=status ) line
           if( status < 0 ) exit ! end of file reached
           line = adjustl( line )
           if( line(1:1) /= "!" ) exit
        end do
        if( status < 0 ) exit ! end of file reached
        ! start a new long_line
        i = i + 1
        long_line = ""
        long_line = line
        ! add continuation lines, if present
        do
           line = adjustr( line )
           if( line(linewidth:linewidth) /= "&" ) exit
           ! remove ampersand
           long_line = adjustr( long_line )
           long_line( len(long_line): ) = " "
           long_line = adjustl( long_line )
           ! skip to next noncomment line
           do
              read( unit=99, fmt="(A)", iostat=status ) line
              if( status < 0 ) exit ! end of file reached
              line = adjustl( line )
              if( line(1:1) /= "!" ) exit
           end do
           if( status < 0 ) exit  !EOF reached
           long_line = trim( long_line ) // " " // trim( line )
        end do
        program_statements(i) = long_line
        if( status < 0 ) exit
     end do
     close(99)
     program_statements_all_files( &
          file_st(fileno) : file_fin(fileno) ) = program_statements(:)
     deallocate( program_statements )
  end do
  !#################################################################
  do i = 1, nlines_all_files
     call remove_extra_spaces( program_statements_all_files(i) )
     call remove_spaces_around_equals( program_statements_all_files(i) )
     call remove_trailing_comments( program_statements_all_files(i) )
  end do
  !#################################################################
  ! find every occurrence of "modname=something"
  no_modnames = 0
  do i = 1, nlines_all_files
     if( index( program_statements_all_files(i), "modname=" ) /= 0 ) &
          no_modnames = no_modnames + 1
  end do
  allocate( modnames( no_modnames ) )
  j = 0
  do i = 1, nlines_all_files
     posn = index( program_statements_all_files(i), "modname=" )
     if( posn == 0 ) cycle
     j = j + 1
     st = posn + 9
     modnames(j) = program_statements_all_files(i)(st:)
     posn = index( modnames(j), '"' )
     modnames(j)(posn:) = ""
  end do
  ijunk = nullify_table( table, err )
  ijunk = tabulate_variable( modnames, table, err )
  ijunk = get_table_values( nonredundant_modnames, table, err )
  !#############################################################
  ! find every occurrence of "subname=something"
  no_subnames = 0
  do i = 1, nlines_all_files
     if( index( program_statements_all_files(i), "subname=" ) /= 0 ) &
          no_subnames = no_subnames + 1
  end do
  allocate( subnames( no_subnames ) )
  j = 0
  do i = 1, nlines_all_files
     posn = index( program_statements_all_files(i), "subname=" )
     if( posn == 0 ) cycle
     j = j + 1
     st = posn + 9
     subnames(j) = program_statements_all_files(i)(st:)
     posn = index( subnames(j), '"' )
     subnames(j)(posn:) = ""
  end do
  ijunk = nullify_table( table, err )
  ijunk = tabulate_variable( subnames, table, err )
  ijunk = get_table_values( nonredundant_subnames, table, err )
  !#############################################################
  ! find every occurrence of "comment=something"
  no_comments = 0
  do i = 1, nlines_all_files
     if( index( program_statements_all_files(i), "comment=" ) /= 0 ) &
          no_comments = no_comments + 1
  end do
  allocate( comments( no_comments ) )
  j = 0
  do i = 1, nlines_all_files
     posn = index( program_statements_all_files(i), "comment=" )
     if( posn == 0 ) cycle
     j = j + 1
     st = posn + 9
     comments(j) = program_statements_all_files(i)(st:)
     posn = index( comments(j), '"' )
     comments(j)(posn:) = ""
  end do
  ijunk = nullify_table( table, err )
  ijunk = tabulate_variable( comments, table, err )
  ijunk = get_table_values( nonredundant_comments, table, err )
  !#############################################################
  open( unit=99, file="icodes_1.h", status="REPLACE", action="WRITE")
  do i = 1, size( nonredundant_modnames )
     write(sInt,"(I12)") i
     sInt = adjustl(sInt)
     write(99,"(A)") '"' // nonredundant_modnames(i) // '",&'
  end do
  close( 99 )
  !#############################################################
  open( unit=99, file="icodes_2.h", status="REPLACE", action="WRITE")
  do i = 1, size( nonredundant_subnames )
     write(sInt,"(I12)") i
     sInt = adjustl(sInt)
     write(99,"(A)") '"' // nonredundant_subnames(i) // '",&'
  end do
  close( 99 )
  !#############################################################
  open( unit=99, file="icodes_3.h", status="REPLACE", action="WRITE")
  do i = 1, size( nonredundant_comments )
     write(sInt,"(I12)") i
     sInt = adjustl(sInt)
     write(99,"(A)") '"' // nonredundant_comments(i) // '",&'
  end do
  close( 99 )
  !#############################################################
  open( unit=99, file="icodes.R", status="REPLACE", action="WRITE")
  write(99,"(A)") "#######################################################"
  write(99,"(A)") "# This R code is automagically generated by errcodes"
  write(99,"(A)") "# When executed within R, it creates the file"
  write(99,"(A)") "# sysdata.rda which should be placed in the"
  write(99,"(A)") "# subdirectory R of package prior to package building."
  write(99,"(A)") "#######################################################"
  write(99,"(A)") "icodesMsg <- list("
  write(sInt,"(I12)") size( nonredundant_modnames )
  sInt = adjustl(sInt)
  write(99,"(A)") "   modnames = character(" // trim(sInt) // "),"
  write(sInt,"(I12)") size( nonredundant_subnames )
  sInt = adjustl(sInt)
  write(99,"(A)") "   subnames = character(" // trim(sInt) // "),"
  write(sInt,"(I12)") size( nonredundant_comments )
  sInt = adjustl(sInt)
  write(99,"(A)") "   comments = character(" // trim(sInt) // ") )"
  do i = 1,  size( nonredundant_modnames )
     write(sInt,"(I12)") i
     sInt = adjustl( sInt )
     write(99,"(A)") "icodesMsg$modnames[" // trim(sInt) &
          // '] <- "' // trim( nonredundant_modnames(i) ) // '"'
  end do
  do i = 1,  size( nonredundant_subnames )
     write(sInt,"(I12)") i
     sInt = adjustl( sInt )
     write(99,"(A)") "icodesMsg$subnames[" // trim(sInt) &
          // '] <- "' // trim( nonredundant_subnames(i) ) // '"'
  end do
  do i = 1,  size( nonredundant_comments )
     write(sInt,"(I12)") i
     sInt = adjustl( sInt )
     write(99,"(A)") "icodesMsg$comments[" // trim(sInt) &
          // '] <- "' // trim( nonredundant_comments(i) ) // '"'
  end do
  write(99,"(A)") 'save("icodesMsg", file="sysdata.rda", compress=TRUE)'
  close(99)
  !#############################################################
  ! revise the source code in error_handler.f90
  mod_filename = "error_handler.f90"
  call count_all_lines( mod_filename, 99, line, nlines )
  allocate( program_statements( nlines ) )
  program_statements(:) = ""
  !### read in all lines
  open( unit=99, file=mod_filename, action="READ", status="OLD" )
  i = 0
  do
     i = i + 1
     read( unit=99, fmt="(A)", iostat=status ) line
     if( status < 0 ) exit ! end of file reached
     program_statements(i) = line
  end do
  close(99)
  do i = 1, nlines
     line = program_statements(i)
     st = index(line, "nonredundant_modnames(" )
     if( st > 0 ) then
        if( index( adjustl(line), "nonredundant_modnames(" ) == 1 ) then
           write(sInt,"(I12)") size( nonredundant_modnames ) + 1
           sInt = adjustl(sInt)
           line(st:) = "nonredundant_modnames( " // trim(sInt) // " ) = (/ &"
           program_statements(i) = line
        end if
     end if
     st = index(line, "nonredundant_subnames(" )
     if( st > 0 ) then
        if( index( adjustl(line), "nonredundant_subnames(" ) == 1 ) then
           write(sInt,"(I12)") size( nonredundant_subnames ) + 1
           sInt = adjustl(sInt)
           line(st:) = "nonredundant_subnames( " // trim(sInt) // " ) = (/ &"
           program_statements(i) = line
        end if
     end if
     st = index(line, "nonredundant_comments(" )
     if( st > 0 ) then
        if( index( adjustl(line), "nonredundant_comments(" ) == 1 ) then
           write(sInt,"(I12)") size( nonredundant_comments ) + 1
           sInt = adjustl(sInt)
           line(st:) = "nonredundant_comments( " // trim(sInt) // " ) = (/ &"
           program_statements(i) = line
        end if
     end if
  end do
  open( unit=99, file=mod_filename, action="WRITE", status="REPLACE" )
  do i = 1, nlines
     write(99,"(A)") trim(program_statements(i))
  end do
  close(99)
  deallocate( program_statements )
  !#############################################################
  deallocate( program_statements_all_files )
  deallocate( file_st )
  deallocate( file_fin )
  deallocate( file_nlines )
  deallocate( modnames )
  deallocate( subnames )
  deallocate( comments )
  ijunk = dyn_dealloc( nonredundant_modnames, err )
  ijunk = dyn_dealloc( nonredundant_subnames, err )
  ijunk = dyn_dealloc( nonredundant_comments, err )
  ijunk = nullify_table( table, err )
  call err_reset(err)
end program errcodes
!#######################################################################
