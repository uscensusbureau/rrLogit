!#####################################################################
! modified on 12/7/2022 to include a dim_int integer kind, and to
! remove old constants that aren't needed within an R package
!#####################################################################
module program_constants
   ! Programming constants used throughout the program.
   ! Unlike most modules, everything here is public.
   implicit none
   public
   ! Define compiler-specific KIND numbers for integers,
   ! single and double-precision reals to help ensure consistency of
   ! performance across platforms:
   integer, parameter :: our_int = selected_int_kind(9), &
        our_sgle = selected_real_kind(6,37), &
        our_dble = selected_real_kind(15,307)
   ! edit next line to define integer type for dimensioning/indexing
   ! arrays that could be large; "9" for 32-bit, "18" for 64-bit
   integer, parameter :: long_int = selected_int_kind(18)
   ! Common integer values returned by all functions to indicate
   ! success or failure:
   integer(kind=our_int), parameter :: RETURN_SUCCESS = 0, &
        RETURN_FAIL = -1
end module program_constants
!#####################################################################
