!
! [Functionality]
! Read wavefunctions from "qe_wfc_char" (file name) in QE (Quantum Espresso)
!
! <namespace io_qe_files is declared in include/io_qe_files.hpp> 
!
! [Other info]
! See Modules/io_base.f90 in QE to know details about QE variables.
!
subroutine read_wfc(nlength, qe_wfc_char, ik, igwx, nbnd, npol, gamma_only, mill, evc) bind(c, name="read_wfc")
  use, intrinsic :: iso_c_binding
  implicit none
  include 'mpif.h'

!! specify the file name to be read (see above)
  integer(c_int), intent(in) :: nlength
  character(c_char), intent(in) :: qe_wfc_char(nlength)

!! QE variables (same names as QE), should be consistent with the xml file of QE
  integer(c_int), intent(in) :: ik     ! k-point index (wfc{ik}.dat will be read)
  integer(c_int), intent(in) :: igwx   ! no. of plane waves for each k-point
  integer(c_int), intent(in) :: nbnd   ! no. of bands
  integer(c_int), intent(in) :: npol   ! (non-collinear) 2 meaning spinor, (otherwise) 1
  logical(c_bool), intent(in) :: gamma_only  ! .true. = memory-saved treatment of wfc in gamma-only calc

!! QE variables (same name as QE), should be extracted from wfc*.dat
  integer(c_int), intent(inout) :: mill(3,igwx)            ! miller index for G vector in evc(:,:)
  complex(c_double), intent(inout) :: evc(npol*igwx,nbnd)  ! wave function

!! QE variables (same name as QE), not used outside this subroutine
  real(8) :: xk(3)       ! k-point coordinate
  integer :: ispin       ! [no-spin or non-collinear] 1 [spin-polarized] 1 = up, 2 = down
  real(8) :: scalef      ! scale factor (usually 1.0)
  integer :: ngw         ! no. of plane waves? not used.
  real(8) :: b1(3), b2(3), b3(3)   ! reciprocal lattice vectors

!! temporary variables in this subroutine
  character(nlength) :: filename_converted
  integer :: iunit, ibnd
  integer(c_int) :: ic1, ic2, ic3
  logical :: ibool
  integer :: ierr

  filename_converted = transfer(qe_wfc_char, filename_converted) ! i.e. filename_converted(i:i) = qe_wfc_char(i)
!  write(6,*) "Read QE wave functions from ", trim(filename_converted)
  iunit = 100
  open(iunit, file = filename_converted, form = 'unformatted', status = 'old')
  
  read(iunit) ic1, xk, ispin, ibool, scalef
  call error_chk_qe_wavefunc_int(ik, ic1, "ik")
  call error_chk_qe_wavefunc_bool(gamma_only, ibool, "gamma_only")
  
  read(iunit) ngw, ic1, ic2, ic3
  call error_chk_qe_wavefunc_int(igwx, ic1, "igwx")
  call error_chk_qe_wavefunc_int(npol, ic2, "npol")
  call error_chk_qe_wavefunc_int(nbnd, ic3, "nbnd")
  
  read(iunit) b1, b2, b3  
  read(iunit) mill(1:3,1:igwx)
  do ibnd = 1, nbnd ! ibnd = band index
     read(iunit) evc(1:npol*igwx,ibnd)
  end do ! ibnd
  evc(:,:) = scalef*evc(:,:) ! multiplied with the scale factor

  close(iunit)       
end subroutine read_wfc

! The same as above but reads non-binary files (for test calculation)
subroutine read_wfc_nonbin(nlength, qe_wfc_char, ik, igwx, nbnd, npol, gamma_only, mill, evc) bind(c, name="read_wfc_nonbin")
  use, intrinsic :: iso_c_binding
  implicit none
  include 'mpif.h'

!! specify the file name to be read (see above)
  integer(c_int), intent(in) :: nlength
  character(c_char), intent(in) :: qe_wfc_char(nlength)

!! QE variables (same names as QE), should be consistent with the xml file of QE
  integer(c_int), intent(in) :: ik     ! k-point index (wfc{ik}.dat will be read)
  integer(c_int), intent(in) :: igwx   ! no. of plane waves for each k-point
  integer(c_int), intent(in) :: nbnd   ! no. of bands
  integer(c_int), intent(in) :: npol   ! (non-collinear) 2 meaning spinor, (otherwise) 1
  logical(c_bool), intent(in) :: gamma_only  ! .true. = memory-saved treatment of wfc in gamma-only calc

!! QE variables (same name as QE), should be extracted from wfc*.dat
  integer(c_int), intent(inout) :: mill(3,igwx)            ! miller index for G vector in evc(:,:)
  complex(c_double), intent(inout) :: evc(npol*igwx,nbnd)  ! wave function

!! QE variables (same name as QE), not used outside this subroutine
  real(8) :: xk(3)       ! k-point coordinate
  integer :: ispin       ! [no-spin or non-collinear] 1 [spin-polarized] 1 = up, 2 = down
  real(8) :: scalef      ! scale factor (usually 1.0)
  integer :: ngw         ! no. of plane waves? not used.
  real(8) :: b1(3), b2(3), b3(3)   ! reciprocal lattice vectors

!! temporary variables in this subroutine
  character(nlength) :: filename_converted
  character(nlength+7) :: filename_converted_nonbin
  character(7) :: nonbin
  integer :: iunit, ibnd
  integer(c_int) :: ic1, ic2, ic3
  logical :: ibool
  integer :: ierr

!! dummy index for do loop
  integer :: i1, i2

  filename_converted = transfer(qe_wfc_char, filename_converted) ! i.e. filename_converted(i:i) = qe_wfc_char(i)
!  write(6,*) "Read QE wave functions from ", trim(filename_converted)
  iunit = 100

  !! non-binary file
  nonbin = "_nonbin"
  filename_converted_nonbin = filename_converted//nonbin
  open(iunit, file = filename_converted_nonbin, status = 'old')

  read(iunit,'(I10,3f20.13,I2,L2,f20.13)') ic1, xk(:), ispin, ibool, scalef  
  call error_chk_qe_wavefunc_int(ik, ic1, "ik")
  call error_chk_qe_wavefunc_bool(gamma_only, ibool, "gamma_only")

  read(iunit, '(4I10)') ngw, ic1, ic2, ic3  
  call error_chk_qe_wavefunc_int(igwx, ic1, "igwx")
  call error_chk_qe_wavefunc_int(npol, ic2, "npol")
  call error_chk_qe_wavefunc_int(nbnd, ic3, "nbnd")
  
  read(iunit, '(9f20.13)') b1(:), b2(:), b3(:)
  do i1 = 1, igwx
     read(iunit, '(3I10)') mill(1:3,i1)
  end do
  do ibnd = 1, nbnd ! ibnd = band index
     do i1 = 1, npol
        do i2 = 1, igwx
           read(iunit, '(2f20.13)') evc((i1-1)*igwx + i2, ibnd) ! complex
        end do
     end do
  end do ! ibnd
  evc(:,:) = scalef*evc(:,:) ! multiplied with the scale factor

  close(iunit)       
end subroutine read_wfc_nonbin

subroutine error_chk_qe_wavefunc_int(ixml, iwfc, varname)
  use, intrinsic :: iso_c_binding
  implicit none
  include 'mpif.h'

  integer(c_int), intent(in) :: ixml   ! input from xml file (QE)
  integer(c_int), intent(in) :: iwfc   ! input from wfc file (QE)
  character(*), intent(in) :: varname  ! name of the variable

  integer :: ierr

  if(ixml.ne.iwfc) then
     write(6,*) " Input from xml: ", ixml, " Input from wfc.dat: ", iwfc
     write(6,*) " Error: ", varname, " is not consistent between them. (read_qe_wfc)"
     call MPI_ABORT(MPI_COMM_WORLD,999,ierr)
  end if
end subroutine error_chk_qe_wavefunc_int

subroutine error_chk_qe_wavefunc_bool(ixml, iwfc, varname)
  use, intrinsic :: iso_c_binding
  implicit none
  include 'mpif.h'

  logical(c_bool), intent(in) :: ixml   ! input from xml file (QE)
  logical, intent(in) :: iwfc   ! input from wfc file (QE)
  character(*), intent(in) :: varname  ! name of the variable

  integer :: ierr
  if((ixml.and.(.not.iwfc)).or.((.not.ixml).and.iwfc)) then
     write(6,*) " Input from xml: ", ixml, " Input from wfc.dat: ", iwfc
     write(6,*) " Error: ", varname, " is not consistent between them. (read_qe_wfc)"
     call MPI_ABORT(MPI_COMM_WORLD,999,ierr)
  end if
end subroutine error_chk_qe_wavefunc_bool
