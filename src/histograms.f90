
module histograms

  use kinds

  implicit none

  type histogram_type
     real(kind=dp) :: hmin,hmax,dh
     integer :: hsize
     real(kind=dp), dimension(:), pointer :: array=>NULL()
     integer :: nupdates
  end type histogram_type

contains

  subroutine create_histogram(histogram,hmin,hmax,dh)
    type(histogram_type), intent(out) :: histogram
    real(kind=dp), intent(in), optional :: hmin
    real(kind=dp), intent(in), optional :: hmax
    real(kind=dp), intent(in), optional :: dh

    real(kind=dp) :: my_hmin=0.0_dp,my_hmax=10.0_dp,my_dh=0.1_dp

    integer :: my_hsize

    if (present(hmin)) my_hmin=hmin
    if (present(hmax)) my_hmax=hmax
    if (present(dh))   my_dh=dh

    if (.not.(my_dh>0)) stop "dh > 0"
    if (.not.(my_hmax>=(my_hmin+my_dh))) stop "hmax >= hmin+dh"

    histogram%hmin=my_hmin
    histogram%hmax=my_hmax
    histogram%dh  =my_dh

    histogram%nupdates=0

    my_hsize=floor((my_hmax-my_hmin)/my_dh)

    histogram%hsize=my_hsize
    allocate(histogram%array(my_hsize))
    histogram%array=0.0_dp

  end subroutine create_histogram

  subroutine update_histogram(histogram,value,amount)
    type(histogram_type), intent(inout) :: histogram
    real(kind=dp), intent(in) :: value
    real(kind=dp), intent(in), optional :: amount

    integer :: index
    real(kind=dp) :: my_amount

    my_amount=1.0_dp
    if (present(amount)) my_amount=amount

    index=floor((value-histogram%hmin)/histogram%dh)+1
    if ((index>=1).and.(index<=histogram%hsize)) then
       histogram%array(index)=histogram%array(index)+my_amount
    end if

    histogram%nupdates=histogram%nupdates+1

  end subroutine update_histogram

  subroutine finalize_histogram(histogram)
    type(histogram_type), intent(inout) :: histogram

    histogram%array=histogram%array/dble(histogram%nupdates)/histogram%dh

  end subroutine finalize_histogram

  subroutine write_histogram(out_unit,histogram,comment)
    integer, intent(in) :: out_unit
    type(histogram_type), intent(in) :: histogram
    character(len=*),intent(in),optional :: comment

    integer :: index
    real(kind=dp) :: h

    if (present(comment)) write(out_unit,'("# ",A)') trim(comment)
    write(out_unit,'("# dh=",F14.6,", hmin=",F14.6,", hmax=",F14.6)') &
         histogram%dh, histogram%hmin,histogram%hmax
    
    do index=1,histogram%hsize
       h=histogram%hmin+(dble(index-1)+0.5_dp)*histogram%dh
       write(out_unit,'(2F14.6)') h,histogram%array(index)
    end do

  end subroutine write_histogram

end module histograms
