
module gofrs

  use kinds
  use selections
  use histograms
  use constants
  use pbc

  implicit none

  type gofr_type
     real(kind=dp) :: limit,dr
     real(kind=dp), dimension(:), pointer :: rdf,integral
     integer :: nrdf
     real(kind=dp) :: nsel1,nsel2
     integer :: nupdates=-1
  end type gofr_type

contains

  subroutine create_gofr(gofr,limit,dr)
    type(gofr_type), intent(out) :: gofr
    real(kind=dp), intent(in), optional :: limit
    real(kind=dp), intent(in), optional :: dr

    real(kind=dp) :: my_limit=10.0_dp,my_dr=0.1_dp

    if (present(limit)) my_limit=limit
    if (present(dr)) my_dr=dr

    gofr%limit=my_limit
    gofr%dr=my_dr

    gofr%nrdf=floor(limit/dr)

    allocate(gofr%rdf(gofr%nrdf))
    allocate(gofr%integral(gofr%nrdf))
    gofr%rdf=0.0_dp
    gofr%integral=0.0_dp

    gofr%nupdates=0
    gofr%nsel1=0.0_dp
    gofr%nsel2=0.0_dp

  end subroutine create_gofr
    
  subroutine update_gofr(gofr,atoms,configuration,selection1,selection2,&
                         box,usepbc)
    type(gofr_type), intent(inout) :: gofr
    type(atom), dimension(:), intent(in) :: atoms
    real(kind=dp), dimension(:,:), intent(in) :: configuration
    logical, dimension(:), intent(in) :: selection1,selection2
    type(box_type), intent(in) :: box
    logical, intent(in), optional :: usepbc

    integer :: iat1, iat2
    real(kind=dp), dimension(3) :: vector
    real(kind=dp) :: dist,dr
    integer :: index

    integer :: natoms
    logical :: my_usepbc=.true.

    dr=gofr%dr

    natoms=size(configuration(1,:))
    if (present(usepbc)) my_usepbc=usepbc

    do iat1=1,natoms
       if (selection1(iat1)) then
          do iat2=1,natoms
             if (selection2(iat2).and.(iat1.ne.iat2)) then
                vector=configuration(:,iat2)-configuration(:,iat1)
                if (my_usepbc) vector=minimum_image(vector,box)
                dist=distance(vector)
                if (dist.lt.0.5) then
                        print*, iat2, configuration(:,iat2)
                        print*, iat1, configuration(:,iat1)
                        print*, configuration(:,iat2)-configuration(:,iat1)
                        print*, vector
                        print*, dist
                        print*, "Error in gofr: distance too shoort"
                        stop
                endif
                index=floor(dist/dr)+1
                if (index<=gofr%nrdf) then
                   gofr%rdf(index)=gofr%rdf(index)+1.0_dp
                end if
             end if
          end do
       end if
    end do

    gofr%nsel1=gofr%nsel1+dble(count_selection(selection1))
    gofr%nsel2=gofr%nsel2+dble(count_selection(selection2))
    
    gofr%nupdates=gofr%nupdates+1

  end subroutine update_gofr

  subroutine finalize_gofr(gofr,box)
    type(gofr_type), intent(inout) :: gofr
    type(box_type), intent(in) :: box

    integer :: index
    real(kind=dp) :: volume,r,fact

    gofr%rdf=gofr%rdf/dble(gofr%nupdates)
    gofr%nsel1=gofr%nsel1/dble(gofr%nupdates)
    gofr%nsel2=gofr%nsel2/dble(gofr%nupdates)

    volume=volume_of_box(box)

    do index=1,gofr%nrdf
       gofr%integral(index)=sum(gofr%rdf(1:index))
    end do
    gofr%integral=gofr%integral/dble(gofr%nsel1)

    do index=1,gofr%nrdf
       r=(dble(index)-0.5_dp)*gofr%dr
       fact=volume/gofr%nsel1/gofr%nsel2/(4.0_dp*pi*r*r*gofr%dr)
       gofr%rdf(index)=gofr%rdf(index)*fact
    end do

  end subroutine finalize_gofr

  subroutine write_gofr(out_unit,gofr,comment)
    integer, intent(in) :: out_unit
    type(gofr_type), intent(in) :: gofr
    character(len=*), intent(in), optional :: comment

    integer :: index
    real(kind=dp) :: r

    if (present(comment)) write(out_unit,'("# ",A)') trim(comment)
    write(out_unit,&
      '("# dr=",F14.6,", limit=",F14.6,", nsel1=",F14.6,", nsel2=",F14.6)') &
      gofr%dr,gofr%limit,gofr%nsel1,gofr%nsel2
    
    do index=1,gofr%nrdf
       r=(dble(index)-0.5_dp)*gofr%dr
       write(out_unit,'(3F14.6)') r,gofr%rdf(index),gofr%integral(index)
    end do

  end subroutine write_gofr

end module gofrs
