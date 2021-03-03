
module conductivity

  use kinds
  use periodic_table
  use selections
  use colvar_utils

  implicit none 

  type conductivity_type

     real(kind=dp) :: tmax, dt
     integer :: nbins
     real(kind=dp), dimension(:), pointer :: displacement
     real(kind=dp) :: nsel
     integer :: nupdates=-1

  end type conductivity_type

contains

  subroutine create_conductivity(conductivity,nbins,dt)
    type(conductivity_type), intent(out) :: conductivity
    real(kind=dp), intent(in) :: dt
    integer, intent(in) :: nbins

    conductivity%dt=dt
    conductivity%nbins=nbins
    
    conductivity%tmax=dble(nbins-1)*dt

    conductivity%nsel=0.0_dp
    conductivity%nupdates=0

    allocate(conductivity%displacement(nbins))
    conductivity%displacement=0.0_dp

  end subroutine create_conductivity

  subroutine update_conductivity(conductivity,trajectory,charge,bin)
    type(conductivity_type), intent(inout) :: conductivity
    real(kind=dp), dimension(:), intent(in) :: charge
    real(kind=dp), dimension(:,:,:), intent(in) :: trajectory
    integer, intent(in) :: bin

    integer :: bin2,dbin,nbins,ntraj
    integer :: iat,natoms

    real(kind=dp), dimension(3) :: vector

    natoms=size(trajectory(1,:,1))
    nbins=conductivity%nbins
    ntraj=nbins
    do dbin=0,nbins-1
       bin2=mod(bin+dbin-1,ntraj)+1
       vector = 0.0
       do iat=1,natoms
             vector=vector+charge(iat)*(trajectory(:,iat,bin)-trajectory(:,iat,bin2))
       end do
       conductivity%displacement(dbin+1)=conductivity%displacement(dbin+1)+&
                  distance(vector)**2

    end do

    conductivity%nupdates=conductivity%nupdates+1
    conductivity%nsel=conductivity%nsel+1

  end subroutine update_conductivity

  subroutine finalize_conductivity(conductivity)
    type(conductivity_type), intent(inout) :: conductivity

    conductivity%displacement=conductivity%displacement/conductivity%nsel

    conductivity%nsel=conductivity%nsel/dble(conductivity%nupdates)

  end subroutine finalize_conductivity

  subroutine write_conductivity(out_unit,conductivity,volume,comment)
    integer, intent(in) :: out_unit
    type(conductivity_type), intent(in) :: conductivity
    real(kind=dp), intent(in) :: volume
    character(len=*), intent(in), optional :: comment

     integer :: index
     real(kind=dp) :: t
 
     if (present(comment)) write(out_unit,'("# ",A)') trim(comment)
     write(out_unit,&
       '("# dt=",F14.6,", tmax=",F14.6,", nsel=",F14.6)') &
       conductivity%dt,conductivity%tmax,conductivity%nsel
     
     do index=2,conductivity%nbins
        t=dble(index-1)*conductivity%dt
        write(out_unit,'(3F14.6)') t,conductivity%displacement(index), volume
     end do


  end subroutine write_conductivity

end module conductivity
