
module msds

  use kinds
  use periodic_table
  use selections
  use colvar_utils

  implicit none 

  type msd_type

     real(kind=dp) :: tmax, dt
     integer :: nbins
     real(kind=dp), dimension(:), pointer :: displacement
     real(kind=dp) :: nsel
     integer :: nupdates=-1

  end type msd_type

contains

  subroutine create_msd(msd,nbins,dt)
    type(msd_type), intent(out) :: msd
    real(kind=dp), intent(in) :: dt
    integer, intent(in) :: nbins

    msd%dt=dt
    msd%nbins=nbins
    
    msd%tmax=dble(nbins-1)*dt

    msd%nsel=0.0_dp
    msd%nupdates=0

    allocate(msd%displacement(nbins))
    msd%displacement=0.0_dp

  end subroutine create_msd

  subroutine update_msd(msd,atoms,trajectory,bin,ntraj,selection)
    type(msd_type), intent(inout) :: msd
    type(atom), dimension(:), intent(in) :: atoms
    real(kind=dp), dimension(:,:,:), intent(in) :: trajectory
    integer, intent(in) :: bin,ntraj
    logical, dimension(:) :: selection

    integer :: bin2,dbin,nbins
    integer :: iat,natoms

    real(kind=dp), dimension(3) :: vector

    natoms=size(atoms)
    nbins=msd%nbins

    do dbin=0,nbins-1
       bin2=mod(bin+dbin-1,ntraj)+1

       do iat=1,natoms
          if (selection(iat)) then
             vector=trajectory(:,iat,bin)-trajectory(:,iat,bin2)
             msd%displacement(dbin+1)=msd%displacement(dbin+1)+&
                  distance(vector)**2
          end if
       end do

    end do

    msd%nupdates=msd%nupdates+1
    msd%nsel=msd%nsel+dble(count_selection(selection))

  end subroutine update_msd

  subroutine finalize_msd(msd)
    type(msd_type), intent(inout) :: msd

    msd%displacement=msd%displacement/msd%nsel

    msd%nsel=msd%nsel/dble(msd%nupdates)

  end subroutine finalize_msd

  subroutine write_msd(out_unit,msd,comment)
    integer, intent(in) :: out_unit
    type(msd_type), intent(in) :: msd
    character(len=*), intent(in), optional :: comment

     integer :: index
     real(kind=dp) :: t
 
     if (present(comment)) write(out_unit,'("# ",A)') trim(comment)
     write(out_unit,&
       '("# dt=",F14.6,", tmax=",F14.6,", nsel=",F14.6)') &
       msd%dt,msd%tmax,msd%nsel
     
     do index=2,msd%nbins
        t=dble(index-1)*msd%dt
        write(out_unit,'(2F14.6)') t,msd%displacement(index)
     end do


  end subroutine write_msd

end module msds
