
program analysis

  use kinds
  use periodic_table
  use atom_methods
  use pbc
  use gofrs
  use msds
  use selections
  use colvar_utils
  use distribution
  use lammps_files

  implicit none 

!                  %----------------%
!                  |   PARAMETERS   |
!                  %----------------%

  !GENERAL PARAMETERS
  type(atom),dimension(:), allocatable :: atoms
  character(len=2), dimension(:), allocatable :: symbol
  integer :: natoms
  real(kind=dp), dimension(:,:), allocatable :: configuration,configuration_old,reference
  real(kind=dp), dimension(:,:), allocatable :: velocities
  real(kind=dp), dimension(:,:,:), allocatable :: traj
  real(kind=dp), dimension(3) :: vector
  logical :: usepbc=.true.  
  logical :: referential ! if T then particle frame otherwise lab frame

  integer :: bin,step,nsteps,skip,nskips,nprint,iatm
  integer :: timestep
  integer :: iostat

  integer, dimension(:), allocatable :: id_type
  real(kind=dp), dimension(:), allocatable :: mass
  real(kind=dp) :: nmass
  real(kind=dp), dimension(3) :: com_box
  integer :: i,j
  logical, dimension(:), allocatable :: selectionT

  !BOX PARAMETERS
  real(kind=dp) :: Lx,Ly,Lz
  type(box_type) :: box

  !XYZ PARAMETERS
  integer :: keytrj,imcon,frame,records
  character(len=200) :: filename
  type(lammpsfile) :: trajfile_lammps
  integer :: trajunit
  character(len=400) :: directory

  !INPUT PARAMETERS
  real(kind=dp) :: limit,dr,dt
  integer :: nbins

  !TYPES
  type(gofr_type) :: gofr
  type(msd_type) :: msd
  type(G_type) :: G

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ directory,filename,nbins,nskips,nprint,referential
  namelist /rdf/ limit,dr
  namelist /msd_info/ nbins, dt
  
  read(5,data)
  read(5,rdf)
  read(5,msd_info)


  write(*,*) "---------------------"
  write(*,*) "DATA INPUT :"
  write(*,*) "file positions =",trim(filename)
  write(*,*) "nskips =",nskips
  write(*,*) "nprint =",nprint
  write(*,*) "nprint =",nbins
  write(*,*) "limit =",limit
  write(*,*) "dr =",dr
  write(*,*) "referential =",referential
  write(*,*) "---------------------"

  trajunit=10
  call open_read_lammps(trajunit,trim(directory)//trim(filename),trajfile_lammps,natoms,timestep,Lx,Ly,Lz,iostat)

  write(*,*) "natoms =",natoms
  write(*,*) "box =",Lx,Ly,Lz
  write(*,*) "---------------------"


!                  %----------------%
!                  |   CREATE BOX   |
!                  %----------------%

  call create_cubic_box(box=box,length=Lx)

!                  %----------------%
!                  |    ALLOCATE    |
!                  %----------------%


  allocate(symbol(natoms))
  allocate(atoms(natoms))
  allocate(configuration(3,natoms))
  allocate(configuration_old(3,natoms))
  allocate(velocities(3,natoms))  

  allocate(id_type(natoms))
  allocate(mass(natoms))
   
  allocate(traj(3,natoms,nbins))
  allocate(reference(3,natoms))
  allocate(selectionT(natoms))
  

 !                 %----------------%
 !                 |  CREATE TOOLS  |
 !                 %----------------%

  selectionT= .true.
  call create_gofr(gofr,limit=limit,dr=dr)
  call create_msd(msd,nbins=nbins,dt=dt)
  call create_G(G,nbins,limit,dr)


 !                 %----------------%
 !                 |    LOOP  !!!   |
 !                 %----------------%


  step=0
  do
  
 !                 %----------------%
 !                 |  READ LAMMPS   |
 !                 %----------------%

     call read_lammps(trajfile_lammps,configuration,velocities,id_type,symbol,iostat)       

     !#############################################################################
     !Unfold trajectory
     if (step.ne.0) then
        do iatm=1,natoms
           vector=configuration(:,iatm)-configuration_old(:,iatm)
           configuration(:,iatm)=configuration_old(:,iatm) +&
                                 minimum_image(vector,box)
        end do
     end if
     configuration_old=configuration

     !Assign mass of particles
     mass = 1.0

     !Center of mass of the box
     com_box=0.0_dp
     do iatm=1,natoms
        com_box(:) = com_box(:) + mass(iatm)*configuration(:,iatm)
     end do
     com_box=com_box/(natoms)

     !Remove center of mass of the box
     do iatm=1,natoms
        configuration(:,iatm) = configuration(:,iatm) - com_box(:)
     end do

     !#############################################################################


     if (iostat.ne.0) exit

     step=step+1
     if (step<=nbins) then
        bin=step
        traj(:,:,bin)=configuration(:,:)
        
     else
        
        bin=mod(step-1,nbins)+1
        reference=traj(:,:,bin)

!                  %----------------%
!                  |  UPDATE TOOLS  |
!                  %----------------%

        call update_gofr(gofr,atoms,configuration,selectionT,&
                      selectionT,box)
        call update_msd(msd,atoms,traj,bin,nbins,selectionT)
        call update_G(G,natoms,reference,configuration,&
                      velocities,selectionT,selectionT,box,usepbc,&
                      referential)        

        traj(:,:,bin)=configuration(:,:)
     end if
     if (mod(step,nprint)==0) write(*,*) step
  end do

  nsteps=step

  write(*,*) "nsteps =",nsteps

  call close_read_lammps(trajfile_lammps)

!                 %----------------%
!                 | FINALIZE TOOLS |
!                 %----------------%

  call finalize_G(G,box)

 !                 %------------------%
 !                 | DUMP OUT OUTPUTS |
 !                 %------------------%

  !Radial distribution function
  open(unit=13,file="gofr.dat")
  call write_gofr(13,gofr,"Ar-Ar")
  close(13)
  !Mean square displacement
  open(unit=14,file="msd.dat")
  call write_msd(14,msd,"Ar-Ar")
  close(14)
  !G barbar
  open(15,file="G.dat")
     call write_G(15,G,"LJ")
  close(15)
  !Density Map
  !call density_map_xy(G,Lx,Ly)
  !!Rot,Div,V_r,V_theta

CONTAINS

  subroutine density_map_xy(G,a1,a2)
    implicit none

    type(G_type), intent(in) :: G
    real(kind=dp),intent(in) :: a1,a2

    real(kind=dp) :: r,x,y
    integer :: index,i,j,L,M
    

    OPEN(unit=16,file="DensityMap_xy.dat",status="unknown")
    WRITE(16,'("# X",F14.6,", Y",F14.6,", Z",F14.6)')
    
    L=int(a1)
    M=int(a2)

    DO i=-L,L
      DO j=-M,M
        x=dble(i)
        y=dble(j)
        r=int(sqrt(x**2+y**2))
        WRITE(16,*)x,y,G%density(int(r))*x/sqrt(x**2+y**2)
      END DO
      WRITE(16,*)
    END DO
    CLOSE(16)
  end subroutine density_map_xy

end program analysis
