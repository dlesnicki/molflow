
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
  character(len=400) :: directory
  character(len=200) :: filename
  integer :: trajunit
  type(atom),dimension(:), allocatable :: atoms
  integer :: natoms
  real(kind=dp), dimension(:,:), allocatable :: configuration,configuration_old,reference
  real(kind=dp), dimension(:,:), allocatable :: velocities
  real(kind=dp), dimension(:,:,:), allocatable :: traj_pos
  real(kind=dp), dimension(:), allocatable :: mass
  real(kind=dp), dimension(3) :: vector
  real(kind=dp), dimension(3) :: com_box
  logical, dimension(:), allocatable :: selectionT
  logical :: usepbc=.true.  
  logical :: referential ! if T then particle frame otherwise lab frame
  integer, dimension(:), allocatable :: id_type
  integer :: bin,step,nsteps,skip,nskips,nprint,iatm
  integer :: iostat
  real(kind=dp) :: nmass
  type(lammpsfile) :: trajfile_lammps
  integer :: i,j
  real(kind=dp), dimension(:,:), allocatable :: com
  real(kind=dp), dimension(:,:), allocatable :: vcm

  !BOX PARAMETERS
  real(kind=dp) :: Lx,Ly,Lz
  type(box_type) :: box

  !DLPOLY PARAMETERS
  integer :: keytrj,imcon,frame,records

  !INPUT PARAMETERS
  real(kind=dp) :: dt
  integer :: nbins,timestep

  !TYPES
  type(gofr_type) :: gofrAr
  type(msd_type) :: msdAr
  type(G_type) :: GAr

  integer :: Nmonomer,Nmolecule

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ directory,filename,nskips,nprint,referential,Nmonomer
  namelist /msd/ nbins,dt
  
  read(5,data)
  read(5,msd)


  write(*,*) "---------------------"
  write(*,*) "DATA INPUT :"
  write(*,*) "nskips =",nskips
  write(*,*) "nprint =",nprint
  write(*,*) "nbins = ",nbins
  write(*,*) "dt =",dt
  write(*,*) "referential =",referential
  write(*,*) "---------------------"

  trajunit=10
  call open_read_lammps(trajunit,trim(directory)//trim(filename),trajfile_lammps,natoms,timestep,Lx,Ly,Lz,iostat)

  write(*,*) "natoms =",natoms
  write(*,*) "box =",Lx,Ly,Lz
  write(*,*) "---------------------"

  Nmolecule=natoms/Nmonomer

!                  %----------------%
!                  |   CREATE BOX   |
!                  %----------------%

  call create_cubic_box(box=box,length=Lx)

!                  %----------------%
!                  |    ALLOCATE    |
!                  %----------------%


  allocate(atoms(natoms))
  allocate(configuration(3,natoms))
  allocate(configuration_old(3,natoms))
  allocate(velocities(3,natoms))  
  
  allocate(id_type(natoms))
  allocate(mass(natoms))
  allocate(com(3,Nmolecule))
  allocate(reference(3,Nmolecule))
  allocate(selectionT(Nmolecule))
  allocate(vcm(3,Nmolecule))
  allocate(traj_pos(3,Nmolecule,nbins))

!                  %----------------%
!                  |  CREATE TOOLS  |
!                  %----------------%

  selectionT= .true.
  call create_msd(msdAr,nbins=nbins,dt=dt)


!                  %----------------%
!                  |    LOOP  !!!   |
!                  %----------------%


  step=0
  do
  
!                  %----------------%
!                  |  READ LAMMPS   |
!                  %----------------%

         
     do skip=1,nskips
        call read_lammps(trajfile_lammps,configuration,velocities,id_type,iostat)
        if (step.ne.0) then
           do iatm=1,natoms
              vector=configuration(:,iatm)-configuration_old(:,iatm)
              configuration(:,iatm)=configuration_old(:,iatm) +&
                                    minimum_image(vector,box)
           end do
        end if
        configuration_old=configuration

        !Assign mass of particles
        do iatm=1,natoms
           if ((id_type(iatm) .eq. 908) .or. (id_type(iatm) .eq. 907)) then
                   mass(iatm) = 12.011
           else if ((id_type(iatm) .eq. 910) .or. (id_type(iatm) .eq. 909)) then
                   mass(iatm) = 1.008
           end if
        end do

        !Mass of one monomer
        nmass = 0.0_dp
        do iatm=1,natoms
           nmass = nmass + mass(iatm)
        end do
        nmass = nmass/Nmolecule 

        com_box = 0.0_dp     
        !Center of mass of the box  
        do iatm=1,natoms 
           com_box(:) = com_box(:) + mass(iatm)*configuration(:,iatm)
        end do
        com_box=com_box/(nmass*Nmolecule)

        do iatm=1,natoms
           configuration(:,iatm) = configuration(:,iatm) - com_box(:)
        end do

        !Center of mass of each monomer
        com=0.0_dp
        !vcm=0.0_dp

        j=0
        do iatm=1,Nmolecule
           do i=1,Nmonomer
           j = j+1
           
           com(:,iatm) = com(:,iatm) + mass(j)*configuration(:,j)
           !vcm(:,iatm) = vcm(:,iatm) + mass(j)*velocities(:,j)
           end do
           com(:,iatm)=com(:,iatm)/nmass
           !vcm(:,iatm)=vcm(:,iatm)/nmass
        end do
        
     end do
      
     if (iostat.ne.0) exit

     step=step+1
     if (step<=nbins) then
        bin=step
        traj_pos(:,:,bin)=com
        
     else
        
        bin=mod(step-1,nbins)+1
        reference=traj_pos(:,:,bin)

!                  %----------------%
!                  |  UPDATE TOOLS  |
!                  %----------------%


        !WRITE(*,*)'UPDATE MSD'
        call update_msd(msdAr,Nmolecule,traj_pos,bin,nbins,selectionT)


        traj_pos(:,:,bin)=com
     end if
     if (mod(step,nprint)==0) write(*,*) step
  end do

  nsteps=step

  write(*,*) "nsteps =",nsteps

  call close_read_lammps(trajfile_lammps)

!                  %----------------%
!                  | FINALIZE TOOLS |
!                  %----------------%

  call finalize_msd(msdAr)

!                  %------------------%
!                  | DUMP OUT OUTPUTS |
!                  %------------------%

  !Mean square displacement
  open(unit=14,file="msd.dat")
  call write_msd(14,msdAr,"Ar-Ar")
  close(14)

CONTAINS

  subroutine density_map_xy(G,a1,a2)
    implicit none

    type(G_type), intent(in) :: G
    real(kind=dp),dimension(3),intent(in) :: a1,a2

    real(kind=dp) :: r,x,y
    integer :: index,i,j,L,M
    

    OPEN(unit=16,file="DensityMap_xy.dat",status="unknown")
    WRITE(16,'("# X",F14.6,", Y",F14.6,", Z",F14.6)')
    
    L=int(a1(1))
    M=int(a2(2))

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
