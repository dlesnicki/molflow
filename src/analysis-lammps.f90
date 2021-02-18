
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
  integer :: natoms
  real(kind=dp), dimension(:,:), allocatable :: configuration,configuration_old,reference
  real(kind=dp), dimension(:,:), allocatable :: velocities
  real(kind=dp), dimension(:,:,:), allocatable :: traj_com
  real(kind=dp), dimension(3) :: vector
  logical, dimension(:), allocatable :: selectionAr
  logical :: usepbc=.true.  
  logical :: referential ! if T then particle frame otherwise lab frame

  integer :: bin,step,nsteps,skip,nskips,nprint,iatm
  integer :: timestep
  integer :: iostat

  integer, dimension(:), allocatable :: id_type
  real(kind=dp), dimension(:), allocatable :: mass
  real(kind=dp) :: nmass
  real(kind=dp), dimension(3) :: com_box
  real(kind=dp), dimension(:,:), allocatable :: com_pos
  real(kind=dp), dimension(:,:), allocatable :: com_vel
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
  integer :: Nmonomer,Nmolecule

  !TYPES
  type(G_type) :: Gcom

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ directory,filename,nbins,nskips,nprint,Nmonomer,referential
  namelist /rdf/ limit,dr
  
  read(5,data)
  read(5,rdf)

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

  Nmolecule=natoms/Nmonomer

  write(*,*) "Nmolecules = ",Nmolecule

 !                 %----------------%
 !                 |   CREATE BOX   |
 !                 %----------------%

  call create_cubic_box(box=box,length=Lx)

 !                 %----------------%
 !                 |    ALLOCATE    |
 !                 %----------------%


  allocate(atoms(natoms))
  allocate(configuration(3,natoms))
  allocate(configuration_old(3,natoms))
  allocate(velocities(3,natoms)) 

  allocate(id_type(natoms))
  allocate(mass(natoms))
  allocate(com_pos(3,Nmolecule))
  allocate(com_vel(3,Nmolecule))
   
  allocate(traj_com(3,Nmolecule,nbins))
  allocate(reference(3,Nmolecule))
  allocate(selectionT(Nmolecule))
  

 !                 %----------------%
 !                 |  CREATE TOOLS  |
 !                 %----------------%

  selectionT= .true.
  call create_G(Gcom,nbins,limit,dr)


 !                 %----------------%
 !                 |    LOOP  !!!   |
 !                 %----------------%


  step=0
  do
  
 !                 %----------------%
 !                 |  READ LAMMPS   |
 !                 %----------------%

     call read_lammps(trajfile_lammps,configuration,velocities,id_type,iostat)       

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
     nmass = nmass/dble(Nmolecule)

     !Center of mass of the box
     com_box=0.0_dp
     do iatm=1,natoms
        com_box(:) = com_box(:) + mass(iatm)*configuration(:,iatm)
     end do
     com_box=com_box/(nmass*dble(Nmolecule))

     !Remove center of mass of the box
     do iatm=1,natoms
        configuration(:,iatm) = configuration(:,iatm) - com_box(:)
     end do

     !Center of mass of each monomer
     com_pos=0.0_dp
     com_vel=0.0_dp
     j=0
     do iatm=1,Nmolecule
        do i=1,Nmonomer
             j = j+1
             com_pos(:,iatm) = com_pos(:,iatm) + mass(j)*configuration(:,j)
             com_vel(:,iatm) = com_vel(:,iatm) + mass(j)*velocities(:,j)
        end do
        com_pos(:,iatm)=com_pos(:,iatm)/nmass
        com_vel(:,iatm)=com_vel(:,iatm)/nmass
     end do
     !#############################################################################


     if (iostat.ne.0) exit

     step=step+1
     if (step<=nbins) then
        bin=step
        traj_com(:,:,bin)=com_pos
        
     else
        
        bin=mod(step-1,nbins)+1
        reference=traj_com(:,:,bin)

!                  %----------------%
!                  |  UPDATE TOOLS  |
!                  %----------------%

        WRITE(*,*)'UPDATE G'
        call update_G(Gcom,Nmolecule,reference,com_pos,&
                      com_vel,selectionT,selectionT,box,usepbc,&
                      referential)        

        traj_com(:,:,bin)=com_pos
     end if
     if (mod(step,nprint)==0) write(*,*) step
  end do

  nsteps=step

  write(*,*) "nsteps =",nsteps

  call close_read_lammps(trajfile_lammps)

 !                 %----------------%
 !                 | FINALIZE TOOLS |
 !                 %----------------%

  call finalize_G(Gcom,box)

 !                 %------------------%
 !                 | DUMP OUT OUTPUTS |
 !                 %------------------%

  !G barbar
  open(15,file="Gcom.dat")
     call write_G(15,Gcom,"com")
  close(15)
  !Density Map
  call density_map_xy(Gcom,Lx,Ly)
  !!Rot,Div,V_r,V_theta
! ! call analysis_tools(GAr,a1,a2)

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
