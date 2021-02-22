
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
  real(kind=dp), dimension(:,:), allocatable :: com
  real(kind=dp), dimension(:,:,:), allocatable :: traj_pos
  real(kind=dp), dimension(:), allocatable :: mass
  real(kind=dp), dimension(:), allocatable :: RG_xx,RG_yy,RG_zz,RG_xy,RG_xz,RG_yz
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
  real(kind=dp), dimension(:,:), allocatable :: RG

  !BOX PARAMETERS
  real(kind=dp) :: Lx,Ly,Lz
  type(box_type) :: box

  !DLPOLY PARAMETERS
  integer :: keytrj,imcon,frame,records

  !INPUT PARAMETERS
  real(kind=dp) :: limit,dr,dt
  integer :: nbins,timestep

  !TYPES
  real(kind=dp) :: sum_RG
  integer :: Nmonomer,Nmolecule
  integer :: Nt

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ directory,filename,nprint,nsteps,Nmonomer
  read(5,data)

  write(*,*) "---------------------"
  write(*,*) "DATA INPUT :"
  write(*,*) "nprint =",nprint
  write(*,*) "---------------------"

  trajunit=10
  call open_read_lammps(trajunit,trim(directory)//trim(filename),trajfile_lammps,natoms,timestep,Lx,Ly,Lz,iostat)

  write(*,*) "natoms =",natoms
  write(*,*) "box =",Lx,Ly,Lz

  Nmolecule=natoms/Nmonomer

  write(*,*) "nmolecules =",Nmolecule
  write(*,*) "---------------------"

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
  allocate(RG(3,Nmolecule))
  allocate(reference(3,Nmolecule))
  allocate(selectionT(Nmolecule))
  allocate(traj_pos(3,Nmolecule,nbins))
  allocate(com(3,Nmolecule))
  allocate(RG_xx(Nmolecule))
  allocate(RG_yy(Nmolecule))
  allocate(RG_zz(Nmolecule))


  open(unit=55,file='RG.dat')
  
!                  %----------------%
!                  |    LOOP  !!!   |
!                  %----------------%
  do step=1,nsteps
  
!                  %----------------%
!                  |  READ LAMMPS   |
!                  %----------------%

        write(*,*)step       
        call read_lammps(trajfile_lammps,configuration,velocities,id_type,iostat)

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
        !write(*,*)"mass monomer = ",nmass

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
        j=0
        do iatm=1,Nmolecule
           do i=1,Nmonomer
           j = j+1
           com(:,iatm) = com(:,iatm) + mass(j)*configuration(:,j)
           end do
           com(:,iatm)=com(:,iatm)/nmass
        end do


        !Radius of gyration
        sum_RG = 0.0_dp

        j=0
        do iatm=1,Nmolecule

           RG_xx(iatm)=0.0_dp
           RG_yy(iatm)=0.0_dp
           RG_zz(iatm)=0.0_dp
           do i=1,Nmonomer
              j = j+1
           
              configuration(:,j) = configuration(:,j) - com(:,iatm)    
              RG_xx(iatm) = RG_xx(iatm) + mass(j)*(configuration(1,j))**2
              RG_yy(iatm) = RG_yy(iatm) + mass(j)*(configuration(2,j))**2
              RG_zz(iatm) = RG_zz(iatm) + mass(j)*(configuration(3,j))**2
           end do

           sum_RG = sum_RG + (RG_xx(iatm)+RG_yy(iatm)+RG_zz(iatm))
        end do
        sum_RG = sum_RG/(Nmolecule*nmass)

!                  %----------------%
!                  |  UPDATE TOOLS  |
!                  %----------------%

        write(55,*)step,sqrt(sum_RG)

  end do

!                  %----------------%
!                  | FINALIZE TOOLS |
!                  %----------------%

  close(55)
  call close_read_lammps(trajfile_lammps)

end program analysis
