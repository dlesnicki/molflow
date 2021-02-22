
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
  logical, dimension(:), allocatable :: selection_C
  logical, dimension(:), allocatable :: selection_C_CH2
  logical, dimension(:), allocatable :: selection_C_CH3
  logical, dimension(:), allocatable :: selection_H_CH2
  logical, dimension(:), allocatable :: selection_H_CH3
  logical :: usepbc=.true.  
  logical :: referential ! if T then particle frame otherwise lab frame
  integer, dimension(:), allocatable :: id_type
  integer :: bin,step,nsteps,skip,nskips,nprint,iatm
  integer :: iostat
  real(kind=dp) :: nmass
  type(lammpsfile) :: trajfile_lammps
  integer :: i,j
  real(kind=dp), dimension(:,:), allocatable :: com

  !BOX PARAMETERS
  real(kind=dp) :: Lx,Ly,Lz
  type(box_type) :: box

  !DLPOLY PARAMETERS
  integer :: keytrj,imcon,frame,records

  !INPUT PARAMETERS
  real(kind=dp) :: limit,dr,dt
  integer :: nbins,timestep

  !TYPES
  type(gofr_type) :: gofrAr
  type(gofr_type) :: gofrC
  type(msd_type) :: msdAr
  type(G_type) :: GAr

  integer :: Nmonomer,Nmolecule
  integer :: Nt

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ directory,filename,nskips,nprint,nsteps,Nmonomer
  namelist /rdf/ limit,dr
  
  read(5,data)
  read(5,rdf)


  write(*,*) "---------------------"
  write(*,*) "DATA INPUT :"
  write(*,*) "nskips =",nskips
  write(*,*) "nsteps =",nsteps
  write(*,*) "nprint =",nprint
  write(*,*) "limit =",limit
  write(*,*) "dr =",dr
  write(*,*) "---------------------"

  trajunit=10
  call open_read_lammps(trajunit,trim(directory)//trim(filename),trajfile_lammps,natoms,timestep,Lx,Ly,Lz,iostat)

  write(*,*) "natoms =",natoms
  write(*,*) "box =",Lx,Ly,Lz
  write(*,*) "---------------------"

  Nmolecule=natoms/Nmonomer

  write(*,*) "Nmolecules = ",Nmolecule

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
  allocate(traj_pos(3,Nmolecule,nbins))

  allocate(selection_C(natoms))
  allocate(selection_C_CH2(natoms))
  allocate(selection_C_CH3(natoms))
  allocate(selection_H_CH2(natoms))
  allocate(selection_H_CH3(natoms))

!                  %----------------%
!                  |  CREATE TOOLS  |
!                  %----------------%

  selectionT= .true.
  call create_gofr(gofrAr,limit=limit,dr=dr)


!                  %----------------%
!                  |    LOOP  !!!   |
!                  %----------------%

  do step=1,nskips
     call read_lammps(trajfile_lammps,configuration,velocities,id_type,iostat)
  end do



  Nt=nsteps-nskips
  do step=1,Nt

        write(*,*)step+nskips

!                  %----------------%
!                  |  READ LAMMPS   |
!                  %----------------%

         
        call read_lammps(trajfile_lammps,configuration,velocities,id_type,iostat)

!        !Assign mass of particles
        do iatm=1,natoms
           if ((id_type(iatm) .eq. 908) .or. (id_type(iatm) .eq. 907)) then
                   mass(iatm) = 12.011
                   selection_C(iatm) = .true.
           else if ((id_type(iatm) .eq. 910) .or. (id_type(iatm) .eq. 909)) then
                   mass(iatm) = 1.008
           end if
        end do
!
!
!!        selection_C_CH2 = .false.
!!        selection_C_CH3 = .false.
!!        selection_H_CH2 = .false.
!!        selection_H_CH3 = .false.
!!        do iatm=1,natoms
!!           if (id_type(iatm) .eq. 908) then !C of CH2
!!                selection_C_CH2(iatm) = .true.           
!!           else if (id_type(iatm) .eq. 907) then !C of CH3 
!!                selection_C_CH3(iatm) = .true.           
!!           else if (id_type(iatm) .eq. 910) then !H of CH2
!!                selection_H_CH2(iatm) = .true.           
!!           else if (id_type(iatm) .eq. 909) then !H of CH3
!!                selection_H_CH3(iatm) = .true.    
!!           end if       
!!        end do

        !Mass of one monomer
        nmass = 0.0_dp
        do iatm=1,natoms
           nmass = nmass + mass(iatm)
        end do
        nmass = nmass/Nmolecule
     
        !Center of mass of the box
        com_box=0.0_dp  
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

!                  %----------------%
!                  |  UPDATE TOOLS  |
!                  %----------------%

        !WRITE(*,*)'UPDATE GOFR'
        call update_gofr_N(gofrAr,Nmolecule,com,selectionT,&
                     selectionT,box)


     !if (mod(step,nprint)==0) write(*,*) step
  end do

  call close_read_lammps(trajfile_lammps)

!                  %----------------%
!                  | FINALIZE TOOLS |
!                  %----------------%

  call finalize_gofr(gofrAr,box)

!                  %------------------%
!                  | DUMP OUT OUTPUTS |
!                  %------------------%

  !Radial distribution function
  open(unit=13,file="gofr-com.dat")
  call write_gofr(13,gofrAr,"Ar-Ar")
  close(13)


end program analysis
