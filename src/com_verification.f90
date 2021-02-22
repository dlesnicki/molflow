
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
character(len=200) :: filename_com
character(len=2) :: symbol
integer :: trajunit
integer :: trajunit_com
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
integer :: bin,step,nsteps,skip,nskips,nprint,iatm,iat
integer :: index
integer :: iostat
real(kind=dp) :: nmass,total_mass
type(lammpsfile) :: trajfile_lammps
integer :: i,j
real(kind=dp), dimension(:,:), allocatable :: com
real(kind=dp), dimension(:,:), allocatable :: com_from_trj

!BOX PARAMETERS
real(kind=dp) :: Lx,Ly,Lz
type(box_type) :: box

!DLPOLY PARAMETERS
integer :: keytrj,imcon,frame,records

!INPUT PARAMETERS
real(kind=dp) :: limit,dr,dt
integer :: nbins,timestep

!TYPES

integer :: Nmonomer,Nmolecule
integer :: time,number,id_mol
real(kind=dp) :: pos_molecule(3)

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ directory,filename,filename_com,nsteps,Nmonomer
  read(5,data)
  
  
  write(*,*) "---------------------"
  write(*,*) "DATA INPUT :"
  write(*,*) "trajectory = ", trim(directory)//trim(filename)
  write(*,*) "com = ", trim(directory)//trim(filename)
  write(*,*) "nsteps = ",nsteps
  write(*,*) "atoms per monomer = ",Nmonomer
  write(*,*) "---------------------"
 

!                  %----------------%
!                  |     LAMMPS     |
!                  %----------------%
  
  trajunit=10
  call open_read_lammps(trajunit,trim(directory)//trim(filename),trajfile_lammps,natoms,timestep,Lx,Ly,Lz,iostat)
  
  write(*,*) "natoms =",natoms
  write(*,*) "box =",Lx,Ly,Lz
  write(*,*) "---------------------"
  
  Nmolecule=natoms/Nmonomer
  
  write(*,*) "Nmolecules = ",Nmolecule

!                  %-------------------------%
!                  |      COM to com.xyz     |
!                  %-------------------------%

  allocate(com(3,Nmolecule))

  ! INPUT
  trajunit_com=12
  open(unit=trajunit_com,file=trim(directory)//trim(filename_com))
  read(trajunit_com,*)
  read(trajunit_com,*)
  read(trajunit_com,*)

  ! OUTPUT
  open(unit=14,file="com.xyz")

  do step=1,nsteps
     read(trajunit_com,*)timestep,Nmolecule
     do iat=1,Nmolecule
        read(trajunit_com,*)index,com(:,iat)
     end do

     write(14,*)Nmolecule
     write(14,*)timestep
     do iat=1,Nmolecule
        write(14,*)"Ar",com(:,iat)
     end do
  end do
  close(14) 
  close(12) 

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
  allocate(com_from_trj(3,Nmolecule))
  allocate(reference(3,Nmolecule))
  allocate(selectionT(Nmolecule))
  allocate(traj_pos(3,Nmolecule,nbins))

  allocate(selection_C(natoms))
  allocate(selection_C_CH2(natoms))
  allocate(selection_C_CH3(natoms))
  allocate(selection_H_CH2(natoms))
  allocate(selection_H_CH3(natoms))


!                  %----------------%
!                  |    LOOP  !!!   |
!                  %----------------%

  open(unit=88,file="com_from_trj.xyz")

  do step=1,nsteps
        
        write(*,*)step

!                  %----------------%
!                  |  READ LAMMPS   |
!                  %----------------%

         
        call read_lammps(trajfile_lammps,configuration,velocities,id_type,iostat)

        !Assign mass of particles
        selection_C = .false.
        do iatm=1,natoms
           if ((id_type(iatm) .eq. 908) .or. (id_type(iatm) .eq. 907)) then
                   mass(iatm) = 12.011_dp
                   selection_C(iatm) = .true.
           else if ((id_type(iatm) .eq. 910) .or. (id_type(iatm) .eq. 909)) then
                   mass(iatm) = 1.008_dp
           end if
        end do


        selection_C_CH2 = .false.
        selection_C_CH3 = .false.
        selection_H_CH2 = .false.
        selection_H_CH3 = .false.
        do iatm=1,natoms
           if (id_type(iatm) .eq. 908) then !C of CH2
                selection_C_CH2(iatm) = .true.           
           else if (id_type(iatm) .eq. 907) then !C of CH3 
                selection_C_CH3(iatm) = .true.           
           else if (id_type(iatm) .eq. 910) then !H of CH2
                selection_H_CH2(iatm) = .true.           
           else if (id_type(iatm) .eq. 909) then !H of CH3
                selection_H_CH3(iatm) = .true.    
           end if       
        end do

        write(88,*)Nmolecule
        write(88,*)
        com_from_trj=0.0_dp
        j=0
        do iatm=1,Nmolecule
           nmass=0.0_dp
           do i=1,Nmonomer
              j = j+1
              nmass=nmass+mass(j)
              com_from_trj(:,iatm) = com_from_trj(:,iatm) + mass(j)*configuration(:,j)
           end do
           com_from_trj(:,iatm)=com_from_trj(:,iatm)/nmass
           write(88,*)"Ar",com_from_trj(:,iatm)
        end do
      
  end do
  call close_read_lammps(trajfile_lammps)
  close(88)

!                  %----------------------%
!                  |  REMOVE COM of BOX   |
!                  %----------------------%

  total_mass=nmass*Nmolecule
  write(*,*)"total_mass",total_mass

  open(unit=90,file="com_from_trj.xyz")
  open(unit=91,file="new_com_from_trj.xyz")
  do i=1,nsteps
     read(90,*)
     read(90,*)
     com_box(:) = 0.0_dp
     do iatm=1,Nmolecule
        read(90,*)symbol,com(:,iatm)
        com_box(:) = com_box(:) + com(:,iatm)*nmass
     end do
     com_box(:) = com_box(:)/total_mass


     write(91,*)Nmolecule
     write(91,*)
     do iatm=1,Nmolecule
        write(91,*)symbol,com(:,iatm)-com_box(:)
     end do
  end do
  close(90)
  close(91)
end program analysis
