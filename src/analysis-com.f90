
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
  integer :: time,number,id_mol
  real(kind=dp) :: pos_molecule(3)

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ directory,filename,nskips,nprint,referential,Nmonomer
  namelist /rdf/ limit,dr
  namelist /msd/ nbins,dt
  
  read(5,data)
  read(5,rdf)
  read(5,msd)


  write(*,*) "---------------------"
  write(*,*) "DATA INPUT :"
  write(*,*) "nskips =",nskips
  write(*,*) "nprint =",nprint
  write(*,*) "limit =",limit
  write(*,*) "dr =",dr
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

!  selectionT= .true.
!  call create_gofr(gofrAr,limit=limit,dr=dr)
!  call create_gofr(gofrC,limit=limit,dr=dr)


!                  %----------------%
!                  |    LOOP  !!!   |
!                  %----------------%

  open(unit=88,file="com.xyz")
  open(unit=77,file="com.out")
  open(unit=75,file="lammps_com.xyz")
  !step=0
  do step=1,11


        read(77,*)time,number
        write(75,*)number
        write(75,*)
        do iatm=1,Nmolecule
                read(77,*)id_mol,pos_molecule(:)
                write(75,*)"C",pos_molecule(:)
        end do

          
!                  %----------------%
!                  |  READ LAMMPS   |
!                  %----------------%

         
        call read_lammps(trajfile_lammps,configuration,velocities,id_type,iostat)
       

        !!if (step.ne.0) then
        !!   do iatm=1,natoms
        !!      vector=configuration(:,iatm)-configuration_old(:,iatm)
        !!      configuration(:,iatm)=configuration_old(:,iatm) +&
        !!                            minimum_image(vector,box)
        !!   end do
        !!end if
        !!configuration_old=configuration

        !Assign mass of particles
        selection_C = .false.
        do iatm=1,natoms
           if ((id_type(iatm) .eq. 908) .or. (id_type(iatm) .eq. 907)) then
                   mass(iatm) = 12.011
                   selection_C(iatm) = .true.
           else if ((id_type(iatm) .eq. 910) .or. (id_type(iatm) .eq. 909)) then
                   mass(iatm) = 1.008
           end if
        end do

        !do iatm=1,natoms
        !   write(88,*)iatm,id_type(iatm),mass(iatm)
        !end do

!        selection_C_CH2 = .false.
!        selection_C_CH3 = .false.
!        selection_H_CH2 = .false.
!        selection_H_CH3 = .false.
!        do iatm=1,natoms
!           if (id_type(iatm) .eq. 908) then !C of CH2
!                selection_C_CH2(iatm) = .true.           
!           else if (id_type(iatm) .eq. 907) then !C of CH3 
!                selection_C_CH3(iatm) = .true.           
!           else if (id_type(iatm) .eq. 910) then !H of CH2
!                selection_H_CH2(iatm) = .true.           
!           else if (id_type(iatm) .eq. 909) then !H of CH3
!                selection_H_CH3(iatm) = .true.    
!           end if       
!        end do

        !Mass of one monomer
        nmass = 0.0_dp
        do iatm=1,natoms
           nmass = nmass + mass(iatm)
        end do
        nmass = nmass/Nmolecule
     
        !!!!Center of mass of the box  
        !!!do iatm=1,natoms 
        !!!   com_box(:) = com_box(:) + mass(iatm)*configuration(:,iatm)
        !!!end do
        !!!com_box=com_box/(nmass*Nmolecule)

        !!!write(88,*)Nmolecule
        !!!write(88,*)
        !!!do iatm=1,natoms
        !!!   configuration(:,iatm) = configuration(:,iatm) - com_box(:)
        !!!end do

        !Center of mass of each monomer

        write(88,*)Nmolecule
        write(88,*)
        com=0.0_dp
        j=0
        do iatm=1,Nmolecule
           do i=1,Nmonomer
           j = j+1
           com(:,iatm) = com(:,iatm) + mass(j)*configuration(:,j)
           end do
           com(:,iatm)=com(:,iatm)/nmass
           write(88,*)"Ar",com(:,iatm)
        end do
      
     !if (iostat.ne.0) exit

     !step=step+1
        

!                  %----------------%
!                  |  UPDATE TOOLS  |
!                  %----------------%

!        !WRITE(*,*)'UPDATE GOFR'
!        call update_gofr_N(gofrAr,Nmolecule,com,selectionT,&
!                      selectionT,box)
!        call update_gofr(gofrC,atoms,configuration,selection_C,&
!                      selection_C,box)
!
!
     if (mod(step,nprint)==0) write(*,*) step
  end do

  close(88)
  !nsteps=step

  !write(*,*) "nsteps =",nsteps
!
!  call close_read_lammps(trajfile_lammps)
!
!!                  %----------------%
!!                  | FINALIZE TOOLS |
!!                  %----------------%
!
!  call finalize_gofr(gofrAr,box)
!  call finalize_gofr(gofrC,box)
!
!!                  %------------------%
!!                  | DUMP OUT OUTPUTS |
!!                  %------------------%
!
!  !Radial distribution function
!  open(unit=13,file="gofr-com.dat")
!  call write_gofr(13,gofrAr,"Ar-Ar")
!  close(13)
!  open(unit=13,file="gofr-C_C.dat")
!  call write_gofr(13,gofrAr,"C-C")
!  close(13)
!
!CONTAINS
!
!  subroutine density_map_xy(G,a1,a2)
!    implicit none
!
!    type(G_type), intent(in) :: G
!    real(kind=dp),dimension(3),intent(in) :: a1,a2
!
!    real(kind=dp) :: r,x,y
!    integer :: index,i,j,L,M
!    
!
!    OPEN(unit=16,file="DensityMap_xy.dat",status="unknown")
!    WRITE(16,'("# X",F14.6,", Y",F14.6,", Z",F14.6)')
!    
!    L=int(a1(1))
!    M=int(a2(2))
!
!    DO i=-L,L
!      DO j=-M,M
!        x=dble(i)
!        y=dble(j)
!        r=int(sqrt(x**2+y**2))
!        WRITE(16,*)x,y,G%density(int(r))*x/sqrt(x**2+y**2)
!      END DO
!      WRITE(16,*)
!    END DO
!    CLOSE(16)
!  end subroutine density_map_xy

end program analysis
