
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
  integer, dimension(:), allocatable :: at_mol
  logical, dimension(:), allocatable :: selectionT
  logical, dimension(:,:), allocatable :: selection_carb
  character(len=2), dimension(:), allocatable :: symbol
  logical :: usepbc=.true.  
  logical :: referential ! if T then particle frame otherwise lab frame
  integer, dimension(:), allocatable :: id_type
  integer :: bin,step,nsteps,skip,nskips,nprint,iatm,iatm2, imol
  integer :: iostat
  type(lammpsfile) :: trajfile_lammps
  integer :: i,j
  real(kind=dp), dimension(:,:), allocatable :: com
  real(kind=dp), dimension(:,:), allocatable :: vcm
  real(kind=dp) :: box_mass, mol_mass

  character(len=2), dimension(4) :: list_labels

  !BOX PARAMETERS
  real(kind=dp) :: Lx,Ly,Lz
  type(box_type) :: box

  !DLPOLY PARAMETERS
  integer :: keytrj,imcon,frame,records

  !INPUT PARAMETERS
  real(kind=dp) :: limit,dr,dt
  integer :: nbins,timestep

  !TYPES
  type(gofr_type),dimension(:,:),allocatable :: gofr
  type(msd_type) :: msdAr
  type(G_type) :: GAr

  integer :: Nmolecule

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ directory,filename,nskips,nprint,referential
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

  Nmolecule=natoms/(2)

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
  allocate(symbol(natoms))
  allocate(com(3,Nmolecule))
  allocate(reference(3,Nmolecule))
  allocate(selectionT(Nmolecule))
  allocate(selection_carb(4,Nmolecule))
  allocate(vcm(3,Nmolecule))
  allocate(traj_pos(3,Nmolecule,nbins))
  allocate(at_mol(Nmolecule))

  allocate(gofr(4,4))

!                  %----------------%
!                  |  CREATE TOOLS  |
!                  %----------------%

  selectionT= .false.
  do iatm=1,4
     do iatm2=iatm,4
        call create_gofr(gofr(iatm,iatm2),limit=limit,dr=dr)
     enddo
  enddo
  call create_msd(msdAr,nbins=nbins,dt=dt)
  call create_G(GAr,nbins,limit,dr)


!                  %----------------%
!                  |    LOOP  !!!   |
!                  %----------------%


  step=0
  do
  
!                  %----------------%
!                  |  READ LAMMPS   |
!                  %----------------%
         
     do skip=1,nskips
        call read_lammps(trajfile_lammps,configuration,velocities,id_type,symbol,iostat)
        if (step.ne.0) then
           do iatm=1,natoms
              vector=configuration(:,iatm)-configuration_old(:,iatm)
              configuration(:,iatm)=configuration_old(:,iatm) +&
                                    minimum_image(vector,box)
           end do
        end if
        configuration_old=configuration

        if (step.eq.0) then
        !Assign mass of particles
        list_labels = [Character(len=2) :: 'C','Li','Na','K']
        imol = 0
        do iatm=1,natoms
           SELECT CASE(symbol(iatm))
           CASE('C')
                   mass(iatm) = 12
                   imol = imol + 1
                   at_mol(imol)=4
                   selection_carb(1,imol)=.true.
           CASE('O')
                   mass(iatm) = 15.999
           CASE('Li') 
                   mass(iatm) = 7
                   imol = imol + 1
                   at_mol(imol)=1
                   selection_carb(2,imol)=.true.
           CASE('Na')
                   mass(iatm) = 39
                   imol = imol + 1
                   at_mol(imol)=1
                   selection_carb(3,imol)=.true.
           CASE('K') 
                   mass(iatm) = 22
                   imol = imol + 1
                   at_mol(imol)=1
                   selection_carb(4,imol)=.true.
           END SELECT
        end do

        !Check number of atoms
        do iatm=1,4
           write(*,*) "There is", count(selection_carb(iatm,:)), "atom of type", list_labels(iatm)
        enddo
        endif
        
        !Center of mass of the box
        box_mass = 0.0  
        do iatm=1,natoms 
           com_box(:) = com_box(:) + mass(iatm)*configuration(:,iatm)
           box_mass   = box_mass + mass(iatm)
        end do
        com_box=com_box/(box_mass)
com_box = 0.0
        do iatm=1,natoms
           configuration(:,iatm) = configuration(:,iatm) - com_box(:)
        end do

        !Center of mass of each monomer
        com=0.0_dp
        vcm=0.0_dp

        j=0
        do iatm=1,Nmolecule
           mol_mass = 0.0
           do i=1,at_mol(iatm)
              j = j+1
              com(:,iatm) = com(:,iatm) + mass(j)*configuration(:,j)
              vcm(:,iatm) = vcm(:,iatm) + mass(j)*velocities(:,j)
              mol_mass = mol_mass + mass(j)
           end do
           com(:,iatm)=com(:,iatm)/mol_mass
           vcm(:,iatm)=vcm(:,iatm)/mol_mass
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

        !WRITE(*,*)'UPDATE GOFR'
        do iatm=1,4
           do iatm2=iatm,4
              call update_gofr(gofr(iatm,iatm2),atoms,com,selectionT,&
                            selectionT,box)
           enddo
        enddo

        !WRITE(*,*)'UPDATE MSD'
        call update_msd(msdAr,atoms,traj_pos,bin,nbins,selectionT)

        WRITE(*,*)'UPDATE G'
        call update_G(GAr,natoms,reference,com,&
                      vcm,selectionT,selectionT,box,usepbc,&
                      referential)        

        traj_pos(:,:,bin)=com
     end if
     if (mod(step,nprint)==0) write(*,*) "step = ",step
  end do

  nsteps=step

  write(*,*) "nsteps =",nsteps

  call close_read_lammps(trajfile_lammps)

!                  %----------------%
!                  | FINALIZE TOOLS |
!                  %----------------%

  do iatm=1,4
     do iatm2=iatm,4
        call finalize_gofr(gofr(iatm,iatm2),box)
     enddo
  enddo
  call finalize_msd(msdAr)
  call finalize_G(GAr,box)

!                  %------------------%
!                  | DUMP OUT OUTPUTS |
!                  %------------------%

  !Radial distribution function
  do iatm=1,4
     do iatm2=iatm,4
        open(unit=13,file="gofr-"//trim(list_labels(iatm))//"-"//trim(list_labels(iatm2))//".dat")
           call write_gofr(13,gofr(iatm,iatm2),trim(list_labels(iatm))//"-"//trim(list_labels(iatm2)))
        close(13)
     enddo
  enddo
  !Mean square displacement
  open(unit=14,file="msd.dat")
  call write_msd(14,msdAr,"Ar-Ar")
  close(14)
  !G barbar
  open(15,file="GAr.dat")
     call write_G(15,GAr,"Ar")
  close(15)
  !Density Map
  !call density_map_xy(GAr,Lx,Ly)
  !!Rot,Div,V_r,V_theta
  !call analysis_tools(GAr,Lx,Ly)

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
