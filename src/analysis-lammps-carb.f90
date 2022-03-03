
program analysis

  use kinds
  use periodic_table
  use atom_methods
  use pbc
  use gofrs
  use msds
  use conductivity
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
  character(len=200) :: title
  integer :: trajunit
  type(atom),dimension(:), allocatable :: atoms
  integer :: natoms
  real(kind=dp), dimension(:,:), allocatable :: configuration,configuration_old,reference
  real(kind=dp), dimension(:),   allocatable  :: charge
  real(kind=dp), dimension(:,:), allocatable :: velocities
  real(kind=dp), dimension(:,:,:), allocatable :: traj_com
  real(kind=dp), dimension(:), allocatable :: mass
  real(kind=dp), dimension(3) :: vector
  real(kind=dp), dimension(3) :: com_box
  integer, dimension(:), allocatable :: at_mol
  logical, dimension(:), allocatable :: selectionT
  logical, dimension(:,:), allocatable :: selection_carb
  logical, dimension(:,:), allocatable :: selection_mol
  character(len=2), dimension(:), allocatable :: symbol
  logical :: usepbc=.true.  
  logical :: referential ! if T then particle frame otherwise lab frame
  integer, dimension(:), allocatable :: id_type
  integer :: bin,step,nsteps,skip,nskips,nprint,iatm,iatm2, imol
  integer :: iostat
  type(lammpsfile) :: trajfile_lammps
  integer :: i,j
  real(kind=dp), dimension(:,:), allocatable :: com
  real(kind=dp), dimension(:), allocatable :: charge_com
  real(kind=dp), dimension(:,:), allocatable :: vcm
  real(kind=dp) :: box_mass, box_charge, mol_mass

  character(len=2), dimension(5) :: list_labels
  character(len=2), dimension(4) :: list_labels_mol
  integer, dimension(5) :: count_labels
  integer, dimension(4) :: count_labels_mol

  !BOX PARAMETERS
  real(kind=dp) :: Lx,Ly,Lz
  type(box_type) :: box

  !DLPOLY PARAMETERS
  integer :: keytrj,imcon,frame,records
  character(len=200) :: orga

  !INPUT PARAMETERS
  real(kind=dp) :: limit,dr,dt
  integer :: nbins,timestep

  !TYPES
  type(gofr_type),dimension(:,:),allocatable :: gofr
  type(msd_type), dimension(:),  allocatable :: msd
  type(conductivity_type)                    :: conduct
  type(G_type), dimension(:,:),  allocatable :: G

  integer :: Nmolecule

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ directory,filename,nskips,nprint,referential
  namelist /lammps/ orga
  namelist /rdf/ limit,dr
  namelist /msd_info/ nbins,dt
 
  
  read(5,data)
  read(5,lammps)
  read(5,rdf)
  read(5,msd_info)
  !read(5,conductivity_info)


  write(*,*) "---------------------"
  write(*,*) "DATA INPUT :"
  write(*,*) "file positions =",trim(filename)
  write(*,*) "nskips =",nskips
  write(*,*) "nprint =",nprint
  write(*,*) "LAMMPS INPUT :"
  write(*,*) "orga =", orga
  write(*,*) "RDF INPUT :"
  write(*,*) "limit =",limit
  write(*,*) "dr =",dr
  write(*,*) "referential =",referential
  write(*,*) "MSD INPUT :"
  write(*,*) "nbins =",nbins
  write(*,*) "dt =",dt
  write(*,*) "---------------------"

  trajunit=10
  call open_read_lammps(trajunit,trim(directory)//trim(filename),trajfile_lammps,natoms,timestep,Lx,Ly,Lz,orga,iostat)

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
  allocate(charge(natoms))
  allocate(configuration_old(3,natoms))
  allocate(velocities(3,natoms))  
  allocate(selection_carb(5,natoms))
  allocate(selection_mol(4,Nmolecule))
  
  allocate(id_type(natoms))
  allocate(mass(natoms))
  allocate(symbol(natoms))
  allocate(com(3,Nmolecule))
  allocate(charge_com(Nmolecule))
  allocate(reference(3,Nmolecule))
  allocate(selectionT(Nmolecule))
  allocate(vcm(3,Nmolecule))
  allocate(traj_com(3,Nmolecule,nbins))
  allocate(at_mol(Nmolecule))

  allocate(gofr(5,5))
  allocate(G(5,5))
  allocate(msd(4))

!                  %----------------%
!                  |  CREATE TOOLS  |
!                  %----------------%

  selectionT= .false.
  do iatm=1,5
     do iatm2=iatm,5
        call create_gofr(gofr(iatm,iatm2),limit=limit,dr=dr)
        call create_G(G(iatm,iatm2),nbins,limit,dr)
     enddo
  enddo
  do imol=1,4
     call create_msd(msd(imol),nbins=nbins,dt=dt)
  enddo
  call create_conductivity(conduct,nbins=nbins,dt=dt)


!                  %----------------%
!                  |    LOOP  !!!   |
!                  %----------------%


  step=0
  do
  
!                  %----------------%
!                  |  READ LAMMPS   |
!                  %----------------%
         
     do skip=1,nskips
        call read_lammps(trajfile_lammps,configuration,velocities,charge,id_type,symbol,iostat)
        if (step.ne.0) then
           do iatm=1,natoms
              vector=configuration(:,iatm)-configuration_old(:,iatm)
              configuration(:,iatm)=configuration_old(:,iatm) +&
                                    minimum_image(vector,box)
           end do
        end if
        configuration_old=configuration

        if ((step.eq.0).AND.(skip.eq.1)) then
        !Assign mass of particles
        list_labels = [Character(len=2) :: 'C','O','Li','Na','K']
        list_labels_mol = [Character(len=2) :: 'C','Li','Na','K']
        imol = 0
        do iatm=1,natoms
           SELECT CASE(symbol(iatm))
           CASE('C')
                   mass(iatm) = 12
                   imol = imol + 1
                   at_mol(imol)=4
                   selection_carb(1,iatm)=.true.
                   selection_mol(1,imol)=.true.
           CASE('O')
                   mass(iatm) = 15.999
                   selection_carb(2,iatm)=.true.
           CASE('Li') 
                   mass(iatm) = 7
                   imol = imol + 1
                   at_mol(imol)=1
                   selection_carb(3,iatm)=.true.
                   selection_mol(2,imol)=.true.
           CASE('Na')
                   mass(iatm) = 22
                   imol = imol + 1
                   at_mol(imol)=1
                   selection_carb(4,iatm)=.true.
                   selection_mol(3,imol)=.true.
           CASE('K') 
                   mass(iatm) = 39
                   imol = imol + 1
                   at_mol(imol)=1
                   selection_carb(5,iatm)=.true.
                   selection_mol(4,imol)=.true.
           END SELECT
        end do

        !Check number of atoms
        do iatm=1,5
           write(*,*) "There is", count(selection_carb(iatm,:)), "atom of type", list_labels(iatm)
           count_labels(iatm)=count(selection_carb(iatm,:))
        enddo
        do imol=1,4
           count_labels_mol(imol) = count(selection_mol(imol,:))
        enddo
        endif

        !Center of mass of the box
        box_mass = 0.0 
        box_charge = 0.0 
        do iatm=1,natoms 
           com_box(:) = com_box(:) + mass(iatm)*configuration(:,iatm)
           box_mass   = box_mass + mass(iatm)
           box_charge   = box_charge + charge(iatm)
        end do
        com_box=com_box/(box_mass)
        do iatm=1,natoms
           configuration(:,iatm) = configuration(:,iatm) - com_box(:)
        end do

        if (box_charge.ne.0) then
                write(*,*) "Total box charge is:", box_charge
                stop
        endif

        !Center of mass of each monomer
        com=0.0_dp
        vcm=0.0_dp
        charge_com=0.0_dp

        j=0
        do iatm=1,Nmolecule
           mol_mass = 0.0
           do i=1,at_mol(iatm)
              j = j+1
              com(:,iatm) = com(:,iatm) + mass(j)*configuration(:,j)
              charge_com(iatm) = charge_com(iatm) + charge(j)
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
        traj_com(:,:,bin)=com
     else
        
        bin=mod(step-1,nbins)+1
        reference=traj_com(:,:,bin)

!                  %----------------%
!                  |  UPDATE TOOLS  |
!                  %----------------%

        do iatm=1,5
           do iatm2=iatm,5
              call update_gofr(gofr(iatm,iatm2),atoms,configuration,selection_carb(iatm,:),&
                            selection_carb(iatm2,:),box)
              !call update_G(G(iatm,iatm2),natoms,reference,com,&
              !        vcm,selection_carb(iatm,:),selection_carb(iatm2,:),box,usepbc,&
              !        referential)        
           enddo
        enddo
        do imol=1,4
           call update_msd(msd(imol),atoms,traj_com,bin,nbins,selection_mol(imol,:))
        enddo
        call update_conductivity(conduct,traj_com,charge_com,bin)



        traj_com(:,:,bin)=com
     end if
     if (mod(step,nprint)==0) write(*,*) "step = ",step
  end do

  nsteps=step

  write(*,*) "nsteps =",nsteps

  call close_read_lammps(trajfile_lammps)
!                  %----------------%
!                  | FINALIZE TOOLS |
!                  %----------------%

  do iatm=1,5
     do iatm2=iatm,5
        call finalize_gofr(gofr(iatm,iatm2),box)
        !call finalize_G(G(iatm,iatm2),box)
     enddo
  enddo
  do imol=1,4
     call finalize_msd(msd(imol))
  enddo
  call finalize_conductivity(conduct)

!                  %------------------%
!                  | DUMP OUT OUTPUTS |
!                  %------------------%

  do iatm=1,5
     if (count_labels(iatm).ne.0) then
     do iatm2=iatm,5
        if (count_labels(iatm2).ne.0) then
        title=trim(list_labels(iatm))//"-"//trim(list_labels(iatm2))
        open(unit=13,file="gofr-"//trim(title)//".dat")
           call write_gofr(13,gofr(iatm,iatm2),trim(title))
        close(13)
        !open(15,file="G-"//trim(title)//".dat")
        !   call write_G(15,G(iatm,iatm2),trim(title))
        !close(15)
        endif
     enddo
     endif
  enddo
  do imol=1,4
     if (count_labels_mol(imol).ne.0) then
        open(unit=14,file="msd-"//trim(list_labels_mol(imol))//".dat")
           call write_msd(14,msd(imol),trim(list_labels_mol(imol)))
        close(14)
     endif
  enddo
  open(unit=14,file="conductivity.dat")
      call write_conductivity(14, conduct, Lx*Ly*Lz,"D")
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
