
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
  real(kind=dp) :: sum_RG_xx,sum_RG_yy,sum_RG_zz,sum_RG_xy,sum_RG_xz,sum_RG_yz
  integer :: Nmonomer,Nmolecule

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
  allocate(RG_xy(Nmolecule))
  allocate(RG_xz(Nmolecule))
  allocate(RG_yz(Nmolecule))

!                  %----------------%
!                  |    LOOP  !!!   |
!                  %----------------%

  open(unit=88,file='gyration.dat')
  open(unit=55,file='RG2.dat')
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
        sum_RG_xx = 0.0_dp
        sum_RG_yy = 0.0_dp
        sum_RG_zz = 0.0_dp
        sum_RG_xy = 0.0_dp
        sum_RG_xz = 0.0_dp
        sum_RG_yz = 0.0_dp

        j=0
        do iatm=1,Nmolecule

           RG_xx(iatm)=0.0_dp
           RG_yy(iatm)=0.0_dp
           RG_zz(iatm)=0.0_dp
           RG_xy(iatm)=0.0_dp
           RG_xz(iatm)=0.0_dp
           RG_yz(iatm)=0.0_dp
           do i=1,Nmonomer
              j = j+1
           
              configuration(:,j) = configuration(:,j) - com(:,iatm)    
              RG_xx(iatm) = RG_xx(iatm) + mass(j)*(configuration(1,j))**2
              RG_yy(iatm) = RG_yy(iatm) + mass(j)*(configuration(2,j))**2
              RG_zz(iatm) = RG_zz(iatm) + mass(j)*(configuration(3,j))**2
              RG_xy(iatm) = RG_xy(iatm) + mass(j)*(configuration(1,j)*configuration(2,j))
              RG_xz(iatm) = RG_xz(iatm) + mass(j)*(configuration(1,j)*configuration(3,j))
              RG_yz(iatm) = RG_yz(iatm) + mass(j)*(configuration(2,j)*configuration(3,j))
           end do
           RG_xx(iatm)=RG_xx(iatm)/nmass
           RG_yy(iatm)=RG_yy(iatm)/nmass
           RG_zz(iatm)=RG_zz(iatm)/nmass
           RG_xy(iatm)=RG_xy(iatm)/nmass
           RG_xz(iatm)=RG_xz(iatm)/nmass
           RG_yz(iatm)=RG_yz(iatm)/nmass
           write(55,*)iatm,RG_xx(iatm),RG_yy(iatm),RG_zz(iatm),RG_xy(iatm),RG_xz(iatm),RG_yz(iatm)

           sum_RG_xx = sum_RG_xx + RG_xx(iatm)
           sum_RG_yy = sum_RG_yy + RG_yy(iatm)
           sum_RG_zz = sum_RG_zz + RG_zz(iatm)
           sum_RG_xy = sum_RG_xy + RG_xy(iatm)
           sum_RG_xz = sum_RG_xz + RG_xz(iatm)
           sum_RG_yz = sum_RG_yz + RG_yz(iatm)
        end do
        sum_RG_xx = sum_RG_xx/Nmolecule
        sum_RG_yy = sum_RG_yy/Nmolecule
        sum_RG_zz = sum_RG_zz/Nmolecule
        sum_RG_xy = sum_RG_xy/Nmolecule
        sum_RG_xz = sum_RG_xz/Nmolecule
        sum_RG_yz = sum_RG_yz/Nmolecule
     end do
      
     if (iostat.ne.0) exit

     step=step+1
!                  %----------------%
!                  |  UPDATE TOOLS  |
!                  %----------------%

     write(88,*)step,sqrt(sum_RG_xx+sum_RG_yy+sum_RG_zz)/3,sum_RG_xx,sum_RG_yy,sum_RG_zz,sum_RG_xy,sum_RG_xz,sum_RG_yz     

     if (mod(step,nprint)==0) write(*,*) step
  end do

  nsteps=step

  write(*,*) "nsteps =",nsteps
  close(55)
  call close_read_lammps(trajfile_lammps)

!                  %----------------%
!                  | FINALIZE TOOLS |
!                  %----------------%

  

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
