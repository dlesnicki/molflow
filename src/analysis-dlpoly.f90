
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
  use xyz_files

  implicit none 

!                  %----------------%
!                  |   PARAMETERS   |
!                  %----------------%

  !GENERAL PARAMETERS
  type(atom),dimension(:), allocatable :: atoms
  integer :: natoms
  real(kind=dp), dimension(:,:), allocatable :: configuration,configuration_old,reference
  real(kind=dp), dimension(:,:), allocatable :: velocities
  real(kind=dp), dimension(:,:,:), allocatable :: traj_pos
  real(kind=dp), dimension(3) :: vector
  logical, dimension(:), allocatable :: selectionAr
  logical :: usepbc=.true.  
  logical :: referential ! if T then particle frame otherwise lab frame

  integer :: bin,step,nsteps,skip,nskips,nprint,iatm
  integer :: iostat

  !BOX PARAMETERS
  real(kind=dp), dimension(3) :: a1,a2,a3
  type(box_type) :: box

  !DLPOLY PARAMETERS
  integer :: keytrj,imcon,frame,records

  !INPUT PARAMETERS
  real(kind=dp) :: limit,dr,dt
  integer :: nbins

  !TYPES
  type(gofr_type) :: gofrAr
  type(msd_type) :: msdAr
  type(G_type) :: GAr

!                  %----------------%
!                  |  READ INPUTS   |
!                  %----------------%

  namelist /data/ nskips,nprint,a1,a2,a3,referential
  namelist /rdf/ limit,dr
  namelist /msd/ nbins,dt
  
  read(5,data)
  read(5,rdf)
  read(5,msd)

  open(unit=11,file='/home/lesnicki/DIFFUSION/src/VdW/PT/HISTORY',iostat=iostat)
  read(unit=11,fmt=*,iostat=iostat)
  read(unit=11,fmt=*,iostat=iostat)keytrj,imcon,natoms,frame,records

  write(*,*) "---------------------"
  write(*,*) "DATA INPUT :"
  write(*,*) "nskips =",nskips
  write(*,*) "nprint =",nprint
  write(*,*) "limit =",limit
  write(*,*) "dr =",dr
  write(*,*) "nbins = ",nbins
  write(*,*) "dt =",dt
  write(*,*) "box =",a1(1),a2(2),a3(3)
  write(*,*) "referential =",referential
  write(*,*) "---------------------"
  write(*,*) "HISTORY :"
  write(*,*) "natoms =",natoms
  write(*,*) "keytrj =",keytrj
  write(*,*) "frame =",frame
  write(*,*) "---------------------"

!                  %----------------%
!                  |   CREATE BOX   |
!                  %----------------%

  call create_general_box(box=box,a1=a1,a2=a2,a3=a3)

!                  %----------------%
!                  |    ALLOCATE    |
!                  %----------------%


  allocate(atoms(natoms))
  allocate(configuration(3,natoms))
  allocate(configuration_old(3,natoms))
  allocate(reference(3,natoms))
  allocate(velocities(3,natoms))  
  allocate(traj_pos(3,natoms,nbins))
  allocate(selectionAr(natoms))
  

!                  %----------------%
!                  |  CREATE TOOLS  |
!                  %----------------%

  selectionAr= select_all(atoms)
!  call create_gofr(gofrAr,limit=limit,dr=dr)
!  call create_msd(msdAr,nbins=nbins,dt=dt)
  call create_G(GAr,nbins,limit,dr)


!                  %----------------%
!                  |    LOOP  !!!   |
!                  %----------------%


  step=0
  do
  
!                  %----------------%
!                  | READ DL_POLY   |
!                  %----------------%

         
     do skip=1,nskips
        call read_DLPOLY(11,natoms,configuration,velocities,iostat)

        if (step.ne.0) then
           do iatm=1,natoms
              vector=configuration(:,iatm)-configuration_old(:,iatm)
              configuration(:,iatm)=configuration_old(:,iatm) +&
                                    minimum_image(vector,box)
           end do
        end if
        configuration_old=configuration
     end do
      
     if (iostat.ne.0) exit

     step=step+1
     if (step<=nbins) then
        bin=step
        traj_pos(:,:,bin)=configuration
        
     else
        
        bin=mod(step-1,nbins)+1
        reference=traj_pos(:,:,bin)

!                  %----------------%
!                  |  UPDATE TOOLS  |
!                  %----------------%

!        WRITE(*,*)'UPDATE GOFR'
!        call update_gofr(gofrAr,atoms,configuration,selectionAr,&
!                      selectionAr,box)

!        WRITE(*,*)'UPDATE MSD'
!        call update_msd(msdAr,atoms,traj_pos,bin,nbins,selectionAr)

        WRITE(*,*)'UPDATE G'
        call update_G(GAr,atoms,reference,configuration,&
                      velocities,selectionAr,selectionAr,box,usepbc,&
                      referential)        

        traj_pos(:,:,bin)=configuration
     end if
     if (mod(step,nprint)==0) write(*,*) step
  end do

  nsteps=step

  write(*,*) "nsteps =",nsteps

  close(11)

!                  %----------------%
!                  | FINALIZE TOOLS |
!                  %----------------%

!  call finalize_gofr(gofrAr,box)
!  call finalize_msd(msdAr)
  call finalize_G(GAr,box)

!                  %------------------%
!                  | DUMP OUT OUTPUTS |
!                  %------------------%

  !Radial distribution function
!  open(unit=13,file="gofr.dat")
!  call write_gofr(13,gofrAr,"Ar-Ar")
!  close(13)
  !Mean square displacement
!  open(unit=14,file="msd.dat")
!  call write_msd(14,msdAr,"Ar-Ar")
!  close(14)
  !G barbar
  open(15,file="GAr.dat")
     call write_G(15,GAr,"Ar")
  close(15)
  !Density Map
  call density_map_xy(GAr,a1,a2)
  !Rot,Div,V_r,V_theta
  call analysis_tools(GAr,a1,a2)

CONTAINS

  subroutine read_DLPOLY(unit,natoms,configuration,velocities,iostat)
     implicit none

     integer, intent(in) :: unit
     integer, intent(in) :: natoms
     real(kind=dp), dimension(:,:), intent(inout) :: configuration
     real(kind=dp), dimension(:,:), intent(inout) :: velocities
     integer, intent(out) :: iostat

     real(kind=dp) :: tstep,time,weight,charge,rsd 
     character*8 :: atmnam,timestep
     real(kind=dp), dimension(3) :: pos,vel
     integer :: i,iatm,nstep,megatm,keytrj,imcon
     real(kind=dp), dimension(3) :: a,b,c
     
     read(unit=unit,fmt=*,iostat=iostat)timestep,nstep,megatm,&
                   keytrj,imcon,tstep,time
     read(unit=unit,fmt=*,iostat=iostat)a(:)
     read(unit=unit,fmt=*,iostat=iostat)b(:)
     read(unit=unit,fmt=*,iostat=iostat)c(:)

     if (iostat.ne.0) return
     do i=1,natoms
       read(unit=unit,fmt=*,iostat=iostat)atmnam,iatm,weight,charge,rsd
       if (iostat.eq.0) then
          read(unit=unit,fmt=*,iostat=iostat)pos(:)
          if (keytrj>0) then
             read(unit=unit,fmt=*,iostat=iostat)vel(:)
          else 
             write(*,*) "WARNING :  NO VELOCITIES"
             STOP
          end if
          configuration(:,i)=pos(:)
          velocities(:,i)=vel(:)
       end if
     end do
     if (iostat.ne.0) stop "Problem reading DLPOLY"
  end subroutine read_DLPOLY

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
