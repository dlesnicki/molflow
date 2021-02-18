
module distribution

  use kinds
  use periodic_table
  use selections
  use colvar_utils
  use constants

  implicit none 

  type G_type

     real(kind=dp) :: tmax, dt
     integer :: nbins
     real(kind=dp) :: limit,dr
     real(kind=dp), dimension(:), pointer :: gofr, density, integral
     real(kind=dp), dimension(:), pointer :: vel_r, vel_theta
     real(kind=dp) :: V1
     integer :: nG
     real(kind=dp) :: nsel1,nsel2
     integer :: nupdates=-1

  end type G_type

contains

  subroutine create_G(G,nbins,limit,dr)
    type(G_type), intent(out) :: G
    real(kind=dp), intent(in), optional :: limit
    real(kind=dp), intent(in), optional :: dr
    integer, intent(in) :: nbins

    G%nbins=nbins
    G%limit=limit
    G%dr=dr
    G%nsel1=0.0_dp
    G%nsel2=0.0_dp
    G%nupdates=0
    G%nG=floor(limit/dr)
    G%nsel1=0.0_dp
    G%nsel2=0.0_dp
    G%V1=0_dp

    allocate(G%gofr(G%nG))
    G%gofr=0.0_dp
    allocate(G%density(G%nG))
    G%density=0.0_dp
    allocate(G%integral(G%nG))
    G%integral=0.0_dp
    allocate(G%vel_r(G%nG))
    G%vel_r=0.0_dp
    allocate(G%vel_theta(G%nG))
    G%vel_theta=0.0_dp

  end subroutine create_G

  subroutine update_G(G,Nmolecule,reference,configuration,&
                      velocities,selection1,selection2,box,usepbc,&
                      referential)

    implicit none
    type(G_type), intent(inout) :: G
    !type(atom), dimension(:), intent(in) :: atoms
    integer, intent(in) :: Nmolecule
    integer :: bin,ntraj,natoms
    real(kind=dp), dimension(:,:), intent(in) :: configuration,reference
    real(kind=dp), dimension(:,:), intent(in) :: velocities
    logical, dimension(:), intent(in) :: selection1,selection2
    type(box_type), intent(in) :: box
    logical, intent(in) :: usepbc
    logical, intent(in) :: referential

    real(kind=dp) :: dist    
    integer :: iat1,iat2
    real(kind=dp), dimension(3) :: U,V,I,U_ortho,V_ortho
    real(kind=dp), dimension(3) :: vector_r,vector_theta
    real(kind=dp) :: produit, dr
    real(kind=dp) :: projRI,projVI
    integer :: index       
    

    dr=G%dr
    U=0_dp
    V=0_dp
    I=0_dp

    U_ortho=0_dp
    V_ortho=0_dp

    projRI=0_dp
    projVI=0_dp

      natoms=Nmolecule!size(atoms) 
      do iat1=1,natoms
        if (selection1(iat1)) then
          U = configuration(:,iat1)-reference(:,iat1)
          G%V1=G%V1+dot_product(U,velocities(:,iat1))
    
          do iat2=1,natoms
            if (selection2(iat2).and.(iat1.ne.iat2)) then

             !Construction Ur
             vector_r = configuration(:,iat2)-configuration(:,iat1)
             if (usepbc) vector_r=minimum_image(vector_r,box)
             dist=distance(vector_r)

             !Projection sur Ur
             I = vector_r/dist
             projRI = dot_product(U,I)
             
             IF (referential) THEN
             V = velocities(:,iat2)-velocities(:,iat1)
             ELSE
             V = velocities(:,iat2)
             END IF
             projVI = dot_product(V,I)

             U_ortho = U - projRI*I
             V_ortho = V - projVI*I

             index=floor(dist/dr)+1 
             if (index<=G%nG) then
                G%gofr(index) = G%gofr(index) + 1.0_dp
                G%density(index) = G%density(index) + projRI
                G%vel_r(index) = G%vel_r(index) + projVI*projRI
                G%vel_theta(index) = G%vel_theta(index) - dot_product(U_ortho,V_ortho)/2_dp
             end if

            end if
          end do
        end if
      end do
    G%nupdates=G%nupdates+1

    G%nsel1=G%nsel1+dble(count_selection(selection1))
    G%nsel2=G%nsel2+dble(count_selection(selection2))   

  end subroutine update_G

  subroutine finalize_G(G,box)

    implicit none
    type(G_type), intent(inout) :: G
    type(box_type), intent(in) :: box    

    integer :: index
    real(kind=dp) :: volume,r,fact

    G%gofr=G%gofr/dble(G%nupdates)    
    G%density=G%density/dble(G%nupdates)
    G%nsel1=G%nsel1/dble(G%nupdates)
    G%nsel2=G%nsel2/dble(G%nupdates)

    G%vel_r=G%vel_r/dble(G%nupdates)
    G%vel_theta=G%vel_theta/dble(G%nupdates)
    G%V1=G%V1/dble(G%nupdates)/G%nsel1
    volume=volume_of_box(box)
    
    do index=1,G%nG
       G%integral(index)=sum(G%density(1:index))
    end do
    G%integral=G%integral/dble(G%nsel1)
    do index=1,G%nG
       r=(dble(index)-0.5_dp)*G%dr
       fact=volume/G%nsel1/G%nsel2/(4.0_dp*pi*r*r*G%dr)

       G%density(index)=G%density(index)*fact
       G%vel_r(index)=G%vel_r(index)*fact
       G%vel_theta(index)=G%vel_theta(index)*fact
       G%gofr(index)= G%gofr(index)*fact
    end do

  end subroutine finalize_G

  subroutine write_G(out_unit,G,comment)

    implicit none
    integer, intent(in) :: out_unit
    type(G_type), intent(in) :: G
    character(len=*), intent(in), optional :: comment

    integer :: index
    real(kind=dp) :: r
 
    if (present(comment)) write(out_unit,'("# ",A)') trim(comment)
      write(out_unit,&
       '("# dr=",F14.6,",limit=",F14.6,", nsel1=",F14.6,", nsel2=",F14.6,&
            "V1=",F14.6)') &
       G%dr,G%limit,G%nsel1,G%nsel2,G%V1
    write(*,*)G%V1
  
    do index=1,G%nG
        r=(dble(index)-0.5_dp)*G%dr
        write(out_unit,'(6F14.6)') r,G%density(index),G%vel_r(index),&
          G%vel_theta(index),G%gofr(index),G%integral(index)
    end do

  end subroutine write_G

  subroutine destroy_G(G)

  implicit none
  type(G_type), intent(inout) :: G

  deallocate(G%gofr)
  deallocate(G%density)
  deallocate(G%integral)
  deallocate(G%vel_r)
  deallocate(G%vel_theta)

  end subroutine destroy_G


  subroutine analysis_tools(G,a1,a2)

  implicit none
  type(G_type), intent(in) :: G
  real(kind=dp), dimension(:), allocatable :: V_r, V_theta, rVtheta
  real(kind=dp) :: D,Rot
  integer :: index,i,j,k,L,M,N
  real(kind=dp) :: r,x,y
  real(kind=dp),dimension(3),intent(in) :: a1,a2 



  N=G%nG
  allocate(V_r(N))
  allocate(V_theta(N))
  allocate(rVtheta(N))

  OPEN(UNIT=55,FILE='Div.dat',STATUS='Unknown')
  !OPEN(UNIT=56,FILE='Velocities.dat',STATUS='Unknown')
!  OPEN(UNIT=57,FILE='Rot.dat',STATUS='Unknown')


  WRITE(*,*) 'Velocity calculations'
  DO index=1,G%nG
     r=(dble(index)-0.5_dp)*G%dr

     IF (G%gofr(index).NE.0.0_dp) THEN

         V_r(index)=G%vel_r(index)!/G%gofr(index)
         V_theta(index)=G%vel_theta(index)!/G%gofr(index)
         rVtheta(index)=G%vel_theta(index)*r!/G%gofr(index)
   !      WRITE(56,'(3F14.6)') r,V_r(index),V_theta(index)
     END IF

  END DO
  !CLOSE(56)
  WRITE(*,*) 'Div & Rot'
  DO index=2,G%nG-1
     r=(dble(index)-0.5_dp)*G%dr

     IF (G%gofr(index).NE.0.0_dp) THEN

         j=index+1
         k=index-1
         D=(V_r(j)-V_r(k))/2.0_dp/G%dr+2.0_dp*(V_r(index)+V_theta(index))/dble(r)
         Rot= (rVtheta(j)-rVtheta(k))/2.0_dp/G%dr/dble(r)+V_r(index)/dble(r)
         WRITE(55,'(2F14.6)')r,D
!         WRITE(57,'(2F14.6)')r,Rot
     END IF
  END DO
  CLOSE(55)
!  CLOSE(56)
!  CLOSE(57)

  WRITE(*,*) 'Flow Lines calculation'
  OPEN(58,file='Flow_lines.dat',status='unknown')
  L=int(a1(1))
  M=int(a2(2))

  DO i=-L,L
    DO j=-M,M
      x=dble(i)
      y=dble(j)
      r=int(sqrt(x**2+y**2))
      index=int(r/G%dr+0.5_dp)
      WRITE(58,*)x,y,- V_r(index)*index*y/x
    END DO
    WRITE(58,*)
  END DO
  CLOSE(58)

  deallocate(V_r)
  deallocate(V_theta)
  deallocate(rVtheta)
  end subroutine analysis_tools


end module distribution
