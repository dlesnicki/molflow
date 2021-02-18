! ******************************************************
!> \file colvar_utils.f90
! ******************************************************

! ******************************************************
!> \brief Defines colvar utils
! ******************************************************

module colvar_utils

  use kinds

  implicit none

contains

  pure function distance(vector)
    real(kind=dp), dimension(3), intent(in) :: vector
    real(kind=dp) :: distance

    distance=dot_product(vector,vector)

    distance=sqrt(distance)

  end function distance

  pure function distance_der(vector)
    real(kind=dp), dimension(3), intent(in) :: vector
    real(kind=dp), dimension(3,2) :: distance_der
    
    distance_der(:,1)=-vector/distance(vector)

    distance_der(:,2)=-distance_der(:,1)

  end function distance_der

  pure function bend(vector1,vector2)
    real(kind=dp), dimension(3), intent(in) :: vector1,vector2
    real(kind=dp) :: bend

    real(kind=dp) :: norm1,norm2,scalar,cosinus
   
    norm1=distance(vector1)
    norm2=distance(vector2)
    scalar=dot_product(vector1,vector2)

    cosinus=scalar/norm1/norm2

    bend=acos(cosinus)
    
  end function bend

  pure function bend_der(vector1,vector2)
    real(kind=dp), dimension(3), intent(in) :: vector1,vector2
    real(kind=dp), dimension(3,3) :: bend_der

    real(kind=dp) :: norm1,norm2,sinus,cosinus,theta
    real(kind=dp), dimension(3) :: u1,u2

    theta=bend(vector1,vector2)

    sinus=sin(theta)
    cosinus=cos(theta)

    norm1=distance(vector1)
    norm2=distance(vector2)

    u1=vector1/norm1
    u2=vector2/norm2

    bend_der(:,1)=(u2*norm2-u1*norm2*cosinus+u1*norm1-u2*norm1*cosinus) &
                  /norm1/norm2/sinus

    bend_der(:,2)=(u1*cosinus-u2)/norm1/sinus
    bend_der(:,3)=(u2*cosinus-u1)/norm2/sinus

  end function bend_der

  pure function vect_product(vector1,vector2)
    real(kind=dp), dimension(3), intent(in) :: vector1,vector2
    real(kind=dp), dimension(3) :: vect_product

    vect_product(1)=vector1(2)*vector2(3)-vector1(3)*vector2(2)
    vect_product(2)=vector1(3)*vector2(1)-vector1(1)*vector2(3)
    vect_product(3)=vector1(1)*vector2(2)-vector1(2)*vector2(1)


  end function vect_product
    
  pure function torsion(vector1,vector2,vector3)
    real(kind=dp), dimension(3), intent(in) :: vector1,vector2,vector3
    real(kind=dp) :: torsion

    real(kind=dp) :: norm1,norm2,norm3,phi2,phi3,tau,cosinus,projection

    real(kind=dp),dimension(3) :: u12,u23,u34,u12u23,u23u34

    norm1=distance(vector1)
    norm2=distance(vector2)
    norm3=distance(vector3)

    u12=vector1/norm1
    u23=vector2/norm2
    u34=vector3/norm3

    u12u23=vect_product(u12,u23)
    u23u34=vect_product(u23,u34)

    phi2=bend(vector1,vector2)
    phi3=bend(vector2,vector3)

    cosinus=dot_product(u12u23,u23u34)/sin(phi2)/sin(phi3)
    tau=acos(cosinus)

    projection=dot_product(u12u23,u34)

    if (projection.gt.0.0_dp) then
       torsion=tau
    else
       torsion=-tau
    end if

  end function torsion

  pure function torsion_der(vector1,vector2,vector3)
    real(kind=dp), dimension(3), intent(in) :: vector1,vector2,vector3
    real(kind=dp), dimension(3,4) :: torsion_der

    real(kind=dp) :: norm1,norm2,norm3,phi2,phi3

    real(kind=dp),dimension(3) :: u12,u23,u34,u12u23,u23u34,u43u32

    norm1=distance(vector1)
    norm2=distance(vector2)
    norm3=distance(vector3)

    u12=vector1/norm1
    u23=vector2/norm2
    u34=vector3/norm3

    u12u23=vect_product(u12,u23)
    u23u34=vect_product(u23,u34)
    u43u32=-u23u34

    phi2=bend(-vector1,vector2)
    phi3=bend(-vector2,vector3)


    torsion_der(:,1)=-u12u23/norm1/sin(phi2)**2

    torsion_der(:,2)=(norm2-norm1*cos(phi2))/ &
                       (norm1*norm2*sin(phi2))*u12u23/sin(phi2) &
                    +cos(phi3)/norm2/sin(phi3)*u43u32/sin(phi3)

    torsion_der(:,3)=(norm2-norm3*cos(phi3))/ &
                       (norm3*norm2*sin(phi3))*u43u32/sin(phi3) &
                    +cos(phi2)/norm2/sin(phi2)*u12u23/sin(phi2)


    torsion_der(:,4)=-u43u32/norm3/sin(phi3)**2

  end function torsion_der

  pure function outplane(vector1,vector2,vector3)
    real(kind=dp), dimension(3), intent(in) :: vector1,vector2,vector3
    real(kind=dp) :: outplane

    real(kind=dp) :: norm1,norm2,norm3,phi1,phi2,phi3,sinus

    real(kind=dp),dimension(3) :: u41,u42,u43,u42u43

    norm1=distance(vector1)
    norm2=distance(vector2)
    norm3=distance(vector3)

    u41=vector1/norm1
    u42=vector2/norm2
    u43=vector3/norm3

    u42u43=vect_product(u42,u43)

    phi1=bend(vector2,vector3)
    phi2=bend(vector1,vector3)
    phi3=bend(vector1,vector2)

    sinus=dot_product(u42u43,u41)/sin(phi1)

    outplane=asin(sinus)

  end function outplane

  pure function outplane_der(vector1,vector2,vector3)
    real(kind=dp), dimension(3), intent(in) :: vector1,vector2,vector3
    real(kind=dp), dimension(3,4) :: outplane_der

    real(kind=dp) :: norm1,norm2,norm3,phi1,phi2,phi3,theta,sinus

    real(kind=dp),dimension(3) :: u41,u42,u43,u42u43,u43u41,u41u42

    norm1=distance(vector1)
    norm2=distance(vector2)
    norm3=distance(vector3)

    u41=vector1/norm1
    u42=vector2/norm2
    u43=vector3/norm3

    u42u43=vect_product(u42,u43)
    u43u41=vect_product(u43,u41)
    u41u42=vect_product(u41,u42)

    phi1=bend(vector2,vector3)
    phi2=bend(vector1,vector3)
    phi3=bend(vector1,vector2)

    sinus=dot_product(u42u43,u41)/sin(phi1)

    theta=asin(sinus)

    outplane_der(:,1)=1.0_dp/norm1*(u42u43/cos(theta)/sin(phi1)-tan(theta)*u41)
    
    outplane_der(:,2)=1.0_dp/norm2*(u43u41/cos(theta)/sin(phi1)-&
                                tan(theta)/sin(phi1)**2*(u42-cos(phi1)*u43))

    outplane_der(:,3)=1.0_dp/norm3*(u41u42/cos(theta)/sin(phi1)-&
                                tan(theta)/sin(phi1)**2*(u43-cos(phi1)*u42))

    outplane_der(:,4)=-outplane_der(:,1)-outplane_der(:,2)-outplane_der(:,3)

  end function outplane_der

end module colvar_utils
