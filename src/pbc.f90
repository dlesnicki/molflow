! *****************************************************
!>\file pbc.f90
! *****************************************************

! *****************************************************
!> \brief Defines box type and periodic boundary conditions methods
! *****************************************************

module pbc

  use kinds
  use colvar_utils

  implicit none

! *****************************************************
!>\brief Generic box type
! *****************************************************

  type box_type
     private
     real(kind=dp), dimension(3,3) ::box_direct
     real(kind=dp), dimension(3,3) ::box_inverse
  end type box_type

  private

  public :: box_type,create_cubic_box,create_general_box
  public :: apply_pbc,minimum_image
  public :: volume_of_box,get_box_cell,get_box_cellinv

contains

! ***************************************************
!>\brief To create a cubic box
!>\param box box to be created
!>\param length size of the cubic box in real space
! ***************************************************

  subroutine create_cubic_box(box,length)
    type(box_type), intent(out) :: box
    real(kind=dp), intent(in) :: length

    box%box_direct=0.0_dp
    box%box_inverse=0.0_dp

    box%box_direct(1,1)=length
    box%box_direct(2,2)=length
    box%box_direct(3,3)=length

    box%box_inverse(1,1)=1.0_dp/length
    box%box_inverse(2,2)=1.0_dp/length
    box%box_inverse(3,3)=1.0_dp/length

  end subroutine create_cubic_box

  subroutine create_general_box(box,a1,a2,a3)
    type(box_type), intent(out) :: box
    real(kind=dp), dimension(3), intent(in) :: a1,a2,a3 

    real(kind=dp), dimension(3,3) :: h_mat,h_inv
    real(kind=dp) :: det,id

    h_mat(:,1) = a1(:)
    h_mat(:,2) = a2(:)
    h_mat(:,3) = a3(:)
    
    det = a1(1)*a2(2)*a3(3) + a2(1)*a3(2)*a1(3) + a3(1)*a1(2)*a2(3) &
         - a1(1)*a3(2)*a2(3) - a2(1)*a1(2)*a3(3) - a3(1)*a2(2)*a1(3)
    id = 1.0_dp / det
    h_inv(1,1) = id * ( a2(2)*a3(3)-a2(3)*a3(2) )
    h_inv(2,1) = id * ( a1(3)*a3(2)-a1(2)*a3(3) )
    h_inv(3,1) = id * ( a1(2)*a2(3)-a1(3)*a2(2) )
    h_inv(1,2) = id * ( a2(3)*a3(1)-a2(1)*a3(3) )
    h_inv(2,2) = id * ( a1(1)*a3(3)-a1(3)*a3(1) )
    h_inv(3,2) = id * ( a1(3)*a2(1)-a1(1)*a2(3) )
    h_inv(1,3) = id * ( a2(1)*a3(2)-a2(2)*a3(1) )
    h_inv(2,3) = id * ( a1(2)*a3(1)-a1(1)*a3(2) )
    h_inv(3,3) = id * ( a1(1)*a2(2)-a1(2)*a2(1) )

    box%box_direct=h_mat
    box%box_inverse=h_inv

  end subroutine create_general_box


! ***********************************************
!>\brief applies pbc on vector
!>\param vector
!>\param box
! ***********************************************

  pure function apply_pbc(vector,box)
    real(kind=dp),dimension(3), intent(in) :: vector
    type(box_type), intent(in) :: box
    real(kind=dp), dimension(3) :: apply_pbc

    apply_pbc=vector-&
              matmul(box%box_direct,floor(matmul(box%box_inverse,vector)))


  end function apply_pbc

! ***********************************************
!>\brief minimum image of vector
!>\param vector
!>\param box
! ***********************************************

  pure function minimum_image(vector,box)
    real(kind=dp),dimension(3), intent(in) :: vector
    type(box_type), intent(in) :: box
    real(kind=dp), dimension(3) :: minimum_image

    minimum_image=vector-&
              matmul(box%box_direct,&
                     floor(matmul(box%box_inverse,vector)+0.5_dp))

  end function minimum_image

  pure function volume_of_box(box) result(vol)
    type(box_type), intent(in) :: box
    real(kind=dp) :: vol

    real(kind=dp), dimension(3,3) :: vectors
    real(kind=dp),dimension(3) :: vector_product

    vectors=box%box_direct

!!$    vector_product=0.0_dp
!!$
!!$    vector_product(1) = vectors(2,1)*vectors(3,2) &
!!$                       -vectors(3,1)*vectors(2,2)
!!$    
!!$    vector_product(2) = vectors(3,1)*vectors(1,2) &
!!$                       -vectors(1,1)*vectors(3,2)
!!$    
!!$    vector_product(3) = vectors(1,1)*vectors(2,2) &
!!$                       -vectors(2,1)*vectors(1,2)
    

    vector_product=vect_product(vectors(:,1),vectors(:,2))

    vol=dot_product(vectors(:,3),vector_product)

  end function volume_of_box

  subroutine get_box_cell(box,cell)
    type(box_type),intent(in) :: box
    real(kind=dp), dimension(3,3), intent(out) :: cell
    
    cell=box%box_direct

  end subroutine get_box_cell

  subroutine get_box_cellinv(box,cellinv)
    type(box_type),intent(in) :: box
    real(kind=dp), dimension(3,3), intent(out) :: cellinv
    
    cellinv=box%box_inverse
  end subroutine get_box_cellinv

end module pbc
