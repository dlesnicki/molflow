

! *************************************************
!> \file atom_methods.f90
! *************************************************

! *************************************************
!> \brief methods for manipulating atoms
! *************************************************

module atom_methods

  use kinds
  use periodic_table

  implicit none

  private

  public :: atom_from_symbol, atom_from_element

contains

! *********************************************
!>\brief Get atom of periodic table from symbol
!>\param symbol
! *********************************************

  function atom_from_symbol(symbol) result(atm)
    character(len=*), intent(in) :: symbol
    type(atom) :: atm

    character(len=2) :: match
    integer :: ielem,it

    call init_periodic_table

    match=""
    it=0
    do ielem=0,nelem
       if (index(trim(symbol),trim(ptable(ielem)%symbol)).eq.1) then
          if ((match.eq."").or.&
               (index(trim(ptable(ielem)%symbol),trim(match)).ne.0)) then
             it=ielem
             match=ptable(ielem)%symbol
          end if
       end if
    end do

    atm=atom_from_element(symbol,ptable(it))

  end function atom_from_symbol


  function atom_from_element(symbol,elmnt) result(atm)
    character(len=*), intent(in) :: symbol
    type(element) :: elmnt
    type(atom) :: atm

    atm%symbol=trim(symbol)
    atm%element=elmnt%symbol
    atm%name=elmnt%name
    atm%number=elmnt%number
    atm%amass=elmnt%amass
    atm%mass=elmnt%mass
    atm%covalent_radius=elmnt%covalent_radius
    atm%vdw_radius=elmnt%vdw_radius
    atm%e_conv=elmnt%e_conv

  end function atom_from_element

end module atom_methods
