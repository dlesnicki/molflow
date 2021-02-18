

module selections

  use kinds
  use periodic_table
  use atom_methods
  use pbc
  use colvar_utils

  implicit none

contains

  function select_all(atoms) result(selection)
    type(atom), dimension(:), intent(in) :: atoms
    
    logical, dimension(size(atoms)) :: selection
    selection=.true.

  end function select_all

  function select_none(atoms) result(selection)
    type(atom), dimension(:), intent(in) :: atoms
    
    logical, dimension(size(atoms)) :: selection
    selection=.false.

  end function select_none

  function select_index(index,atoms,op) result(selection)
    integer, intent(in) :: index
    type(atom), dimension(:), intent(in) :: atoms
    character(len=*), intent(in),optional :: op
   
    logical, dimension(size(atoms)) :: selection

    character(len=2) :: oper
    integer natoms

    oper="eq"
    if(present(op)) then
       oper=op
       oper(len(op)+1:2)=" "
    end if

    natoms=size(atoms)
    selection=.false.

    select case(trim(oper))
       
    case("eq","==")
       selection(index)=.true.
       
    case("lt","<")
       selection(1:index-1)=.true.
       
    case("le","<=")
       selection(1:index)=.true.
       
    case("gt",">")
       selection(index+1:natoms)=.true.
       
    case("ge",">=")
       selection(index:natoms)=.true.
       
    case default
       stop "Operator not defined for selection by index"
       
    end select

  end function select_index

  function select_indexes(indexes,atoms) result(selection)
    integer, dimension(:),intent(in) :: indexes
    type(atom), dimension(:), intent(in) :: atoms
    
    logical, dimension(size(atoms)) :: selection

    integer :: i,index

    selection=.false.

    do i=1,size(indexes)
       index=indexes(i)
       selection(index)=.true.
    end do

  end function select_indexes

  function select_atom(atoms, &
                       symbol, &
                       name, &
                       element, &
                       number, &
                       amass, &
                       mass, &
                       covalent_radius, &
                       vdw_radius, &
                       op) &
                                           result(selection)

    type(atom), dimension(:), intent(in) :: atoms

    CHARACTER ( LEN = * ), intent(in), optional :: symbol
    CHARACTER ( LEN = * ), intent(in), optional :: name
    CHARACTER ( LEN = * ), intent(in), optional :: element
    INTEGER, intent(in), optional :: number
    REAL (KIND=dp), intent(in), optional :: amass
    REAL (KIND=dp), intent(in), optional :: mass 
    REAL (KIND=dp), intent(in), optional :: covalent_radius
    REAL (KIND=dp), intent(in), optional :: vdw_radius
    
    character(len=*), intent(in), optional :: op
    
    logical, dimension(size(atoms)) :: selection


    character(len=2) :: oper
    integer :: natoms,iatom

    oper="eq"
    if(present(op)) then
       oper=op
       oper(len(op)+1:2)=" "
    end if


    selection=.false.

    natoms=size(atoms)


    if (present(symbol)) then
       do iatom=1,natoms
          if(trim(atoms(iatom)%symbol).eq.trim(symbol)) selection(iatom)=.true.
       end do

    else if (present(name)) then
       do iatom=1,natoms
          if(trim(atoms(iatom)%name).eq.trim(name)) selection(iatom)=.true.
       end do

    else if (present(element)) then
       do iatom=1,natoms
          if(trim(atoms(iatom)%element).eq.trim(element)) &
               selection(iatom)=.true.
       end do

    else if (present(number)) then
       select case(trim(oper))
          
       case("eq","==")
          do iatom=1,natoms
             if(atoms(iatom)%number==number) selection(iatom)=.true.
          end do
          
       case("lt","<")
          do iatom=1,natoms
             if(atoms(iatom)%number<number) selection(iatom)=.true.
          end do
          
       case("le","<=")
          do iatom=1,natoms
             if(atoms(iatom)%number<=number) selection(iatom)=.true.
          end do
       
       case("gt",">") 
          do iatom=1,natoms
             if(atoms(iatom)%number>number) selection(iatom)=.true.
          end do
      
       case("ge",">=")
          do iatom=1,natoms
             if(atoms(iatom)%number>=number) selection(iatom)=.true.
          end do
       
       case default
          stop "Operator not defined for selection by index"
       
       end select
    
    else if (present(amass)) then
       select case(trim(oper))
          
       case("eq","==")
          do iatom=1,natoms
             if(atoms(iatom)%amass==amass) selection(iatom)=.true.
          end do
          
       case("lt","<")
          do iatom=1,natoms
             if(atoms(iatom)%amass<amass) selection(iatom)=.true.
          end do
          
       case("le","<=")
          do iatom=1,natoms
             if(atoms(iatom)%amass<=amass) selection(iatom)=.true.
          end do
       
       case("gt",">") 
          do iatom=1,natoms
             if(atoms(iatom)%amass>amass) selection(iatom)=.true.
          end do
      
       case("ge",">=")
          do iatom=1,natoms
             if(atoms(iatom)%amass>=amass) selection(iatom)=.true.
          end do
       
       case default
          stop "Operator not defined for selection by index"
       
       end select
    

    else if (present(mass)) then
       select case(trim(oper))
          
       case("eq","==")
          do iatom=1,natoms
             if(atoms(iatom)%mass==mass) selection(iatom)=.true.
          end do
          
       case("lt","<")
          do iatom=1,natoms
             if(atoms(iatom)%mass<mass) selection(iatom)=.true.
          end do
          
       case("le","<=")
          do iatom=1,natoms
             if(atoms(iatom)%mass<=mass) selection(iatom)=.true.
          end do
       
       case("gt",">") 
          do iatom=1,natoms
             if(atoms(iatom)%mass>mass) selection(iatom)=.true.
          end do
      
       case("ge",">=")
          do iatom=1,natoms
             if(atoms(iatom)%mass>=mass) selection(iatom)=.true.
          end do
       
       case default
          stop "Operator not defined for selection by index"
       
       end select
    

    else if (present(covalent_radius)) then
       select case(trim(oper))
          
       case("eq","==")
          do iatom=1,natoms
             if(atoms(iatom)%covalent_radius==covalent_radius) &
                  selection(iatom)=.true.
          end do
          
       case("lt","<")
          do iatom=1,natoms
             if(atoms(iatom)%covalent_radius<covalent_radius) &
                  selection(iatom)=.true.
          end do
          
       case("le","<=")
          do iatom=1,natoms
             if(atoms(iatom)%covalent_radius<=covalent_radius) &
                  selection(iatom)=.true.
          end do
       
       case("gt",">") 
          do iatom=1,natoms
             if(atoms(iatom)%covalent_radius>covalent_radius) &
                  selection(iatom)=.true.
          end do
      
       case("ge",">=")
          do iatom=1,natoms
             if(atoms(iatom)%covalent_radius>=covalent_radius) &
                  selection(iatom)=.true.
          end do
       
       case default
          stop "Operator not defined for selection by index"
       
       end select
    

    else if (present(vdw_radius)) then
       select case(trim(oper))
          
       case("eq","==")
          do iatom=1,natoms
             if(atoms(iatom)%vdw_radius==vdw_radius) selection(iatom)=.true.
          end do
          
       case("lt","<")
          do iatom=1,natoms
             if(atoms(iatom)%vdw_radius<vdw_radius) selection(iatom)=.true.
          end do
          
       case("le","<=")
          do iatom=1,natoms
             if(atoms(iatom)%vdw_radius<=vdw_radius) selection(iatom)=.true.
          end do
       
       case("gt",">") 
          do iatom=1,natoms
             if(atoms(iatom)%vdw_radius>vdw_radius) selection(iatom)=.true.
          end do
      
       case("ge",">=")
          do iatom=1,natoms
             if(atoms(iatom)%vdw_radius>=vdw_radius) selection(iatom)=.true.
          end do
       
       case default
          stop "Operator not defined for selection by index"
       
       end select
    

    end if


  end function select_atom



  function select_symbol(symbol,atoms) result(selection)
    character(len=*), intent(in) :: symbol
    type(atom), dimension(:), intent(in) :: atoms
    
    logical, dimension(size(atoms)) :: selection

    integer :: natoms,iatom

    selection=.false.

    natoms=size(atoms)
    do iatom=1,natoms
       if(trim(atoms(iatom)%symbol).eq.trim(symbol)) selection(iatom)=.true.
    end do

  end function select_symbol

  function select_within(cutoff,selection1,atoms,configuration,box) &
    result(selection)

    real(kind=dp), intent(in) :: cutoff
    logical, dimension(:),intent(in) :: selection1
    type(atom), dimension(:), intent(in) :: atoms
    real(kind=dp), dimension(:,:), intent(in) :: configuration
    type(box_type), intent(in), optional :: box
    
    logical, dimension(size(atoms)) :: selection

    integer ::natoms,iatom,iatom1
    real(kind=dp),dimension(3) :: vector

    selection=.false.

    natoms=size(atoms)
    do iatom=1,natoms
       do iatom1=1,natoms
          if (selection1(iatom1)) then
             vector=configuration(:,iatom1)-configuration(:,iatom)
             if (present(box)) vector=minimum_image(vector,box)
             if (distance(vector)<=cutoff) selection(iatom)=.true.
          end if
       end do
    end do

  end function select_within

  function select_within_box(sup,inf,selection1,atoms,configuration) &
    result(selection)

    real(kind=dp), intent(in) :: sup
    real(kind=dp), intent(in) :: inf
    logical, dimension(:),intent(in) :: selection1
    type(atom), dimension(:), intent(in) :: atoms
    real(kind=dp), dimension(:,:), intent(in) :: configuration

    logical, dimension(size(atoms)) :: selection

    integer ::natoms,iatom
    
    selection=.false.    

    natoms=size(atoms)
    
    do iatom=1,natoms
          if (selection1(iatom)) then
             if ((configuration(1,iatom)<=sup) .and. &
                 & (configuration(1,iatom)>=inf) .and. &
                 & (configuration(2,iatom)<=sup) .and. &
                 & (configuration(2,iatom)>=inf)) then  
                   selection(iatom)=.true.
             end if
          end if
    end do

  end function select_within_box

  function select_within_box_bis(sup_x,sup_y,inf_x,inf_y,selection1,atoms,configuration) &
    result(selection)

    real(kind=dp), intent(in) :: sup_x
    real(kind=dp), intent(in) :: inf_x
    real(kind=dp), intent(in) :: sup_y
    real(kind=dp), intent(in) :: inf_y
    logical, dimension(:),intent(in) :: selection1
    type(atom), dimension(:), intent(in) :: atoms
    real(kind=dp), dimension(:,:), intent(in) :: configuration

    logical, dimension(size(atoms)) :: selection

    integer ::natoms,iatom

    selection=.false.

    natoms=size(atoms)

    do iatom=1,natoms
          if (selection1(iatom)) then
             if ((configuration(1,iatom)<=sup_x) .and. &
                 & (configuration(1,iatom)>=inf_x) .and. &
                 & (configuration(2,iatom)<=sup_y) .and. &
                 & (configuration(2,iatom)>=inf_y)) then
                   selection(iatom)=.true.
             end if
          end if
    end do

  end function select_within_box_bis


  function select_closest(n,selection1,selection2,atoms,configuration,box) & 
       result(selection)

    integer, intent(in) :: n
    logical, dimension(:),intent(in) :: selection1,selection2
    type(atom), dimension(:), intent(in) :: atoms
    real(kind=dp), dimension(:,:), intent(in) :: configuration
    type(box_type), intent(in), optional :: box
    
    logical, dimension(size(atoms)) :: selection

    integer ::natoms,iatom,iatom1,imax,i
    real(kind=dp),dimension(3) :: vector

    integer, dimension(:), allocatable :: indexes
    real(kind=dp), dimension(:), allocatable :: distances

    real(kind=dp) :: mindist,dmax

    selection=.false.

    allocate(indexes(n))
    allocate(distances(n))
    indexes=0
    distances=huge(1.0_dp)

    natoms=size(atoms)

    do iatom=1,natoms
       if (selection1(iatom).and.(.not.selection2(iatom))) then
          mindist=huge(1.0_dp)
          do iatom1=1,natoms
             if (selection2(iatom1)) then
                vector=configuration(:,iatom1)-configuration(:,iatom)
                if (present(box)) vector=minimum_image(vector,box)
                if (distance(vector)<=mindist) mindist=distance(vector)
             end if
          end do
          imax=maxloc(distances,dim=1)
          dmax=maxval(distances,dim=1)
          if (mindist<dmax) then
             indexes(imax)=iatom
             distances(imax)=mindist
          end if
       end if
    end do

    do i=1,n
       selection(indexes(i))=.true.
    end do

    deallocate(indexes)
    deallocate(distances)

  end function select_closest


  function select_coordinate(axe,value,atoms,configuration,op) &
       result(selection)

    character(len=1), intent(in) :: axe
    real(kind=dp), intent(in) :: value
    type(atom), dimension(:), intent(in) :: atoms
    real(kind=dp), dimension(:,:), intent(in) :: configuration
    character(len=*), intent(in), optional :: op
    
    logical, dimension(size(atoms)) :: selection

    character(len=2) :: oper
    integer natoms,iaxe,iat

    oper=">="
    if(present(op)) then
       oper=op
       oper(len(op)+1:2)=" "
    end if

    natoms=size(atoms)
    selection=.false.

    select case(axe)

    case("x","X")
       iaxe=1
       
    case("y","Y")
       iaxe=2
       
    case("z","Z")
       iaxe=3
       
    case default
       stop "axe not defined"
       
    end select
    
    select case(trim(oper))
       
    case("eq","==")
       do iat=1,natoms
          if(configuration(iaxe,iat)==value) selection(iat)=.true.
       end do
       
    case("lt","<")
       do iat=1,natoms
          if(configuration(iaxe,iat)<value) selection(iat)=.true.
       end do
       
    case("le","<=")
       do iat=1,natoms
          if(configuration(iaxe,iat)<=value) selection(iat)=.true.
       end do
       
    case("gt",">") 
       do iat=1,natoms
          if(configuration(iaxe,iat)>value) selection(iat)=.true.
       end do
      
    case("ge",">=")
       do iat=1,natoms
          if(configuration(iaxe,iat)>=value) selection(iat)=.true.
       end do
       
    case default
       stop "Operator not defined for selection by index"
       
    end select


  end function select_coordinate

  function count_selection(selection) result(count)
    logical,dimension(:), intent(in) :: selection
    integer :: count

    integer :: index

    count=0
    do index=1,size(selection)
       if (selection(index)) count=count+1
    end do

  end function count_selection

  subroutine selection_to_indexes(selection,indexes)
    logical, dimension(:), intent(in) :: selection
    integer, dimension(:), intent(inout) :: indexes

    integer index,count

    if (size(indexes).lt.count_selection(selection)) then
       stop "size of indexes too small" 
    end if
    
    count=0
    do index=1,size(selection)
       if (selection(index)) then
          count=count+1
          indexes(count)=index
       end if
    end do

  end subroutine selection_to_indexes

end module selections
