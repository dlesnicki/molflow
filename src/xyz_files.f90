
! ******************************************
!>\file xyz_files.f90
! ******************************************

! ******************************************
!>\brief tools to read and write xyz files
! ******************************************

module xyz_files

  use kinds
  use periodic_table
  use atom_methods

  implicit none

  type xyzfile
     integer :: unit,iostat=0
     character(len=80) :: filename='NULL'
     integer :: natoms=0
     type(atom), dimension(:), pointer :: atomlist => NULL()
  end type xyzfile

  public

contains

! **********************************************
!>\brief get number of atoms from an xyz file.
!>On exit get_natoms_xyz=-1 if a problem occured 
!>\param unit file unit of the xyz file
! ***********************************************

  subroutine open_read_xyz(unit,filename,file,natoms,iostat)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename
    integer, intent(inout) :: natoms
    type(xyzfile), intent(out) :: file
    integer, intent(out),optional :: iostat

    character(len=200) :: line

    if(present(iostat))file%iostat=iostat
    open(unit=unit,file=trim(filename),status='old',iostat=file%iostat)
    if(present(iostat))iostat=file%iostat

    if(file%iostat>0)then
      write(0,'(a,i3,a,a)')'Error ',file%iostat,'! Opening file: ',trim(filename)
      return
    endif

    read(unit=unit,fmt='(a)',iostat=file%iostat) line
    if(present(iostat))iostat=file%iostat

    if (file%iostat.ne.0) then
      write(0,'(a,i3,a,a)')'Error ',file%iostat,'! Reading line 1 of ',trim(filename)
      return
    endif

    read(line,*) natoms

    rewind(unit)
    file%unit=unit
    file%filename=trim(filename)
    file%natoms=natoms
    nullify(file%atomlist)

  end subroutine open_read_xyz


  subroutine close_read_xyz(file)
    type(xyzfile), intent(inout) :: file
    close(file%unit)
    file%unit=0
    file%natoms=0
    nullify(file%atomlist)
  end subroutine close_read_xyz

! ********************************************
!>\brief read an xyz file configuration
!>\param unit file to be read
!>\param iostat
!>\param natoms (can be determined from get_natoms_xyz())
!>atoms and configuration should not be present
!>\param atoms array containing atom types as read from file
!>\param configuration 
! ********************************************

  subroutine read_xyz_atoms(file,atoms,iostat)
    type(xyzfile), intent(in) :: file 
    integer, intent(out) :: iostat
    type(atom), dimension(:) :: atoms


    character(len=200) :: line
    character(len=4) :: symbol
    real(kind=dp), dimension(3) :: vector
    
    integer :: iat,natoms,unit

    unit=file%unit
    natoms=file%natoms

    read(unit=unit,fmt='(a)',iostat=iostat) line
    
    read(unit=unit,fmt='(a)',iostat=iostat) line
    
    do iat=1,natoms
       read(unit,fmt='(a)',iostat=iostat) line

       if (iostat.eq.0) then
          read(line,*) symbol, vector(:)
          
          atoms(iat)=atom_from_symbol(symbol)
       end if
    end do
      
    rewind(file%unit)
      
  end subroutine read_xyz_atoms

  subroutine read_xyz_comment(file,comment,iostat)
    type(xyzfile), intent(in) :: file 
    integer, intent(out) :: iostat
    character(len=200) :: comment

    integer :: unit
    character(len=200) :: line

    unit=file%unit

    read(unit=unit,fmt='(a)',iostat=iostat) line
    read(unit=unit,fmt='(a)',iostat=iostat) comment

    rewind(unit)

  end subroutine read_xyz_comment


  subroutine read_xyz_config(file,configuration,iostat)
    type(xyzfile), intent(in) :: file
    integer, intent(out) :: iostat
    real(kind=dp), dimension(:,:) :: configuration


    character(len=200) :: line
    character(len=2) :: symbol
    real(kind=dp), dimension(3) :: vector
    
    integer :: iat,natoms,unit

    unit=file%unit
    natoms=file%natoms

    read(unit=unit,fmt='(a)',iostat=iostat) line
    
    if (iostat.ne.0) return

    read(unit=unit,fmt='(a)',iostat=iostat) line
    
    do iat=1,natoms
       read(unit,fmt='(a)',iostat=iostat) line

       if (iostat.eq.0) then
          read(line,*) symbol, vector(:)

          configuration(:,iat)=vector
       end if

    end do
      
    if (iostat.ne.0) stop "Problem reading xyz file"
      
  end subroutine read_xyz_config

  subroutine open_write_xyz(unit,filename,file,natoms,iostat)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename
    integer, intent(in) :: natoms
    type(xyzfile), intent(out) :: file
    integer, intent(out) :: iostat

    open(unit=unit,file=trim(filename),iostat=iostat)

    file%unit=unit
    file%filename=trim(filename)
    file%natoms=natoms
    nullify(file%atomlist)

  end subroutine open_write_xyz

  subroutine close_write_xyz(file)
    type(xyzfile), intent(inout) :: file
    close(file%unit)
    file%unit=0
    file%natoms=0
    deallocate(file%atomlist)
    nullify(file%atomlist)
  end subroutine close_write_xyz

  subroutine write_xyz_atoms(file,atoms,iostat)
    type(xyzfile), intent(inout) :: file 
    integer, intent(out) :: iostat
    type(atom), dimension(:) :: atoms

    integer :: natoms

    natoms=file%natoms

    allocate(file%atomlist(natoms))
    file%atomlist=atoms
    
    iostat=0

  end subroutine write_xyz_atoms

  subroutine write_xyz_config(file,configuration,comment,iostat)
    type(xyzfile), intent(in) :: file
    integer, intent(out) :: iostat
    real(kind=dp), dimension(:,:) :: configuration
    character(len=*) :: comment

    character(len=4) :: symbol
    real(kind=dp), dimension(3) :: vector
    
    integer :: iat,natoms,unit

    unit=file%unit
    natoms=file%natoms

    write(unit=unit,fmt=*,iostat=iostat) natoms
    
    write(unit=unit,fmt='(a)',iostat=iostat) trim(comment)
    
    do iat=1,natoms

       vector=configuration(:,iat)
       symbol=file%atomlist(iat)%symbol

       write(unit=unit,fmt='(a4,3(2x,F14.6))',iostat=iostat) symbol,vector

    end do
      
    if (iostat.ne.0) stop "Problem writing xyz file"
      
  end subroutine write_xyz_config

end module xyz_files
