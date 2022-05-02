
! ******************************************
!>\file lammps_files.f90
! ******************************************

! ******************************************
!>\brief tools to read and write lammps files
! ******************************************

module lammps_files

  use kinds
  use periodic_table
  use atom_methods

  implicit none

  type lammpsfile
     integer :: unit,iostat=0
     character(len=80) :: filename='NULL'
     integer :: natoms=0
     integer :: timestep=0
     real(kind=dp) :: a,b,c
     character(len=200) :: orga
     type(atom), dimension(:), pointer :: atomlist => NULL()
  end type lammpsfile

  public

contains

! **********************************************
!>\brief get number of atoms from an lammps file.
!>On exit get_natoms_lammps=-1 if a problem occured 
!>\param unit file unit of the lammps file
! ***********************************************

  subroutine open_read_lammps(unit,filename,file,natoms,timestep,a,b,c,orga,iostat)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename
    integer, intent(inout) :: natoms
    integer, intent(inout) :: timestep
    real(kind=dp), intent(inout) :: a,b,c
    character(len=200), intent(inout) :: orga
    type(lammpsfile), intent(out) :: file
    integer, intent(out),optional :: iostat

    character(len=200) :: line1
    character(len=200) :: line2
    character(len=200) :: line3
    character(len=200) :: line4
    character(len=200) :: line5
    character(len=200) :: line6
    character(len=200) :: line7
    character(len=200) :: line8

    real(kind=dp) :: Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max

    if(present(iostat))file%iostat=iostat
    open(unit=unit,file=trim(filename),status='old',iostat=file%iostat)
    if(present(iostat))iostat=file%iostat

    if(file%iostat>0)then
      write(0,'(a,i3,a,a)')'Error ',file%iostat,'! Opening file: ',trim(filename)
      return
    endif

    read(unit=unit,fmt='(a)',iostat=file%iostat) line1
    read(unit=unit,fmt='(a)',iostat=file%iostat) line2
    read(unit=unit,fmt='(a)',iostat=file%iostat) line3
    read(unit=unit,fmt='(a)',iostat=file%iostat) line4
    read(unit=unit,fmt='(a)',iostat=file%iostat) line5
    read(unit=unit,fmt='(a)',iostat=file%iostat) line6
    read(unit=unit,fmt='(a)',iostat=file%iostat) line7
    read(unit=unit,fmt='(a)',iostat=file%iostat) line8
    if(present(iostat))iostat=file%iostat

    if (file%iostat.ne.0) then
      write(0,'(a,i3,a,a)')'Error ',file%iostat,'! Reading line 1 of ',trim(filename)
      return
    endif

    read(line1,*) 
    read(line2,*) timestep
    read(line3,*) 
    read(line4,*) natoms
    read(line5,*) 
    read(line6,*)Lx_min,Lx_max 
    read(line7,*)Ly_min,Ly_max
    read(line8,*)Lz_min,Lz_max

    a = Lx_max - Lx_min
    b = Ly_max - Ly_min
    c = Lz_max - Lz_min

    rewind(unit)
    file%unit=unit
    file%filename=trim(filename)
    file%natoms=natoms
    file%timestep=timestep
    file%orga = orga
    nullify(file%atomlist)

  end subroutine open_read_lammps


  subroutine close_read_lammps(file)
    type(lammpsfile), intent(inout) :: file
    close(file%unit)
    file%unit=0
    file%natoms=0
    nullify(file%atomlist)
  end subroutine close_read_lammps

! ********************************************
!>\brief read an lammps file configuration
!>\param unit file to be read
!>\param iostat
!>\param natoms (can be determined from get_natoms_xyz())
!>atoms and configuration should not be present
!>\param atoms array containing atom types as read from file
!>\param configuration 
! ********************************************

  subroutine read_lammps_atoms(file,atoms,iostat)
    type(lammpsfile), intent(in) :: file 
    integer, intent(out) :: iostat
    type(atom), dimension(:) :: atoms


    character(len=200) :: line
    character(len=200) :: line1,line2,line3,line4,line5
    character(len=200) :: line6,line7,line8,line9
    character(len=4) :: symbol
   
    integer :: id_at, id_mol, id_type 
    integer :: iat,natoms,unit

    unit=file%unit
    natoms=file%natoms

    read(unit=unit,fmt=*,iostat=iostat) line1
    read(unit=unit,fmt=*,iostat=iostat) line2
    read(unit=unit,fmt=*,iostat=iostat) line3
    read(unit=unit,fmt=*,iostat=iostat) line4
    read(unit=unit,fmt=*,iostat=iostat) line5
    read(unit=unit,fmt=*,iostat=iostat) line6
    read(unit=unit,fmt=*,iostat=iostat) line7
    read(unit=unit,fmt=*,iostat=iostat) line8
    read(unit=unit,fmt=*,iostat=iostat) line9
     
    do iat=1,natoms
       read(unit=unit,fmt='(a)',iostat=iostat) line
        
       if (iostat.eq.0) then
          read(line,*)id_at,id_mol,id_type,symbol
          atoms(iat)=atom_from_symbol(symbol)
       end if
    end do
      
    rewind(file%unit)
      
  end subroutine read_lammps_atoms

  subroutine read_lammps(file,configuration,velocities,charge,type,symbol,iostat)
    type(lammpsfile), intent(in) :: file
    integer, intent(out) :: iostat
    real(kind=dp), dimension(:,:) :: configuration
    real(kind=dp), dimension(:,:) :: velocities
    real(kind=dp), dimension(:) :: charge
    integer, dimension(:) :: type 
    character(len=2),dimension(:) :: symbol

    character(len=500) :: line
    character(len=200) :: line1,line2,line3,line4
    character(len=200) :: line5,line6,line7,line8,line9
    real(kind=dp), dimension(3) :: vector_pos,vector_vel
    real(kind=dp), dimension(3) :: vector_frc
    real(kind=dp) :: q
    integer :: id_at,id_mol,id_type 
    integer :: iat,natoms,unit
    real(kind=dp) :: mass
    character(len=2) :: symbol_at

    unit=file%unit
    natoms=file%natoms

    read(unit=unit,fmt='(a)',iostat=iostat) line1
    read(unit=unit,fmt='(a)',iostat=iostat) line2
    read(unit=unit,fmt='(a)',iostat=iostat) line3
    read(unit=unit,fmt='(a)',iostat=iostat) line4
    read(unit=unit,fmt='(a)',iostat=iostat) line5
    read(unit=unit,fmt='(a)',iostat=iostat) line6
    read(unit=unit,fmt='(a)',iostat=iostat) line7
    read(unit=unit,fmt='(a)',iostat=iostat) line8
    read(unit=unit,fmt='(a)',iostat=iostat) line9
    
    if (iostat.ne.0) return
    
    do iat=1,natoms
       read(unit,fmt='(a)',iostat=iostat) line

       if (iostat.eq.0) then
          SELECT CASE(file%orga)
          CASE("carbonate")
              read(line,*) id_at,id_mol,id_type,symbol_at,q,vector_pos(:),vector_vel(:), vector_frc(:)
          CASE("lj")
              read(line,*) id_at,id_type,symbol_at,vector_pos(:),vector_vel(:),vector_frc(:)
              q = 0
          CASE DEFAULT
              read(line,*) id_at,id_type,symbol_at,vector_pos(:),vector_vel(:),vector_frc(:)
          END SELECT
          charge(iat)=q
          configuration(:,iat)=vector_pos
          velocities(:,iat)=vector_vel
          type(iat)=id_type
          symbol(iat)=symbol_at
       end if

    end do
      
    if (iostat.ne.0) stop "Problem reading xyz file"
      
  end subroutine read_lammps

end module lammps_files
