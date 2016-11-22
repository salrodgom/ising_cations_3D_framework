program main
 implicit none
 integer                    :: op1,i,j,k,ii,configuration_file=456
 integer                    :: err_apertura,nop1
 integer,parameter          :: n_atoms = 60
 integer,parameter          :: n_configurations = 5
 character(len=80)          :: file_name,line
 integer, parameter         :: NOPMAX=10000
 real,dimension(NOPMAX,3,3) :: mgroup1
 real,dimension(NOPMAX,3)   :: vgroup1
 real,dimension(n_configurations,n_atoms,3)             :: cryst_coor
 character(len=4),dimension(n_configurations,n_atoms,2) :: label
!
 open(12,file="SGO")
 read(12,*)
 read(12,*)op1
 do while(op1>0)
  do i=1,3
    read (12,*)(mgroup1(op1,i,j),j=1,3),vgroup1(op1,i)
  enddo
  nop1=op1
  READ (12,*) op1
 enddo
!
 do ii=1,n_configurations
  ! c0004.gin
  !write(6,*) LEN(char(ii)) 
  if(LEN(char(ii))==1) write (file_name, '( "gulp-1-subs/c000", I1, ".gin" )' ) ii
  write(6,*)'Configuration',ii
  write(6,*)'==============='
  open(unit=configuration_file,file=file_name,status='old',iostat=err_apertura)
  if(err_apertura/=0) stop "[ERROR] reading file // file_name // "
  reading_file: do
   READ (configuration_file,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit reading_file
   !write(6,'(a)') line
   if(line(1:5)==' frac') then
    do i=1,n_atoms
     !Si   core    0.0298000    0.8771000    0.4669700
     read (configuration_file,'(A)',iostat=err_apertura) line
     if ( err_apertura /= 0 ) exit reading_file
     if ( line(6:9)=='core') then
      read(line,*)(label(ii,i,j),j=1,2),(cryst_coor(ii,i,j),j=1,3)
      write(6,'(a4,1x,a4)') (label(ii,i,j),j=1,2)
      write(6,*) (cryst_coor(ii,i,j),j=1,3)
     end if
    end do
   end if
  end do reading_file
 end do 
 stop
end program
