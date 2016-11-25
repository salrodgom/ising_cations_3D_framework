program main
 implicit none
 integer                    :: op1,ii,deg,configuration_file=456
 integer                    :: i,j,k,l,m,n,jj,kk,ll,mm,nn
 integer                    :: err_apertura,nop1
 integer,parameter          :: n_atoms =   60 ! 60 T + 120 Os + 120 Oc + 12 F
 integer,parameter          :: n_T_atoms = 60
 integer,parameter          :: n_Ge = 24
 integer,parameter          :: MC_steps = 100
 integer,parameter          :: n_configurations = 0
 character(len=80)          :: file_name,line
 real,parameter             :: temperature = 100.0
 integer, parameter         :: NOPMAX=10000
 integer                    :: delta_1(n_atoms,n_atoms),k_max_1,k_max_2
 real                       :: ener_0 = -7745.86721305,epsilon_
 real           :: ener_1(n_atoms) = 0.0
 real           :: deg_1(n_atoms) = 0.0
 real           :: ener_2(n_atoms,n_atoms) = 0.0
 real           :: deg_2(n_atoms,n_atoms) = 0.0
 real           :: cell_0(1:6)
 integer                    :: delta_2(n_atoms,n_atoms,n_atoms,n_atoms)
 real,dimension(NOPMAX,3,3) :: mgroup1
 real,dimension(NOPMAX,3)   :: vgroup1
 real,dimension(0:n_configurations,n_atoms,3)             :: cryst_coor
 character(len=4),dimension(0:n_configurations,n_atoms,2) :: label
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
 do ii=0,n_configurations
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
   if(line(1:5)=='cell ') then
    write(6,*)'Cell parameters'
    read(configuration_file,'(A)',iostat=err_apertura) line
    if ( err_apertura /= 0 ) exit reading_file
    read(line,*) ( cell_0(j),j=1,6 )
    write(6,*) ( cell_0(j),j=1,6 )
   end if
   if(line(1:5)==' frac') then
    write(6,*)'Atom coordinates'
    do i=1,n_atoms
     !Si   core    0.0298000    0.8771000    0.4669700
     read (configuration_file,'(A)',iostat=err_apertura) line
     if ( err_apertura /= 0 ) exit reading_file
     !if ( line(6:9)=='core') then
     if(line(1:7)=='species') exit reading_file
      read(line,*)(label(ii,i,j),j=1,2),(cryst_coor(ii,i,j),j=1,3)
      write(6,'(a4,1x,a4)') (label(ii,i,j),j=1,2)
      write(6,*) (cryst_coor(ii,i,j),j=1,3)
     !end if
    end do
   end if
  end do reading_file
  close(configuration_file)
 end do 
!
 write(6,*)'Read relation between unique and dependent configurations'
! {{{ 
 write(6,*)'One substitution:'
 delta_1(1:n_atoms,1:n_atoms) = 0
 deg_1(1:n_atoms) = 0.0
 k_max_1 = 0
 ener_1(1:n_atoms)  = 0.0
 open(unit=111,file='gulp-1-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=112,file='gulp-1-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=113,file='gulp-1-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"
 read_matrix_1: do
  READ(111,'(A)',IOSTAT=err_apertura) line
  IF( err_apertura /= 0 ) exit read_matrix_1
  read(line,*)k,deg,l
  delta_1(l,l)=1
  deg_1(l) = real(deg)
  READ(113,*)epsilon_
  ener_1(l)=epsilon_ - ener_0
  write(6,*)l,l,ener_1(l),deg_1(l)
  do j=2,deg
   READ (112,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_1
   read(line,*)k,m,i,n
   delta_1(l,i)=1
   delta_1(i,l)=1
   write(6,*)l,i,ener_1(i),deg_1(i)
  end do
 end do read_matrix_1
! final 
 close(111)
 close(112)
 close(113)
! }}}
 write(6,*)'Two substitutions:'
 delta_2(1:n_atoms,1:n_atoms,1:n_atoms,1:n_atoms) = 0
 deg_2(1:n_atoms,1:n_atoms) = 0.0
 deg_2 = 0.0
 k_max_2 = 0
 ener_2(1:n_atoms,1:n_atoms) = 0.0
 open(unit=121,file='gulp-2-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=122,file='gulp-2-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=123,file='gulp-2-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"
 read_matrix_2: do
  READ(121,'(A)',IOSTAT=err_apertura) line
  IF( err_apertura /= 0 ) exit read_matrix_2
  read(line,*)k,deg,i,j
  delta_2(i,j,i,j)=1
  delta_2(i,j,j,i)=delta_2(i,j,i,j)
  delta_2(j,i,i,j)=delta_2(i,j,i,j)
  delta_2(j,i,j,i)=delta_2(i,j,i,j)
  deg_2(i,j) = real(deg)
  deg_2(j,i) = deg_2(i,j)
  READ(123,*)epsilon_
  ener_2(i,j)=epsilon_ - ener_0 -ener_1(i)-ener_1(j)
  ener_2(j,i)=ener_2(i,j)
  write(6,*)i,j,i,j,ener_2(i,j),deg_2(i,j)
  do k=2,deg
   READ (122,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_2
   read(line,*)jj,kk,l,m,n
   delta_2(i,j,l,m)=1
   delta_2(j,i,l,m)=delta_2(i,j,l,m)
   delta_2(i,j,m,l)=delta_2(i,j,l,m)
   delta_2(j,i,m,l)=delta_2(i,j,l,m)
   delta_2(l,m,i,j)=delta_2(i,j,l,m)
   delta_2(m,l,i,j)=delta_2(i,j,l,m)
   delta_2(l,m,j,i)=delta_2(i,j,l,m)
   delta_2(m,l,j,i)=delta_2(i,j,l,m)
   write(6,*)i,j,l,m,ener_2(l,m),deg_2(l,m)
  end do
 end do read_matrix_2
! final 
 close(121)
 close(122)
 close(123)
! {{{
call farthrest_nodes_subtitutions(n_atoms,n_T_atoms,n_Ge,delta_1,delta_2,ener_0,&
     ener_1,ener_2,deg_1,deg_2,cell_0,cryst_coor,0,label,MC_steps,temperature)
 stop
end program main
