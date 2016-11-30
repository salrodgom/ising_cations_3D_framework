program main
 implicit none
 integer                    :: op1,ii,deg,configuration_file=456
 integer                    :: i,j,k,l,m,n,jj,kk,ll,mm,nn
 integer                    :: err_apertura,nop1
 integer,parameter          :: n_atoms =   60 ! 60 T + 120 Os + 120 Oc + 12 F
 integer,parameter          :: n_T_atoms = 60
 integer,parameter          :: n_Ge = 4
 integer,parameter          :: MC_steps = 1000
 integer,parameter          :: n_configurations = 4
 character(len=80)          :: file_name,line
 real,parameter             :: temperature = 5000.0
 integer, parameter         :: NOPMAX=10000
 integer                    :: k_max_1,k_max_2,k_max_3
 real                       :: ener_0 = -7745.86721305,epsilon_
 !
 real                       :: ener_1(n_atoms) = 0.0
 real                       :: deg_1(n_atoms) = 0.0
 real                       :: ener_2(n_atoms,n_atoms) = 0.0
 real                       :: deg_2(n_atoms,n_atoms) = 0.0
 real                       :: ener_3(n_atoms,n_atoms,n_atoms) = 0.0
 real                       :: deg_3(n_atoms,n_atoms,n_atoms) = 0.0
 real                       :: ener_4(n_atoms,n_atoms,n_atoms,n_atoms) = 0.0
 real                       :: deg_4(n_atoms,n_atoms,n_atoms,n_atoms) = 0.0
 real                       :: cell_0(1:6)
 logical                    :: MC_flag = .false.
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
 do ii=1,n_configurations
  !write(6,'(i4.4)') ii
  write (file_name, '( "gulp-",i1,"-subs/c", I5.5, ".gin" )' )n_Ge,ii
  write(6,*)'Configuration',ii, file_name
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
 k_max_1 = 0
 open(unit=111,file='gulp-1-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=112,file='gulp-1-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=113,file='gulp-1-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"
 read(111,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k !<-N.configuraciones unicas
!
 read(112,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k_max_1 !<-N.configuraciones totales
 !allocate( ener_1(1:k_max_1) )
 !allocate( deg_1(1:k_max_1))
 ener_1(1:n_atoms)  = 0.0
 deg_1(1:n_atoms)   = 0.0
 if( err_apertura /=0 ) stop
!
 read_matrix_1: do
  read(111,'(a)',iostat=err_apertura) line
  IF( err_apertura /= 0 ) exit read_matrix_1
  read(line,*)i,deg,k
  read(113,*,iostat=err_apertura) epsilon_
  IF( err_apertura /= 0 ) exit read_matrix_1
  ener_1(k)=epsilon_-ener_0
  deg_1(k) =deg
  !write(6,*)k,k,ener_1(k),deg_1(k)
  do j=2,deg
   READ (112,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_1
   read(line,*)l,m,i,n
   ener_1(i)=ener_1(k)
   deg_1(i) =deg_1(k)
   !write(6,*)k,i,ener_1(i),deg_1(i)
  end do
 end do read_matrix_1
! final 
 close(111)
 close(112)
 close(113)
! }}}
 write(6,*)'Two substitutions:'
 k_max_2 = 0
 open(unit=121,file='gulp-2-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=122,file='gulp-2-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=123,file='gulp-2-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"

 read(121,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k !<-N.configuraciones unicas
 read(122,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k_max_2 !<-N.configuraciones totales
 !allocate( ener_2(1:k_max_2,1:k_max_2))
 !allocate( deg_2(1:k_max_2,1:k_max_2))
 ener_2(1:n_atoms,1:n_atoms)  = 0.0
 deg_2(1:n_atoms,1:n_atoms)   = 0.0
 if( err_apertura /=0 ) stop
 ii=0
 read_matrix_2: do
  read(121,'(A)',IOSTAT=err_apertura) line
  IF( err_apertura /= 0 ) exit read_matrix_2
  read(line,*)k,deg,i,j
  READ(123,*,IOSTAT=err_apertura) epsilon_
  IF( err_apertura /= 0 ) exit read_matrix_2
  ii=ii+1
  epsilon_ =   epsilon_-ener_0-ener_1(i)-ener_1(j)
  ener_2(i,j)= epsilon_
  ener_2(j,i)= ener_2(i,j)
  deg_2(i,j) = deg
  deg_2(j,i) = deg_2(i,j)
  !write(6,*)i,j,i,j,ener_2(i,j),deg_2(i,j)
  do k=2,deg
   read(122,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_2
   ii=ii+1
   read(line,*)jj,kk,l,m,n
   ener_2(l,m) = epsilon_
   ener_2(m,l) = epsilon_
   deg_2(l,m) = deg
   deg_2(m,l) = deg
   !write(6,*)i,j,l,m,ener_2(l,m),deg_2(l,m),ii
  end do
 end do read_matrix_2
!!!
 write(6,*)'3 substitutions:'
 k_max_3 = 0
 open(unit=131,file='gulp-3-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=132,file='gulp-3-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=133,file='gulp-3-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"

 read(131,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k !<-N.configuraciones unicas
 read(132,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k_max_3 !<-N.configuraciones totales
 !allocate( ener_3(1:k_max_3,1:k_max_3,1:k_max_3))
 !allocate( deg_3 (1:k_max_3,1:k_max_3,1:k_max_3))
 ener_3(1:n_atoms,1:n_atoms,1:n_atoms)  = 0.0
 deg_3 (1:n_atoms,1:n_atoms,1:n_atoms)  = 0.0
 if( err_apertura /=0 ) stop 
 ii=0
 read_matrix_3: do
  read(131,'(A)',IOSTAT=err_apertura) line
  IF( err_apertura /= 0 ) exit read_matrix_3
  read(line,*)k,deg,i,j,k
  READ(133,*,IOSTAT=err_apertura)epsilon_
  IF( err_apertura /= 0 ) exit read_matrix_3
  ii=ii+1
  epsilon_ = epsilon_-ener_0-ener_1(i)-ener_1(k)-ener_1(j)-ener_2(i,j)-ener_2(j,k)-ener_2(i,k)
  ener_3(i,j,k) = epsilon_
  ener_3(j,i,k) = epsilon_
  ener_3(i,k,j) = epsilon_
  ener_3(k,j,i) = epsilon_
  ener_3(k,i,j) = epsilon_
  ener_3(j,k,i) = epsilon_
  deg_3(i,j,k) = deg
  deg_3(j,i,k) = deg
  deg_3(i,k,j) = deg
  deg_3(k,j,i) = deg
  deg_3(k,i,j) = deg
  deg_3(j,k,i) = deg
  !write(6,*)i,j,k,i,j,k,ener_3(i,j,k),deg_3(i,j,k)
  do k=2,deg
   READ (132,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_3
   ii=ii+1
   read(line,*)jj,kk,l,m,n,ii
   ener_3(l,m,n) = epsilon_
   ener_3(m,l,n) = epsilon_
   ener_3(l,n,m) = epsilon_
   ener_3(n,m,l) = epsilon_
   ener_3(m,n,l) = epsilon_
   ener_3(n,l,m) = epsilon_
   deg_3(l,m,n) = deg
   deg_3(m,l,n) = deg
   deg_3(l,n,m) = deg
   deg_3(n,m,l) = deg
   deg_3(m,n,l) = deg
   deg_3(n,l,m) = deg
   !write(6,*)i,j,k,l,m,n,ener_3(l,m,n),deg_3(l,m,n),ii
  end do
 end do read_matrix_3
! final 
 close(121)
 close(122)
 close(123)
! {{{
 do ii=1,n_configurations
  if(MC_flag)then
  call farthrest_nodes_subtitutions(n_atoms,n_T_atoms,n_Ge,ener_0,ener_1,&
     ener_2,ener_3,deg_1,deg_2,deg_3,cell_0,cryst_coor,n_configurations,label,&
     MC_steps,temperature)
  else
  call geometrical_properties(n_atoms,n_T_atoms,n_Ge,cell_0,cryst_coor,n_configurations,&
       label,deg_1,deg_2,deg_3,ener_0,ener_1,ener_2,ener_3,ii)
  end if
 end do
 stop
end program main

