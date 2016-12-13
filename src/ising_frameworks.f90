program main
 implicit none
 integer                    :: op1,ii,deg,configuration_file=456
 integer                    :: i,j,k,ijk,l,m,n,lmn,jj,kk,ll,mm,nn,o
 integer                    :: err_apertura,nop1
 integer,parameter          :: n_atoms =   60 ! 60 T + 120 Os + 120 Oc + 12 F
 integer,parameter          :: n_T_atoms = 60
 integer,parameter          :: n_Ge = 4
 integer,parameter          :: MC_steps = 100
 integer,parameter          :: n_configurations = 0
 character(len=120)          :: file_name,line
 real,parameter             :: temperature = 5000.0
 integer, parameter         :: NOPMAX=10000
 integer                    :: k_max_1,k_max_2,k_max_3,k_max_4
 real                       :: ener_0 = -7745.86721305,epsilon_,energy_123
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
 logical                    :: MC_flag = .true.,no_presente=.true.
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
  write(6,*)k,k,ener_1(k),deg_1(k)
  do j=2,deg
   READ (112,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_1
   read(line,*)l,m,i,n
   ener_1(i)=ener_1(k)
   deg_1(i) =deg_1(k)
   write(6,*)k,i,ener_1(i),deg_1(i)
  end do
 end do read_matrix_1
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
  write(6,*)i,j,i,j,ener_2(i,j),deg_2(i,j)
  do k=2,deg
   read(122,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_2
   ii=ii+1
   read(line,*)jj,kk,l,m,n
   ener_2(l,m) = epsilon_
   ener_2(m,l) = epsilon_
   deg_2(l,m) = deg
   deg_2(m,l) = deg
   write(6,*)i,j,l,m,ener_2(l,m),deg_2(l,m),ii
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
 ener_3(1:n_atoms,1:n_atoms,1:n_atoms)  = 0.0
 deg_3 (1:n_atoms,1:n_atoms,1:n_atoms)  = 0.0
 if( err_apertura /=0 ) stop 
 read_matrix_3: do
  read(131,'(A)',IOSTAT=err_apertura) line
  IF( err_apertura /= 0 ) exit read_matrix_3
  read(line,*)ijk,deg,i,j,k
  read(133,*,IOSTAT=err_apertura) epsilon_
  IF( err_apertura /= 0 ) exit read_matrix_3
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
 close(131)
 close(132)
 close(133)
 write(6,*)'4 substitutions:'
 k_max_4 = 0
 open(unit=141,file='gulp-4-subs/OUTSOD',status='old',iostat=err_apertura)
 open(unit=142,file='gulp-4-subs/OUTSOD_modified',status='old',iostat=err_apertura)
 open(unit=143,file='gulp-4-subs/ENERGIES',status='old',iostat=err_apertura)
 if(err_apertura/=0) stop "[ERROR] reading file with matrix 1"
 read(141,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k !<-N.configuraciones unicas
 read(142,'(a)',iostat=err_apertura) line
 if( err_apertura /=0 ) stop
 read(line,*) k_max_4 !<-N.configuraciones totales
 ener_4(1:n_atoms,1:n_atoms,1:n_atoms,1:n_atoms)  = 0.0
 deg_4(1:n_atoms,1:n_atoms,1:n_atoms,1:n_atoms)  = 0.0
 if( err_apertura /=0 ) stop
 read_matrix_4: do
  read(141,'(A)',IOSTAT=err_apertura) line
  if( err_apertura /= 0 ) exit read_matrix_4
  read(line,*)ijk,deg,i,j,k,l
  no_presente=.true.
  check_energy_4: do
   read(143,*,iostat=err_apertura) epsilon_,jj
   if( err_apertura /= 0 )then
    rewind(143)
    no_presente=.true.
    exit check_energy_4
   end if
   if(jj==ijk)then
    no_presente = .false.
    epsilon_ = epsilon_-ener_0-ener_1(i)-ener_1(k)-ener_1(j)-ener_1(l)-ener_2(i,j)-ener_2(j,k)-&
     ener_2(i,k)-ener_2(i,l)-ener_2(j,l)-ener_2(k,l)-&
     ener_3(i,j,k)-ener_3(i,k,j)-ener_3(j,i,k)-ener_3(j,k,i)-ener_3(k,i,j)-ener_3(k,j,i)-&
     ener_3(i,j,l)-ener_3(i,l,j)-ener_3(j,i,l)-ener_3(j,l,i)-ener_3(l,i,j)-ener_3(l,j,i)-&
     ener_3(i,l,k)-ener_3(i,k,l)-ener_3(l,i,k)-ener_3(l,k,i)-ener_3(l,k,i)-ener_3(k,i,l)-&
     ener_3(k,l,i)-ener_3(l,j,k)-ener_3(l,k,j)-ener_3(j,l,k)-ener_3(j,k,l)-ener_3(k,l,j)-&
     ener_3(j,l,k)
    exit check_energy_4
   end if
  end do check_energy_4 
  if(no_presente)then
   label(0,i,1)='Ge  '
   label(0,j,1)='Ge  '
   label(0,k,1)='Ge  '
   label(0,l,1)='Ge  '
   epsilon_=energy_123(n_atoms,n_T_atoms,label(0,1:1:n_T_atoms,1),ener_0,ener_1,&
    ener_2,ener_3,deg_1,deg_2,deg_3)
   label(0,i,1)='Ge  '
   label(0,j,1)='Ge  '
   label(0,k,1)='Ge  '
   label(0,l,1)='Ge  '
  end if
  ener_4(i,j,k,l) = epsilon_
  ener_4(i,j,l,k) = epsilon_
  ener_4(i,k,l,j) = epsilon_
  ener_4(i,k,j,l) = epsilon_
  ener_4(i,l,j,k) = epsilon_
  ener_4(i,l,k,j) = epsilon_
  ener_4(j,i,k,l) = epsilon_
  ener_4(j,i,l,k) = epsilon_
  ener_4(j,k,i,l) = epsilon_
  ener_4(j,k,l,i) = epsilon_
  ener_4(j,l,k,i) = epsilon_
  ener_4(j,l,i,k) = epsilon_
  ener_4(k,i,j,l) = epsilon_
  ener_4(k,i,l,j) = epsilon_
  ener_4(k,j,i,l) = epsilon_
  ener_4(k,j,l,i) = epsilon_
  ener_4(k,l,i,j) = epsilon_
  ener_4(k,l,j,i) = epsilon_
  ener_4(l,i,j,k) = epsilon_
  ener_4(l,i,k,j) = epsilon_
  ener_4(l,j,i,k) = epsilon_
  ener_4(l,j,k,i) = epsilon_
  ener_4(l,k,i,j) = epsilon_
  ener_4(l,k,j,i) = epsilon_
  deg_4(i,j,k,l) = deg
  deg_4(i,j,l,k) = deg
  deg_4(i,k,l,j) = deg
  deg_4(i,k,j,l) = deg
  deg_4(i,l,j,k) = deg
  deg_4(i,l,k,j) = deg
  deg_4(j,i,k,l) = deg
  deg_4(j,i,l,k) = deg
  deg_4(j,k,i,l) = deg
  deg_4(j,k,l,i) = deg
  deg_4(j,l,k,i) = deg
  deg_4(j,l,i,k) = deg
  deg_4(k,i,j,l) = deg
  deg_4(k,i,l,j) = deg
  deg_4(k,j,i,l) = deg
  deg_4(k,j,l,i) = deg
  deg_4(k,l,i,j) = deg
  deg_4(k,l,j,i) = deg
  deg_4(l,i,j,k) = deg
  deg_4(l,i,k,j) = deg
  deg_4(l,j,i,k) = deg
  deg_4(l,j,k,i) = deg
  deg_4(l,k,i,j) = deg
  deg_4(l,k,j,i) = deg
  !write(6,*)i,j,k,l,i,j,k,l,ener_4(i,j,k,l),deg_4(i,j,k,l)
  ll=l
  do kk=2,deg
   READ (142,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_4
   read(line,*)jj,lmn,l,m,n,o,ii
    ener_4(l,m,n,o) = epsilon_
    ener_4(l,m,o,n) = epsilon_
    ener_4(l,n,o,m) = epsilon_
    ener_4(l,n,m,o) = epsilon_
    ener_4(l,o,m,n) = epsilon_
    ener_4(l,o,n,m) = epsilon_
    ener_4(m,l,n,o) = epsilon_
    ener_4(m,l,o,n) = epsilon_
    ener_4(m,n,l,o) = epsilon_
    ener_4(m,n,o,l) = epsilon_
    ener_4(m,o,n,l) = epsilon_
    ener_4(m,o,l,n) = epsilon_
    ener_4(n,l,m,o) = epsilon_
    ener_4(n,l,o,m) = epsilon_
    ener_4(n,m,l,o) = epsilon_
    ener_4(n,m,o,l) = epsilon_
    ener_4(n,o,l,m) = epsilon_
    ener_4(n,o,m,l) = epsilon_
    ener_4(o,l,m,n) = epsilon_
    ener_4(o,l,n,m) = epsilon_
    ener_4(o,m,l,n) = epsilon_
    ener_4(o,m,n,l) = epsilon_
    ener_4(o,n,l,m) = epsilon_
    ener_4(o,n,m,l) = epsilon_
    deg_4(l,m,n,o) = deg
    deg_4(l,m,o,n) = deg
    deg_4(l,n,o,m) = deg
    deg_4(l,n,m,o) = deg
    deg_4(l,o,m,n) = deg
    deg_4(l,o,n,m) = deg
    deg_4(m,l,n,o) = deg
    deg_4(m,l,o,n) = deg
    deg_4(m,n,l,o) = deg
    deg_4(m,n,o,l) = deg
    deg_4(m,o,n,l) = deg
    deg_4(m,o,l,n) = deg
    deg_4(n,l,m,o) = deg
    deg_4(n,l,o,m) = deg
    deg_4(n,m,l,o) = deg
    deg_4(n,m,o,l) = deg
    deg_4(n,o,l,m) = deg
    deg_4(n,o,m,l) = deg
    deg_4(o,l,m,n) = deg
    deg_4(o,l,n,m) = deg
    deg_4(o,m,l,n) = deg
    deg_4(o,m,n,l) = deg
    deg_4(o,n,l,m) = deg
    deg_4(o,n,m,l) = deg
    !write(6,*)i,j,k,ll,l,m,n,o,ener_4(l,m,n,o),deg_4(l,m,n,o)
  end do
 end do read_matrix_4
! final 
 close(141)
 close(142)
 close(143)

! {{{
 if(MC_flag)then
    call farthrest_nodes_subtitutions(n_atoms,n_T_atoms,n_Ge,ener_0,ener_1,&
     ener_2,ener_3,ener_4,deg_1,deg_2,deg_3,deg_4,cell_0,cryst_coor,&
     n_configurations,label,MC_steps,temperature)
 else
  do ii=1,n_configurations
   call geometrical_properties(n_atoms,n_T_atoms,n_Ge,cell_0,cryst_coor,n_configurations,&
       label,deg_1,deg_2,deg_3,deg_4,ener_0,ener_1,ener_2,ener_3,ener_4,ii)
  end do
 end if
 stop
end program main

real FUNCTION energy_123(n,n_T,lab,ener_0,ener_1,ener_2,ener_3,deg_1,deg_2,deg_3)
  implicit none
  integer, intent(in) :: n,n_T
  character (LEN=4),intent(in)   :: lab(1:n)
  integer             :: i,j,k,l
  real,intent(in)     :: ener_0,ener_1(n),ener_2(n,n),ener_3(n,n,n)
  real,intent(in)     :: deg_1(n),deg_2(n,n),deg_3(n,n,n)
  integer             :: spin(n),ss
  ss = 0
  spin(1:n) = 0
  energy_123 = ener_0
!
  get_spin: do i=1,n
   if(lab(i)=='Ge  ') then
    spin(i) = 1
    ss=ss+1
   end if
  end do get_spin
  if ( ss==0 ) stop 'no hay Ge, bro'
!
  do1subs: do i=1,n
   energy_123 = energy_123 + ener_1(i)*spin(i)
   do2subs: do j=i+1,n
    energy_123 = energy_123 + ener_2(i,j)*spin(i)*spin(j)
    do3subs: do k=j+1,n
     energy_123 = energy_123 + ener_3(i,j,k)*spin(i)*spin(j)*spin(k)
    end do do3subs
   end do do2subs
  end do do1subs
  return
END FUNCTION energy_123
