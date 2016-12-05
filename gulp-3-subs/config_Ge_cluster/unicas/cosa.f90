program main
 implicit none
 integer                    :: ii,configuration_file=456
 integer                    :: i,j,k,l,m,n,jj,kk,ll,mm,nn
 integer                    :: err_apertura,nop1
 integer,parameter          :: n_atoms =   60 ! 60 T + 120 Os + 120 Oc + 12 F
 integer,parameter          :: n_T_atoms = 60
 integer,parameter          :: n_Ge = 4
 integer,parameter          :: n_configurations = 41793
 character(len=80)          :: file_name(n_configurations),line
 integer                    :: k_max_1,k_max_2,k_max_3
 character(len=64)          :: word_confg(n_configurations)
 character(len=64)          :: scratch1,scratch2
 !doubleprecision            :: scratch1,scratch2,scr(n_configurations)
 integer                    :: unicas(n_configurations),repetidas(n_configurations,n_configurations)
 !
 real                       :: cell_0(1:6)
 real,dimension(0:n_configurations,n_atoms,3)             :: cryst_coor
 character(len=4),dimension(0:n_configurations,n_atoms,2) :: label
 logical                    :: flag
!
 do i=1,n_atoms
 word_confg(1:n_configurations)(i:i)="."
 end do
 init_search: do ii=1,n_configurations
  if(ii<=9999)then
   write (file_name(ii), '( "gulp-4/c", I4.4, ".gin" )' ) ii
  else if (ii > 9999.and.ii<=40890) then
   write (file_name(ii), '( "gulp-4/c", I5.5, ".gin" )' ) ii
  else
   write (file_name(ii), '( "u", I4.4, ".gin" )' ) ii-40890
  end if
  inquire(file=file_name(ii),exist=flag)
  !call system ("wc -l " // file_name(ii) )
  if(flag.eqv..false.) cycle init_search
  open(unit=configuration_file,file=file_name(ii),status='old',iostat=err_apertura)
  if(err_apertura/=0)  cycle init_search !stop "[ERROR] reading file // file_name // "
  mm=0
  reading_file: do
   READ (configuration_file,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) exit reading_file
   if(line(1:5)=='cell ') then
    read(configuration_file,'(A)',iostat=err_apertura) line
    if ( err_apertura /= 0 ) exit reading_file
    read(line,*) ( cell_0(j),j=1,6 )
   end if
   if(line(1:4)=='frac') then
    do i=1,n_atoms
     read (configuration_file,'(A)',iostat=err_apertura) line
     if ( err_apertura /= 0 ) exit reading_file
     if(line(1:7)=='species') exit reading_file
      read(line,*)(label(ii,i,j),j=1,2),(cryst_coor(ii,i,j),j=1,3)
      if(label(ii,i,1)(1:2)=="Ge")then
        mm=mm+1
        word_confg(ii)(i:i)="1"
      else
        word_confg(ii)(i:i)="0"
      end if
    end do
   end if
  end do reading_file
  if(mm/=n_Ge) stop
  do i=61,64
   word_confg(ii)(i:i)="0"
  end do
  !read(word_confg(ii)(1:64),'(b64.64)') scr(ii)
  !write(6,'(a64,1x,a20,1x,ES14.7)')word_confg(ii),file_name(ii),scr(ii)
  close(configuration_file)
 end do init_search
 check_fake: do ii=40891,n_configurations
  inquire(file=file_name(ii),exist=flag)
  if(flag)then
  check_fake_inner: do jj=1,40890
   !if(len_trim(word_confg(ii))==len_trim(word_confg(jj)))then
    !read(word_confg(jj)(1:64),'(b64.64)') scratch1
    scratch1=word_confg(ii)
    scratch2=word_confg(jj)
    !if(scr(ii).eq.scratch1)then
    if(scratch1==scratch2)then
     write(6,'(a64,1x,a20,1x,a64,1x,a20)')word_confg(ii),file_name(ii),word_confg(jj),file_name(jj)
     cycle check_fake
    end if
   !end if
  end do check_fake_inner
  write(6,'(a64,1x,a20,1x,a)')word_confg(ii),file_name(ii),'problem!'
  end if
 end do check_fake
 stop
end program main

