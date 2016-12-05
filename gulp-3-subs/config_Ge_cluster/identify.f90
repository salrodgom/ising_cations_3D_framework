program main
 implicit none
 integer                    :: ii,configuration_file=456
 integer                    :: i,j,k,l,m,n,jj,kk,ll,mm,nn
 integer                    :: err_apertura,nop1
 integer,parameter          :: n_atoms =   60 ! 60 T + 120 Os + 120 Oc + 12 F
 integer,parameter          :: n_T_atoms = 60
 integer,parameter          :: n_Ge = 4
 integer,parameter          :: n_configurations = 904
 character(len=80)          :: file_name(n_configurations),line
 integer                    :: k_max_1,k_max_2,k_max_3
 character(len=n_T_atoms)   :: word_confg(n_configurations)
 integer                    :: unicas(n_configurations),repetidas(n_configurations,n_configurations)
 !
 real                       :: cell_0(1:6)
 real,dimension(0:n_configurations,n_atoms,3)             :: cryst_coor
 character(len=4),dimension(0:n_configurations,n_atoms,2) :: label
!
 do ii=1,n_configurations
  write (file_name(ii), '( "u", I4.4, ".gin" )' ) ii
  open(unit=configuration_file,file=file_name(ii),status='old',iostat=err_apertura)
  if(err_apertura/=0) stop "[ERROR] reading file // file_name // "
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
        word_confg(ii)(i:i)="#"
      else
        word_confg(ii)(i:i)="."
      end if
    end do
   end if
  end do reading_file
  write(6,'(a60,1x,a10)')word_confg(ii),file_name(ii)
  close(configuration_file)
 end do 
 do ii=1,n_configurations
  unicas(ii)=ii
  repetidas(ii,1:n_configurations)=0
 end do
 check_fake: do ii=1,n_configurations
  check_fake_inner: do jj=ii+1,n_configurations
   if(len_trim(word_confg(ii))==len_trim(word_confg(jj)))then
    if(word_confg(ii)==word_confg(jj))then
     write(6,'(a60,1x,a10,1x,a60,1x,a10)')word_confg(ii),file_name(ii),word_confg(jj),file_name(jj)
     repetidas(ii,jj)=jj
    end if
   end if
  end do check_fake_inner
  !write(6,'(i4,1x,904(i4,1x))')ii,(repetidas(ii,jj),jj=1,n_configurations)
 end do check_fake
 
 unicas(1)=1
 do ii=2,n_configurations
  do jj=1,n_configurations
   if(ii/=jj)then
    do kk=1,n_configurations
     if(ii==repetidas(jj,kk))then
      unicas(ii)=0 
     end if
    end do
   end if
  end do
 end do
 do ii=1,n_configurations
  if(unicas(ii)/=0) write(6,'(a20)') file_name(unicas(ii))
 end do
 stop
end program main

