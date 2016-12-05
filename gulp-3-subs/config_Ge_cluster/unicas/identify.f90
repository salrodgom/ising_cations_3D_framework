program main
 implicit none
 integer                    :: configuration_file=456
 integer                    :: ijk,i,j,k,l,m,n,ii,jj,kk,ll,mm,nn
 integer                    :: ei,ej,ek
 integer                    :: err_apertura,nop1
 integer,parameter          :: n_atoms =   60 ! 60 T + 120 Os + 120 Oc + 12 F
 integer,parameter          :: n_T_atoms = 60
 integer,parameter          :: n_Ge = 4
 integer,parameter          :: n_configurations = 904
 character(len=80)          :: file_name(n_configurations),line
 integer                    :: k_max_1,k_max_2,k_max_3,k_max_4
 character(len=64)   :: word_confg(n_configurations)
 character(len=64)   :: word_confg_sod,tmp_sod2,tmp_sod,make_word
 integer                    :: unicas(n_configurations),repetidas(n_configurations,n_configurations)
 !
 real                       :: cell_0(1:6)
 doubleprecision            :: scr(n_configurations),scratch1,scratch2
 real,dimension(0:n_configurations,n_atoms,3)             :: cryst_coor
 character(len=4),dimension(0:n_configurations,n_atoms,2) :: label
 logical                    :: flag,compare_two_strings
!
 init_search: do ii=1,n_configurations
  write (file_name(ii), '( "u", I4.4, ".gin" )' ) ii
  inquire(file=file_name(ii),exist=flag)
  if(flag.eqv..false.) cycle init_search
  open(unit=configuration_file,file=file_name(ii),status='old',iostat=err_apertura)
  if(err_apertura==6)  cycle init_search
  if(err_apertura/=0)  cycle init_search !stop "[ERROR] reading file // file_name // "
  reading_file: do
   nn=0
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
      if(label(ii,i,1)(1:4)=="Ge  ")then
        word_confg(ii)(i:i)="1"
        nn=nn+1
      else
        word_confg(ii)(i:i)="0"
      end if
    end do
    if(nn/=n_Ge)write(6,*)'Diferente número de Ge'
    do i=61,64
     word_confg(ii)(i:i)="0"
    end do
    read(word_confg(ii)(1:64),'(b64.64)') scr(ii)
   end if
  end do reading_file
  write(6,'(a60,1x,a20)')word_confg(ii),file_name(ii)
  close(configuration_file)
 end do init_search
 ii=1
 nop1=0
 check_fake: do !ii=1,n_configurations
  if(ii>=n_configurations) exit check_fake
  inquire(file=file_name(ii),exist=flag)
  open(unit=configuration_file,file=file_name(ii),status='old',iostat=err_apertura)
  !write(6,*)file_name(ii),flag,err_apertura
  if(flag.eqv..false..or.err_apertura==6) then
   ii=ii+1
   close(configuration_file)
   cycle check_fake
  end if
  open(unit=111,file='gulp-4/OUTSOD',status='old',iostat=err_apertura)
  open(unit=112,file='gulp-4/OUTSOD_modified',status='old',iostat=err_apertura)
  read(111,'(a)',iostat=err_apertura) line
  if( err_apertura /=0 ) stop
  read(line,*) k       ! <- N.configuraciones unicas
  read(112,'(a)',iostat=err_apertura) line
  if( err_apertura /=0 ) stop
  read(line,*) k_max_4 ! <- N.configuraciones totales
  if( err_apertura /=0 ) stop
  read_matrix_3: do
   read(111,'(A)',iostat=err_apertura) line
   IF( err_apertura /= 0 ) exit read_matrix_3
   read(line,*)kk,ll,i,j,k,l
   word_confg_sod=make_word(64,i,j,k,l)
   !read(word_confg_sod(1:64),'(b64.64)') scratch1
   !scratch2=scr(ii)
   if(word_confg_sod==word_confg(ii))then
   !flag=compare_two_strings(n_atoms,word_confg_sod,tmp_sod2,ijk)
   !if( flag .or. ijk==0 ) then
    if(kk<=9999) write (line, '( "c", I4.4, ".gin" )' ) kk
    if(kk> 9999) write (line, '( "c", I5.5, ".gin" )' ) kk
    write(6,'(i9,1x,i9,1x,a64,1x,a64,1x,a15,1x,a15)')ii,kk,word_confg(ii),word_confg_sod,file_name(ii),line
    ii=ii+1
    nop1=0
    cycle check_fake
   else
    tmp_sod=word_confg_sod
    do k=2,ll
     read (112,*,iostat=err_apertura) ei,ej,m,n,mm,nn,ek
     IF( err_apertura /= 0 ) exit read_matrix_3
     word_confg_sod=make_word(64,m,n,mm,nn)
     !read(word_confg_sod(1:64),'(b64.64)') scratch1
     if(scratch1==scratch2)then
     !flag= compare_two_strings(n_atoms,word_confg_sod,tmp_sod2,ijk)
     !if(flag .or. ijk==0) then
      if(kk<=9999) write (line, '( "c", I4.4, ".gin" )' ) kk
      if(kk> 9999) write (line, '( "c", I5.5, ".gin" )' ) kk
      write(6,'(i9,1x,i9,1x,a64,1x,a64,1x,a20,1x,a20)')ii,kk,word_confg(ii),tmp_sod,file_name(ii),line
      write(6,'(9x,1x,9x,1x,64x,1x,a64)')word_confg_sod
      ii=ii+1
      nop1=0
      cycle check_fake
     end if
    end do
   end if
  end do read_matrix_3
  close(111)
  close(112)
  !write(6,'(i9,1x,9x,1x,a64,1x,a64,1x,a20)')ii,word_confg(ii),'No encontrado',file_name(ii)
  !ii=ii+1
  nop1=nop1+1
  if(nop1>=10) then
   write(6,'(i9,1x,9x,1x,a64,1x,a64,1x,a20)')ii,word_confg(ii),'No encontrado',file_name(ii)
   ii=ii+1
   nop1=0
  end if
 end do check_fake
 stop
end program main

logical function compare_two_strings(n,word1,word2,comp)
 implicit none
 integer,intent(in)    :: n
 character(len=n)      :: word1,word2
 integer               :: i
 integer,intent(out)   :: comp
 !compare_two_strings = .true.
 compare_strings: do i=1,n
  comp=ichar(word1(i:i))-ichar(word2(i:i))
  if(comp==0)then
   compare_two_strings=.true.
  else 
   compare_two_strings=.false.
   exit compare_strings
  end if
 end do compare_strings
 return
end function compare_two_strings

subroutine variations_ijkl(n,i,j,k,l,word)
 implicit none
 integer,intent(in)            :: n
 integer,intent(in)            :: i,j,k,l
 integer                       :: ii,jj,kk,ll,id(4),ijk,var
 character(len=n)              :: word(24),make_word
 id(1)=i
 id(2)=j
 id(3)=k
 id(4)=l
 var=0
 do ii=1,4
  do jj=1,4
   do kk=1,4
    do ll=1,4
     if(ii/=jj.and.ii/=kk.and.ii/=ll.and.jj/=kk.and.jj/=ll.and.kk/=ll) then
      var=var+1
      word(var)=make_word(n,id(ii),id(jj),id(kk),id(ll))
      write(6,*)var,word(var),ii,jj,kk,ll,id(ii),id(jj),id(kk),id(ll)
     end if
    enddo
   end do
  end do
 enddo
 return
end subroutine

character(len=*) function make_word(n,i,j,k,l)
 implicit none
 integer,intent(in)             :: n,i,j,k,l
 integer                        :: ijk
 do ijk=1,n
  make_word(ijk:ijk)="0"
 end do
 make_word(i:i)="1"
 make_word(j:j)="1"
 make_word(k:k)="1"
 make_word(l:l)="1"
 return
end function make_word
