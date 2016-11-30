!call farthrest_nodes_subtitutions(n_atoms,n_T_atoms,n_Ge,ener_0,ener_1,&
!     ener_2,ener_3,deg_1,deg_2,deg_3,cell_0,cryst_coor,0,label,&
!     MC_steps,temperature)
subroutine farthrest_nodes_subtitutions(n_atoms,n_T,n_Al,ener_0,ener_1,&
      ener_2,ener_3,deg_1,deg_2,deg_3,cell_0,cryst_coor,n_configurations,&
      labels,MC_cycles,T)
 implicit none
 integer,intent(in) :: n_atoms,n_T
 real,intent(in)    :: deg_1(n_atoms),deg_2(n_atoms,n_atoms),deg_3(n_atoms)
 real,intent(in)    :: ener_0, ener_1(n_atoms),ener_2(n_atoms,n_atoms)
 real,intent(in)    :: ener_3(n_atoms,n_atoms,n_atoms)
 integer,intent(in) :: n_configurations,n_Al,MC_cycles
 real,intent(in)    :: T
 real,intent(in)    :: cell_0(6),cryst_coor(0:n_configurations,n_atoms,3)
 character(len=4),intent(inout) :: labels(0:n_configurations,n_atoms,2)
!
 real               :: dist_matrix(n_atoms,n_atoms)!,energy_matrix(n_T,n_T)
 real               :: energy,repulsive_potential_Lowenstein,distance
 real               :: r4_uniform
 real               :: xcryst(0:3,1:n_atoms)
 character(len=4)   :: label(n_atoms), id(1:n_atoms)
 character(len=1)   :: adj_char(1:n_atoms,1:n_atoms)
 integer            :: pivots(1:n_Al), adj(n_atoms,n_atoms)
!
 INTEGER            :: n_Si = 0
 INTEGER            :: walker_position,walker_drunked
 REAL               :: r1  = 1.2
 REAL               :: r2  = 2.0
 REAL,    PARAMETER :: infinite = 9999999.999999
 REAL,    PARAMETER :: solera   = 1e-10
 INTEGER            :: IERR,SEED
 INTEGER            :: spacegroup = 1
 INTEGER            :: i,j,k,l,m
 REAL               :: atom(3),ouratom(3),pot_dist
 REAL               :: cost,p,q,r,s,rv(3,3),vr(3,3)
 
 CHARACTER (LEN=6)  :: potential,atomc
 CHARACTER (LEN=80) :: line,string
 CHARACTER (LEN=4)  :: mol
 LOGICAL            :: FLAG      = .false.
! {{
 do i=1,n_atoms
  label(i)=labels(n_configurations,i,1)
  do j=1,3
  xcryst(j,i) = cryst_coor(0,i,j)
  enddo
 enddo
!
 CALL get_seed(SEED)
 WRITE(6,'(a,2x,i10)')'Random Seed:',SEED
 WRITE(6,'(a,2x,i4)')'Number of atoms:',n_atoms
 WRITE(6,'(a,2x,i3)')'Space Group:',spacegroup
 CALL cell(rv,vr,cell_0)
! {{ creamos el grafo de atomos }}
 make_graph: IF( FLAG.EQV..true. ) THEN
    pairs: DO i=1,n_atoms
      DO j=i,n_atoms
      adj(i,j)=0
      adj_CHAR(i,j)=' '
      if(labels(0,i,2)=='shel') cycle pairs
      if(labels(0,j,2)=='shel') cycle pairs
      IF(j/=i)THEN
       FORALL ( k=1:3 )
          atom(k)= xcryst(k,j)
          ouratom(k) = xcryst(k,i)
       END FORALL
       call make_distances(cell_0,ouratom,atom,rv,r)
       IF(r>=r1.AND.r<r2)THEN
         adj(i,j)=1
         adj_CHAR(I,J)='@'
       ENDIF
      ENDIF
      adj(j,i)=adj(i,j)
      adj_CHAR(J,I)=adj_CHAR(I,J)
     ENDDO
    ENDDO pairs
    conectivity: DO i=1,n_atoms ! conectividad de los ATOM  ! k=4 y k=2    en zeolitas.
     k=0
     DO j=1,n_atoms
        k=k+adj(i,j)
     ENDDO
     xcryst(0,i)=k
    ENDDO conectivity
    WRITE(6,'(A)')'[loops] Connectivity between nodes [Degree]:'
    WRITE(6,'(1000(I1))')(int(xcryst(0,k)),k=1,n_atoms)
 END IF make_graph
! make dist_matrix
 DO i=1,n_atoms
    dist_matrix(i,i)=0.0
    DO j=i+1,n_atoms
       forall ( k=1:3 )
        ouratom(k)=xcryst(k,i)
        atom(k)   =xcryst(k,j)
       end forall
       call make_distances(cell_0,ouratom,atom,rv,s)
       dist_matrix(i,j)=s
       dist_matrix(j,i)=dist_matrix(i,j)
    END DO
 END DO
! {{ sustituimos los aluminios
!    {{ el primer Al es aleatorio
 pivots(1)=INT(r4_uniform(1.0,n_T+1.0,SEED))
 DO WHILE ( label(pivots(1))/='Si' ) 
  pivots(1)=INT(r4_uniform(1.0,n_T+1.0,SEED))
 ENDDO
 m = pivots(1)
 label(m)='Ge'
 write(6,*)m
 q=3.40
 WRITE(6,'(a)')'Repulsive substitution with cost-function:'
 WRITE(6,'(a,f5.2,a)')'cost = 1/(s - p0),  p0 =',q,' 10^-10 m, if s >= p0'
 WRITE(6,'(a,f14.3,a)')'cost = ',infinite,' if s <  p0'
 WRITE(6,'(a)')'     Crystalographic positions                  Cost     Label'
 write(6,'(60a2)')( label(j), j=1,n_atoms )
 r = energy(n_atoms,n_T,label,ener_0,ener_1,ener_2,ener_3,deg_1,deg_2,deg_3)
 write(6,*)m,label(m),r
 if(r==0.0) STOP 'energy 0???'
 !WRITE(6,'(3(f14.6,1x),f20.10,1x,i3,1x,a2)')xcryst(1,pivots(1)),xcryst(2,pivots(1)),xcryst(3,pivots(1)),&
 !      0.0,pivots(1),label( pivots(1) )
! }}
! {{ colocamos los siguienets Al
 walker_position = pivots(1)
 walker_drunked = 0
 choose_Al: DO k=2,n_Al
  m = SEED
  r = infinite
  jump: DO walker_position=1,n_T
     IF (label(walker_position)=='Si') THEN
       pot_dist=0.0
       p=energy(n_atoms,n_T,label,ener_0,ener_1,ener_2,ener_3,deg_1,deg_2,deg_3)
       DO j=1,k-1
          walker_drunked = pivots(j)
          s = dist_matrix(walker_position,walker_drunked)
          q = repulsive_potential_Lowenstein(s,3.40,infinite,potential)
          pot_dist = pot_dist + q/real(k)
       END DO
       pot_dist = pot_dist + p
       IF ( pot_dist < r ) THEN
          r = pot_dist
          m = walker_position
       END IF
     END IF
  END DO jump
  IF (m==SEED) THEN
     m=INT(r4_uniform(1.0,n_T+1.0,SEED))
     DO WHILE ( label(m)/='Si' )
       m=INT(r4_uniform(1.0,n_T+1.0,SEED))
     ENDDO
  END IF
  pivots(k) = m
  label(pivots(k))='Ge'
  write(6,'(60a2)')( label(j), j=1,n_atoms )
  r = energy(n_atoms,n_T,label,ener_0,ener_1,ener_2,ener_3,deg_1,deg_2,deg_3)
  write(6,*)pivots(k),label(pivots(k)),r
 END DO choose_Al                       ! }}i
 q=infinite
 k=0
 WRITE(6,'(a)')'=============================================='
 WRITE(6,'(a)')'MonteCarlo: Metropolis with identity changes'
 WRITE(6,'(a,f10.5)')'Si > Ge, Ge > Si. STOP when occurrence < ',solera
 WRITE(6,'(a)')'=============================================='
 WRITE(6,'(a,5x,a,5x,a)')'Cost_function','Occurrence %','Scan-Cost'
 !DO WHILE ( q >= solera .and. k<=MC_cycles )
 DO WHILE ( k<=MC_cycles )
    ! Metropolis:
    call MonteCarlo(n_atoms,n_T,dist_matrix,n_Al,SEED,xcryst,label,cell_0,rv,q,&
         ener_0,ener_1,ener_2,ener_3,deg_1,deg_2,deg_3,cost,T)
    k=k+1
 END DO

 do i=1,n_atoms
  labels(n_configurations,i,1)=label(i)
 enddo

 call write_gin(cell_0,xcryst,n_atoms,labels)
 WRITE(6,'(a)')'=============================================='
 WRITE(6,'(a,i3,a)')'geometrical properties of ',n_Al,' Ge:'
 q=0.0
 r=0.0
 m=0
 do i=1,n_atoms
   if(label(i)=='Ge')then
    do j=1,n_atoms
     if(i/=j)then
      if(label(j)=='Ge')then
       s = dist_matrix(i,j)
       m=m+1
       r=r+s
       q=q+s*s
      end if
     endif
    end do
   end if
 end do
 WRITE(6,*)'distancia media y dispersion'
 WRITE(6,*)r/real(m),sqrt(q/real(m)-(r*r)/real(m*m))
 pot_dist=r/real(m)
 q=0.0
 r=0.0
 m=0
 do i=1,n_atoms
   if(label(i)=='Ge')then
    p=infinite
    do j=1,n_atoms
     if(i/=j)then
      if(label(j)=='Ge')then
       s=dist_matrix(i,j)
       if( s <= p ) p=s
      end if
     endif
    end do
    m=m+1
    r=r+p
    q=q+p*p
   endif
 end do
 WRITE(6,*)'distancia minima media y dispersion'
 WRITE(6,*)r/real(m),sqrt(q/real(m)-(r*r)/real(m*m))
 WRITE(6,*)'STOP',pot_dist,r/real(m),cost
end subroutine
!
 REAL FUNCTION repulsive_potential_Lowenstein(r,p0,p1,p2)
  IMPLICIT NONE
  REAL              :: r,p0,p1
  CHARACTER (LEN=6) :: p2
  if(p2/='coulom'.and.p2/='london') p2='zerooo'
  IF (r <  p0) repulsive_potential_Lowenstein = p1
  if (r >= p0.and.p2=='coulom') repulsive_potential_Lowenstein = (r-p0)**(-1)
  if (r >= p0.and.p2=='london') repulsive_potential_Lowenstein = (r-p0)**(-6)
  if (r >= p0.and.p2=='zerooo') repulsive_potential_Lowenstein = 0.00
 END FUNCTION repulsive_potential_Lowenstein
!

 SUBROUTINE MonteCarlo(n_atoms,n_T,dist_matrix,n_Al,SEED,xcryst,label,cell_0,rv,exito,&
  ener_0,ener_1,ener_2,ener_3,deg_1,deg_2,deg_3,coste,temperature)
  IMPLICIT NONE
  INTEGER, intent(in) :: n_atoms,n_Al,SEED
  REAL,    intent(in) :: xcryst(0:3,1:n_atoms),dist_matrix(n_atoms,n_atoms)
  REAL,    intent(in) :: cell_0(1:6),rv(1:3,1:3)
  integer,intent(in)  :: n_T
  real,intent(in)     :: ener_1(n_atoms),ener_2(n_atoms,n_atoms),ener_0
  real,intent(in)     :: ener_3(n_atoms,n_atoms,n_atoms)
  real,intent(in)     :: deg_1(n_atoms),deg_2(n_atoms,n_atoms),deg_3(n_atoms,n_atoms,n_atoms)
  real                :: R4_UNIFORM,energy
  integer,parameter   :: max_ident_numer=10
  CHARACTER (LEN=4)   :: label(1:n_atoms),cation_distribution(0:max_ident_numer,n_T)
  INTEGER             :: i,j,k,l,m,iii,jjj
  REAL                :: energia(0:2) = 0.0
  REAL                :: eta = 0.0
  REAL,intent(in)     :: temperature
  real                :: k_B = 8.61734E-5
  REAL                :: delta = 0.0
  REAL,    PARAMETER  :: infinite = 9999999.999999
  REAL,    intent(out):: exito,coste
  exito = 0.0
  MC_step: DO k=0,n_T-1
    if(k==0)then
      energia(1) = energy(n_atoms,n_T,label,ener_0,ener_1,ener_2,ener_3,deg_1,deg_2,deg_3)
      eta = energia(1)
    end if
    cation_distribution(0,1:n_T) = label(1:n_T)
    m = INT(R4_UNIFORM(1.0,real(max_ident_numer+1),SEED))
    do l=1,m
     scan_repetitions: do while(1>0)
      call identity_change_move(label,n_T,n_atoms,seed)
      cation_distribution(l,1:n_T) = label(1:n_T)
      jjj=0
      do iii=1,n_T
       if(cation_distribution(l,iii)=='Ge  ') jjj=jjj+1
      end do
      if(jjj==n_Al)  exit scan_repetitions
     end do scan_repetitions
    end do 
    energia(2) = energy(n_atoms,n_T,label,ener_0,ener_1,ener_2,ener_3,deg_1,deg_2,deg_3)
    IF ( energia(2) > energia(1) ) THEN
       eta = R4_UNIFORM(0.0,1.0,SEED)
       IF ( eta > exp(-(energia(2)-energia(1))/(k_B*temperature)) ) THEN
          do l=1,n_T
           label(l) = cation_distribution(0,l)
          end do
          energia(0) = energia(1)
       ELSE ! se aprueba aun siendo la energia mayor
          exito = exito + 1.0
          energia(0) = energia(1)
          energia(1) = energia(2)
       ENDIF
    ELSE ! si es menor, se aprueba
       energia(0) = energia(1)
       energia(1) = energia(2)
       exito = exito + 1.0
    END IF
    write(6,'(f20.10,1x,60a2)')energia(1),(label(j)(1:2),j=1,n_atoms)
  END DO MC_step
  exito   = exito/REAL(n_atoms*n_atoms)
  coste = energia(1)
  eta = coste - eta
! coste, ocurrencia, exito
  !write(6,*) coste,exito,eta
  RETURN
 END SUBROUTINE MonteCarlo
!
 SUBROUTINE identity_change_move(label,n_T,n_atoms,seed)
  IMPLICIT NONE
  integer,intent(in)              :: n_atoms,n_T,seed
  CHARACTER (LEN=4),intent(inout) :: label(1:n_atoms)
  INTEGER             :: i,j
  real                            :: R4_UNIFORM
!
  i = INT(R4_UNIFORM(1.0,real(n_T)+1.0,SEED))
  DO WHILE ( label(i) /= 'Si' )
    i = INT(R4_UNIFORM(1.0,real(n_T)+1.0,SEED))
  END DO
  j = INT(R4_UNIFORM(1.0,real(n_T)+1.0,SEED))
  DO WHILE ( label(j) /= 'Ge' )
    j = INT(R4_UNIFORM(1.0,real(n_T)+1.0,SEED))
  END DO
! identity_change
  label(i) = 'Ge'
  label(j) = 'Si'
  RETURN
 END SUBROUTINE identity_change_move
!
 REAL FUNCTION energy(n,n_T,lab,ener_0,ener_1,ener_2,ener_3,deg_1,deg_2,deg_3)
  IMPLICIT NONE
  INTEGER, intent(in) :: n,n_T
  CHARACTER (LEN=4),intent(in)   :: lab(1:n)
  INTEGER             :: i,j,k,l
  real,intent(in)     :: ener_0,ener_1(n),ener_2(n,n),ener_3(n,n,n)
  real,intent(in)     :: deg_1(n),deg_2(n,n),deg_3(n,n,n)
  integer             :: spin(n),ss
  ss = 0
  spin(1:n) = 0
  energy = ener_0
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
   energy = energy + ener_1(i)*spin(i)
   do2subs: do j=i+1,n
    energy = energy + ener_2(i,j)*spin(i)*spin(j)
    do3subs: do k=j+1,n
     energy = energy + ener_3(i,j,k)*spin(i)*spin(j)*spin(k)
    end do do3subs
   end do do2subs
  end do do1subs
  !do i =1,n
  ! do j=1,n
  !  if(deg_1(j)/=0.0.and.delta_1(j,i)*spin(i)==1)then
  !   energy = energy + (ener_1(j))*delta_1(j,i)*spin(i)
  !  end if
  ! end do
  !end do
  !if(ss>=2)then
  !do i =1, n
  ! do j =1+i,n
  !  do k=1,n
  !   do l=1+k,n
  !    if(deg_2(k,l)/=0.0.and.delta_2(k,l,i,j)*spin(i)*spin(j)==1)then
  !     energy = energy + (ener_2(k,l))*delta_2(k,l,i,j)*spin(i)*spin(j)
  !    end if
  !   end do
  !  end do
  ! end do
  !end do
  !end if
  return
 END FUNCTION energy
!
 SUBROUTINE cell(rv,vr,cell_0)
! ======================
! GULP
! ======================
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = ACOS(-1.0)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: DEGTORAD
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
 call inverse(rv,vr,3)
 WRITE(6,'(a)') 'Cell:'
 WRITE(6,'(6F14.7)')( cell_0(j), j=1,6 )
 WRITE(6,'(a)')'Box:'
 DO i=1,3
    WRITE(6,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 WRITE(6,'(a)')'Inverse Box:'
 DO i=1,3
    WRITE(6,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 RETURN
 END SUBROUTINE cell
!
 SUBROUTINE uncell(rv,cell_0)
! ======================
! GULP
! ======================
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0) 
  radtodeg=180.0/PI
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
!
 SUBROUTINE inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!==========================================================
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 L=0.0
 U=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  L(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 RETURN
 END SUBROUTINE inverse
!
 SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
! {{ prepara para el calculo de distancias en una celda triclinica }}
 IMPLICIT NONE
 REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist
 real                 :: distance
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
! REAL                 :: distance
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
! {{ calculamos la matriz de distancias }}
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
   enddo
  enddo
  dist=MINVAL(d_image)
  RETURN
 END SUBROUTINE
!
 REAL FUNCTION DISTANCE(atom,ouratom,rv)
  IMPLICIT NONE
  INTEGER :: j
  REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
  REAL :: rv(3,3)
  FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
  END FORALL
  DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
 END FUNCTION
!
 subroutine get_seed(seed)
  implicit none
  integer day,hour,i4_huge,seed,milli,minute,month,second,year
  parameter (i4_huge=2147483647)
  double precision temp
  character*(10) time
  character*(8) date
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  temp=0.0D+00
  temp=temp+dble(month-1)/11.0D+00
  temp=temp+dble(day-1)/30.0D+00
  temp=temp+dble(hour)/23.0D+00
  temp=temp+dble(minute)/59.0D+00
  temp=temp+dble(second)/59.0D+00
  temp=temp+dble(milli)/999.0D+00
  temp=temp/6.0D+00
!  Force 0 < TEMP <= 1.
  DO
    IF(temp<=0.0D+00 )then
       temp=temp+1.0D+00
       CYCLE
    ELSE
       EXIT
    ENDIF
  ENDDO
  DO
    IF(1.0D+00<temp)then
       temp=temp-1.0D+00
       CYCLE
    else
       EXIT
    ENDIF
  ENDDO
  seed=int(dble(i4_huge)*temp)
  if(seed==0)       seed = 1
  if(seed==i4_huge) seed = seed-1
  RETURN
 END subroutine get_seed

 REAL function r4_uniform(b1,b2,seed)
  implicit none
  real b1,b2
  integer i4_huge,k,seed
  parameter (i4_huge=2147483647)
  if(seed==0)then
       write(*,'(b1)')' '
       write(*,'(b1)')'R4_UNIFORM - Fatal error!'
       write(*,'(b1)')'Input value of SEED = 0.'
      stop '[ERROR] r4_uniform'
  end if
  k=seed/127773
  seed=16807*(seed-k*17773)-k*2836
  if(seed<0) then
    seed=seed+i4_huge
  endif
  r4_uniform=b1+(b2-b1)*real(dble(seed)* 4.656612875D-10)
  RETURN
 end function r4_uniform
!
 SUBROUTINE write_grafo(A,N,U)
  integer, intent(in) :: N,A(N,N),U
  integer :: X,Y
  nodos: DO X=1,N
   DO Y=X+1,N
      IF(A(X,Y)>0.0) THEN
        IF(X>Y)THEN
         WRITE(U,'(I5,1x,I5,1x,A)') X,Y,'+'
        ELSE
         WRITE(U,'(I5,1x,I5,1x,A)') Y,X,'+'
        ENDIF
      ENDIF
   ENDDO
  ENDDO nodos
  RETURN
 END SUBROUTINE write_grafo
!
 SUBROUTINE ESCRITURA_GRAFO(A,N,U)
!--------------------------------------------
! ADAPTAMOS LA MATRIS DE ADYACENCIA PARA QUE EL PROGRAMA DE
! REPRESENTACION DE GRAFOS LO ENTIENDA. LA SALIDA LA TIENEN UN ARCHIVO
! LLAMADO: "GRAFO.DOT"
!--------------------------------------------
  INTEGER N,X,Y,A(N,N),U
  OPEN(U,file='grafo.dot')
  write(U,*)'graph G {'
  WRITE(U,*)'node[label="",height=0.1,width=0.1,fontsize=1]'
  do X=1,N
    do Y=X+1,N
       if(A(X,Y)>0.0)then
         write(U,'(A,I4,A,I4,A)')'   ',X,' -- ',Y,' ;'
       endif
    enddo
  enddo
  write(U,*)'}'
  close(u)
  RETURN
 END SUBROUTINE
!
 SUBROUTINE write_gin(cell_0,xcryst,n_atoms,labels)
  ! Catlow potential
  IMPLICIT NONE
  INTEGER,          intent(in) :: n_atoms
  character (len=4),intent(in) :: labels(0:0,n_atoms,2)
  REAL,             intent(in) :: xcryst(0:3,1:n_atoms),cell_0(1:6)
  integer                      :: i,j
  OPEN(999,file='out.gin')
  WRITE(999,'(A)')'single conv nosym qok'
  WRITE(999,'(A)')'cell'
  WRITE(999,'(6(f9.5,1x))') (cell_0(j) , j=1,6)
  WRITE(999,'(A)')'fractional'
  do i=1,n_atoms
    WRITE(999,*)labels(0,i,1),labels(0,i,2),xcryst(1,i),xcryst(2,i),xcryst(3,i)
  enddo
  WRITE(999,'(A)')'species 6'
  WRITE(999,'(A)')'Na     core    1.000000'
  WRITE(999,'(A)')'Sr     core    2.000000'
  WRITE(999,'(A)')'Al     core    3.000000'
  WRITE(999,'(A)')'Si     core    4.000000'
  WRITE(999,'(A)')'O2     core    0.869020'
  WRITE(999,'(A)')'O2     shel   -2.869020'
  WRITE(999,'(A)')'buck Na core O2 shel  1226.84000  0.306500  0.000000 0.00 10.00'
  WRITE(999,'(A)')'buck Na core Na core  7895.40000  0.170900  0.000000 0.00 10.00'
  WRITE(999,'(A)')'buck Sr core O2 shel  1952.39000  0.336850 19.22000  0.00 10.00'
  WRITE(999,'(A)')'buck Al core O2 shel  1460.30000  0.299120  0.000000 0.00 10.00'
  WRITE(999,'(A)')'buck Si core O2 shel  1283.90700  0.320520 10.66158  0.00 10.00'
  WRITE(999,'(A)')'buck O2 shel O2 shel 22764.0000   0.149000 27.87900  0.00 12.00'
  WRITE(999,'(A)')'spring'
  WRITE(999,'(A)')'O2 74.920000'
  WRITE(999,'(A)')"three regular Si core O2 shel O2 shel 2.0972 109.47000 &"
  WRITE(999,'(A)')'  0.000  1.800  0.000  1.800  0.000  3.200'
  WRITE(999,'(A)')"three regular Al core O2 shel O2 shel 2.0972 109.47000 &"
  WRITE(999,'(A)')'  0.000  1.800  0.000  1.800  0.000  3.200'
  WRITE(999,'(A)')'cutp     16.00000    1.00000'
  WRITE(999,'(A)')'switch_min rfo   gnorm     0.001000'
  WRITE(999,'(A)')'dump every 1     out.res'
  WRITE(999,'(A)')'output cif out.cif'
  CLOSE(999)
  RETURN
 END SUBROUTINE write_gin
