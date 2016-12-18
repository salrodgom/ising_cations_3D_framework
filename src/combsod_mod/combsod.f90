!*******************************************************************************
!    Copyright (c) 2014 Ricardo Grau-Crespo, Said Hamad
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SOD.  If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

    PROGRAM combsod 
       IMPLICIT NONE

       INTEGER,PARAMETER :: NSPMAX=10 , NATMAX=10000, NOPMAX=5000, NCELLMAX=1000
       INTEGER,PARAMETER :: NLINEAMAX=200 
       REAL,PARAMETER :: tol0=0.0001, kB=8.61734E-5
      
       INTEGER :: i,j,k,l,xx,yy,suma,t,ina,inb,inc,elei, factorial, sub
       INTEGER :: op1, nop1, op, nop, opsc, nopsc,ifound,op1new,nop1new,aux,opc,nopc
       INTEGER :: sp, nsp, spmap, nspmap, cumnatsp, cumnatspmap, sptarget, sptargetmap,ssp
       INTEGER :: at0, nat0, at1, nat1, nat1r, at, nat , atmap, natmap, at1r, at1i,attmp,att,mapno,FILER, MAPPER
       INTEGER :: na,nb,nc,nsubs,nsubsmap,atini,atfin,atinimap,atfinmap
       INTEGER :: pos,npos,count,ntc,ntcmax,nic,equivcount,iequiv,indcount
       LOGICAL :: found,foundnoind,mores
       INTEGER,DIMENSION(:), ALLOCATABLE:: newconf
       INTEGER,DIMENSION(2)   :: newshell, newshellmap
       INTEGER,DIMENSION(NOPMAX)   :: op1good
       INTEGER,DIMENSION(:), ALLOCATABLE:: degen
       INTEGER,DIMENSION(NSPMAX) :: natsp0, natsp1, natsp, natspmap, snatsp, ishell, ishellmap
       INTEGER,DIMENSION(NATMAX) :: spat0, spat1, spat, spatmap, spat1r
       REAL,DIMENSION(NATMAX,3) :: coords0, coords1, coords, coords1r, result, coordsmap
       REAL,DIMENSION(3) :: coordstemp
       REAL,DIMENSION(3,3) :: cellvector
       REAL,DIMENSION(NOPMAX,3,3) :: mgroup1, mgroup, mgroup1new
       REAL,DIMENSION(NOPMAX,3)   :: vgroup1, vgroup,vgroup1sca,vgroup1scb,vgroup1scc,vgroup1new
       REAL,DIMENSION(NCELLMAX,3)   :: vt
       REAL,DIMENSION(3) :: x, vcorr
       INTEGER,DIMENSION(NOPMAX,NATMAX)   :: fulleqmatrix,eqmatrixtarget
       INTEGER,DIMENSION(:,:), ALLOCATABLE:: conf,indconf,indconfmap,equivconf
       INTEGER,DIMENSION(:,:), ALLOCATABLE:: mapping
       INTEGER,DIMENSION(NATMAX)   :: as,test
       REAL :: a1,b1,c1,alpha,beta,gamma,a,b,c,cc, amap, bmap, cmap, alphamap, betamap, gammamap
       REAL :: prod, xs, ientropyguess, ientropy, perc
       CHARACTER,DIMENSION(NSPMAX) :: symbol*3, symbolmap*3, ssymbol*3
       CHARACTER,DIMENSION(2) :: newsymbol*3, newsymbolmap*3
       CHARACTER :: linea*85, title*15, tmptitle*15, runtitle*40
       CHARACTER :: sindcount*5,snumber*4


!!!!!! Input files

       OPEN (UNIT= 9,FILE="INSOD")
       OPEN (UNIT=12,FILE="SGO")

!!!!!! Output files

       OPEN (UNIT=25,FILE="coordinates.xyz")
       OPEN (UNIT=26,FILE="matrix")
       OPEN (UNIT=27,FILE="conf")
       OPEN (UNIT=30,FILE="OUTSOD")
       OPEN (UNIT=31,FILE="supercell")
       OPEN (UNIT=43,FILE="filer")
       OPEN (UNIT=46,FILE="operators")
       OPEN (UNIT=47,FILE="cSGO")
	

!
! DEFINITION OF VARIABLES:
!
! op                  Index for the operators in the supercell
! op1                 Index for the operators in the unit cell
! nop                 Total number of operators in the supercell
! nop1                Total number of operators in the unit cell
! sp                  Index for the species 
! nsp                 Total number of species
! sptarget            Number of the species to be substituted
! at0,at1,at          Indexes for the atoms in the assymetric unit, unit cell and supercell
! nat0,nat1,nat       Total numbers of atoms in the assymetric unit, unit cell and supercell
! at1r,nat1r	      Idem for the atoms in the redundant cell (with repeated positions)
! atini,atfin	      Initial and final atom indexes of the species to be substituted
! npos		      Number of atoms of the target species 
! conf		      List of all configurations (each configuration is a list of the substituted positions)
! count		      Index for the configurations (conf)     
! ntc                 Total number of configurations in conf (count=1,ntc)
! indconf             List of independent configurations 
! indcount            Index for the independent configurations
! nic                 Total number of independent configurations in indconf (indcount=1,nic)
! equivconf           Temporary list containing the equivalent configurations at every step of the algorithm
! equivcount          Index for the equivalent configurations (equivconf)     
! tol0 		      General tolerance
! tol1 		      Tolerance used for correcting the x-FLOOR(x) function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the input file with the uc space group information: SGO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       READ (12,*)
       READ (12,*) op1
       do while(op1>0)
	do i=1,3
          READ (12,*)(mgroup1(op1,i,j),j=1,3),vgroup1(op1,i)
	enddo
        nop1=op1
        READ (12,*) op1
       enddo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the INSOD file 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       READ (9,*)

       READ (9,*) runtitle
       READ (9,*)
       READ (9,*)
       READ (9,*) a1,b1,c1,alpha,beta,gamma
       READ (9,*)
       READ (9,*)
       READ (9,*) nsp
       READ (9,*)
       READ (9,*)
       READ (9,*)(symbol(sp),sp=1,nsp)
       READ (9,*)
       READ (9,*)
       READ (9,*)(natsp0(sp),sp=1,nsp)
       READ (9,*)
       READ (9,*)

       nat0=0
       do sp=1,nsp
         nat0 = nat0 + natsp0(sp)
       enddo

       do at0=1,nat0
         READ (9,*)(coords0(at0,i),i=1,3)
       enddo
       READ (9,*)
       READ (9,*)
       READ (9,*) na,nb,nc
       READ (9,*)
       READ (9,*)
       READ (9,*) sptarget
       READ (9,*)
       READ (9,*)
       READ (9,*) nsubs
       READ (9,*)
       READ (9,*)
       READ (9,*)
       READ (9,*) (newsymbol(i),i=1,2)
       READ (9,*)
       READ (9,*)
       READ (9,*)
       READ (9,*)
       READ (9,*)
       READ (9,*) FILER, MAPPER
       READ (9,*)

       if (FILER<10) then
       	READ (9,*)
       	READ (9,*)
       	READ (9,*) (ishell(sp),sp=1,nsp)
       	READ (9,*)
       	READ (9,*) (newshell(i),i=1,2)
       endif


!cccccccccccccccccccccccccccccccccccc
! Generating spat0 array
!cccccccccccccccccccccccccccccccccccc

       do at0=1, natsp0(1)
          spat0(at0)=1
          cumnatsp=natsp0(1)
       enddo 
       do sp=2, nsp
          do at0=cumnatsp+1, cumnatsp+natsp0(sp)
             spat0(at0)=sp
          enddo
          cumnatsp=cumnatsp+natsp0(sp) 
       enddo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      First create the redundant unit cell from the asymmetric unit 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
       WRITE (*,*) "" 
       WRITE (*,*) "Generating the supercell..." 
       WRITE (*,*) "" 

       at1r=0
       do at0=1,nat0
          do op1=1,nop1
    	     at1r=at1r+1
             coords1r(at1r,1:3)= MATMUL(mgroup1(op1,1:3,1:3),coords0(at0,1:3))+vgroup1(op1,1:3)


!  For some reason the matrix multiplication above has not worked in some machines. 
!  In that case it can be expanded like this: 
!                                 
!             coords1r(at1r,1)= mgroup1(op1,1,1)*coords0(at0,1)&
!                              +mgroup1(op1,1,2)*coords0(at0,2)&
!                              +mgroup1(op1,1,3)*coords0(at0,3)&
!                              +vgroup1(op1,1)
!             coords1r(at1r,2)= mgroup1(op1,2,1)*coords0(at0,1)&
!                              +mgroup1(op1,2,2)*coords0(at0,2)&
!                              +mgroup1(op1,2,3)*coords0(at0,3)&
!                              +vgroup1(op1,2)
!             coords1r(at1r,3)= mgroup1(op1,3,1)*coords0(at0,1)&
!                              +mgroup1(op1,3,2)*coords0(at0,2)&
!                              +mgroup1(op1,3,3)*coords0(at0,3)&
!                              +vgroup1(op1,3)

                 coords1r(at1r,1)=cc(coords1r(at1r,1))
	         coords1r(at1r,2)=cc(coords1r(at1r,2))
	         coords1r(at1r,3)=cc(coords1r(at1r,3))
             spat1r(at1r)=spat0(at0)
          enddo
       enddo
       nat1r=at1r

 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Get rid of the redundant atoms, to create the unit cell 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	coords1(1,:)=coords1r(1,:)
	at1r=1
	at1=1
	coords1(1,:)=coords1r(1,:)
	spat1(1)=spat1r(1)
	do at1r=2,nat1r
           found=.FALSE.
	   do at1i=1,at1
              prod=DOT_PRODUCT(coords1r(at1r,:)-coords1(at1i,:),&
                            coords1r(at1r,:)-coords1(at1i,:)) 
	      if (prod.le.tol0) found=.TRUE.
	   enddo
           if (found.eqv..FALSE.) then    
	       at1=at1+1
	       coords1(at1,:)=coords1r(at1r,:)
		spat1(at1)=spat1r(at1r)
	   endif
        enddo  
        nat1=at1

	do sp=1,nsp
	   natsp1(sp)=0
	enddo
	do at1=1,nat1
	   natsp1(spat1(at1))=natsp1(spat1(at1))+1
	enddo

	natsp(:)=na*nb*nc*natsp1(:)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Create the traslation vectors of the supercell
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	nopsc=nop1
	t=0
	do ina=0,na-1
	  do inb=0,nb-1
	    do inc=0,nc-1
		t=t+1
		vt(t,1)=real(ina)/real(na)
		vt(t,2)=real(inb)/real(nb)
		vt(t,3)=real(inc)/real(nc)
	    enddo 
	  enddo
	enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Create the traslation vectors of the supercell
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do op1=1,nop1
	   op1good(op1)=1
	enddo

	if (na.eq.nb) then 
	   if (na.eq.nc) then 
	      goto 255
	      else
              do op1=1,nop1
                 do i=1,3
                    if ((mgroup1(op1,1,3).ne.0).OR.(mgroup1(op1,3,1).ne.0).OR.&
		        (mgroup1(op1,2,3).ne.0).OR.(mgroup1(op1,3,2).ne.0)) then
	            op1good(op1)=0
		    endif
                 enddo
              enddo
	   endif
	else
	   if (na.eq.nc) then 
              do op1=1,nop1
                 do i=1,3
                    if ((mgroup1(op1,1,2).ne.0).OR.(mgroup1(op1,2,1).ne.0).OR.&
		        (mgroup1(op1,3,2).ne.0).OR.(mgroup1(op1,2,3).ne.0)) then
	            op1good(op1)=0
	            endif
                 enddo
              enddo
	   endif
	   if (nb.eq.nc) then 
              do op1=1,nop1
                 do i=1,3
                    if ((mgroup1(op1,2,1).ne.0).OR.(mgroup1(op1,1,2).ne.0).OR.&
		        (mgroup1(op1,3,1).ne.0).OR.(mgroup1(op1,1,3).ne.0)) then
	            op1good(op1)=0
	            endif
                 enddo
              enddo
	   endif
	   if ((nb.ne.nc).AND.(na.ne.nc)) then 
              do op1=1,nop1
                 do i=1,3
                    if ((mgroup1(op1,2,1).ne.0).OR.(mgroup1(op1,1,2).ne.0).OR.&
		        (mgroup1(op1,3,1).ne.0).OR.(mgroup1(op1,1,3).ne.0).OR.&
		        (mgroup1(op1,3,2).ne.0).OR.(mgroup1(op1,2,3).ne.0)) then
	            op1good(op1)=0
	            endif
                 enddo
              enddo
	   endif
	endif      
 255       continue

	op1new=0
        do op1=1,nop1
	   if (op1good(op1).eq.1) then
 	       op1new=op1new+1
	       mgroup1new(op1new,:,:)=mgroup1(op1,:,:)
	       vgroup1new(op1new,:)=vgroup1(op1,:)
	   endif
	enddo
	nop1new=op1new

	nop1=nop1new
	do op1=1,nop1
           mgroup1(op1,:,:)=mgroup1new(op1,:,:)
           vgroup1(op1,:)=vgroup1new(op1,:)
	enddo



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Generate the supercell coordinates, by applying these vectors
!      to the unit cell. Also calculate the supercell parameters.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	a=na*a1
	b=nb*b1
	c=nc*c1


	at=0
	do at1=1,nat1
	   do t=1,na*nb*nc
		at=at+1
		coords(at,1)=vt(t,1)+coords1(at1,1)/na
		coords(at,2)=vt(t,2)+coords1(at1,2)/nb
		coords(at,3)=vt(t,3)+coords1(at1,3)/nc
                spat(at)=spat1(at1)
	   enddo
	enddo
	nat=at

	at=0
	do at1=1,nat1
	   do t=1,na*nb*nc
		at=at+1
 100    format(3(f11.7,2x))
	   enddo
	enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Write a file with the fractional coordinates of the supercell
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        WRITE(31,111) a,b,c,alpha,beta,gamma 
 111    format(6(f10.4,2x))
        WRITE(31,*) natsp(1:nsp)
        do at=1,nat
                 WRITE(31,110) symbol(spat(at)),coords(at,1),coords(at,2),coords(at,3)
 110    format(a3,2x,3(f11.7,2x))
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the supercell operators
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	op=0
	do op1=1,nop1
	   do t=1,na*nb*nc
		op=op+1
		mgroup(op,:,:)=mgroup1(op1,:,:)		
		vgroup(op,1)=vgroup1(op1,1)/na+vt(t,1)
		vgroup(op,2)=vgroup1(op1,2)/nb+vt(t,2)
		vgroup(op,3)=vgroup1(op1,3)/nc+vt(t,3)
!		WRITE(*,*) op,op1,t
	   enddo
	enddo

	nop=op

	WRITE(46,*) nop
	do op=1,nop
	   WRITE(46,*) "Operator number ",op
	   do i=1,3
	   WRITE(46,*) mgroup(op,i,1:3),vgroup(op,i)
 163   format(4(f10.8))
	   enddo
	enddo
       
       WRITE (*,*) "       Composition of the original supercell (parent structure):" 
       do sp=1,nsp
          WRITE (*,*) "                                                         ", symbol(sp),natsp(sp)
       enddo     
       WRITE (*,*) " "
       WRITE (*,*) "       Number of symmetry operators in the supercell:       ", nop 
       WRITE (*,*) " " 
       WRITE (*,*) "       Composition of the substituted supercell:" 
       do sp=1,nsp
          if (sp==sptarget) then
             WRITE (*,*) "                                                         ", newsymbol(1),nsubs
             WRITE (*,*) "                                                         ", newsymbol(2),natsp(sp)-nsubs
          else   
             WRITE (*,*) "                                                         ", symbol(sp),natsp(sp)
          endif
       enddo

       WRITE (*,*) ""
       WRITE (*,*) ""
       WRITE (*,*) "Generating the complete configurational space..." 
! 161   format(a48)
       WRITE (*,*) " " 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Create Full Equivalence Matrix
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	coordstemp(:)=MATMUL(mgroup(1,:,:),coords(1,:))+vgroup(1,:)

	do op=1,nop
	   do at=1,nat
	        coordstemp(:)=MATMUL(mgroup(op,:,:),coords(at,:))+vgroup(op,:)
                coordstemp(1)=cc(coordstemp(1))
                coordstemp(2)=cc(coordstemp(2))
                coordstemp(3)=cc(coordstemp(3))
	        found=.FALSE.	
		attmp=0
		do while((found.eqv..FALSE.).AND.(attmp.lt.nat))
		    attmp=attmp+1
                    if((abs(coordstemp(1)-coords(attmp,1)).lt.tol0).AND.&
                       (abs(coordstemp(2)-coords(attmp,2)).lt.tol0).AND.&
                       (abs(coordstemp(3)-coords(attmp,3)).lt.tol0))& 
                       found=.TRUE.
	    enddo
		If( (attmp.eq.nat).And.(.Not.found) ) Then 
             Write(*,*) "Error!!! Operator",& 
             	        op,"applied on atom",at,"does not produce another atom in the list!!!"
             attmp=0
        Endif

        fulleqmatrix(op,at)=attmp

	   enddo
	enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Calculate the initial and final nat of the target species
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	if (sptarget.eq.1) then 
	atini=1

	else
	   atini=1
	   do sp=1,sptarget-1
	      atini=atini+natsp(sp)
	   enddo
	endif      

	   atfin=atini+natsp(sptarget)-1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Obtain the Equivalence Matrix for the target species from the Full Equivalence Matrix
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	npos=natsp(sptarget)
	do op=1,nop
	att=0
	do at=atini,atfin
	   att=att+1
           eqmatrixtarget(op,att)=fulleqmatrix(op,at)-atini+1
		WRITE(26,*)op,att,eqmatrixtarget(op,att)
	enddo
 162   format(50(i3))
	enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Generate the list of configurations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
        xs=REAL(nsubs)/REAL(npos)
        ientropyguess=-1*npos*kB*(xs*log(xs)+(1-xs)*log(1-xs))
!              that was using Stirling's formula,
!              real ientropy= factorial(npos)/(factorial(npos-nsubs)*factorial(nsubs))	
        ntcmax=abs(exp(ientropyguess/kB))
        WRITE(*,*) "       Maximum entropy for this composition:              ", ientropyguess, " eV/K"
        write(*,*) ntcmax
        WRITE(*,*) " " 
        

!!!!!!!!Allocating array sizes

        ALLOCATE(degen(1:ntcmax))
        ALLOCATE(newconf(1:nsubs))
        ALLOCATE(conf(1:ntcmax,1:nsubs))
        ALLOCATE(indconf(1:ntcmax,1:nsubs))
        ALLOCATE(equivconf(1:ntcmax,1:nsubs))

   	mores=.false.
	as(:)=1
	    
	count=1
	   call ksubset(npos,nsubs,as,mores)
           nsubs=abs(nsubs)
	   conf(1,1:nsubs)=as(1:nsubs)
	do while(mores)
	   call ksubset(npos,nsubs,as,mores)
           count=count+1
	   conf(count,1:nsubs)=as(1:nsubs)
	enddo

        ntc=count
        ientropy=kB*log(ntc*1.0)

        WRITE(*,*) "       Maximum entropy for this composition and supercell:", ientropy, " eV/K"
        WRITE(*,*) " " 
        WRITE(*,*) "       Total number of configurations in the supercell:     ", ntc 
        WRITE(*,*) " " 
        deallocate( degen )
        deallocate( newconf )
        deallocate( conf )
        deallocate( indconf )
        deallocate( equivconf )
    END PROGRAM combsod 
