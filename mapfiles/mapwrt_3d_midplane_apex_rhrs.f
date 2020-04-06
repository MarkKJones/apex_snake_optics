      program mapwrt_apex
      integer N,T
      PARAMETER (N=130)
        real posx(n,n,n),bfx(n,n,n)
        real posy(n,n,n),bfy(n,n,n)
        real posz(n,n,n),bfz(n,n,n)
      PARAMETER (T=3000000)
      REAL BOUTX(T),XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN
      REAL BOUTY(T),BOUTZ(T)
      INTEGER NDX,NDY,NDZ
      CHARACTER*1 A(3)
	character*80 FILE1,FILE2
      DATA INDFI/5/
      DATA XMAX/-10000./,XMIN/10000./,ZMIN/10000./
      DATA YMIN/10000./,YMAX/-10000./,ZMAX/-10000./
      DATA A/'N','N','N'/
      WRITE(6,*)' TYPE IN INPUT FILE NAME'
      READ(5,'(A80)')FILE1
      WRITE(6,*)' TYPE IN OUTPUT FILE NAME'
      READ(5,'(A80)')FILE2
	PRINT*,FILE1
      CALL READ_TOS_MAP(FILE1,NDX,NDY,ndz,posx,posy,posz
     >    ,bfx,bfy,bfz,ierr)
         write(*,*) 'finished reading tosca file'  
c       YMIN=posX(1,1,1)
       XMIN=posX(1,1,1)
       XMax=posX(ndx,ndy,ndz)
c
       XMIN=-posx(ndx,ndy,ndz)
       XMAX=posX(1,1,1)
c
       YMIN=-posy(ndx,ndy,ndz)
       YMax=posy(ndx,ndy,ndz)
       ZMAX=-posz(ndx,ndy,ndz)
       ZMIN=-posz(1,1,1)

      L=0
C     
c x = 0 to 40cm in 0.5cm steps (81)
c y= 0 to 15cm in 0.5cm steps  (31)
c z= -175 to 75cm in 2cm steps (126)
c  
      ndy = ndy+30
c      ndx = ndx+80
      INDEX=NDX*ndy*ndz
       write(*,*) ' index = ',index
      L=0
       DO 41 ix=ndx,1,-1   ! horizontal zsnake             
       DO 431 iz=1,ndz ! along beam ysnake
       DO 420 iy=1,ndy ! vertical xsnake
       L = L + 1
       iytemp=iy-30
       if (iy .le. 30) iytemp = 32 - iy
c       ixtemp=ix-80
c       if (ix .le. 80) ixtemp = 82 - iy
c       iytemp=iy
       ixtemp=ix
       BOUTX(L)=BfX(ixtemp,iytemp,iz)
       BOUTY(L)=BfY(ixtemp,iytemp,iz)
       BOUTZ(L)=BfZ(ixtemp,iytemp,iz)
c when vertical position is negative, then change horizontal field sign.
       if (iy .le. 30) BOUTZ(L)=-BOUTZ(L)
c       write(*,'(i10,6(f10.5,1x))') l,posx(ixtemp,iytemp,iz)
c     >,posy(ixtemp,iytemp,iz)
c     >,posz(ixtemp,iytemp,iz)
c     > ,bouty(l),boutz(l),boutx(l)
       IF(L.GT.INDEX)STOP 'TOO MANY POINTS'
 420   CONTINUE
 431   CONTINUE
 41   CONTINUE
      OPEN(4,FILE=FILE2,status='new')
      WRITE(4,99)INDFI
  99  FORMAT(I2)
      WRITE(4,100)(A(I),I=1,3)
 100  FORMAT(3A1)
      WRITE(4,*)NDY,YMIN,YMAX
      WRITE(4,*)NDZ,ZMIN,ZMAX
      WRITE(4,*)NDX,XMIN,XMAX
      WRITE(4,*)'(5G15.8)'
      WRITE(4,110)(BOUTY(I),BOUTZ(I),BOUTX(I),I=1,INDEX)
  110 FORMAT(5G15.8)
      Write(6,*)' Finished!'
      write(6,*) ' Index = ', INDEX,' L = ', L
      STOP
      END
C

	subroutine read_tos_map(file,nx,ny,nz
     >     ,posx,posy,posz,bfx,bfy,bfz,ierr)
C---------------------------------------------------------
C...	read a TOSCA MAP file in ascii formated file
C---------------------------------------------------------
C	input:
C		file C* file name (default type .out)
C	output:
C		nx,x	int,real	number of x values in array x
C		ny,y	int,real	number of y values in array y
C		nz,z	int,real	number of z values in array z
C		bx	real		x-comp of field (3D array)		
C		by	real		y-comp of field (3D array)		
C		bz	real		z-comp of field (3D array)		
C		ierr	int		.ne.0	fail
C---------------------------------------------------------
        integer N
	PARAMETER (N=130)
	REAL pot
        real posx(n,n,n),bfx(n,n,n)
        real posy(n,n,n),bfy(n,n,n)
        real posz(n,n,n),bfz(n,n,n)
	character*(*) file
        character*20 junk

	lun=99
	close(unit=lun)
	ierr=0
	open(unit=lun,file=file,status='old',form='formatted',
     &  	err=88)
	goto 80
88	ierr=1
	WRITE(6,100)FILE
 100  FORMAT('F - READ_TOS_MAP: Error openning file - ',80A1)
	return

80	READ(lun,*) NZ,Ny,Nx
	print*,nx,ny,nz
c
        do i=1,7
           READ(lun,*) junk
           enddo
c
	nn=0
	DO IX=1,Nx
	DO IY=1,Ny                 
	DO IZ=1,Nz
	read(lun,*) posx(ix,iy,iz),posy(ix,iy,iz),posz(ix,iy,iz)
     > ,bfx(ix,iy,iz),bfy(ix,iy,iz),bfz(ix,iy,iz)	
c   convert cm,gauss to mm,tesla. 	
		posx(ix,iy,iz)=posx(ix,iy,iz)*10.
		posy(ix,iy,iz)=posy(ix,iy,iz)*10.
		posz(ix,iy,iz)=posz(ix,iy,iz)*10.
		bfx(ix,iy,iz)=bfx(ix,iy,iz)*1.e-4
		bfy(ix,iy,iz)=bfy(ix,iy,iz)*1.e-4	
		bfz(ix,iy,iz)=bfz(ix,iy,iz)*1.e-4
	
c        write(*,*) posx(ix,iy,iz),posy(ix,iy,iz),posz(ix,iy,iz)
c     >,bfx(ix,iy,iz),bfy(ix,iy,iz),bfz(ix,iy,iz)
	ENDDO
	ENDDO
	ENDDO

	close(unit=lun)
	return
	END
