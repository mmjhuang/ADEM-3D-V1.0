!These works are writtten in Hong Kong Polytechnic University from Oct.2008-Oct.2009
!In order to make the flow chart clearly and add the fluid phase, it is rewrittn in NTNU
!Author: Yrjo^:) Jun HUANG	
!version Swendish Meatball (V.3.22)    
!11st Nov 2013. After coming back from Chengdu

!30 June, 2014, Canada
!10 Sept.2014, Montreal, Canada, V4.0 for ellipse particle
!20 Nov.2014, Shanghai, "timepara" is added. simulations can be restared from the previous data saved
!   subroutine delepart is divied into two subroutines, delepart and numfull
!01 Jan 2015, Shanghai, added subroutine "ghostwall", one wall only. The rest walls would be added by DONG SU
!10 Feb 2015, Change the nb from table index to chain index, 40% faster than before
!20 Sept.2015, Y.J aaded the rest four walls using subroutine ghostwall. 
!26 Sept.2015, Y.J added subroutine coat to avoid the particle overlap with neighboring particles. 
!in the previous version 20150210beta, there is a bug:
!if the number of particles created by particlecreate is fewer or more, sometime shows an error. 
!this version solved this problem, in subrotine: "for i =1 to 3000" 3000 was changed to numfull.
! 2016 to write the ellipse codes to spherical particle
! 23,Dec.2015 finished this version, versio:Bei
! 12, Dec 2016 Add surface friciton to shpere particle
! 13 Dec 2016 Add baby and reset are finished for complex 3D particles
! 15 April 2018 change int to anint
! 10 OCT 2018, ADD A NEW MATRIX ABOUT ALL PARTICLES, smtrx
! 16 OCT Reported by Fei that tangenrial force due to rotation is not added in TDM_model subroutine.
! 24 OCT 2018, replace the radius of small particle using a reference radius. 
! 30 OCT 2018, add the subroutins for plates and cylinder. 
! 30 Nov 2018, Try OpenMP   using command gfortran -gopenmp SphereMP.f90
! If it shows Segmentation fault (core dumped), try type "ulimit -s unlimited ;export KMP_STACKSIZE=20480000000" in terminal  for numpfull=300000
! Change 20480000 to a big enough number.
! 5 DEC 2018  try to expand all matrix matul functions to speed


	module face		 
	integer, allocatable:: nb(:,:)
	integer, allocatable:: mv(:,:,:,:)  	   
        end module	
!	=================== 
	program ApenDEM
!	=================== 
!The main code	
!input:   nump--maximum number of particles allowed; 
!         numcol--number of columns in array stma;
!         stma--main array;  xlower,xupper, ylower, yupper--compuational domain;  
!         tp0--initial time; tpfinal--finial time; rhop density;
!         tpout--output time gap;  fricp--surface roughness;f
!         young--Young's modulus;  pois-- poision's ratio;
!         e--coefficient of resitution; 
!         youngr--relative young's modulus; 
!internal:istep--index for time steps between two screen outputs
!         iip--index for time steps between two files saved to harddisk 
!         kedm--time step of EDM / time step of TDM
!         iiv--time step for free particle / time step of TDM
!         numfull--number of rows in arrayment stma is not zero
!         timepara--parameters about time, including tp0, tpout and dtp0
!         nw -- particles close to each wall
!               (eg. nw(5000,2) means for 2 walls & each wall has 5000 neighboring particles)
       
        use face 	
	implicit double precision (a-h, o-z)                             
	parameter (nump=30000)                  
	parameter (numcol=35)                 
	dimension stma(numcol, nump) 
        dimension smtrx(3,3, nump)  		
	dimension lsg(3,nump)	  
        dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)	
	dimension stma2(numcol, nump)
        dimension stma0(numcol, nump)
        dimension timepara(5)
        dimension nwS(18000,1)
        dimension nwp(18000,2)
        dimension nwC(18000,1)
        dimension walls(4,1)
        dimension wallp(4,2)
        dimension wallc(4,1)
        nb2d=nump*600
	allocate (nb(2, nb2d))  

	tpfinal=5.
	tpout=1.e-2  
	pi=3.1415927   
	istep=1    ! index for screen         
	                    
	open (10, file='file000398.dat')
	read (10,110) stma2
 110    format(35e16.8)
	close(10)
        stma=stma2
                           
	open (13, file='time000398.dat')
  	read (13,113) timepara
 113    format(5e18.10)
	close(13) 
         
        tp0=timepara(3)   
        print *,'numfull 001=', timepara(5)        

        smtrx=0.                                                                                   
                                                    
        if (tp0<1.d-20)  then
         iip=1
          
!         call delepart (stma, numcol, nump, numfull)              
!         call addbaby (stma, numcol, nump, numfull, stma0)
           call addbaby_dice (stma, numcol, nump, numfull, stma0)
!         call addcoat (stma, numcol, nump, numfull, ratio)
                
          call getnumfull (stma, numcol, nump, numfull)                          
          call intstat(nump, stma, numcol)   
                   
          call field (stma, numcol, nump, tp0) 
                                
          call dtpint (stma, youngr, vmax, numcol, nump, radmax, dtp0, pmassminr, radmin, iedm, dtfree, numfull)	 

          call xyminmax (stma, vmax, numcol, nump, xlower, xupper, ylower, yupper, zlower,  & 
 zupper, cellp, radmax, numfull)           
                
       	call searchsize (xlower, xupper, ylower, yupper, zlower, zupper, radmax,  mvx, mvy, mvz, cellp)    
             
	call search (stma, numcol, nump, numfull, xlower, xupper, ylower, yupper, zlower, zupper, cellp, lsg)  

 !       call BoundaryCellFix_sphere (stma, xlower, xupper, ylower, yupper, zlower, zupper, &
 ! radmax,  mvx, mvy, mvz, cellp, nwS, walls, numcol, nump, tp0)  

 !       CALL BoundaryCellFix_cylinder (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax, &
 ! mvx, mvy, mvz, cellp, nwc, wallc, numcol, nump, tp0) 

       call Search_sphere (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nws, walls, numcol, nump, tp0)

       CALL Search_cylinder (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax, &
   mvx, mvy, mvz, cellp, nwc, wallc, numcol, nump, tp0)

!!       call BoundaryCellFix_plate (stma, xlower, xupper, ylower, yupper, zlower, zupper, &
!  radmax,  mvx, mvy, mvz, cellp, nwp, wallp, numcol, nump, tp0)  

      call  Search_Plate (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nwp, wallp, numcol, nump, tp0)
              
          iip=1
          iiv=1

        endif 
                            
        if (tp0>=1.d-20) then 
        stma2=stma  

	open (18, file='file_init.dat')
  	read (18,118) stma0
 118    format(35e16.8)
	close(18)
                
        call reset (stma, numcol, nump, tp0, dtp0, stma0, smtrx) 
          
 !       call addcoat (stma, numcol, nump, numfull, ratio)
        iip=anint(timepara(4))+1
        tpout=timepara(2) 
        print *, 'iip=', iip
        numfull = anint(timepara(5))

        call dtpint (stma, youngr, vmax, numcol, nump, radmax, dtp0, pmassminr, radmin, iedm, dtfree,numfull)
        call xyminmax (stma, vmax, numcol, nump, xlower, xupper, ylower, yupper, zlower,  & 
 zupper, cellp, radmax, numfull)      
   	call searchsize (xlower, xupper, ylower, yupper, zlower, zupper, radmax,  mvx, mvy, mvz, cellp)   
        call search (stma, numcol, nump, numfull, xlower, xupper, ylower, yupper, zlower, zupper, cellp, lsg)
       
 !       CALL BoundaryCellFix_sphere (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax, &
 ! mvx, mvy, mvz, cellp, nwS, walls, numcol, nump, tp0)

      !  CALL BoundaryCellFix_cylinder (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax, &
  ! mvx, mvy, mvz, cellp, nwc, wallc, numcol, nump, tp0)

        call Search_sphere (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nws, walls, numcol, nump, tp0)

        CALL Search_cylinder (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax, &
   mvx, mvy, mvz, cellp, nwc, wallc, numcol, nump, tp0)

!        call BoundaryCellFix_plate (stma, xlower, xupper, ylower, yupper, zlower, zupper, &
! radmax,  mvx, mvy, mvz, cellp, nwp, wallp, numcol, nump, tp0) 

      call  Search_Plate (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nwp, wallp, numcol, nump, tp0)

        iiv=1000000 
        endif     

100	do while (tp0<tpfinal) 
                
	 if (iiv>=200) then	                        ! change from 500 to 400 to 20

 !       call xyminmax (stma, vmax, numcol, nump, xlower, xupper, ylower, yupper, zlower,  & 
 !zupper, cellp, radmax, numfull)      

 !  	call searchsize (xlower, xupper, ylower, yupper, zlower, zupper, radmax,  mvx, mvy, mvz, cellp) 

	 call search (stma, numcol, nump, numfull, xlower, xupper, ylower, yupper, zlower, zupper, cellp, lsg)	               
                         
 !       call BoundaryCellFix_sphere (stma, xlower, xupper, ylower, yupper, zlower, zupper, &
 ! radmax,  mvx, mvy, mvz, cellp, nwS, walls, numcol, nump, tp0)  

  !       CALL BoundaryCellFix_cylinder (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax, &
  ! mvx, mvy, mvz, cellp, nwc, wallc, numcol, nump, tp0)

        call Search_sphere (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nws, walls, numcol, nump, tp0)

        CALL Search_cylinder (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax, &
   mvx, mvy, mvz, cellp, nwc, wallc, numcol, nump, tp0)

        call  Search_Plate (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nwp, wallp, numcol, nump, tp0) 
 
!        call BoundaryCellFix_plate (stma, xlower, xupper, ylower, yupper, zlower, zupper, &
! radmax,  mvx, mvy, mvz, cellp, nwp, wallp, numcol, nump, tp0) 
                  
	 call dtp1 (stma, youngr, vmax, numcol, nump, dtp0, pmassminr, radmin, iedm, kedm, dtedm, tp0, tpminus1, dtfree, numfull)   
       
	   iiv=0
	   else
	   iiv=iiv+1
        endif     
             
	call reset (stma, numcol, nump, tp0, dtp0, stma0, smtrx) 
           
        call field (stma, numcol, nump, tp0)

        call binary (stma, numcol, nump, pai, paj,  alpha,  youngr, dist, pi,  dtp0, kedm, dtedm, deltamax, ratio) 

!   	call ghostwall_sphere (stma, numcol, nump, pai, paj, pak, pal,  alpha, youngr, dist, pi, & 
! dtp0, kedm, dtedm, nwS, walls, tp0)    

 	call ghostwall_cylinder (stma, numcol, nump, pai, paj, pak, pal,  alpha, youngr, dist, pi, & 
  dtp0, kedm, dtedm, nwc, wallc, tp0)    

       call ghostwall_plate (stma, numcol, nump, pai, paj, pak, pal, alpha,  youngr,  &
  dist, pi,  dtp0, kedm, dtedm, nwp, wallp, tp0)      
            	   
       call velocityPE(stma, numcol, nump, dtp0, kedm, dtedm, dtfree, iiv, numfull, smtrx, pi)           
	
            tp0=tp0+dtp0            
 	call outputp (stma, nump, numcol, dtp0, tpout, istep, tp0, vmax, iip, dudt, pressF, partF, FricF, &
numfull, deltamax, ratio, stma2, pi, wallc)           	
	enddo    
 
        print *, 'finish'  
	end	


!!=====================================
    subroutine GetNumfull (stma, numcol, nump, numfull)
!=====================================  
!     If the code is not start from the inital state (baby particles were added)
!     using this subroutine to instead delpartII
!     Empty rows should be removed from the array stma
!     in the subroutine delepart already.

      use face
      implicit double precision (a-h, o-z)        
      dimension stma(numcol, nump)

      IPN=0    ! REPORTED BY Fei Wang 20180404, that the number of row of STMA should > nump. We addded ipn

      do i=1, nump        
        if   (anint(stma(28,i))==0) then
        numfull = i-1 
        ipn=1                      ! REPORTED BY Fei Wang 20180404
        goto 5030 
        endif
      enddo 

      if (ipn==0) numfull=nump     ! REPORTED BY Fei Wang 20180404 

5030  continue 
        
      end

!!=====================================
    subroutine addbaby (stma, numcol, nump, numfull,stma0)
!=====================================  
!     make the mass of all ghost partile infinite
!     move all empty rows to the end of array stma
!     input:stma; 
!     output:stma

      use face
      implicit double precision (a-h, o-z)        
      dimension stma(numcol, nump), stma0(numcol, nump)
          
 !    real :: a(3)=(/0.9,0.75,0.6/)
 !     real :: b(3)=(/0.98,1.88,2.60/)
      real :: a(2)=(/0.8, 0.6/)
      real :: b(2)=(/0.4, 0.7/) 

      do i=1, nump        
        if   (anint(stma(18,i))/=0 .or. anint(stma(28,i))==0 ) then
        numfull = i-1 
        goto 5040 
        endif
      enddo 
5040  continue  

      numfull2=numfull 

        do j=1, numfull  
       if (anint(stma(28,j))==3) then 
        stma(3,j)=stma(3,j)  *1.4   
        stma(29,j)=stma(29,j)  *4. !  added 2018-06-30 
        stma(30,j)=stma(29,j)  *4.
        stma(31,j)=stma(29,j)  *4.
        stma(19,j)=stma(19,j)* 1.1
        stma(20,j)=stma(3,j)

       goto 12312
       do 5050 m=1,2,1
        do i=-1,1,2 
     	stma(1,numfull2+1) = a(m)*stma(1,j) 
     	stma(2,numfull2+1) = stma(2,j)     	
	stma(3,numfull2+1) = stma(3,j)             
     	stma(4,numfull2+1) = stma(4,j)
        stma(5,numfull2+1) = stma(5,j) + i* b(m)* stma(1,j) * cos(stma(32,j)) 

        stma(6,numfull2+1) = stma(6,j)  +  i* b(m)* stma(1,j)   * sin(stma(32,j))
     	stma(7,numfull2+1) = stma(7,j)
     	stma(8,numfull2+1) = stma(8,j)
     	stma(9,numfull2+1) = stma(9,j)
     	stma(10,numfull2+1) = stma(10,j)
     	stma(11,numfull2+1) = stma(11,j)
     	stma(12,numfull2+1) = stma(12,j)         
     	stma(13,numfull2+1) = abs(i* b(m) * stma(1,j))
        stma(14,numfull2+1) = stma(14,j)    ! added 2018-06-30
     	 do k=15,17
       	 stma(k,numfull2+1) = stma(k,j)
       	 enddo
     	stma(18,numfull2+1)=j
        stma(19,numfull2+1)=stma(5,numfull2+1)-stma(5,j)
        stma(20,numfull2+1)=stma(6,numfull2+1)-stma(6,j)
        stma(21,numfull2+1)=stma(7,numfull2+1)-stma(7,j)

        if (stma(19,numfull2+1)>0. .and. stma(19,numfull2+1)<10d-15)    stma(19,numfull2+1)=10d-15
        if (stma(20,numfull2+1)>0. .and. stma(20,numfull2+1)<10d-15)    stma(20,numfull2+1)=10d-15
        if (stma(21,numfull2+1)>0. .and. stma(21,numfull2+1)<10d-15)    stma(21,numfull2+1)=10d-15
        if (stma(19,numfull2+1)<0. .and. stma(19,numfull2+1)>-10d-15)   stma(19,numfull2+1)=-10d-15
        if (stma(20,numfull2+1)<0. .and. stma(20,numfull2+1)>-10d-15)   stma(20,numfull2+1)=-10d-15
        if (stma(21,numfull2+1)<0. .and. stma(21,numfull2+1)>-10d-15)   stma(21,numfull2+1)=-10d-15

    	 do k=22,27
       	 stma(k,numfull2+1) = stma(k,j)
       	 enddo
     	stma(28,numfull2+1)=3     

     	stma(29,numfull2+1) = stma(29,j)
     	stma(30,numfull2+1) = stma(30,j)
        stma(31,numfull2+1) = stma(31,j)  


    	 do k=29,34
       	 stma(k,numfull2+1) = stma(k,j)
       	 enddo
       
     	stma(28,numfull2+1)=3 
 
    	numfull2=numfull2+1
       enddo    	
5050       enddo       
12312 continue
       endif   
   
      enddo 
 
          print *, 'total number', numfull2, stma(5,141)
          stma0=stma
        
	 open(77, file='file_init.dat')
	 write (77, 177) stma0
	 177 format(35e16.8)
	 close(77)

      end

!
!!=====================================
    subroutine addcoat (stma, numcol, nump, numfull, ratio)
!=====================================  
!     the particle diameters are times by a ratio, ratio
!     ratio -- is this ratio

      use face
      implicit double precision (a-h, o-z)        
      dimension stma(numcol, nump)

      ratio =1.005
      do i=1, nump        
       if (anint(stma(28,i))==0) goto 1488
       stma(1,i)=stma(1,i) * ratio
      enddo    
1488  continue    
      end

!
!=====================================
      subroutine reset (stma, numcol, nump, tp0, dtp0, stma0, smtrx)
!=====================================
!     before each time step, make the accleration zero
!     input:stma
!     output:stma

      use face
      implicit double precision (a-h, o-z)  
      dimension stma(numcol, nump), stma0(numcol, nump), smtrx(3,3,nump)
        real:: pmx(3,3), fmx(3,3), trm(3,3), trmt(3,3) 
        real:: omx1(3,1), omx2(3,1), ammx(3,3)  


       ! call omp_set_num_threads(8)
       ! !$OMP  PARALLEL DO   &
       ! !$OMP  DEFAULT(SHARED)   &
       ! !$OMP  private(  i, j, pmx,fmx,trm, trmt,  omx1, omx2, ammx, d2h, kk) &
       ! !$OMP  private( cosphi, cospsi, costhet, sinphi, sinpsi, sinthet, psi, phi, thet)  
           
      do 1425 i=1,nump	 
        if (anint(stma(28,i))/=0) then
	stma(15,i)=0. 	   
        stma(16,i)=0.
	stma(17,i)=0.
	stma(25,i)=0. 	   
        stma(26,i)=0. 
	stma(27,i)=0. 

        if (anint(stma(18,i))==0) then   
          psi=stma(34,i)    
          phi=stma(33,i)    
	  thet=stma(32,i)  

        cosphi=cos(phi)
        cospsi=cos(psi)
        costhet=cos(thet)
        sinphi=sin(phi)
        sinpsi=sin(psi) 
        sinthet=sin(thet)       
 
	!trm(1,1)=cos(phi)*cos(psi);   trm(1,2)=-cos(phi)*sin(psi); trm(1,3)=sin(phi)
	!trm(2,1)=cos(thet)*sin(psi)+sin(thet)*sin(phi)*cos(psi)
	!trm(2,2)=cos(thet)*cos(psi)-sin(thet)*sin(phi)*sin(psi)
	!trm(2,3)=-sin(thet)*cos(phi)
	!trm(3,1)=sin(thet)*sin(psi)-cos(thet)*sin(phi)*cos(psi)
	!trm(3,2)=sin(thet)*cos(psi)+cos(thet)*sin(phi)*sin(psi)
	!trm(3,3)=cos(thet)*cos(phi)

	trm(1,1)=cosphi*cospsi;   
        trm(1,2)=-cosphi*sinpsi; 
        trm(1,3)=sinphi
	trm(2,1)=costhet*sinpsi+sinthet*sinphi*cospsi
	trm(2,2)=costhet*cospsi-sinthet*sinphi*sinpsi
	trm(2,3)=-sinthet*cosphi
	trm(3,1)=sinthet*sinpsi-costhet*sinphi*cospsi
	trm(3,2)=sinthet*cospsi+costhet*sinphi*sinpsi
	trm(3,3)=costhet*cosphi

	!trmt(1,1)=trm(1,1);  trmt(1,2)=trm(2,1);   trmt(1,3)=trm(3,1) 
	!trmt(2,1)=trm(1,2);  trmt(2,2)=trm(2,2);   trmt(2,3)=trm(3,2) 
	!trmt(3,1)=trm(1,3);  trmt(3,2)=trm(2,3);   trmt(3,3)=trm(3,3)  

	smtrx(1,1,i)=trm(1,1);  smtrx(1,2,i)=trm(1,2);   smtrx(1,3,i)=trm(1,3) 
	smtrx(2,1,i)=trm(2,1);  smtrx(2,2,i)=trm(2,2);   smtrx(2,3,i)=trm(2,3) 
	smtrx(3,1,i)=trm(3,1);  smtrx(3,2,i)=trm(3,2);   smtrx(3,3,i)=trm(3,3) 
        endif

        if (anint(stma(18,i))>0) then
        j=anint(stma(18,i)) 

	omx1(1,1)=stma0(19,i);  omx1(2,1)=stma0(20,i); omx1(3,1)=stma0(21,i)      ! corrected 2018-06-30 

	ammx(1,1)=stma0(29,i);  ammx(1,2)=0;   ammx(1,3)=0;
	ammx(2,1)=0;  ammx(2,2)=stma0(30,i); ammx(2,3)=0;
	ammx(3,1)=0; ammx(3,2)=0; ammx(3,3)=stma0(31,i)

	smtrx(1,1,i)=smtrx(1,1,j);  smtrx(1,2,i)=smtrx(1,2,j);   smtrx(1,3,i)=smtrx(1,3,j) 
        smtrx(2,1,i)=smtrx(2,1,j);  smtrx(2,2,i)=smtrx(2,2,j);   smtrx(2,3,i)=smtrx(2,3,j) 
        smtrx(3,1,i)=smtrx(3,1,j);  smtrx(3,2,i)=smtrx(3,2,j);   smtrx(3,3,i)=smtrx(3,3,j) 
 
        trm(1,1)=smtrx(1,1,i);   trm(1,2)=smtrx(1,2,i);   trm(1,3)=smtrx(1,3,i)
        trm(2,1)=smtrx(2,1,i);   trm(2,2)=smtrx(2,2,i);   trm(2,3)=smtrx(2,3,i)
        trm(3,1)=smtrx(3,1,i);   trm(3,2)=smtrx(3,2,i);   trm(3,3)=smtrx(3,3,i)

	trmt(1,1)=trm(1,1);  trmt(1,2)=trm(2,1);   trmt(1,3)=trm(3,1) 
	trmt(2,1)=trm(1,2);  trmt(2,2)=trm(2,2);   trmt(2,3)=trm(3,2) 
	trmt(3,1)=trm(1,3);  trmt(3,2)=trm(2,3);   trmt(3,3)=trm(3,3) 

	!omx2=matmul(trm, omx1)   
        omx2(1,1)=trm(1,1)*omx1(1,1)+trm(1,2)*omx1(2,1)+trm(1,3)*omx1(3,1)        
        omx2(2,1)=trm(2,1)*omx1(1,1)+trm(2,2)*omx1(2,1)+trm(2,3)*omx1(3,1)  
        omx2(3,1)=trm(3,1)*omx1(1,1)+trm(3,2)*omx1(2,1)+trm(3,3)*omx1(3,1)        

	!pmx=matmul(ammx,trmt)
        pmx(1,1)=ammx(1,1)*trmt(1,1)+ammx(1,2)*trmt(2,1)+ammx(1,3)*trmt(3,1)
        pmx(1,2)=ammx(1,1)*trmt(1,2)+ammx(1,2)*trmt(2,2)+ammx(1,3)*trmt(3,2)
        pmx(1,3)=ammx(1,1)*trmt(1,3)+ammx(1,2)*trmt(2,3)+ammx(1,3)*trmt(3,3)

        pmx(2,1)=ammx(2,1)*trmt(1,1)+ammx(2,2)*trmt(2,1)+ammx(2,3)*trmt(3,1)
        pmx(2,2)=ammx(2,1)*trmt(1,2)+ammx(2,2)*trmt(2,2)+ammx(2,3)*trmt(3,2)
        pmx(2,3)=ammx(2,1)*trmt(1,3)+ammx(2,2)*trmt(2,3)+ammx(2,3)*trmt(3,3)

        pmx(3,1)=ammx(3,1)*trmt(1,1)+ammx(3,2)*trmt(2,1)+ammx(3,3)*trmt(3,1)
        pmx(3,2)=ammx(3,1)*trmt(1,2)+ammx(3,2)*trmt(2,2)+ammx(3,3)*trmt(3,2)
        pmx(3,3)=ammx(3,1)*trmt(1,3)+ammx(3,2)*trmt(2,3)+ammx(3,3)*trmt(3,3)

	!fmx=matmul(trm,pmx)
        fmx(1,1)=trm(1,1)*pmx(1,1)+trm(1,2)*pmx(2,1)+trm(1,3)*pmx(3,1)  
        fmx(1,2)=trm(1,1)*pmx(1,2)+trm(1,2)*pmx(2,2)+trm(1,3)*pmx(3,2)  
        fmx(1,3)=trm(1,1)*pmx(1,3)+trm(1,2)*pmx(2,3)+trm(1,3)*pmx(3,3)  

        fmx(2,1)=trm(2,1)*pmx(1,1)+trm(2,2)*pmx(2,1)+trm(2,3)*pmx(3,1)  
        fmx(2,2)=trm(2,1)*pmx(1,2)+trm(2,2)*pmx(2,2)+trm(2,3)*pmx(3,2)  
        fmx(2,3)=trm(2,1)*pmx(1,3)+trm(2,2)*pmx(2,3)+trm(2,3)*pmx(3,3)  

        fmx(3,1)=trm(3,1)*pmx(1,1)+trm(3,2)*pmx(2,1)+trm(3,3)*pmx(3,1)  
        fmx(3,2)=trm(3,1)*pmx(1,2)+trm(3,2)*pmx(2,2)+trm(3,3)*pmx(3,2)  
        fmx(3,3)=trm(3,1)*pmx(1,3)+trm(3,2)*pmx(2,3)+trm(3,3)*pmx(3,3)  

	stma(19,i)=omx2(1,1);  stma(20,i)=omx2(2,1); stma(21,i)=omx2(3,1)   ! distance to host particle 
	stma(29,j)=fmx(1,1);    stma(30,j)=fmx(2,2); stma(31,j)=fmx(3,3)    !20180630    space momtuntum  

          stma(5,i)=stma(5,j) + stma(19,i) 
          stma(6,i)=stma(6,j) + stma(20,i) 
          stma(7,i)=stma(7,j) + stma(21,i)
 
          stma(10,i)=stma(10,j)  + stma(21,i)*stma(23,j) - stma(20,i)*stma(24,j) 
          stma(11,i)=stma(11,j)  - stma(21,i)*stma(22,j) + stma(19,i)*stma(24,j) 
          stma(12,i)=stma(12,j)  + stma(20,i)*stma(22,j) - stma(19,i)*stma(23,j)  
            
          stma(22,i) =  stma(22,j)  
          stma(23,i) =  stma(23,j)    
          stma(24,i) =  stma(24,j)  

          stma(25,i) =  stma(25,j)  
          stma(26,i) =  stma(26,j)    
          stma(27,i) =  stma(27,j) 

          do kk=29,34
          stma(kk,i) =  stma(kk,j)  
          enddo 

	!if (anint(stma(18,i)) /= anint(stma(18,i-1)))   then 
	!if (omx2(1,1)>0. .and. omx2(1,1)<10d-15)   omx2(1,1)=10d-15
	!if (omx2(2,1)>0. .and. omx2(2,1)<10d-15)   omx2(2,1)=10d-15
	!if (omx2(3,1)>0. .and. omx2(3,1)<10d-15)   omx2(3,1)=10d-15
	!if (omx2(1,1)<0. .and. omx2(1,1)>-10d-15)   omx2(1,1)=-10d-15
	!if (omx2(2,1)<0. .and. omx2(2,1)>-10d-15)   omx2(2,1)=-10d-15
	!if (omx2(3,1)<0. .and. omx2(3,1)>-10d-15)   omx2(3,1)=-10d-15

	!if (stma(19,i)>0. .and. stma(19,i)<10d-15)   stma(19,i)=10d-15
	!if (stma(20,i)>0. .and. stma(20,i)<10d-15)   stma(20,i)=10d-15
	!if (stma(21,i)>0. .and. stma(21,i)<10d-15)   stma(21,i)=10d-15
	!if (stma(19,i)<0. .and. stma(19,i)>-10d-15)   stma(19,i)=-10d-15
	!if (stma(20,i)<0. .and. stma(20,i)>-10d-15)   stma(20,i)=-10d-15
	!if (stma(21,i)<0. .and. stma(21,i)>-10d-15)   stma(21,i)=-10d-15
	!endif
                               
      endif    ! for if (anint(stma(18,i)) ...
      endif    ! for if (anint(stma(28,i))==0) ...

1425  enddo 

   ! !$OMP  END PARALLEL DO 
    
 

	end	
!=====================================
	subroutine field (stma, numcol, nump, tp0)
!=====================================
!        add accelration due to body force to each particle. 
!        Here, only gravitation is added
!input:  stma
!output: stma

       use face
       implicit double precision (a-h, o-z)  
       dimension stma(numcol, nump)

        gravity=-9.81  
        if (tp0<0.05) gravity=0.0

call omp_set_num_threads(8)
!$OMP  PARALLEL DO   &
!$OMP  DEFAULT(SHARED) &
!$OMP  private(i) 

	do 1605 i=1,nump	 
         if (anint(stma(28,i))==0)  then
         continue
         else
            if (anint(stma(28,i))==3)  then
            stma(15,i)= 0.
            stma(16,i)= 0.
            stma(17,i)= gravity
            endif
         endif 
                                      
1605	enddo 

!$OMP  END PARALLEL DO 

       end	

!
!=====================================
	subroutine search(stma, numcol, nump, numfull, xlower, xupper, ylower, yupper, zlower, zupper, cellp, lsg)
!=====================================
!       search potential neighbors for each particle
! input:   stma;
!          cellp--size of background grid;
!          xlower,xupper, ylower, yupper--compuational domain;  
! output:  mv--list of all particles in each background grid;
!          nb--list of number to each particle;
! internal:lsg--the location of all particles in the background
!          mvxv--the x-location for a particle
!          mvyv--the y-location for a particle    
!          iib--the first element in a row of nb     
!          ix, iy -- index for the neighboring cells for HA cell search, this idea is from SU Dong
!                    !! the orginial codes for HA-cell is from Dong Su, Y.J  and Fei polished it.  

	use face		
	implicit double precision (a-h, o-z)   
	dimension lsg(3,nump) 
	dimension stma(numcol, nump)             
	                                              
        nb=0  
	mv=0	                       
        inb = size(nb,2)   ! size of nb
                     
        nb(1,inb)=1


        do 1305 i=1,numfull	
        if (anint(stma(28,i))==0) then
        continue
        else	
                                 			 
	mvxv=int((stma(5,i)-xlower)/cellp)+3
	mvyv=int((stma(6,i)-ylower)/cellp)+3	
        mvzv=int((stma(7,i)-zlower)/cellp)+3 
	                               	                        
	lsg(1,i)=mvxv
	lsg(2,i)=mvyv   
        lsg(3,i)=mvzv
		
        jmv=0	      
	mv(mvxv,mvyv,mvzv,1)=mv(mvxv,mvyv,mvzv,1)+1          
        jmv=mv(mvxv,mvyv,mvzv,1)+1                           
        mv(mvxv,mvyv,mvzv, jmv)=i   
        endif          
1305	enddo	
          
       	do 1320 i=1,numfull
        if (anint(stma(28,i))==0) then
        continue
        else 
            
	 do 1322 ix=-1,1	     
	 do 1325 iy=-1,1
         do 1327 iz=-1,1
            
           if ( (ix+iy+iz)<0 .or. ((ix+iy+iz)==0 .and. (ix>iy))   .or. stma(28, i)==0 ) then  ! 2nd case, ix\=iy\=iz ( equals 0,-1,1 respectively) 
          ! if ( (ix+iy+iz)<0  .or. stma(28, i)==0 ) then  ! 2nd case, ix\=iy\=iz
           continue
           else
         
              if  (mv((lsg(1,i)+ix), (lsg(2,i)+iy),(lsg(3,i)+iz),1)==0)  then
              continue
              else
                            
	        do 1328 jjp=2, mv((lsg(1,i)+ix),(lsg(2,i)+iy),(lsg(3,i)+iz), 1)+1	
                k=  mv(lsg(1,i)+ix, lsg(2,i)+iy, lsg(3,i)+iz, jjp)

                dx=stma(5,i)-stma(5,k)
                dy=stma(6,i)-stma(6,k) 
                dz=stma(7,i)-stma(7,k)

                dist2=dx**2+dy**2+dz**2

               if (dist2<cellp**2) then  
                 if (ix==0 .and. iy==0 .and. iz==0 .and. (mv(lsg(1,i)+ix, lsg(2,i)+iy, lsg(3,i)+iz, jjp)>=i)) then
                ! if ((ix+iy+iz==0) .and. (mv(lsg(1,i)+ix, lsg(2,i)+iy, lsg(3,i)+iz, jjp)>=i)) then  
                  continue
                  else
                      
                   nb(1,nb(1,inb))= k   
                   nb(2,nb(1,inb))= i  

                   nb(1,inb)=nb(1,inb)+1                     
                      if ( nb(1,inb)>nump*200-10)   print *, 'NB is full.' 
                  endif                                            
               endif
                      
1328         enddo
           endif	  
         endif 

1327   enddo  					
1325   enddo   
1322   enddo  

      endif	
1320  enddo 
      
       end
!
!
!=====================================
	subroutine dtpint (stma, youngr, vmax, numcol, nump, radmax, dtp0, pmassminr, radmin, iedm, dtfree, numfull)		 
!=====================================
!       calculate the time step
! input:   stma
! output:  dtp0--time step for TDM algorithm
!          dtedm--time step for EDM algorithm
!          dtfree--time step for free flight particle
!          pmassmin--the minmum mass of all particle 
!          radmin--the radius of the smallest particle 
!          radmax--the radius of the biggest particle!         
!          pmassminr--ralative mass of two particles with mass of pmassmin
! internal:vmax--the velocity of the fastest particle !          
!          vsp--speed of each particle

        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump)  

	vmax=0.d0
	pmassmin=1.d2	   ! a big enough number
	radmin=stma(1,1) *100.		   ! a big enough number
	radmax=stma(1,1) * 0.0001		   ! a small enough number
        youngr=1.e9
	do 1110 i=1, numfull
          if (anint(stma(28, i))<=0) goto 1112
	  vsp=dsqrt(stma(10,i)*stma(10,i)+stma(11,i)*stma(11,i)+ stma(12,i)*stma(12,i) )     !vsp: velocity of single particle
	  if (vsp>vmax) vmax=vsp                                             
	  if (stma(3,i)<pmassmin) pmassmin=stma(3,i)	 	                                     
	  if (stma(1,i)<radmin) radmin=stma(1,i)  
	  if (stma(1,i)>radmax) radmax=stma(1,i)
          if (stma(8,i)>youngr) youngr=stma(8,i)
     
1110    enddo	
1112    continue		
		
	if (vmax<1.d0) vmax=1.d0	
	pmassminr=pmassmin/2.  
        youngr=youngr/2.0  

        dtp0=0.5e-2 *2.87d0*(pmassminr*pmassminr/4/(0.5d0*radmin*youngr *youngr*vmax))**0.2   
	if (dtp0>1.d-2*radmin/vmax) dtp0=1.d-2*radmin/vmax

	end
!
!=====================================
	subroutine intstat(nump, stma, numcol)
!=====================================
! define/change the initial state of each particle in array stma
! input: stma
! putput: stma

	use face
        implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump) 

call omp_set_num_threads(8)
!$OMP  PARALLEL DO   &
!$OMP  DEFAULT(SHARED)  
			     				   
	do 1010 i=1, nump	
        if (anint(stma(28,i))==0)  then 
        continue
        else
	stma(10,i)=0.0d0
	  if (anint(stma(28,i))==3)     then
          stma(10,i)=sin(i*1.00)
          stma(11,i)=sin(i*2.0)
          stma(12,i)=sin(i*3.0)+0.1  ! initial velocity
          stma(22,i)=50.*sin(i*1.0)
          stma(23,i)=40.*sin(i*1.0)
          stma(24,i)=60.*sin(i*1.0)
          endif
       endif     
1010  enddo 

!$OMP  END PARALLEL DO        
      end
!
!
!
!=====================================
	subroutine dtp1 (stma, youngr, vmax, numcol, nump, dtp0, pmassminr, radmin, iedm, kedm, dtedm, tp0, tpminus1, dtfree, numfull)
!=====================================
! This subroutine is the same as subroutine dtpint, but simpiler
! cause  radmin, radmax, pmassminr and pmassmin are calcualted in dtpinf
! input:   stma
! output:  dtp0--time step for TDM algorithm
!          dtedm--time step for EDM algorithm
!          dtfree--time step for free flight particle        
! internal:vmax--the velocity of the fastest particle          
!          vsp--speed of each particle

        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump)  
                      
	vmax=1.d-2	  
	do 1910 i=1, numfull
        if (anint(stma(28,i))==0)  goto 1918                                                       	
	vsp=dsqrt(stma(10,i)*stma(10,i)+stma(11,i)*stma(11,i)+stma(12,i)*stma(12,i)) 
	if (vsp>vmax) vmax=vsp 
1910    enddo	
1918    continue
	 	
        dtp0=(5.d-3)*2.87d0*(pmassminr*pmassminr/4/(0.5d0*radmin*youngr*youngr*vmax))**0.2   !contact time
	if (dtp0>0.5d-2*radmin/vmax) dtp0=0.5d-2*radmin/vmax  

	end
!
!
!=========================
	subroutine binary (stma, numcol, nump, pai, paj,  alpha, youngr, dist, pi,  dtp0, kedm, dtedm, deltamax, ratio) 
!========================
!  This subroutine calculates the accelration (call TDM)
!  or velocity (call EDM) due to contact.
!  input: nb--list of neighbor
!         e, youngr, fricp, kedm, pi -- defined in the main program
!         dtedm-- time step for EDM 
!  outpur: stma
! internal: pai,paj--information for two colliding particle, to save memory        
!           alpha--angle of two particles to the cooridnates
!           dist--distance between two particles
!           coat-- thick of imagine coat
!    contactX & contactY -- cooridantes of the contact point

        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump)  	
	dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)
	dimension s1(nump),  s2(nump),s3(nump),s4(nump),s5(nump),s6(nump)
        real:: emx(3,3),fmx(3,3), gmx(3,3), hmx(3,3), for1(3,1), for2(3,1), tor1(3,1), tor2(3,1)
        real::  fori(3,1), forj(3,1),tmx(3,3)
                      
        deltamax=0.d-4     
        inb = size(nb,2)  

	do i=1, nump
	       s1(i)= stma(15,i)
	       s2(i)= stma(16,i)	
	       s3(i)= stma(17,i)
	       s4(i)= stma(25,i)
	       s5(i)= stma(26,i)
	       s6(i)= stma(27,i)
	end do
         call omp_set_num_threads(8)
         !$OMP  PARALLEL DO schedule(static) &
         !$OMP  DEFAULT(SHARED)&
         !$OMP  private(k0p, i,j,e,disqr,dist,kjp,ljp,logi,logj,jp,pai,paj,pak,pal,alpha,xy,gama,contactX)&
         !$OMP  private(contactY,contactZ,deltamax,pmassr,radr,delX,delY,delZ,tmx,fori,forj,tor1,tor2   )& 
         !$OMP  REDUCTION (+: s1,s2,s3,s4,s5,s6 )           
	        
        do 1530 k0p=1, nb(1,inb)               
        i=nb(1,k0p)
        j=nb(2,k0p)  
      
          if  (i==0 .and. j==0) goto 15532  

         e=0.125*(stma(14,i)+stma(14,j))    
         e=0.05
   
    	  disqr=((stma(5,i)-stma(5,j))**2 + (stma(6,i)-stma(6,j))**2 +(stma(7,i)-stma(7,j))**2)
	  dist=dsqrt(disqr) 
  

         if (anint(stma(18,i))==j .or. anint(stma(18,j))==i .or.   & 
          (anint(stma(18,i))==anint(stma(18,j)) .and. anint(stma(18,i))/=0))  dist=1.e10
     
           if (dist<(stma(1,i)+stma(1,j)) ) then
 
              kjp=anint(stma(18,i))
              ljp=anint(stma(18,j))   
           logi=0
           logj=0
           if (kjp/=0) logi=1
           if (ljp/=0) logj=1  
                                       
	   do 1511 jp=1,14    
	    pai(jp)=stma(jp,i)
	    paj(jp)=stma(jp,j)
              if (logi==1)  pak(jp)=stma(jp, kjp) 
              if (logj==1)  pal(jp)=stma(jp, ljp)
1511	   end do

	   do 1512 jp=18,24
	      pai(jp)=stma(jp,i)
	      paj(jp)=stma(jp,j)
              if (logi==1) pak(jp)=stma(jp,kjp) 
              if (logj==1) pal(jp)=stma(jp,ljp)
1512	   enddo
	 	
	   do 1514 jp=15,17
	      pai(jp)=0.
	      paj(jp)=0.
              if (logi==1)  pak(jp)=0.
              if (logj==1)  pal(jp)=0.
1514	   end do

	   do 1515 jp=25,27
	      pai(jp)=0.
	      paj(jp)=0.
              if (logi==1)  pak(jp)=0.
              if (logj==1)  pal(jp)=0.
1515	   end do

	   do 1516 jp=29,34
	    pai(jp)=stma(jp,i)
	    paj(jp)=stma(jp,j)
              if (logi==1)  pak(jp)=stma(jp, kjp) 
              if (logj==1)  pal(jp)=stma(jp, ljp)
1516	   end do

           pai(35)=1
           paj(35)=1

           if (abs(pai(5)-paj(5))<1.e-15)   pai(5)=pai(5)+2.d-15  
           alpha=datan((paj(6)-pai(6))/(paj(5)-pai(5))) 
           if (isnan(alpha)) print *, 'alpha is no a number in binary'
           xy=dsqrt( (paj(6)-pai(6))**2 +  (paj(5)-pai(5))**2 )
           if (xy<1.e-20)   xy=1.e-20
     
           if ((paj(5)-pai(5))<0.) then
           ialpha=1
           alpha=alpha+pi
           else

           ialpha=0
           endif
           if (alpha<0.) alpha=alpha+2.*pi  
           
           gama=datan( (paj(7)-pai(7))/xy ) 
          
           contactX=stma(5,i)+ (stma(5,j)-stma(5,i)) * stma(1,i) / (stma(1,i)+stma(1,j))     
           contactY=stma(6,i)+ (stma(6,j)-stma(6,i)) * stma(1,i) / (stma(1,i)+stma(1,j))   
           contactZ=stma(7,i)+ (stma(7,j)-stma(7,i)) * stma(1,i) / (stma(1,i)+stma(1,j))
  
           if ((pai(1)+paj(1)-dist)>deltamax) deltamax=pai(1)+paj(1)-dist   
           if (logi==0 .and. logj==0) then
            pmassr=(stma(20,i)*stma(20,j))/ (stma(20,i)+stma(20,j)) 
	    radr=(stma(19,i)*stma(19,j))/ (stma(19,i)+stma(19,j)) 
           endif 
    
           if (logi==1 .and. logj==0) then
            pmassr=(stma(20,kjp)*stma(20,j))/ (stma(20,kjp)+stma(20,j)) 
	    radr=(stma(19,kjp)*stma(19,j))/ (stma(19,kjp)+stma(19,j)) 
           endif 

           if (logi==0 .and. logj==1) then
            pmassr=(stma(20,i)*stma(20,ljp))/ (stma(20,i)+stma(20,ljp)) 
	    radr=(stma(19,i)*stma(19,ljp))/ (stma(19,i)+stma(19,ljp)) 
           endif 

           if (logi==1 .and. logj==1) then
            pmassr=(stma(20,kjp)*stma(20,ljp))/ (stma(20,kjp)+stma(20,ljp)) 
	    radr=(stma(19,kjp)*stma(19,ljp))/ (stma(19,kjp)+stma(19,ljp)) 
           endif 

            
           call tdm_model (pai, paj, alpha, gama, e, youngr, dist, pi, numcol, contactX, &
contactY, contactz, ialpha, i, j, fori, forj,  pmassr, radr)

           pai(15)=fori(1,1)
           pai(16)=fori(2,1)   ! Here pai(15-17), paj(15-17) are forces 
           pai(17)=fori(3,1)
           paj(15)=forj(1,1)
           paj(16)=forj(2,1)
           paj(17)=forj(3,1)

!	   do 1542 jp=15,17   
!	   if (logi==0) stma(jp,i)=stma(jp,i) + pai(jp)/pai(3) 
!	   if (logj==0) stma(jp,j)=stma(jp,j) + paj(jp)/paj(3) 
!	   if (logi==1) stma(jp,kjp)=stma(jp,kjp) + pai(jp)/pai(3)
!	   if (logj==1) stma(jp,ljp)=stma(jp,ljp) + paj(jp)/paj(3)  
!1542       enddo 

	 if (logi==0) then
	 s1(i)=s1(i)+pai(15)/pai(3)
	 s2(i)=s2(i)+pai(16)/pai(3)
	 s3(i)=s3(i)+pai(17)/pai(3)
	 endif

	 if (logj==0) then
	 s1(j)=s1(j)+paj(15)/paj(3)
	 s2(j)=s2(j)+paj(16)/paj(3)
	 s3(j)=s3(j)+paj(17)/paj(3)
	 endif

	 if (logi==1) then
	 s1(kjp)=s1(kjp)+pai(15)/pai(3)
	 s2(kjp)=s2(kjp)+pai(16)/pai(3)
	 s3(kjp)=s3(kjp)+pai(17)/pai(3)
	 endif

	if (logj==1) then
	s1(ljp)=s1(ljp)+paj(15)/paj(3)
	s2(ljp)=s2(ljp)+paj(16)/paj(3)
	s3(ljp)=s3(ljp)+paj(17)/paj(3)
	endif

              if (logi==1)  then           !&&&&              
              delX=contactX-stma(5,kjp)
              delY=contactY-stma(6,kjp)
              delZ=contactZ-stma(7,kjp)

               if (delX>0. .and. delX<1.e-10)  delX=1.e-10; if (delX<0. .and. delX>-1.e-10)  delX=-1.e-10
               if (delY>0. .and. delY<1.e-10)  delY=1.e-10; if (delY<0. .and. delY>-1.e-10)  delY=-1.e-10
               if (delZ>0. .and. delZ<1.e-10)  delZ=1.e-10; if (delZ<0. .and. delZ>-1.e-10)  delZ=-1.e-10

              tmx(1,1)=0.;     tmx(1,2)=-delZ;           tmx(1,3)=delY
              tmx(2,1)=delZ;  tmx(2,2)=0.;             tmx(2,3)=-delX
              tmx(3,1)=-delY;   tmx(3,2)=delX;          tmx(3,3)=0.
        
             ! Tor1=matmul(tmx,fori)
              tor1(1,1)= tmx(1,2)*fori(2,1) + tmx(1,3)*fori(3,1)
              tor1(2,1)=tmx(2,1)*fori(1,1) + tmx(2,3)*fori(3,1)
              tor1(3,1)=tmx(3,1)*fori(1,1) + tmx(3,2)*fori(2,1) 
             
              s4(kjp)=s4(kjp)+tor1(1,1)/pak(29)
 	      s5(kjp)=s5(kjp)+tor1(2,1)/pak(30)
	      s6(kjp)=s6(kjp)+tor1(3,1)/pak(31)
              !stma(25,kjp)=stma(25,kjp)+ tor1(1,1)/pak(29)      !20171111 
              !stma(26,kjp)=stma(26,kjp)+ tor1(2,1)/pak(30)   !20171111 , 20180630
              !stma(27,kjp)=stma(27,kjp)+ tor1(3,1)/pak(31)      !20171111 , 20180630
              endif                 !&&&&

            		  if (logi==0)  then           !&&&&              
           		  delX=contactX-stma(5,i)
            		  delY=contactY-stma(6,i)
            		  delZ=contactZ-stma(7,i)

               if (delX>0. .and. delX<1.e-10)  delX=1.e-10; if (delX<0. .and. delX>-1.e-10)  delX=-1.e-10
               if (delY>0. .and. delY<1.e-10)  delY=1.e-10; if (delY<0. .and. delY>-1.e-10)  delY=-1.e-10
               if (delZ>0. .and. delZ<1.e-10)  delZ=1.e-10; if (delZ<0. .and. delZ>-1.e-10)  delZ=-1.e-10


           		  tmx(1,1)=0.;     tmx(1,2)=-delZ;    tmx(1,3)=delY
           		  tmx(2,1)=delZ;  tmx(2,2)=0.;        tmx(2,3)=-delX
           		  tmx(3,1)=-delY;   tmx(3,2)=delX;    tmx(3,3)=0.        

             		 ! Tor1=matmul(tmx, fori)
              		tor1(1,1)= tmx(1,2)*fori(2,1) + tmx(1,3)*fori(3,1)
                        tor1(2,1)=tmx(2,1)*fori(1,1) + tmx(2,3)*fori(3,1)
                        tor1(3,1)=tmx(3,1)*fori(1,1) + tmx(3,2)*fori(2,1) 

             	          s4(i)=s4(i)+tor1(1,1)/pai(29)
 	    	          s5(i)=s5(i)+tor1(2,1)/pai(30)
	    	          s6(i)=s6(i)+tor1(3,1)/pai(31)            		
                          !stma(25,i)=stma(25,i) + tor1(1,1)/pai(29)    ! 20171111
                    	  !stma(26,i)=stma(26,i) + tor1(2,1)/pai(30)    !  20171111, 20180630
              		  !stma(27,i)=stma(27,i) + tor1(3,1)/pai(31)   ! 20171111, 20180630
                     
              		  endif                 !&&&&

             		  if (logj==0)  then           !&&&&              
           		  delX=contactX-stma(5,j)
            		  delY=contactY-stma(6,j)
            		  delZ=contactZ-stma(7,j)

               if (delX>0. .and. delX<1.e-10)  delX=1.e-10; if (delX<0. .and. delX>-1.e-10)  delX=-1.e-10
               if (delY>0. .and. delY<1.e-10)  delY=1.e-10; if (delY<0. .and. delY>-1.e-10)  delY=-1.e-10
               if (delZ>0. .and. delZ<1.e-10)  delZ=1.e-10; if (delZ<0. .and. delZ>-1.e-10)  delZ=-1.e-10

           		  tmx(1,1)=0.;     tmx(1,2)=-delZ;           tmx(1,3)=delY
           		  tmx(2,1)=delZ;  tmx(2,2)=0.;               tmx(2,3)=-delX
           		  tmx(3,1)=-delY;   tmx(3,2)=delX;           tmx(3,3)=0.     

             		  !Tor2=matmul(tmx, forj)
              		tor2(1,1)= tmx(1,2)*forj(2,1) + tmx(1,3)*forj(3,1)
                        tor2(2,1)=tmx(2,1)*forj(1,1) + tmx(2,3)*forj(3,1)
                        tor2(3,1)=tmx(3,1)*forj(1,1) + tmx(3,2)*forj(2,1) 

                          s4(j)=s4(j)+tor2(1,1)/paj(29)
 	    	          s5(j)=s5(j)+tor2(2,1)/paj(30)
	    	          s6(j)=s6(j)+tor2(3,1)/paj(31)
                          !stma(25,j)=stma(25,j)+ tor2(1,1)/paj(29)   ! 20171111
                    	  !stma(26,j)=stma(26,j)+ tor2(2,1)/paj(30)   !  20171111 , 20180630
              	          !stma(27,j)=stma(27,j)+ tor2(3,1)/paj(31)   ! 20171111  , 20180630                        
              		  endif                 !&&&&

	      if (logj==1)  then          
              delX=contactX-stma(5,ljp)
              delY=contactY-stma(6,ljp)
              delZ=contactZ-stma(7,ljp)  

               if (delX>0. .and. delX<1.e-10)  delX=1.e-10; if (delX<0. .and. delX>-1.e-10)  delX=-1.e-10
               if (delY>0. .and. delY<1.e-10)  delY=1.e-10; if (delY<0. .and. delY>-1.e-10)  delY=-1.e-10
               if (delZ>0. .and. delZ<1.e-10)  delZ=1.e-10; if (delZ<0. .and. delZ>-1.e-10)  delZ=-1.e-10

              tmx(1,1)=0.;     tmx(1,2)=-delZ;           tmx(1,3)=delY
              tmx(2,1)=delZ;  tmx(2,2)=0.;             tmx(2,3)=-delX
              tmx(3,1)=-delY;   tmx(3,2)=delX;          tmx(3,3)=0.
            
             ! Tor2=matmul(tmx, forj)
              	tor2(1,1)= tmx(1,2)*forj(2,1) + tmx(1,3)*forj(3,1)
                tor2(2,1)=tmx(2,1)*forj(1,1) + tmx(2,3)*forj(3,1)
                tor2(3,1)=tmx(3,1)*forj(1,1) + tmx(3,2)*forj(2,1) 

               s4(ljp)=s4(ljp)+tor2(1,1)/pal(29)
 	       s5(ljp)=s5(ljp)+tor2(2,1)/pal(30)
 	       s6(ljp)=s6(ljp)+tor2(3,1)/pal(31)            
              !stma(25,ljp)=stma(25,ljp)+ tor2(1,1)/pal(29)    ! 20171111
              !stma(26,ljp)=stma(26,ljp)+ tor2(2,1)/pal(30)   !  20171111 , 20180630
              !stma(27,ljp)=stma(27,ljp)+ tor2(3,1)/pal(31)     ! 20171111         , 20180630   

              endif                 !&&&&	 
           	
          endif

15532     continue 
   
1530     enddo

        !$OMP  END PARALLEL DO 

	do i=1, nump
	        stma(15,i)= s1(i)  
		stma(16,i)= s2(i)	
                stma(17,i)= s3(i)  
		stma(25,i)= s4(i)
                stma(26,i)= s5(i)  
		stma(27,i)= s6(i)
	end do

      end	
!
!
!================================
	subroutine  ghostwall_plate (stma, numcol, nump, pai, paj, pak, pal, alpha,  youngr,  &
dist, pi,  dtp0, kedm, dtedm, nwp, wallp, tp0) 
!!================================!
!  judge wheather the particles contact with the ghostwall, if it is ture 
        use face
        use omp_lib
	implicit double precision (a-h, o-z)   

	dimension stma(numcol, nump) 
        dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)
        dimension nwp(18000,2) 
        dimension wallp(4,2)
        real::  fori(3,1), forj(3,1), tmx(3,3), tor2(3,1)

        numpar=size(nwp,1)
        numwall=size(nwp,2) 
            
       !if (tp0<0.15  .and. wallp(4,2)>-0.06) then 
       ! wallp(4,2)=-0.003-tp0*0.3
       ! else
       ! wallp(4,2)=-0.003-0.15*0.3
       ! endif 

        ftan=0.
        fnor=0.

        call omp_set_num_threads(8)
            !$OMP  PARALLEL DO schedule(static) &
            !$OMP  DEFAULT(SHARED)&
            !$OMP  private( ni,ii, i, j, kjp, e, dist,pai,paj, alpha,xy, ialpha, gama )&
            !$OMP  private( radr, pmassr, contactX, contactY, contactZ, delX, delY, delZ, fori, tmx, tor2 ) 

        do ni=1, numwall

        if  (nwp(numpar,ni)==0)  goto 1950
        do 1955 ii=1, nwp(numpar,ni)
             
        i=nwp(ii, ni)     
             
        if (anint(STMA(28,I))==3)  then      !it is a real solid particle   1906     
    
           do 1940 j=1,4
	   pai(j)=stma(j,i)  
           paj(j)=pai(j)
1940	   enddo 

          pai(5)=stma(5,i) 
          pai(6)=stma(6,i) 
          pai(7)=stma(7,i) 

	   pai(8)=stma(8,i)  
           paj(8)=pai(8)
	   pai(9)=stma(9,i)  
           paj(9)=pai(9)

           do 1944 j=10,12	   
           paj(j)=0.
1944	   enddo 
       
	   pai(13)=stma(13,i)  
           paj(13)=0.
           pai(14)=stma(14,i)
           paj(14)=pai(14)
           e=0.25*pai(14) 
           e=0.05     
    
           do 1952 j=15,17  !accelartion 
           pai(j)=0. 
           paj(j)=0.
1952       enddo
 
           paj(18)=0  

           do 1947 j=22,24   !augular velocity
           paj(j)=0.
1947      enddo         

           do 1945 j=25,27   !augular accelartion
           pai(j)=0.
           paj(j)=0.        ! 20181022
1945      enddo

           do 1941 j=29,31
	   pai(j)=stma(j,i)  
           paj(j)=pai(j)
1941	   enddo 
    
        pai(35)=3.
        paj(35)=3.5
         
        paj(5)=pai(5)-2.*wallp(1,ni)*(wallp(1,ni)*pai(5)+wallp(2,ni)*pai(6)+wallp(3,ni)*pai(7)+wallp(4,ni))      &
               / (wallp(1,ni)*wallp(1,ni)+wallp(2,ni)*wallp(2,ni)+ wallp(3,ni)*wallp(3,ni))

        paj(6)=pai(6)-2.*wallp(2,ni)*(wallp(1,ni)*pai(5)+wallp(2,ni)*pai(6)+wallp(3,ni)*pai(7)+wallp(4,ni))      &
               / (wallp(1,ni)*wallp(1,ni)+wallp(2,ni)*wallp(2,ni)+ wallp(3,ni)*wallp(3,ni))

        paj(7)=pai(7)-2.*wallp(3,ni)*(wallp(1,ni)*pai(5)+wallp(2,ni)*pai(6)+wallp(3,ni)*pai(7)+wallp(4,ni))     &
               / (wallp(1,ni)*wallp(1,ni)+wallp(2,ni)*wallp(2,ni)+ wallp(3,ni)*wallp(3,ni))
    
          dist=( (pai(6)-paj(6))**2 + (pai(5)-paj(5))**2 + (pai(7)-paj(7))**2)**0.5  
 
           if (dist<2.*pai(1)) then   !1904
             
           if (abs(pai(5)-paj(5))<1.e-15)   pai(5)=pai(5)-2.d-15   
           alpha=datan((paj(6)-pai(6))/(paj(5)-pai(5))) 
           if (isnan(alpha)) print *,'error alpha is Not a number, plz check'
           xy=dsqrt( (paj(6)-pai(6))**2 +  (paj(5)-pai(5))**2 )
           if (xy<1.e-20)   xy=1.e-20
    
           if ((paj(5)-pai(5))<0.) then
           ialpha=1
           alpha=alpha+pi
           else
           ialpha=0
           endif

           if (alpha<0.) alpha=alpha+2*pi  
           
           gama=datan( (paj(7)-pai(7))/xy ) 
           
                contactX=0.5* (pai(5)+ paj(5))   ! contact point
                contactY=0.5* (pai(6)+ paj(6)) 
                contactZ=0.5* (pai(7)+ paj(7))  

           kjp=0         
           kjp=anint(stma(18,i)) 


           if (kjp==0) then
            pmassr= stma(20,i) *0.5
	    radr= stma(19,i) *0.5
           endif
    
           if (kjp/=0) then
            pmassr= stma(20,kjp) *0.5!
	    radr= stma(19,kjp) *0.5
!              
           endif           
           
           call tdm_model (pai, paj,  alpha, gama, e, youngr, dist, pi,  numcol,  &
            contactX, contactY, contactz, ialpha, i,j, fori, forj,  pmassr, radr)

           if (kjp==0) then

           delX=contactX-stma(5,i)  ! 20171212
           delY=contactY-stma(6,i)  ! 20171212
           delZ=contactZ-stma(7,i)  ! 20171212  

           if (delX>0. .and. delX<1.e-6)  delX=1.e-6; if (delX<0. .and. delX>-1.e-6)  delX=-1.e-6
           if (delY>0. .and. delY<1.e-6)  delY=1.e-6; if (delY<0. .and. delY>-1.e-6)  delY=-1.e-6
           if (delZ>0. .and. delZ<1.e-6)  delZ=1.e-6; if (delZ<0. .and. delZ>-1.e-6)  delZ=-1.e-6

           pai(15)=fori(1,1)/pai(3)
           pai(16)=fori(2,1)/pai(3)
           pai(17)=fori(3,1)/pai(3)   

           tmx(1,1)=0.;      tmx(1,2)=-delZ;         tmx(1,3)=delY
           tmx(2,1)=delZ;    tmx(2,2)=0.;           tmx(2,3)=-delX
           tmx(3,1)=-delY;   tmx(3,2)=delX;         tmx(3,3)=0.

         !  Tor2=matmul(tmx, fori)
            tor2(1,1)= tmx(1,2)*fori(2,1) + tmx(1,3)*fori(3,1)
            tor2(2,1)= tmx(2,1)*fori(1,1) + tmx(2,3)*fori(3,1)
            tor2(3,1)=tmx(3,1)*fori(1,1) + tmx(3,2)*fori(2,1) 

           pai(25)=tor2(1,1)/pai(29)    
           pai(26)=tor2(2,1)/pai(30)    
           pai(27)=tor2(3,1)/pai(31) 

                             
	   do 1946 j=15,17
	   stma(j,i)=pai(j)+stma(j,i)                                     
1946       enddo 

	   do 1948 j=25,27
             stma(j,i)=stma(j,i) + pai(J) 
1948       enddo   
       endif  ! for  (kjp==0) 

       if (kjp/=0) then
          delX=contactX-stma(5,kjp)  ! 20171212
          delY=contactY-stma(6,kjp)  ! 20171212
          delZ=contactZ-stma(7,kjp)  ! 20171212 

               if (delX>0. .and. delX<1.e-10)  delX=1.e-10; if (delX<0. .and. delX>-1.e-10)  delX=-1.e-10
               if (delY>0. .and. delY<1.e-10)  delY=1.e-10; if (delY<0. .and. delY>-1.e-10)  delY=-1.e-10
               if (delZ>0. .and. delZ<1.e-10)  delZ=1.e-10; if (delZ<0. .and. delZ>-1.e-10)  delZ=-1.e-10
  
          pai(15)=fori(1,1)/pai(3)    ! here pai is the host (main) particle 
          pai(16)=fori(2,1)/pai(3)
          pai(17)=fori(3,1)/pai(3) 

          tmx(1,1)=0.;      tmx(1,2)=-delZ;          tmx(1,3)=delY
          tmx(2,1)=delZ;    tmx(2,2)=0.;             tmx(2,3)=-delX
          tmx(3,1)=-delY;   tmx(3,2)=delX;           tmx(3,3)=0.

          !Tor2=matmul(tmx, fori)
            tor2(1,1)= tmx(1,2)*fori(2,1) + tmx(1,3)*fori(3,1)
            tor2(2,1)= tmx(2,1)*fori(1,1) + tmx(2,3)*fori(3,1)
            tor2(3,1)= tmx(3,1)*fori(1,1) + tmx(3,2)*fori(2,1) 

          pai(25)=tor2(1,1)/pai(29)    
          pai(26)=tor2(2,1)/pai(30)    
          pai(27)=tor2(3,1)/pai(31)     

	   do 1956 j=15,17
	   stma(j,kjp)=pai(j)+stma(j,kjp)                                     
1956       enddo 

	   do 1958 j=25,27
             stma(j,kjp)=stma(j,kjp)+pai(J) 
1958       enddo   

       endif   ! for  (kjp/=0)                
           
        endif            !1904
        endif             !1906

1955   enddo
1950  continue
       enddo    

       !$OMP  END PARALLEL DO                   
       end

!=====================================
subroutine outputp (stma, nump, numcol, dtp0, tpout, istep, tp0, vmax, iip, dudt, pressF, &
                    partF, FricF, numfull, deltamax, ratio, stma2, pi, wallc)	   	
!===============================
	use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump) 
	dimension stma2(numcol, nump)
        dimension outfull(numcol, nump)
        dimension outhalf(numcol, nump) 
        dimension wallc(4,1)
			 	 
	character * 17 fileout
	character * 17 file2out
        character * 17 file3out
        character * 17 file4out
        character * 17 time
	character * 6 nout0	   

	! save to computer 	
        iip=anint(tp0/tpout)   ! 20180415 by Fei wang, we rewrote the experssion of iip.
	if ((iip*tpout>=tp0-dtp0) .and. (iip*tpout<=tp0 + dtp0)) then   !20180416 reported by Feiwang and we changed the range
         
         stma2=stma  

         do i=1, nump        
         if (anint(stma(28,i))==0) goto 1477
         enddo    
1477     continue 

         if (iip/=0) then

         outfull=stma2    !all abt outfull and ourhalf was added in 20181121: from here
         outhalf=stma2
         do j=1, nump

         if (anint(stma(18,j))==0) then 
         if  (abs(stma2(5,j)-wallc(1,1))<1.e-12) then
         theta=pi*0.5
         else
         theta=atan((stma2(6,j)-wallc(2,1))/(stma2(5,j)-wallc(1,1)))
         ENDIF
 
         outfull(20,j)=theta
         outfull(2,j)= sin(theta) * stma2(10,j) + cos(theta) * stma2(11,j)
         outhalf(2,j)=outfull(2,j)
         outfull(3,j)= cos(theta) * stma2(10,j) - sin(theta) * stma2(11,j)
         outhalf(3,j)=outfull(3,j)
         outfull(4,j)= sqrt(stma2(10,j)*stma2(10,j) + stma2(11,j)*stma2(11,j) + stma2(12,j)*stma2(12,j))
         outhalf(4,j)=outfull(4,j)
      
        outfull(34,j) = abs((stma2(34,j)-theta) - (0.5*pi)*anint((stma2(34,j)-theta)/(0.5*pi))) 

         outfull(8,j)= abs (stma2(32,j) - (0.5*pi)*anint((stma2(32,j))/(0.5*pi))) 
          outfull(8,j)=abs(outfull(8,j)-0.25*pi)
         outhalf(8,j)=outfull(8,j)

         outfull(9,j)= abs (stma2(33,j) - (0.5*pi)* anint((stma2(33,j))/(0.5*pi))) 
         outfull(9,j)=abs(outfull(9,j)-0.25*pi)
         outhalf(9,j)=outfull(9,j)

         outfull(13,j)= abs ((stma2(34,j)-theta) - (0.5*pi)* anint((stma2(34,j)-theta)/(0.5*pi)))    
         outfull(13,j)=abs(outfull(13,j)-0.25*pi)
         outhalf(13,j)=outfull(13,j)

         outfull(14,j)=sqrt(outfull(8,j)**2+outfull(9,j)**2+outfull(13,j)**2)
         outhalf(14,j)=outfull(14,j)

         else
         jk=anint(stma(18,j))
         outfull(2,j)=outfull(2,jk)
         outfull(3,j)=outfull(3,jk)
         outfull(4,j)=outfull(4,jk)
         outfull(8,j)=outfull(8,jk)
         outfull(9,j)=outfull(9,jk)
         outfull(13,j)=outfull(13,jk)
         outfull(14,j)=outfull(14,jk)
         outfull(34,j)=outfull(34,jk)
         outfull(20,j)=outfull(20,jk)
 
         outhalf(2,j)=outfull(2,jk)
         outhalf(3,j)=outfull(3,jk)
         outhalf(4,j)=outfull(4,jk)
        outhalf(8,j)=outfull(8,jk)
         outhalf(9,j)=outfull(9,jk)
         outhalf(13,j)=outfull(13,jk)
         outhalf(14,j)=outfull(14,jk)
         outhalf(34,j)=outfull(34,jk)
         endif

          l=anint(stma2(18,j))
          if ((stma2(5,j)<0.06 .and. anint(stma2(18,j))==0) .or. (stma2(5,l)<0.06) ) then           
              do k=1, numcol
              outhalf(k,j)=0.00
              enddo
           endif
         enddo

	write (nout0, '(i6.6)') iip
	fileout='file'//nout0//'.dat'
	time='time'//nout0//'.dat'
        file3out='full'//nout0//'.dat'	   
        file4out='half'//nout0//'.dat'	

	open(71, file=file3out)
	write (71, 171)  outfull
	 171 format(35e16.8)
	close(71)

	open(72, file=file4out)
	write (72, 172)  outhalf
	 172 format(35e16.8)
	close(72)                          

	open(70, file=fileout)
	write (70, 170)  stma2
	 170 format(35e16.8)
	close(70)

    	open(73, file=time)
	write (73, 173) dtp0, tpout, tp0, real(iip), real(numfull)
	173 format(5e18.10)
        close(73)
        endif

	endif 

	!! this is for screen
	nsrn=100000   	! very nsrn steps, show the result in sreeen	 
	if (istep==nsrn) then
	istep=0
	tke=0.d0
        pte=0.d0
	omax=0.d0
        rote=0.0
	do i=1, nump	 
	tke=tke+0.5*stma(3,i)*(stma(10,i)*stma(10,i)+stma(11,i)*stma(11,i)+stma(12,i)*stma(12,i))
        pte=pte + stma(3,i)* stma(7,i) * 9.81 
        rote= rote+ 0.5*     &
     (stma(29,i)* stma(22,i)*stma(22,i)+stma(30,i)* stma(23,i)*stma(23,i) +stma(31,i)*stma(24,i)*stma(24,i) )
	if (dabs(stma(22,i))>omax) omax=dabs(stma(22,i))
         	
	enddo
	print*,'time=', tp0, 's    dt=', dtp0,'s'
	print*,'vmax=', vmax, 'ke=', tke,'rote',rote,'total e=', tke+pte+rote,'omax=',omax,'deltamax=',deltamax	   
        Print *, '    '        
	endif
	istep=istep+1 
	 
	end
!
!
!=====================================
      subroutine VelocityPE (stma, numcol, nump,  dtp0, kedm, dtedm, dtfree, iiv, numfull, smtrx, pi)
!===============================
        use face
	implicit double precision (a-h, o-z)  
	dimension stma(numcol, nump) , smtrx(3,3,nump), tmx(3,3), umx(3,1), vmx(3,1)	

        call omp_set_num_threads(8)
        !$OMP  PARALLEL DO schedule(static) &
        !$OMP  DEFAULT(SHARED)&
        !$OMP  private(  i,tmx,vmx,umx)           

	do i=1, numfull	  
        if (anint(stma(18,i))==0)    then ! goto 1420	 
        
	stma(10,i) = stma(10,i)+stma(15,i)*dtp0
	stma(11,i) = stma(11,i)+stma(16,i)*dtp0
	stma(12,i) = stma(12,i)+stma(17,i)*dtp0

 
	stma(5,i) = stma(5,i)+stma(10,i)*dtp0 
        stma(6,i) = stma(6,i)+stma(11,i)*dtp0
        stma(7,i) = stma(7,i)+stma(12,i)*dtp0

	stma(22,i) = stma(22,i)+stma(25,i)*dtp0   
	stma(23,i) = stma(23,i)+stma(26,i)*dtp0 
	stma(24,i) = stma(24,i)+stma(27,i)*dtp0 

        tmx(1,1)=smtrx(1,1,i);    tmx(1,2)=smtrx(1,2,i);   tmx(1,3)=smtrx(1,3,i)
        tmx(2,1)=smtrx(2,1,i);    tmx(2,2)=smtrx(2,2,i);   tmx(2,3)=smtrx(2,3,i)
        tmx(3,1)=smtrx(3,1,i);    tmx(3,2)=smtrx(3,2,i);   tmx(3,3)=smtrx(3,3,i)

      !  tmx(1,1)=1.;    tmx(1,2)=0.;   tmx(1,3)=0.
      !  tmx(2,1)=0.;    tmx(2,2)=1.;   tmx(2,3)=0.
      !  tmx(3,1)=0.;    tmx(3,2)=0.;   tmx(3,3)=1.

        vmx(1,1)=stma(22,i) ;  vmx(2,1)=stma(23,i);  vmx(3,1)=stma(24,i)

        !umx=matmul(tmx, vmx) 
        umx(1,1)=tmx(1,1)*vmx(1,1) + tmx(1,2)*vmx(2,1) + tmx(1,3)*vmx(3,1)
        umx(2,1)=tmx(2,1)*vmx(1,1) + tmx(2,2)*vmx(2,1) + tmx(2,3)*vmx(3,1)
        umx(3,1)=tmx(3,1)*vmx(1,1) + tmx(3,2)*vmx(2,1) + tmx(3,3)*vmx(3,1)

	stma(22,i) = umx(1,1)
	stma(23,i) = umx(2,1)
	stma(24,i) = umx(3,1)

	stma(32,i) = stma(32,i)+umx(1,1)*dtp0          ! corrected 2018-06-30 
	stma(33,i) = stma(33,i)+umx(2,1)*dtp0          ! corrected 2018-06-30 
	stma(34,i) = stma(34,i)+umx(3,1)*dtp0          ! corrected 2018-06-30 

        if (stma(32,i)<0.)  stma(32,i)=stma(32,i)+2.*pi
        if (stma(32,i)>(2.*pi))  stma(32,i)=stma(32,i)-2.*pi

        if (stma(33,i)<0.)  stma(33,i)=stma(33,i)+2.*pi
        if (stma(33,i)>(2.*pi))  stma(33,i)=stma(33,i)-2.*pi

        if (stma(34,i)<0.)  stma(34,i)=stma(34,i)+2.*pi
        if (stma(34,i)>(2.*pi))  stma(34,i)=stma(34,i)-2.*pi

        endif     
        enddo
 
        !$OMP  END PARALLEL DO 
	end
!!
!=====================================
	subroutine searchsize(xlower, xupper, ylower, yupper, zlower, zupper, radmax,  mvx, mvy, mvz, cellp)
!=====================================
	use face
	! integer, allocatable:: mv (:,:,:,:)	
	! mv is the grid for locate the particle for find their neighbours
	! this subroutine give the size of mv
	implicit double precision (a-h, o-z)  
  
	mvx=int((xupper-xlower)/cellp)+4
	mvy=int((yupper-ylower)/cellp)+4 
        mvz=int((zupper-zlower)/cellp)+4  
	allocate (mv(mvx,mvy,mvz,200))
	end

!
!================================
   subroutine tdm_model (pai, paj, alpha, gama, e, youngr, dist, pi, numcol, contactX,   &
contactY, contactz, ialpha, i,j, fori, forj, pmassr, radr)	   	
!===============================f
! TDM collison model 
! input: pai,paj--information of two colliding particles
!        dist--distance between two particles
!
	implicit double precision (a-h, o-z)  
	dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)
        real:: amx(3,3),bmx(3,3),cmx(3,3), dmx(3,3), rmxi(3,3), rmxj(3,3)
        real:: posi(3,1),posj(3,1),veli(3,1), velj(3,1), fori(3,1), forj(3,1), ppp(3,1), rotvi(3,1), rotvj(3,1)
        real:: posix(3,1),posjx(3,1),velix(3,1), veljx(3,1), forix(3,1), forjx(3,1), pppx(3,1), contact(3,1)
        real:: posiy(3,1),posjy(3,1)
        real:: armi(3,1), armj(3,1), roti(3,1), rotj(3,1), rotix(3,1), rotjx(3,1)

        !alpha=-alpha ;         gama=-gama
        caph=cos(alpha)
        saph=-sin(alpha)
        cgm=cos(gama)
        sgm=sin(-gama)

        posi(1,1)=pai(5);    posi(2,1)=pai(6);    posi(3,1)=pai(7)
        posj(1,1)=paj(5);   posj(2,1)=paj(6);     posj(3,1)=paj(7)   

        veli(1,1)=pai(10);        veli(2,1)=pai(11);        veli(3,1)=pai(12)
        velj(1,1)=paj(10);        velj(2,1)=paj(11);        velj(3,1)=paj(12)
        
        roti(1,1)=pai(22);   roti(2,1)=pai(23);   roti(3,1)=pai(24)     !20181020  angular velocities abt i in old system
        rotj(1,1)=paj(22);   rotj(2,1)=paj(23);   rotj(3,1)=paj(24)    !20181020

        fori=0.;        forj=0.;         forix=0.;          forjx=0.;       

         amx(1,1)=caph; amx(1,2)=-saph;  amx(1,3)=0.
         amx(2,1)=saph; amx(2,2)=caph;   amx(2,3)=0.
         amx(3,1)=0.;         amx(3,2)=0.;           amx(3,3)=1.

         bmx(1,1)=cgm;  bmx(1,2)=0.;           bmx(1,3)=-sgm 
         bmx(2,1)=0.;         bmx(2,2)=1.;           bmx(2,3)=0.
         bmx(3,1)=sgm;  bmx(3,2)=0.;           bmx(3,3)=cgm

        ! cmx=matmul(bmx, amx)
         cmx(1,1)=bmx(1,1)*amx(1,1)  
         cmx(1,2)=bmx(1,1)*amx(1,2)  
         cmx(1,3)=  bmx(1,3) 
         cmx(2,1)=  amx(2,1) 
         cmx(2,2)=  amx(2,2)  
         cmx(2,3)=  0.
         cmx(3,1)=bmx(3,1)*amx(1,1)   
         cmx(3,2)=bmx(3,1)*amx(1,2)  
         cmx(3,3)=  bmx(3,3) 
              

       ! velix=matmul(cmx,veli)
        velix(1,1)= cmx(1,1)*veli(1,1) +  cmx(1,2)*veli(2,1) +  cmx(1,3)*veli(3,1)
        velix(2,1)= cmx(2,1)*veli(1,1) +  cmx(2,2)*veli(2,1)  
        velix(3,1)= cmx(3,1)*veli(1,1) +  cmx(3,2)*veli(2,1) +  cmx(3,3)*veli(3,1)

       ! veljx=matmul(cmx,velj)
        veljx(1,1)= cmx(1,1)*velj(1,1) +  cmx(1,2)*velj(2,1) +  cmx(1,3)*velj(3,1)
        veljx(2,1)= cmx(2,1)*velj(1,1) +  cmx(2,2)*velj(2,1)  
        veljx(3,1)= cmx(3,1)*velj(1,1) +  cmx(3,2)*velj(2,1) +  cmx(3,3)*velj(3,1)

       ! rotix=matmul(cmx,roti)       !20181020    angular velocities abt i in new system
        rotix(1,1)= cmx(1,1)*roti(1,1) +  cmx(1,2)*roti(2,1) +  cmx(1,3)*roti(3,1)
        rotix(2,1)= cmx(2,1)*roti(1,1) +  cmx(2,2)*roti(2,1)  
        rotix(3,1)= cmx(3,1)*roti(1,1) +  cmx(3,2)*roti(2,1) +  cmx(3,3)*roti(3,1)

       !  rotjx=matmul(cmx,rotj)       !20181020
        rotjx(1,1)= cmx(1,1)*rotj(1,1) +  cmx(1,2)*rotj(2,1) +  cmx(1,3)*rotj(3,1)
        rotjx(2,1)= cmx(2,1)*rotj(1,1) +  cmx(2,2)*rotj(2,1)  
        rotjx(3,1)= cmx(3,1)*rotj(1,1) +  cmx(3,2)*rotj(2,1) +  cmx(3,3)*rotj(3,1)

        contact=posjx* pai(1)/(pai(1)+paj(1))    !20181020   contact point in arrayment 
        armi=contact-posix   !20181020
        armj=contact-posjx   !20181020

        rmxi(1,1)=0;            rmxi(1,2)=-rotix(3,1);   rmxi(1,3)=rotix(2,1)       !20181020
        rmxi(2,1)=rotix(3,1);   rmxi(2,2)=0.;            rmxi(2,3)=-rotix(1,1)      !20181020
        rmxi(3,1)=-rotix(2,1);  rmxi(3,2)=rotix(1,1) ;   rmxi(3,3)=0.               !20181020        

        rmxj(1,1)=0;            rmxj(1,2)=-rotjx(3,1);   rmxj(1,3)=rotjx(2,1)       !20181020
        rmxj(2,1)=rotjx(3,1);   rmxj(2,2)=0.;            rmxj(2,3)=-rotjx(1,1)      !20181020
        rmxj(3,1)=-rotjx(2,1);  rmxj(3,2)=rotjx(1,1) ;   rmxj(3,3)=0.               !20181020
        
       ! rotvi=matmul(rmxi,armi)       !20181020   ! velocity due rotation 
         rotvi(1,1)=rmxi(1,2)*armi(2,1) + rmxi(1,3)*armi(3,1)
         rotvi(2,1)=rmxi(2,1)*armi(1,1) + rmxi(2,3)*armi(3,1)
         rotvi(3,1)=rmxi(3,1)*armi(1,1)+ rmxi(3,2)*armi(2,1) 

       ! rotvj=matmul(rmxj,armj)       !20181020
         rotvj(1,1)=rmxj(1,2)*armj(2,1) + rmxj(1,3)*armj(3,1)
         rotvj(2,1)=rmxj(2,1)*armj(1,1) + rmxj(2,3)*armj(3,1)
         rotvj(3,1)=rmxj(3,1)*armj(1,1) + rmxj(3,2)*armj(2,1)
    
	pmassr=(pai(3)*paj(3))/ (pai(3)+paj(3))  
	radr=(pai(1)*paj(1))/ (pai(1)+paj(1)) 
     
        youngrr= 1/ (((1-pai(9)*pai(9))/pai(8)) + ((1-paj(9)*paj(9))/paj(8)))
        fricp=(pai(4)+paj(4))*0.5  
       
	coe1=4/3*youngrr*dsqrt(radr)
	coe2=-2*0.9129*(dlog(e)/dsqrt(dlog(e)*dlog(e)+pi*pi)) *dsqrt(pmassr*2*youngrr*radr**0.5)
		
	fnor=coe1*(pai(1)+paj(1)-dist)**1.5   + coe2*(velix(1,1)-veljx(1,1))*(pai(1)+paj(1)-dist)**0.5 
                   
	if (fnor<0) then 
        fori(1,1)=0. ;  fori(2,1)=0.; fori(3,1)=0.      
        forj(1,1)=0. ;  forj(2,1)=0.; forj(3,1)=0.   
        else  
         delv2=veljx(2,1)+rotvj(2,1) -velix(2,1)- rotvi(2,1)     !20181020
         delv3=veljx(3,1)+rotvj(3,1) -velix(3,1)- rotvi(3,1)     !20181020
        delvs=delv2*delv2+ delv3*delv3
        delv=sqrt(delvs)  
            
         if (delvs>1.) then
        ftan=fnor*fricp    
           else
         ftan=fnor * fricp  
         endif
  
        if  (delv>1.e-5) then
        scal2=abs(delv2/delv) 
        scal3=abs(delv3/delv) 
        else
        scal2=sqrt(2.0)/2
        scal3=scal2
        endif        

        forix(1,1)=-fnor 
        forjx(1,1)= fnor 

        if (delv2>0.) then   !20181022
        forix(2,1)=ftan*scal2 
        forjx(2,1)=-ftan*scal2 
        else 
        forix(2,1)=-ftan*scal2
        forjx(2,1)=ftan*scal2
        endif
        
         if (delv3>0.) then  !20181022
        forix(3,1)=ftan*scal3
        forjx(3,1)=-ftan*scal3
        else 
        forix(3,1)=-ftan*scal3
        forjx(3,1)=ftan*scal3
        endif
               
         amx(1,2)=-amx(1,2); amx(2,1)=-amx(2,1);
         bmx(3,1)=-bmx(3,1); bmx(1,3)=-bmx(1,3)

         !dmx=matmul(amx, bmx)
         dmx(1,1)=amx(1,1)*bmx(1,1) 
         dmx(1,2)=amx(1,2) 
         dmx(1,3)=amx(1,1)*bmx(1,3)  
         dmx(2,1)=amx(2,1)*bmx(1,1)  
         dmx(2,2)=amx(2,2) 
         dmx(2,3)=amx(2,1)*bmx(1,3) 
         dmx(3,1)=bmx(3,1)
         dmx(3,2)=0.
         dmx(3,3)=bmx(3,3)

        ! fori=matmul(dmx,forix)
         fori(1,1)=dmx(1,1)*forix(1,1) + dmx(1,2)*forix(2,1) + dmx(1,3)*forix(3,1)
         fori(2,1)=dmx(2,1)*forix(1,1) + dmx(2,2)*forix(2,1) + dmx(2,3)*forix(3,1)
         fori(3,1)=dmx(3,1)*forix(1,1) + dmx(3,3)*forix(3,1)

         !forj=matmul(dmx,forjx)
         forj(1,1)=dmx(1,1)*forjx(1,1) + dmx(1,2)*forjx(2,1) + dmx(1,3)*forjx(3,1)
         forj(2,1)=dmx(2,1)*forjx(1,1) + dmx(2,2)*forjx(2,1) + dmx(2,3)*forjx(3,1)
         forj(3,1)=dmx(3,1)*forjx(1,1) + dmx(3,3)*forjx(3,1)
          
        endif

	end	

!=====================================
	subroutine BoundaryCellFix_sphere (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nwS, walls, numcol, nump, tp0)
!=====================================
!       give the neighobring particles to any fixed boundary (such as wall) 
	! integer, allocatable:: mv (:,:,:,:)	
	! mv is the grid for locate the particle for find their neighbours
	! this subroutine give the size of mv
        ! k is the id of each   wall
        ! numwell-- number of walls in the codes
        ! (x-x_o)^2 + (y-y_o)^2 + (z-z_(o))^2=R2  
        ! there are 4 elements in each row of array walls (wall data), which means
        ! x_o, y_o, z_o, and R

	use face
	implicit double precision (a-h, o-z)  

        dimension stma(numcol, nump)
        dimension dis(8)
        dimension nwS(18000,1)
        dimension walls(4,1)
        numpar=size(nwS,1)
        numwall=size(nwS,2)                                                       !        print *, 'enter BCF'    
        walls(:,1)=(/0.06, 0.06, 0.06, 0.058/) 
       
        do ni=1,numwall
        nwS(numpar, ni)=0
        enddo		 
              
        do i=1, mvx                
        do j=1, mvy
        do k=1, mvz
                   
        do ni=1,numwall 
              
        x= (i+0.5 -3)* cellp + xlower
        y= (j+0.5 -3)* cellp + ylower
        z= (k+0.5 -3)* cellp + zlower
        node=1
                     
           do ii=-1,1,2
           do jj=-1,1,2 
           do kk=-1,1,2
                   
            xcorner=x+ii*cellp*2.
            ycorner=y+jj*cellp*2.
            zcorner=z+kk*cellp*2.     

             d= sqrt((walls(1,ni)-xcorner)**2 + (walls(2,ni)-ycorner)**2 +(walls(3,ni)-zcorner)**2)-walls(4,ni)
             dis(node)=d     
            node=node+1     
           enddo
           enddo 
           enddo
          
        dtd = minval(dis)*maxval(dis)
                                       
        if (dtd<=0.   .or. (abs(maxval(dis))<3.*cellp .and.  abs(minval(dis))<1.*cellp) ) then    
         if (mv(i,j,k,1) ==0)  goto 2530                           
         do ia=2, (mv(i,j,k,1)+1)    
         if (mv(i,j,k,ia)==0) goto 2532     
         nwS(numpar,ni)=nwS(numpar,ni)+1 
         nwS(nwS(numpar,ni),ni)=mv(i,j,k,ia)    

           if (nwS(numpar, ni)>(numpar-20)) then
  print *, 'arrangement NW_S is full. please change the size of NW_S or check the codes!',  nwS(numpar,ni), numpar  
           stop
           endif

2532    continue
         enddo
2530    continue
        endif  

        enddo
        enddo
        enddo

        enddo
	end
!
!================================
	subroutine  ghostwall_sphere (stma, numcol, nump, pai, paj, pak, pal, alpha,  youngr,  &
dist, pi,  dtp0, kedm, dtedm, nwS, walls, tp0) 
!!================================!
!  judge wheather the particles contact with the spherical ghostwall, if it is ture 

        use face
	implicit double precision (a-h, o-z)   

	dimension stma(numcol, nump) 
        dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)
        dimension nwS(18000,1) 
        dimension walls(4,1)
        real::  fori(3,1), forj(3,1), tmx(3,3), tor2(3,1)

        numpar=size(nwS,1)
        numwall=size(nwS,2) 
            
        ftan=0.
        fnor=0.

        call omp_set_num_threads(8)
        !$OMP  PARALLEL DO schedule(static) &
        !$OMP  DEFAULT(SHARED)&
        !$OMP  private(  ni ,i,j,kjp,e,dx,dy,dz, dc ,dist,pai,paj,xr,yr,zr,alpha,xy,ialpha,gama )&
        !$OMP  private(  radr,pmassr,contactX,contactY,contactZ,delX,delY,delZ, fori, tmx ,tor2 )

        do ni=1, numwall

        if  (nwS(numpar,ni)==0)  goto 1900
        do 1905 ii=1, nwS(numpar,ni)
             
        i=nwS(ii, ni)     
             
        if (anint(STMA(28,I))==3)  then      !it is a real solid particle   1906     
    
           do 1920 j=1,4
	   pai(j)=stma(j,i)  
           paj(j)=pai(j)
1920	   enddo 

          pai(5)=stma(5,i) 
          pai(6)=stma(6,i) 
          pai(7)=stma(7,i) 

	   pai(8)=stma(8,i)  
           paj(8)=pai(8)
	   pai(9)=stma(9,i)  
           paj(9)=pai(9)


           do 1924 j=10,12	   
           paj(j)=0.
1924	   enddo 
       
	   pai(13)=stma(13,i)  
           paj(13)=0.
           pai(14)=stma(14,i)
           paj(14)=pai(14)
           e=0.25*pai(14) 
           e=0.05     
    
           do 1922 j=15,17  !accelartion 
           pai(j)=0. 
           paj(j)=0.
1922       enddo
 
           paj(18)=0  

           do 1927 j=22,24   !augular velocity
           paj(j)=0.
1927      enddo
         

           do 1925 j=25,27   !augular accelartion
           pai(j)=0.
           paj(j)=0.        ! 20181022
1925      enddo

           do 1921 j=29,31
	   pai(j)=stma(j,i)  
           paj(j)=pai(j)
1921	   enddo 
         
          dc=sqrt ((pai(5)-walls(1,ni))**2+( pai(6)-walls(2,ni) )**2+(  pai(7)-walls(3,ni) )**2)
      
          XR= (pai(5)-walls(1,ni))/DC
          YR= (pai(6)-walls(2,ni))/DC
          ZR= (pai(7)-walls(3,ni))/DC

          paj(5)=  XR*(walls(4,ni) + PAI(1))   + walls(1,ni)
          paj(6)=  YR*(walls(4,ni) + PAI(1))   + walls(2,ni)
          paj(7)=  ZR*(walls(4,ni) + PAI(1))   + walls(3,ni)           
    
          dist=( (pai(6)-paj(6))**2 + (pai(5)-paj(5))**2 + (pai(7)-paj(7))**2  )**0.5  
                       
           if (dist<2.*pai(1)) then   !1904
             
           if (abs(pai(5)-paj(5))<1.e-15)   pai(5)=pai(5)-2.d-15   
           alpha=datan((paj(6)-pai(6))/(paj(5)-pai(5))) 
           if (isnan(alpha)) print *,'error alpha is Not a number, plz check'
           xy=dsqrt( (paj(6)-pai(6))**2 +  (paj(5)-pai(5))**2 )
           if (xy<1.e-20)   xy=1.e-20
    
           if ((paj(5)-pai(5))<0.) then
           ialpha=1
           alpha=alpha+pi
           else
           ialpha=0
           endif

           if (alpha<0.) alpha=alpha+2*pi  
           
           gama=datan( (paj(7)-pai(7))/xy ) 
           
                contactX=0.5* (pai(5)+ paj(5))   ! contact point
                contactY=0.5* (pai(6)+ paj(6)) 
                contactZ=0.5* (pai(7)+ paj(7))  

           kjp=0         
           kjp=anint(stma(18,i)) 


           if (kjp==0) then
            pmassr= stma(20,i) *0.5
	    radr= stma(19,i) *0.5
           endif
    
           if (kjp/=0) then
            pmassr= stma(20,kjp) *0.5 
	    radr= stma(19,kjp) *0.5 
           endif
                      
           call tdm_model (pai, paj,  alpha, gama, e, youngr, dist, pi,  numcol,  &
            contactX, contactY, contactz, ialpha, i,j, fori, forj, pmassr, radr)


           if (kjp==0) then

           delX=contactX-stma(5,i)  ! 20171212
           delY=contactY-stma(6,i)  ! 20171212
           delZ=contactZ-stma(7,i)  ! 20171212  

           if (delX>0. .and. delX<1.e-6)  delX=1.e-6; if (delX<0. .and. delX>-1.e-6)  delX=-1.e-6
           if (delY>0. .and. delY<1.e-6)  delY=1.e-6; if (delY<0. .and. delY>-1.e-6)  delY=-1.e-6
           if (delZ>0. .and. delZ<1.e-6)  delZ=1.e-6; if (delZ<0. .and. delZ>-1.e-6)  delZ=-1.e-6

           pai(15)=fori(1,1)/pai(3)
           pai(16)=fori(2,1)/pai(3)
           pai(17)=fori(3,1)/pai(3)   

           tmx(1,1)=0.;      tmx(1,2)=-delZ;         tmx(1,3)=delY
           tmx(2,1)=delZ;    tmx(2,2)=0.;           tmx(2,3)=-delX
           tmx(3,1)=-delY;   tmx(3,2)=delX;         tmx(3,3)=0.

           !Tor2=matmul(tmx, fori)
            tor2(1,1)= tmx(1,2)*fori(2,1) + tmx(1,3)*fori(3,1)
            tor2(2,1)= tmx(2,1)*fori(1,1) + tmx(2,3)*fori(3,1)
            tor2(3,1)=tmx(3,1)*fori(1,1) + tmx(3,2)*fori(2,1) 

           pai(25)=tor2(1,1)/pai(29)    
           pai(26)=tor2(2,1)/pai(30)    
           pai(27)=tor2(3,1)/pai(31) 
                             
	   do 1926 j=15,17
	   stma(j,i)=pai(j)+stma(j,i)                                     
1926       enddo 

	   do 1928 j=25,27
             stma(j,i)=stma(j,i) + pai(J)
1928       enddo   
       endif  ! for  (kjp==0) 

       if (kjp/=0) then
          delX=contactX-stma(5,kjp)  ! 20171212
          delY=contactY-stma(6,kjp)  ! 20171212
          delZ=contactZ-stma(7,kjp)  ! 20171212 

               if (delX>0. .and. delX<1.e-10)  delX=1.e-10; if (delX<0. .and. delX>-1.e-10)  delX=-1.e-10
               if (delY>0. .and. delY<1.e-10)  delY=1.e-10; if (delY<0. .and. delY>-1.e-10)  delY=-1.e-10
               if (delZ>0. .and. delZ<1.e-10)  delZ=1.e-10; if (delZ<0. .and. delZ>-1.e-10)  delZ=-1.e-10
  
          pai(15)=fori(1,1)/pai(3)    ! here pai is the host (main) particle 
          pai(16)=fori(2,1)/pai(3)
          pai(17)=fori(3,1)/pai(3) 

          tmx(1,1)=0.;      tmx(1,2)=-delZ;          tmx(1,3)=delY
          tmx(2,1)=delZ;    tmx(2,2)=0.;             tmx(2,3)=-delX
          tmx(3,1)=-delY;   tmx(3,2)=delX;           tmx(3,3)=0.

         ! Tor2=matmul(tmx, fori)
            tor2(1,1)= tmx(1,2)*fori(2,1) + tmx(1,3)*fori(3,1)
            tor2(2,1)= tmx(2,1)*fori(1,1) + tmx(2,3)*fori(3,1)
            tor2(3,1)=tmx(3,1)*fori(1,1) + tmx(3,2)*fori(2,1) 

          pai(25)=tor2(1,1)/pai(29)    
          pai(26)=tor2(2,1)/pai(30)    
          pai(27)=tor2(3,1)/pai(31)     

	   do 1936 j=15,17
	   stma(j,kjp)=pai(j)+stma(j,kjp)                                     
1936       enddo 

	   do 1938 j=25,27
             stma(j,kjp)=stma(j,kjp)+pai(J) 
1938       enddo   

       endif   ! for  (kjp/=0) 
               
           
        endif            !1904
        endif             !1906

1905   enddo

 1900  continue
       enddo  
    !$OMP  END PARALLEL DO
                   
       end

!!=====================================
    subroutine delepart (stma, numcol, nump, numfull)
!=====================================  
!     make the mass of all ghost partile infinite
!     move all empty rows and rows for baby particles to the end of array stma
!     input:stma; 
!     output:stma

      use face
      implicit double precision (a-h, o-z)        
      dimension stma(numcol, nump)

      do i=1, nump
        if (anint(stma(28,i))==0) then 
 
         do j=i+1, nump
	   do k=1, numcol
             stma(k,i)=stma(k,j)
           enddo
         enddo 
            do k=1, numcol
             stma(k,nump)=0.
            enddo   
       endif
       if (anint(stma(28,i))==-3)  stma(3,i)=1.d6
      enddo 

      goto 5030

      do i=1, nump        
        if   (anint(stma(28,i))==0) then
        numfull = i-1 
        goto 5030 
        endif
      enddo 
5030  continue   

     print  *, 'numfull', numfull          
     end

!=====================================
	subroutine xyminmax (stma, vmax, numcol, nump, xlower, xupper, ylower, yupper, zlower, zupper, cellp, radmax, numfull)		 
!=====================================
        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump)  

        
	xlower=1.e10
        xupper=-1.e10
        ylower=1.e10
        yupper=-1.e10
        zlower=1.e10
        zupper=-1.e10
 
	do 11103 i=1, numfull
          if (anint(stma(28, i))<=0) print *, 'error in xyminmax'	   
	  if (stma(5,i) < xlower)      xlower=stma(5,i)   
          if (stma(5,i) > xupper)      xupper=stma(5,i)  
	  if (stma(6,i) < ylower)      ylower=stma(6,i)   
          if (stma(6,i) > yupper)      yupper=stma(6,i) 
	  if (stma(7,i) < zlower)      zlower=stma(7,i)   
          if (stma(7,i) > zupper)      zupper=stma(7,i) 
11103   enddo		
  
       cellp=2.5d0*radmax    

	xlower= xlower-4.*cellp
        xupper= xupper+4.*cellp
        ylower= ylower-4.*cellp
        yupper= yupper+4.*cellp
        zlower= zlower-4.*cellp
        zupper= zupper+4.*cellp

   
        xlower= -3.*cellp
        xupper= 0.12
        ylower=-3.*cellp
        yupper= 0.12
        zlower= -3.*cellp
        zupper= 0.7
	 
         end 

!=====================================
	subroutine BoundaryCellFix_plate (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nwp, wallp, numcol, nump, tp0)
!=====================================
!       give the neighobring particles to any fixed boundary (such as wall) 
	! integer, allocatable:: mv (:,:,:,:)	
	! mv is the grid for locate the particle for find their neighbours
	! this subroutine give the size of mv
        ! k is the id of each   wall
        ! numwall-- number of plate walls in the codes
        ! (x-x_o)^2 + (y-y_o)^2 + (z-z_(o))^2=R2  
        ! there are 4 elements in each row of array wallp (wall data), which means
        ! Ax+By+Cz+D=0, and the four elements are A,B,C and D

	use face
	implicit double precision (a-h, o-z)  

        dimension stma(numcol, nump)
        dimension dis(8)
        dimension nwp(18000,2)
        dimension wallp(4,2)
        numpar=size(nwp,1)
        numwall=size(nwp,2)                                                      
        wallp(:,1)=(/0.0, 0.0, 1., -0.6/) 
        wallp(:,2)=(/0.0, 0.0, 1., -0.002/) 
           
        do ni=1,numwall
        nwP(numpar, ni)=0
        enddo		 
              
        do i=1, mvx                
        do j=1, mvy
        do k=1, mvz
                   
        do ni=1,numwall 
              
        x= (i+0.5 -3)* cellp + xlower    ! center of each cell
        y= (j+0.5 -3)* cellp + ylower
        z= (k+0.5 -3)* cellp + zlower
        node=1
                     
           do ii=-1,1,2
           do jj=-1,1,2 
           do kk=-1,1,2
                   
            xcorner=x+ii*cellp*2.
            ycorner=y+jj*cellp*2.
            zcorner=z+kk*cellp*2.     
         
             d= (wallp(1,ni)*xcorner+wallp(2,ni)*ycorner+wallp(3,ni)*zcorner+wallp(4,ni))  &
                / sqrt(wallp(1,ni)*wallp(1,ni)+wallp(2,ni)*wallp(2,ni)+wallp(3,ni)*wallp(3,ni))
             dis(node)=d     
            node=node+1     
           enddo
           enddo 
           enddo
          
        dtd = minval(dis)*maxval(dis)
                                       
        if (dtd<=0. ) then   
         if (mv(i,j,k,1) ==0)  goto 2511                           
         do ia=2, (mv(i,j,k,1)+1)    
         if (mv(i,j,k,ia)==0) goto 1551     
         nwP(numpar,ni)=nwP(numpar,ni)+1 
         nwP(nwP(numpar,ni),ni)=mv(i,j,k,ia)    

           if (nwP(numpar, ni)>(numpar-20)) then
           print *, 'NW_PLATE is full. please change the size of NW_PLATE or check the codes!',  nwP(numpar,ni), numpar  
           stop
           endif

1551   continue
         enddo
2511    continue
        endif 

        enddo

        enddo
        enddo
        enddo
	end
!

!!=====================================
    subroutine addbaby_dice (stma, numcol, nump, numfull, stma0)
!=====================================  
!     make the mass of all ghost partile infinite
!     move all empty rows to the end of array stma
!     input:stma; 
!     output:stma

      use face
      implicit double precision (a-h, o-z)        
      dimension stma(numcol, nump), stma0(numcol, nump)

      real :: a = 0.5 ! sqrt(3.)-1 
      real :: b=  0.5  !1/sqrt(3.) 
  
      do i=1, nump        
        if   (anint(stma(18,i))/=0 .or. anint(stma(28,i))==0 ) then
        numfull = i-1 
        goto 5040 
        endif
      enddo 
5040  continue  

      numfull2=numfull

        do j=1, numfull  
       if (anint(stma(28,j))==3) then 

       stma(1,j)=stma(1,j)*0.5 
       stma(3,j)=stma(3,j)  *1.4   
        stma(29,j)=stma(29,j)  *2.3 !  added 2018-06-30 
        stma(30,j)=stma(29,j)  *2.3
        stma(31,j)=stma(29,j)  *2.3

       do 5050 m=-1,1,1
       do 5051 n=-1,1,1
       do 5052 l=-1,1,1  
       if  (abs(m)==1  .or. abs(n)==1 .or. abs(l)==1) then  !EXCLUDING PARTICLES IN THE CENTER OF BODY
            
     	stma(1,numfull2+1) = stma(1,j)   
     	stma(2,numfull2+1) = stma(2,j)     	
	stma(3,numfull2+1) = stma(3,j)             
     	stma(4,numfull2+1) = stma(4,j)
     
        stma(5,numfull2+1) = stma(5,j)  +  m*b* stma(1,j)  
        stma(6,numfull2+1) = stma(6,j)  +  n*b* stma(1,j)  
        stma(7,numfull2+1) = stma(7,j)  +  l*b* stma(1,j)  
    
     	stma(8,numfull2+1) = stma(8,j)
     	stma(9,numfull2+1) = stma(9,j)
     	stma(10,numfull2+1) = stma(10,j)
     	stma(11,numfull2+1) = stma(11,j)
     	stma(12,numfull2+1) = stma(12,j)         
     	stma(13,numfull2+1) = stma(1,j)
        stma(14,numfull2+1) = stma(14,j)

     	do k=15,17
       	 stma(k,numfull2+1) = stma(k,j)
        enddo

     	stma(18,numfull2+1)=j
        stma(19,numfull2+1)=stma(5,numfull2+1)-stma(5,j)
        stma(20,numfull2+1)=stma(6,numfull2+1)-stma(6,j)
        stma(21,numfull2+1)=stma(7,numfull2+1)-stma(7,j)

        if (stma(19,numfull2+1)>0. .and. stma(19,numfull2+1)<10d-15)   stma(19,numfull2+1)=10d-15
        if (stma(20,numfull2+1)>0. .and. stma(20,numfull2+1)<10d-15)   stma(20,numfull2+1)=10d-15
        if (stma(21,numfull2+1)>0. .and. stma(21,numfull2+1)<10d-15)   stma(21,numfull2+1)=10d-15
        if (stma(19,numfull2+1)<0. .and. stma(19,numfull2+1)>-10d-15)   stma(19,numfull2+1)=-10d-15
        if (stma(20,numfull2+1)<0. .and. stma(20,numfull2+1)>-10d-15)   stma(20,numfull2+1)=-10d-15
        if (stma(21,numfull2+1)<0. .and. stma(21,numfull2+1)>-10d-15)   stma(21,numfull2+1)=-10d-15

    	 do k=22,27
       	 stma(k,numfull2+1) = stma(k,j)
       	 enddo
     	stma(28,numfull2+1)=3. 
      
  	 do k=29,35
       	 stma(k,numfull2+1) = stma(k,j)
       	 enddo

    	numfull2=numfull2+1
    !    ENDIF
        
        endif

5052       enddo
5051       enddo    	
5050       enddo       
12312 continue
       endif   
   
      enddo 
          print *, 'total number', numfull2, stma(5,141)
          stma0=stma

	 open(77, file='file_init.dat')
	 write (77, 177) stma0
	 177 format(35e16.8)
	 close(77)

      end


!=====================================
	subroutine BoundaryCellFix_cylinder (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nwc, wallc, numcol, nump, tp0)
!=====================================
!       give the neighobring particles to any fixed boundary (such as wall) 
	!  in this subroutine the axis of the cyliner should be parllal to one of the axix, 
        !  and the coresponding element in wallc is zero.
        ! 4 elements in wallc are x0, y0, z0 and radius
        
	use face
	implicit double precision (a-h, o-z)  

        dimension stma(numcol, nump)
        dimension dis(8)
        dimension nwc(18000,1)
        dimension wallc(4,1), wallc2(4,1)
        numpar=size(nwc,1)
        numwall=size(nwc,2)                                                       !        print *, 'enter BCF'    
        wallc(:,1)=(/0.06, 0.06, 2.e10, 0.058/)  

        wallc(4,1)=0.058-tp0*0.3
        if (wallc(4,1)<0.05)    wallc(4,1)=0.05      

         wallc2=wallc

       
        do ni=1,numwall
        nwc(numpar, ni)=0
        enddo		 
              
        do i=1, mvx                
        do j=1, mvy
        do k=1, mvz
                   
        do ni=1,numwall 
              
        x= (i+0.5 -3)* cellp + xlower
        y= (j+0.5 -3)* cellp + ylower
        z= (k+0.5 -3)* cellp + zlower

            if (wallc(1,ni)>1.e10)  wallc2(1,ni)=x
            if (wallc(2,ni)>1.e10)  wallc2(2,ni)=y
            if (wallc(3,ni)>1.e10)  wallc2(3,ni)=z
     
        d= sqrt((x-wallc2(1,ni))*(x-wallc2(1,ni))+(y-wallc2(2,ni))*(y-wallc2(2,ni))+(z-wallc2(3,ni))*(z-wallc2(3,ni)) )
        d= sqrt((x-wallc2(1,ni))*(x-wallc2(1,ni))+(y-wallc2(2,ni))*(y-wallc2(2,ni)))
        if ( abs(d-wallc2(4,ni))<4.*cellp   ) then
         if (mv(i,j,k,1) ==0)  goto 2510                           
         do ia=2, (mv(i,j,k,1)+1)    
         if (mv(i,j,k,ia)==0) goto 1520     
         nwc(numpar,ni)=nwc(numpar,ni)+1 
         nwc(nwc(numpar,ni),ni)=mv(i,j,k,ia)   
 

           if (nwc(numpar, ni)>(numpar-20)) then
           print *, 'arrangement NW_Cylinder is full. please change the size of NW_C!',  nwc(numpar,ni), numpar  
           stop
           endif

1520    continue
         enddo
2510    continue
        endif  

        enddo
        enddo
        enddo

        enddo
	end
!
!================================
	subroutine  ghostwall_cylinder (stma, numcol, nump, pai, paj, pak, pal, alpha,  youngr,  &
dist, pi,  dtp0, kedm, dtedm, nwc, wallc, tp0) 
!!================================!
!  judge wheather the particles contact with the spherical ghostwall, if it is ture 

        use face
	implicit double precision (a-h, o-z)   

	dimension stma(numcol, nump) 
        dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)
        dimension nwC(18000,1) 
        dimension wallc(4,1), wallc2(4,1)
        real::  fori(3,1), forj(3,1),  tmx(3,3), tor2(3,1)

        numpar=size(nwC,1)
        numwall=size(nwC,2) 
        wallc(:,1)=(/0.06, 0.06, 2.e10, 0.058/)  

        wallc(4,1)=0.058-tp0*0.3
        if (wallc(4,1)<0.05)        wallc(4,1)=0.05  

        wallc2=wallc  
            
        ftan=0.
        fnor=0.

        call omp_set_num_threads(8)
        !$OMP  PARALLEL DO schedule(static) &
        !$OMP  DEFAULT(SHARED)&
        !$OMP  private(  ni ,i,j,kjp,e,dx,dy,dz, dc,zoom,dist,pai,paj,xr,yr,zr,ialpha, gama)&         
        !$OMP  private(  radr,pmassr,contactX,contactY,contactZ,delX,delY,delZ, fori, tmx ,tor2 )  

        do ni=1, numwall

        if  (nwC(numpar,ni)==0)  goto 1970
        do 1975 ii=1, nwC(numpar,ni)
             
        i=nwC(ii, ni)  

           kjp=0         
           kjp=anint(stma(18,i))   
             
        if (anint(STMA(28,I))==3)  then      !it is a real solid particle       
    
           do 1960 j=1,4
	   pai(j)=stma(j,i)  
           paj(j)=pai(j)
1960	   enddo 

           pai(5)=stma(5,i) 
           pai(6)=stma(6,i) 
           pai(7)=stma(7,i) 

	   pai(8)=stma(8,i)  
           paj(8)=pai(8)
	   pai(9)=stma(9,i)  
           paj(9)=pai(9)

           if (abs(pai(5)-wallc(2,1)) < 1.e-14)  pai(5)=pai(5)+1.e-13
           angle=atan ((pai(6)-wallc(1,1))/(pai(5)-wallc(2,1)) )
           if ((pai(5)-wallc(2,1))<0.)   angle = angle + pi
           vela = 1.* sin(20.*tp0)
           paj(10)=vela * sin(angle)
           paj(11)=vela * cos(angle)
           paj(12)=0.

!           do 1964 j=10,12	   
!           paj(j)=0.
!1964	   enddo 
       
	   pai(13)=stma(13,i)  
           paj(13)=0.
           pai(14)=stma(14,i)
           paj(14)=pai(14)
           e=0.25*pai(14)  
    
           do 1962 j=15,17  !accelartion 
           pai(j)=0. 
           paj(j)=0.
1962       enddo
 
           paj(18)=0  

           do 1967 j=22,24   !augular velocity
           paj(j)=0.
1967      enddo
         

           do 1965 j=25,27   !augular accelation
           pai(j)=0.
           paj(j)=0.        ! 20181022
1965      enddo

           do 1961 j=29,31
	   pai(j)=stma(j,i)  
           paj(j)=pai(j)
1961	   enddo 

          pai(35)=2.
          paj(35)=2.5

            if (wallc(1,ni)>1.e10)  wallc2(1,ni)=pai(5) 
            if (wallc(2,ni)>1.e10)  wallc2(2,ni)=pai(6)
            if (wallc(3,ni)>1.e10)  wallc2(3,ni)=pai(7)  

          dx= pai(5)-wallc2(1,ni)    
          dy= pai(6)-wallc2(2,ni)
          dz= pai(7)-wallc2(3,ni)
 
         
          dc=sqrt(dx*dx+dy*dy+dz*dz)
              
           if (dc>wallc2(4,ni)) then
            print *, 'particle leaves the cylinder !',i, dc, wallc2(4,ni), pai(5),pai(6), pai(7), dz        
           endif

          XR= dx/DC
          YR= dy/DC
          ZR= dz/DC

           zoom=(wallc2(4,ni)+pai(1))/dc  

           paj(5)=dx*zoom+wallc2(1,ni)
           paj(6)=dy*zoom+wallc2(2,ni) 
           paj(7)=dz*zoom+wallc2(3,ni)         
              
           dist=sqrt((pai(5)-paj(5))**2 + (pai(6)-paj(6))**2 + (pai(7)-paj(7))**2  )  
                       
           if (dist<(2.*pai(1))  ) then   
             if (abs(pai(5)-paj(5))<1.e-13)   pai(5)=pai(5)-2.e-13   
             alpha=datan((paj(6)-pai(6))/(paj(5)-pai(5))) 
             if (isnan(alpha)) print *,'error alpha is Not a number, plz check cylinder'
             xy=dsqrt( (paj(6)-pai(6))**2 +  (paj(5)-pai(5))**2 )
             if (xy<1.e-13)   xy=1.e-13
    
             if ((paj(5)-pai(5))<0.) then
             ialpha=1
             alpha=alpha+pi
             else
             ialpha=0
             endif

             if (alpha<0.) alpha=alpha+2*pi  
           
             gama=datan( (paj(7)-pai(7))/xy) 
           
                contactX=0.5* (pai(5)+ paj(5))   ! contact point
                contactY=0.5* (pai(6)+ paj(6)) 
                contactZ=0.5* (pai(7)+ paj(7)) 

           kjp=0         
           kjp=anint(stma(18,i)) 

           if (kjp==0) then
            pmassr= stma(20,i) *0.5
	    radr= stma(19,i) *0.5
           endif
    
           if (kjp/=0) then
            pmassr= stma(20,kjp) *0.5!
	    radr= stma(19,kjp) *0.5
           endif
                      
           call tdm_model (pai, paj,  alpha, gama, e, youngr, dist, pi,  numcol,  &
            contactX, contactY, contactz, ialpha, i,j, fori, forj, pmassr, radr)

           if (kjp==0) then

           delX=contactX-stma(5,i)  ! 20171212
           delY=contactY-stma(6,i)  ! 20171212
           delZ=contactZ-stma(7,i)  ! 20171212  

           if (delX>0. .and. delX<1.e-6)  delX=1.e-6; if (delX<0. .and. delX>-1.e-6)  delX=-1.e-6
           if (delY>0. .and. delY<1.e-6)  delY=1.e-6; if (delY<0. .and. delY>-1.e-6)  delY=-1.e-6
           if (delZ>0. .and. delZ<1.e-6)  delZ=1.e-6; if (delZ<0. .and. delZ>-1.e-6)  delZ=-1.e-6

           pai(15)=fori(1,1)/pai(3)
           pai(16)=fori(2,1)/pai(3)
           pai(17)=fori(3,1)/pai(3)   

           tmx(1,1)=0.;      tmx(1,2)=-delZ;         tmx(1,3)=delY
           tmx(2,1)=delZ;    tmx(2,2)=0.;           tmx(2,3)=-delX
           tmx(3,1)=-delY;   tmx(3,2)=delX;         tmx(3,3)=0.

          ! Tor2=matmul(tmx, fori)
            tor2(1,1)= tmx(1,2)*fori(2,1) + tmx(1,3)*fori(3,1)
            tor2(2,1)= tmx(2,1)*fori(1,1) + tmx(2,3)*fori(3,1)
            tor2(3,1)= tmx(3,1)*fori(1,1) + tmx(3,2)*fori(2,1) 

           pai(25)=tor2(1,1)/pai(29)    
           pai(26)=tor2(2,1)/pai(30)    
           pai(27)=tor2(3,1)/pai(31) 

                             
	   do 1966 j=15,17
	   stma(j,i)=pai(j)+stma(j,i)                                     
1966       enddo 

	   do 1968 j=25,27
             stma(j,i)=stma(j,i) + pai(J)
1968       enddo   
       endif  ! for  (kjp==0) 

       if (kjp/=0) then
          delX=contactX-stma(5,kjp)  ! 20171212
          delY=contactY-stma(6,kjp)  ! 20171212
          delZ=contactZ-stma(7,kjp)  ! 20171212 

               if (delX>0. .and. delX<1.e-10)  delX=1.e-10; if (delX<0. .and. delX>-1.e-10)  delX=-1.e-10
               if (delY>0. .and. delY<1.e-10)  delY=1.e-10; if (delY<0. .and. delY>-1.e-10)  delY=-1.e-10
               if (delZ>0. .and. delZ<1.e-10)  delZ=1.e-10; if (delZ<0. .and. delZ>-1.e-10)  delZ=-1.e-10
  
          pai(15)=fori(1,1)/pai(3)    ! here pai is the host (main) particle 
          pai(16)=fori(2,1)/pai(3)
          pai(17)=fori(3,1)/pai(3) 

          tmx(1,1)=0.;      tmx(1,2)=-delZ;          tmx(1,3)=delY
          tmx(2,1)=delZ;    tmx(2,2)=0.;             tmx(2,3)=-delX
          tmx(3,1)=-delY;   tmx(3,2)=delX;           tmx(3,3)=0.

         ! Tor2=matmul(tmx, fori)
            tor2(1,1)= tmx(1,2)*fori(2,1) + tmx(1,3)*fori(3,1)
            tor2(2,1)= tmx(2,1)*fori(1,1) + tmx(2,3)*fori(3,1)
            tor2(3,1)=tmx(3,1)*fori(1,1)  + tmx(3,2)*fori(2,1)
 
          pai(25)=tor2(1,1)/pai(29)    
          pai(26)=tor2(2,1)/pai(30)    
          pai(27)=tor2(3,1)/pai(31)     

	   do 1976 j=15,17
	   stma(j,kjp)=pai(j)+stma(j,kjp)                                     
1976       enddo 

	   do 1978 j=25,27
           stma(j,kjp)=stma(j,kjp)+pai(J) 
1978       enddo   

       endif   ! for  (kjp/=0)                
           
        endif            !1904
        endif             !1906

1975   enddo

1970  continue
       enddo      

      !$OMP  END PARALLEL DO               
       end
!
!=====================================
	subroutine Search_cylinder (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nwc, wallc, numcol, nump, tp0)
!=====================================
!       give the neighobring particles to any fixed boundary (such as wall) 
	!  in this subroutine the axis of the cyliner should be parllal to one of the axix, 
        !  and the coresponding element in wallc is zero.
        ! 4 elements in wallc are x0, y0, z0 and radius
        
	use face
	implicit double precision (a-h, o-z)  

        dimension stma(numcol, nump)
        dimension dis(8)
        dimension nwc(18000,1)
        dimension wallc(4,1), wallc2(4,1)
        numpar=size(nwc,1)
        numwall=size(nwc,2)                                                       !        print *, 'enter BCF'    
        wallc(:,1)=(/0.06, 0.06, 2.e10, 0.058/)  

        wallc(4,1)=0.058-tp0*0.3
        if (wallc(4,1)<0.05)         wallc(4,1)=0.05    

        wallc2=wallc
       
        do ni=1,numwall
        nwc(numpar, ni)=0        	

       do i=1, nump

            if (wallc(1,ni)>1.e10)  wallc2(1,ni)=stma(5,i)
            if (wallc(2,ni)>1.e10)  wallc2(2,ni)=stma(6,i)
            if (wallc(3,ni)>1.e10)  wallc2(3,ni)=stma(7,i)

       ds = (stma(5,i)-wallc2(1,ni))*(stma(5,i)-wallc2(1,ni)) +(stma(6,i)-wallc2(2,ni)) &
          *(stma(6,i)-wallc2(2,ni))+(stma(7,i)-wallc2(3,ni))*(stma(7,i)-wallc2(3,ni))  

         if (abs(sqrt(ds)-wallc(4,ni))<4.*cellp  .and.  (stma(1,i)>1.e-20)) then
         nwc(numpar,ni)=nwc(numpar,ni)+1 
         nwc(nwc(numpar,ni),ni)=i
         endif
       enddo
      enddo
     end
!
!
!=====================================
	subroutine Search_Plate (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nwp, wallp, numcol, nump, tp0)
!=====================================
!       give the neighobring particles to any fixed boundary (such as wall) 
	!  in this subroutine the axis of the cyliner should be parllal to one of the axix, 
        !  and the coresponding element in wallc is zero.
        ! 4 elements in wallc are x0, y0, z0 and radius
        
	use face
	implicit double precision (a-h, o-z)  

        dimension stma(numcol, nump)
        dimension nwp(18000,2)
        dimension wallp(4,2)
        numpar=size(nwp,1)
        numwall=size(nwp,2)                                                      
        wallp(:,1)=(/0.0, 0.0, 1., -0.5/) 
        wallp(:,2)=(/0.0, 0.0, 1., -0.002/) 
            
        do ni=1,numwall
        nwp(numpar, ni)=0        	

       do i=1, nump 

        if (  abs(stma(7,i)+wallp(4,ni)) <(4.*cellp) ) then
         nwp(numpar,ni)=nwp(numpar,ni)+1 
         nwp(nwP(numpar,ni),ni)=i
       endif
       enddo

      enddo

     end
! 
!=====================================
	subroutine Search_sphere (stma, xlower, xupper, ylower, yupper, zlower, zupper, radmax,  &
 mvx, mvy, mvz, cellp, nws, walls, numcol, nump, tp0)
!=====================================
!       give the neighobring particles to any fixed boundary (such as wall) 
	!  in this subroutine the axis of the cyliner should be parllal to one of the axix, 
        !  and the coresponding element in wallc is zero.
        ! 4 elements in wallc are x0, y0, z0 and radiuswwxueshenxuemaomao
        
	use face
	implicit double precision (a-h, o-z)  

        dimension stma(numcol, nump) 
        dimension nws(18000,1)
        dimension walls(4,1)
        numpar=size(nws,1)
        numwall=size(nws,2)                                    !        print *, 'enter BCF'    
        walls(:,1)=(/0.06, 0.06, 0.06, 0.058/)  
         
      do ni=1,numwall
        nws(numpar, ni)=0        	

       do i=1, nump 
       ds = (stma(5,i)-walls(1,ni))*(stma(5,i)-walls(1,ni)) +(stma(6,i)-walls(2,ni)) &
          *(stma(6,i)-walls(2,ni))+(stma(7,i)-walls(3,ni))*(stma(7,i)-walls(3,ni))  

       if (abs(sqrt(ds)-walls(4,ni))<4.*cellp  .and.  (stma(1,i)>1.e-20)) then
         nws(numpar,ni)=nws(numpar,ni)+1 
         nws(nws(numpar,ni),ni)=i
       endif
      enddo

     enddo
    end
!


