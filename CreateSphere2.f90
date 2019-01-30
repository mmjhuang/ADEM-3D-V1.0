     integer iseed, i, j, n, k, m, l
     real    rnd, r, rho, pi
     real data1(35,30000) 

     iseed = 425001

     r=5.e-3
     rho=2524
     din= r *1.0  
     dil= r *1.0  
     n=1
     pi=2.*acos(0.)

do k=1,42
do j=1,14
do i=1,14

     data1(5,n)=i*2.8*r   + R*2.8
     data1(6,n)=din  + R*2.8
     data1(7,n)=dil  + R*2.8

 !  if ( ((data1(5,n)-0.06)**2 + (data1(6,n)-0.06)**2 + (data1(7,n)-0.06)**2 )< 0.04**2) then  ! for sphere
    if (((data1(5,n)-0.06)**2+(data1(6,n)-0.06)**2<0.043**2) .and. (data1(7,n)>0.11) .and. (data1(7,n)<0.62)) then  !for cylinder

     data1(1,n)=r
     data1(2,n)=rho    
     data1(3,n)=4*pi* data1(1,n)**3 *rho/3
     data1(4,n)=0.1   ! fricp(
     data1(5,n)=i* 2.8*r    + R*2.8
     data1(6,n)=din  + R*2.8	
    data1(7,n)=dil  + R*2.8

       data1(8,n)=6.3e10
       data1(9,n)=0.25   

       call random_number(rnd) 
       data1(10,n)=-1. + 2.*rnd
       data1(11,n)=0. 
       data1(12,n)=0.
       data1(13,n)=0.
       data1(14,n)=0.1

       do l=15,35
       data1(l,n)=0.
       enddo

       data1(19,n)=data1(1,n)
       data1(20,n)=data1(3,n)

       data1(28,n)=3.
       data1(29,n)=0.4*data1(3,n)*r*r   ! angular momentum 
       data1(30,n)= data1(29,n)
       data1(31,n)= data1(29,n)           

       data1(22,n)= (rnd-0.5)    !   angular velocity
       data1(23,n)= (rnd-0.5)
       data1(24,n)= (rnd-0.5)
       data1(35,n)=n*1.0000
        n=n+1
     endif 

 !    if ((data1(5,n)-0.006)**2+(data1(6,n)-0.006)**2+(data1(7,n)-0.006)**2>=0.005**2) then
 !    data1(5,n)=0.
 !    data1(6,n)=0.
 !    data1(7,n)=0.
 !    endif

end do
      din=din+2.8*r
end do
      dil=dil+2.8*r
      din= r *2.8
enddo

print *, "REAL" , N-1


open(20, file='file000000.dat')
 write (20, 110) data1
 110 format (35E16.8)
close(20)

PRINT *, "total", n-1	 
end
