      program main
      implicit real*8 (a-h,o-z)
      integer,parameter::imax=5,jmax=5,nxy=25,nb=1
      dimension bxnb(nb),bynb(nb),txnb(nb),tynb(nb),qnxy(nxy),qnb(nb)
      real*8::dx=1.0,dy=1.0

      integer i,j
      do i=1,nxy
      qnxy(i) = i/10.
      enddo
      
      do i=1,nb
      bxnb(i) = 1.1
      bynb(i) = 1.1
      txnb(i) = -0.1
      tynb(i) = 2.0
      qnb(i) = 10.4
      enddo


      call solve(imax,jmax,dx,dy,qnxy,nxy,nb,bxnb,bynb,
     .    txnb,tynb,qnb)
      end program

      subroutine solve(imax,jmax,dx,dy,qnxy,nxy,nb,bxnb,bynb,
     .    txnb,tynb,qnb)
cc    imax: the num of nodes in x direction
cc    jmax: the num of nodes in y direction
cc    dxnxy: dimension in x direction of a cell 
cc    dynxy: dimension in y direction of a cell 
cc    qnxy: value stored in nodes
cc    nxy: total num of nodes 

      implicit real*8 (a-h,o-z)
      integer imax,jmax,nxy,nb
      dimension bxnb(nb),bynb(nb),txnb(nb),tynb(nb),qnxy(nxy),qnb(nb)
      real*8 dx,dy


      integer i,j,bxmin,bxmax,bymin,bymax,p2,p3
cc    closetp1 must be dimensioned for using minloc
      integer,dimension(1)::closetp1
      integer,dimension(4)::surrps1,surrps2
      dimension csurrpsx(4),csurrpsy(4),distob(4)
      real*8 bx,by,cosval,distobtem

      do i=1, nb
      bx = bxnb(i)
      by = bynb(i)

      bxmin = floor(bx/dx)
      bxmax = bxmin + 1
      bymin = floor(by/dy)
      bymax = bymin + 1

      surrps1(1) = bymin*imax + bxmin + 1
      surrps1(2) = bymin*imax + bxmax + 1
      surrps1(3) = bymax*imax + bxmin + 1
      surrps1(4) = bymax*imax + bxmax + 1

      csurrpsx(1) = bxmin * dx
      csurrpsx(2) = bxmax * dx
      csurrpsx(3) = bxmin * dx
      csurrpsx(4) = bxmax * dx

      csurrpsy(1) = bymin * dy
      csurrpsy(2) = bymin * dy
      csurrpsy(3) = bymax * dy
      csurrpsy(4) = bymax * dy
cc    find the closet point  
      do j=1, 4
cc    caculate angle between boundary point and every nodes
      tx2 = csurrpsx(j) - bx
      ty2 = csurrpsy(j) - by
      distobtem = sqrt(tx2**2 + ty2**2)

      cosval = (txnb(i)*tx2 + tynb(i)*ty2)/
     . (distobtem*sqrt(txnb(i)**2+tynb(i)**2))

      write(*,*) 'cos is:',cosval

      if(cosval .lt. 0) then
cc    if node locate inner solid then set dis to maximum          
          distob(j) = sqrt(dx**2+dy**2) + 1
      else 
          distob(j) = distobtem
      endif
      enddo
      closetp1 = minloc(distob)

      write(*,100) closetp1
      enddo

          
  100 format('result is',I10)   
      end subroutine

